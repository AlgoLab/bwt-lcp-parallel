#include <zlib.h>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <cmath>
#include <cstring>
#include <cassert>
#include <unistd.h>
#include <sys/mman.h>
#include "kseq.h"

#include <vector>
#include <algorithm>

#include <thread>
#include <atomic>
#include <mutex>

#include "base_types.hpp"
#include "common_types.hpp"

#include "encoding.hpp"
#include "io.hpp"


KSEQ_INIT(gzFile, gzread)


static unsigned concurrency_level;


template <typename T, unsigned int N>
inline
void fill_with_0(T* const v) {
  for (unsigned i = 0; i < N; ++i) v[i] = static_cast<T>(0);
}

struct base_env_t {
  const strlen_t k;
  const strnum_t m;
  const fnames_t& IX_names;
  lcp_writer_t& LCP;
  qp_readers_t Qp;
  strlen_t p;
  luint_t* const midpoints;
  std::vector<bool> E;
  std::mutex Emutex;

  base_env_t(
      const strlen_t _k,
      const strnum_t _m,
      const fnames_t& _IX_names,
      lcp_writer_t& _LCP
  )
    : k(_k),
      m(_m),
      IX_names(_IX_names),
      LCP(_LCP),
      p(0),
      midpoints{new luint_t[concurrency_level]},
      E(static_cast<luint_t>(_m)*(_k+1), false)
  {
    const luint_t tot = static_cast<luint_t>(m)*(k+1);
    const luint_t step = tot/concurrency_level + ((tot%concurrency_level)>(concurrency_level/2) ? 1 : 0);
    for (luint_t i = 0, pos = step; i<concurrency_level-1; ++i, pos += step) {
      midpoints[i] = pos;
    }
    midpoints[concurrency_level-1] = tot;

    E[tot-1] = true;
  }

  base_env_t(const base_env_t&) = delete;

  ~base_env_t() {
    delete [] midpoints;
  }

};

struct bfv_env_t {
  lcp_writer_t LCP;
  const strlen_t k;
  bool has_pseg_wider_than_1;

  bfv_env_t(const base_env_t& base_env)
    : LCP(base_env.LCP.clone()),
      k{base_env.k},
      has_pseg_wider_than_1(false)
  {}

  bfv_env_t(const bfv_env_t&) = delete;

  bfv_env_t(bfv_env_t&& other)
    : LCP(std::move(other.LCP)),
      k{other.k},
      has_pseg_wider_than_1(other.has_pseg_wider_than_1)
  {
  }

  ~bfv_env_t() {
  }

};

void compute_partial_BWT(
    const strlen_t k, const strnum_t m,
    const fnames_t& S_tpl,
    const fnames_t& B_tpl
) {
  const fnames_t N_tpl[] = {
    template_generator("supportfile-N0_%d.tmp", 0, ALPHSIZE),
    template_generator("supportfile-N1_%d.tmp", 0, ALPHSIZE),
  };
  { // B0 (and N0)
    copy_file(S_tpl[0], B_tpl[0]);
    multi_strnum_writer_t N0(N_tpl[0]);
    for (strnum_t i = 0; i < m; ++i) {
      N0.write(0, i);
    }
  }
  { // B1...Bk
    for (strlen_t l = 1; l <= k; ++l) {
      {
        enc_fasta_reader_t Bl1(B_tpl[l-1]);
        multi_strnum_reader_t Nl1(N_tpl[(l + 1) % 2]);
        multi_strnum_writer_t Nl(N_tpl[l % 2]);
        for (strnum_t i = 0; i < m; ++i) {
          Nl.write(Bl1.read(), Nl1.read());
        }
      }
      {
        enc_fasta_rnd_reader_t Sl(S_tpl[l], m);
        multi_strnum_reader_t Nl(N_tpl[l % 2]);
        enc_fasta_writer_t Bl(B_tpl[l]);
        for (strnum_t i = 0; i < m; ++i) {
          const strnum_t Nli = Nl.read();
          const enc_nucl_t SlNli = Sl.read(Nli);
          Bl.write(SlNli);
        }
      }
    }
  }
  // Cleaning
  remove_files(N_tpl[0]);
  remove_files(N_tpl[1]);
}

inline
bool compute_Lpos(
  luint_t* const Lpos,
  const luint_t b,
  const luint_t* Lcount
) {
  Lpos[0] = b;
  for (enc_nucl_t sigma = 1; sigma<ALPHSIZE; ++sigma) {
    Lpos[sigma] = Lpos[sigma-1] + Lcount[sigma-1];
  }
  bool this_has_pseg_wider_than_1 = false;
  for (enc_nucl_t sigma = 1; sigma<ALPHSIZE; ++sigma) {
    this_has_pseg_wider_than_1 |= (Lcount[sigma] > 1);
  }
  return this_has_pseg_wider_than_1;
}

void update_LCP(
  lcp_writer_t& LCP,
  const strlen_t p,
  const luint_t base_pos,
  const luint_t* Lcount
) {
  LCP.seek(base_pos);
  for (enc_nucl_t sigma = 1; sigma<ALPHSIZE; ++sigma) {
    if (Lcount[sigma]==0) continue;
    LCP.skip();
    for (luint_t rh = 1; rh < Lcount[sigma]; ++rh) {
      LCP.write(p);
    }
  }
}

void compute_Q1_l(
    const strlen_t k,
    const occ_countv_t& occv,
    const fnames_t& Q_tpls
) {
  for (strlen_t l = 1; l <= k; ++l) {
    plain_fasta_writer_t Ql(Q_tpls[l]);
    for (enc_nucl_t sigma = 0; sigma < ALPHSIZE; ++sigma) {
      for (luint_t i = 0; i < occv[l][sigma]; ++i) {
        Ql.write(sigma);
      }
    }
  }
}

#define MARGIN (64) //(__CHAR_BIT__*sizeof(unsigned long long))

void compute_interleave_LCP_segment(
                                    base_env_t& base_env,
                                    bfv_env_t& bfv_env,
                                    const luint_t segm_b, const luint_t segm_e,
                                    const strlen_t p,
                                    const std::vector<strnum_t>& Qppos,
                                    std::vector<luint_t>& end_pos,
                                    std::vector<std::vector<strnum_t>>& end_pos_Qp,
                                    unsigned idx
) {
  std::cout << ".";
  const luint_t tot = segm_e - segm_b;
  const strlen_t k = base_env.k;
  auto& E = base_env.E;

  auto& LCP = bfv_env.LCP;
  std::vector<strnum_t> Qppos_int = Qppos;

  mmap_reader_t<strlen_t> IX0(base_env.IX_names[(p-1)%2]);
  mmap_writer_t<strlen_t, ALPHSIZE> IX1(base_env.IX_names[p%2]);
  qp_readers_t Qp, Qpbis;
  for (strlen_t l = p+1; l <= k+1; ++l) {
    Qp.emplace_back(base_env.Qp[l-p-1].cloneAndSeek(Qppos[l-1]));
    Qpbis.emplace_back(base_env.Qp[l-p-1].cloneAndSeek(Qppos[l-1]));
  }
  luint_t Lcount[ALPHSIZE];
  fill_with_0<luint_t, ALPHSIZE>(Lcount);
  luint_t Lpos[ALPHSIZE];
  luint_t b = segm_b;
  bool has_pseg_wider_than_1 = false;
  IX0.seek(segm_b);
  IX1.seek(ALPHSIZE-1, segm_b);
  for (luint_t i = segm_b; i < segm_e; ++i) {
    const strlen_t IX0i = IX0.read();
    const enc_nucl_t c = (IX0i < p) ? 0 : Qp[IX0i-p].read();
    ///printf("p=%d i=%lu IX=%d c=%x L%x=%d\n", p, i, IX0i, c, c, IX0i);
    ++Qppos_int[IX0i];
    ++Lcount[c];
    if (E[i]) {
      ///printf("finished segment at pos %lu\n", i);
      if ((i+1 >= base_env.midpoints[idx-1]) & (i+1 < end_pos[idx] || end_pos[idx] == 0)) {
        end_pos[idx] = i+1;
        std::copy_n(Qppos_int.begin(), k+1, end_pos_Qp[idx].begin());
        ++idx;
      }
      if (b==i) {
        // Special case of a single-element segment
        if (IX0i >= p) Qpbis[IX0i-p].skip();  // advance Qpbis as well
        IX1.write(ALPHSIZE-1, IX0i);
        Lcount[c] = 0;
        b = i + 1;
      } else {
        const bool this_has_pseg_wider_than_1 = compute_Lpos(Lpos, b, Lcount);
        has_pseg_wider_than_1 |= this_has_pseg_wider_than_1;
        for (enc_nucl_t sigma = 0; sigma<ALPHSIZE; ++sigma) {
          if (Lcount[sigma] > 0) {
            const luint_t new_e = Lpos[sigma]+Lcount[sigma]-1-segm_b;
            ///printf("ending pos for sigma=%d in idx=%lu\n", sigma, Lpos[sigma]+Lcount[sigma]-1);
            if ((new_e < MARGIN) | (tot-new_e < MARGIN)) {
              std::lock_guard<std::mutex> lg(base_env.Emutex);
              E[segm_b+new_e] = true;
            } else {
              E[segm_b+new_e] = true;
            }
          }
        }
        if (this_has_pseg_wider_than_1) {
          update_LCP(LCP, p, Lpos[1], Lcount);
        }
        IX0.seek(b);
        IX1.seek(Lpos);
        for (luint_t r = b; r <= i; ++r) {
          const strlen_t IX0r = IX0.read();
          const enc_nucl_t cr = (IX0r < p) ? 0 : Qpbis[IX0r-p].read();
          ///printf("*p=%d i=%lu IX=%d c=%x L%x=%d\n", p, r, IX0r, cr, cr, IX0r);
          IX1.write(cr, IX0r);
        }
        fill_with_0<luint_t, ALPHSIZE>(Lcount);
        b = i + 1;
      }
    }
  }
  bfv_env.has_pseg_wider_than_1 = has_pseg_wider_than_1;
}

void worker_IX(const unsigned idx,
               base_env_t& base_env,
               bfv_env_t& bfv_env,
               const std::vector<luint_t>& this_end_pos,
               const std::vector<std::vector<strnum_t>>& this_end_pos_Qp,
               std::vector<luint_t>& end_pos,
               std::vector<std::vector<strnum_t>>& end_pos_Qp
               ) {
  const luint_t segm_b = this_end_pos[idx-1];
  const luint_t segm_e = this_end_pos[idx];
  compute_interleave_LCP_segment(base_env, bfv_env, segm_b, segm_e, base_env.p,
                                 this_end_pos_Qp[idx-1], end_pos, end_pos_Qp, idx);
}

void worker_Qp(std::atomic<luint_t>& shared_l,
               const strlen_t k, const strnum_t m,
               const occ_countv_t& occv,
               const fnames_t& B_tpl,
               const fnames_t& prevQ_tpl,
               const fnames_t& currQ_tpl
               ) {
  luint_t l;
  while ((l = shared_l.fetch_add(1)) <= k) {
    assert(l <= k);
    enc_fasta_reader_t Bl1(B_tpl[l]);
    plain_fasta_reader_t Ql1(prevQ_tpl[l-1]);
    create_empty_file<enc_nucl_t>(currQ_tpl[l], m);
    mmap_writer_t<enc_nucl_t, ALPHSIZE> Ql(currQ_tpl[l]);
    luint_t cumsum = 0;
    for (enc_nucl_t sigma = 0; sigma < ALPHSIZE; ++sigma) {
      Ql.seek(sigma, cumsum);
      cumsum += occv[l][sigma];
    }
    for (strnum_t i = 0; i < m; ++i) {
      const enc_nucl_t Ql1i = Ql1.read();
      const enc_nucl_t Bl1i = Bl1.read();
      Ql.write(Bl1i, Ql1i);
    }
    remove_file(prevQ_tpl[l-1]);
  }
}

void compute_interleave_LCP(
    const strlen_t k, const strnum_t m,
    const occ_countv_t& occv,
    const fnames_t& B_tpl,
    const fname_t& IXk_name
) {
  const fnames_t IX_names = {
    "supportfile-IX0.tmp",
    "supportfile-IX1.tmp",
  };
  const fnames_t Q_tpls[] = {
    template_generator("supportfile-Q0_%.3d.tmp", 0, k+1),
    template_generator("supportfile-Q1_%.3d.tmp", 0, k+1),
  };

  const luint_t tot = static_cast<luint_t>(m)*(k+1);
  const fnames_t L_tpl = template_generator("supportfile-L_%d.tmp", 0, ALPHSIZE);

  const fname_t LCP_name = "outfile-LCP.bin";
  { // IX0/1 and LCP
    strlen_writer_t IX0(IX_names[0]);
    for (strlen_t l = 0; l <= k; ++l) {
      for (strnum_t i = 0; i < m; ++i) {
        IX0.write(l);
      }
    }
    create_empty_file<strlen_t>(LCP_name, tot);
    create_empty_file<strlen_t>(IX_names[1], tot);
  }

  compute_Q1_l(k, occv, Q_tpls[0]);

  lcp_writer_t LCP(LCP_name);
  LCP.seek(0);
  LCP.write(static_cast<strlen_t>(-1));

  bool has_pseg_wider_than_1 = true;
  base_env_t base_env{k, m, IX_names, LCP};

  std::vector<luint_t> end_pos(concurrency_level+1, 0);
  end_pos[0] = 0;
  end_pos[1] = tot;
  std::vector<std::vector<strnum_t>> end_pos_Qp(concurrency_level+1, std::vector<strnum_t>(k+1, 0));
  std::vector<bfv_env_t> bfv_envs;
  for (unsigned idx = 0; idx < concurrency_level; ++idx)
    bfv_envs.emplace_back(base_env);

  while (has_pseg_wider_than_1) {
    has_pseg_wider_than_1 = false;
    ++base_env.p;
    std::cout << " - iteration " << static_cast<fuint_t>(base_env.p) << std::endl;
    const fnames_t& prevQ_tpl = Q_tpls[(base_env.p-1)%2];
    const fnames_t& currQ_tpl = Q_tpls[base_env.p%2];
    std::vector<luint_t> this_end_pos = end_pos;
    std::vector<std::vector<strnum_t>> this_end_pos_Qp = end_pos_Qp;
    {
      std::cout << "   - compute IX";
      for (strlen_t l = base_env.p+1; l <= k+1; ++l) {
        base_env.Qp.emplace_back(prevQ_tpl[l-1]);
      }
      std::vector<std::thread> threads;
      unsigned idx;
      for (idx = 1; idx <= concurrency_level && this_end_pos[idx-1] < tot; ++idx) {
        threads
          .emplace_back(worker_IX,
                        idx, std::ref(base_env), std::ref(bfv_envs[idx-1]),
                        std::cref(this_end_pos), std::cref(this_end_pos_Qp),
                        std::ref(end_pos), std::ref(end_pos_Qp)
                        );
      }
      for (auto& t: threads) {
        t.join();
      }
      for (unsigned j = 0; j < idx-1; ++j) {
        has_pseg_wider_than_1 |= bfv_envs[j].has_pseg_wider_than_1;
      }
      base_env.Qp.clear();
      std::cout << std::endl;
    }
    if (has_pseg_wider_than_1) {
      std::cout << "   - compute Qp" << std::endl;
      remove_file(currQ_tpl[base_env.p-1]);
      remove_file(prevQ_tpl[k]);
      std::vector<std::thread> threads;
      std::atomic<luint_t> shared_l(base_env.p+1);
      for (unsigned idx = 0; idx < concurrency_level; ++idx) {
        threads
          .emplace_back(worker_Qp,
                        std::ref(shared_l), k, m, std::cref(occv),
                        std::cref(B_tpl), std::cref(prevQ_tpl), std::cref(currQ_tpl)
                        );
      }
      for (auto& t: threads) {
        t.join();
      }
    } else {
      remove_file(currQ_tpl[base_env.p-1]);
      remove_file(prevQ_tpl[k]);
      for (strlen_t l = base_env.p+1; l <= k; ++l) {
        remove_file(prevQ_tpl[l-1]);
      }
    }
  }
  rename_file(IX_names[base_env.p%2], IXk_name);
}

void compute_BWT(
    const strlen_t k, const strnum_t m,
    const fnames_t& B_tpl,
    const fname_t& IXk_name
) {
  const strlen_t k1 = k+1;
  const luint_t tot = static_cast<luint_t>(m)*k1;
  strlen_reader_t IXk(IXk_name);
  std::vector<enc_fasta_reader_t> Bl;
  for (const fname_t& fname: B_tpl) {
    Bl.emplace_back(fname);
  }
  enc_fasta_writer_t B("outfile-BWT.bin");

  for (luint_t i = 0; i < tot; ++i) {
    const strlen_t l = IXk.read();
    const enc_nucl_t c = Bl[(l+1)%k1].read();
    B.write(c);
  }
}

int main(int argc, char *argv[]) {
  // Initialization
  gzFile fp;
  kseq_t *seq;

  // Error if no file is given as input
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <in.fasta> [<nthreads>]" << std::endl;
    return 1;
  }

  if (argc == 3) {
    int nthreads = std::atoi(argv[2]);
    if (nthreads <= 0) {
      std::cerr << "Usage: " << argv[0] << " <in.fasta> [<nthreads>]" << std::endl;
      std::cerr << "nthreads cannot be less than 1. Given: " << nthreads << std::endl;
      return 1;
    }
    concurrency_level = static_cast<unsigned>(nthreads);
  } else {
    unsigned ncores = std::thread::hardware_concurrency();
    concurrency_level = ncores == 0 ? 1 : std::min(ncores, 4u);
  }
  std::cout << "Using " << concurrency_level << " threads." << std::endl;

  fp = gzopen(argv[1], "r");
  seq = kseq_init(fp);

  strnum_t read_num = 0;
  int maxl = 0;

  {
    // int is required by library kseq
    int l = 0;

    // Calculate the length of the longest read
    while ((l = kseq_read(seq)) >= 0) {
      ++read_num;
      if (l > maxl) maxl = l;
    }

    if (maxl > 254) {
      std::cerr << "Input reads are too long! Current = "
                << maxl << ", max = 254"
                << std::endl;
      kseq_destroy(seq);
      gzclose(fp);
      return 1;
    }
  }

  const strlen_t read_max_length = static_cast<strlen_t>(maxl);

  std::cout << "longest read length = " << static_cast<fuint_t>(read_max_length) << std::endl;
  std::cout << "number of reads = " << read_num << std::endl;

  const fnames_t S_tpl= template_generator("supportfile-S%.3d.tmp", 0, read_max_length+1);
  const fnames_t B_tpl= template_generator("supportfile-B%.3d.tmp", 0, read_max_length+1);
  const fname_t IXk_name("supportfile-IXk.tmp");

  std::cout << "Preparing arrays S_l..." << std::endl;
  occ_countv_t occv(read_max_length+1);
  {
    multi_enc_fasta_writer_t S(S_tpl);
    kseq_rewind(seq);
    gzrewind(fp);

    int l = 0;
    while ((l = kseq_read(seq)) >= 0) {
      const strlen_t len = static_cast<strlen_t>(l);
      for (strlen_t i = 0; i <= len; ++i) {
        const enc_nucl_t c = encode(seq->seq.s[len-i]);
        S.write(i, c);
        ++occv[i][c];
      }
      for (strlen_t i = len+1; i <= read_max_length; ++i) {
        S.write(i, ALPHSIZE-1);
        ++occv[i][ALPHSIZE-1];
      }
    }

  }

  kseq_destroy(seq);
  gzclose(fp);

  std::cout << "Computing partial BWT..." << std::endl;
  compute_partial_BWT(read_max_length, read_num, S_tpl, B_tpl);
  remove_files(S_tpl);

  std::cout << "Computing interleave and LCP..." << std::endl;
  compute_interleave_LCP(read_max_length, read_num, occv, B_tpl, IXk_name);

  std::cout << "Computing BWT..." << std::endl;
  compute_BWT(read_max_length, read_num, B_tpl, IXk_name);

  std::cout << "Done." << std::endl;

  return 0;
}
