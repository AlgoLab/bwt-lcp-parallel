#ifndef IO_HPP
#define IO_HPP

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#include <algorithm>
#include <vector>
#include <string>

#include "base_types.hpp"

template <typename T>
class seq_reader_t {
 public:
  typedef T t;

 private:
  typedef seq_reader_t<T> my_t;

  FILE* f;
  T el;

 public:
  seq_reader_t(const fname_t& fname, const char* mode = "rb") {
    f = fopen(fname.c_str(), mode);
    assert(f != nullptr);
    int fno = fileno(f);
    posix_fadvise(fno, 0, 0, POSIX_FADV_SEQUENTIAL);
    posix_fadvise(fno, 0, 0, POSIX_FADV_WILLNEED);
    read();
  }

  seq_reader_t(my_t&& other)
      : f(other.f), el(other.el)
  {
    other.f = nullptr;
  }

  my_t& operator=(my_t&& other) {
    if (f != nullptr) fclose(f);
    f = other.f;
    el = other.el;
    other.f = nullptr;
    return *this;
  }

  ~seq_reader_t() {
    if (f != nullptr) fclose(f);
  }

  T read() {
    T res = el;
    if (sizeof(T) == sizeof(unsigned char)) {
      el = static_cast<T>(getc_unlocked(f));
    } else {
      fread_unlocked(&el, sizeof(T), 1, f);
    }
    return res;
  }

  bool has_next() {
    return !feof_unlocked(f);
  }
};

template <typename T>
class seq_writer_t {
 public:
  typedef T t;

 private:
  FILE* f;

  seq_writer_t(const seq_writer_t<T>&);
  seq_writer_t& operator=(const seq_writer_t<T>&);

 public:
  seq_writer_t(const fname_t& fname, const char* mode = "wb") {
    f = fopen(fname.c_str(), mode);
  }

  seq_writer_t(seq_writer_t<T>&& other)
      : f(other.f)
  {
    other.f = nullptr;
  }

  ~seq_writer_t() {
    if (f != nullptr) fclose(f);
  }

  void write(const T el) {
    if (sizeof(T) == sizeof(unsigned char)) {
      putc_unlocked(static_cast<int>(el), f);
    } else {
      fwrite_unlocked(&el, sizeof(T), 1, f);
    }
  }
};

// L **must** be a multiple of sysconf(_SC_PAGESIZE) (usually 4096)

template <typename T, suint_t N, fuint_t L=PAGE_SIZE*1024/sizeof(T)>
class mmap_writer_t {
 public:
  typedef T t;

 private:
  int fno;
  luint_t pos[N];
  luint_t off[N];
  t* bases[N];

 public:
  mmap_writer_t(const fname_t& fname, const luint_t* new_pos = nullptr) {
    fno = open(fname.c_str(), O_RDWR, 0);
    std::fill_n(pos, N, 0);
    std::fill_n(off, N, 0);
    std::fill_n(bases, N, nullptr);
    seek(new_pos != nullptr ? new_pos : pos);
  }

  ~mmap_writer_t() {
    for (suint_t i= 0; i < N; ++i) {
      if (bases[i] != nullptr) munmap(bases[i], L);
      bases[i] = nullptr;
    }
    if (fno >= 0) close(fno);
  }

  void seek(const suint_t i, const luint_t new_pos) {
    if (bases[i] == nullptr || new_pos < pos[i] || new_pos >= pos[i]+L) {
      if (bases[i] != nullptr) {
        munmap(bases[i], L);
      }
      pos[i] = new_pos - (new_pos % L);
      bases[i] = static_cast<t*>(
          mmap(nullptr, sizeof(t)*L, PROT_READ | PROT_WRITE, MAP_SHARED, fno,
               static_cast<__off_t>(pos[i]))
      );
      assert(bases[i] != MAP_FAILED);
    }
    off[i] = new_pos - pos[i];
  }

  void seek(const luint_t* new_pos) {
    for (suint_t i= 0; i < N; ++i) {
      seek(i, new_pos[i]);
    }
  }

  void write(const suint_t i, const T el) {
    assert(i<N);
    assert(bases[i] != nullptr);
    assert(off[i] < L);
    if (bases[i][off[i]] != el) bases[i][off[i]] = el;
    ++off[i];
    if (off[i]>=L) seek(i, pos[i]+off[i]);
  }

  void skip(const suint_t i) {
    ++off[i];
    if (off[i]>=L) seek(i, pos[i]+off[i]);
  }

};


template <typename T, fuint_t L=PAGE_SIZE*1024/sizeof(T)>
class mmap_swriter_t {

 private:
  int fno;
  luint_t pos;
  luint_t off;
  T* base;
  bool is_clone;

  mmap_swriter_t(const mmap_swriter_t<T, L>& other)
    : fno(other.fno), pos(0), off(0), base(nullptr), is_clone(true)
  {
    if (other.base != nullptr) seek(other.pos+other.off);
  }

public:
  mmap_swriter_t(const fname_t& fname)
    : fno(open(fname.c_str(), O_RDWR, 0)), pos(0), off(0), base(nullptr), is_clone(false)
  { }

  mmap_swriter_t(mmap_swriter_t<T, L>&& other)
    : fno(other.fno), pos(other.pos), off(other.off), base(other.base), is_clone(other.is_clone)
  {
    other.fno = -1;
    other.base = nullptr;
  }

  ~mmap_swriter_t() {
    if (base != nullptr) munmap(base, L);
    if ((fno >= 0) & !is_clone) close(fno);
  }

  void seek(const luint_t new_pos) {
    if ((base == nullptr) | (new_pos < pos) | (new_pos >= pos+L)) {
      if (base != nullptr) munmap(base, L);
      pos = new_pos - (new_pos % L);
      base = static_cast<T*>(
          mmap(nullptr, sizeof(T)*L, PROT_READ | PROT_WRITE, MAP_SHARED, fno,
               static_cast<__off_t>(pos))
      );
      assert(base != MAP_FAILED);
    }
    off = new_pos - pos;
  }

  void write(const T el) {
    assert(base != nullptr);
    assert(off < L);
    if (base[off] != el) base[off] = el;
    ++off;
    if (off>=L) seek(pos+off);
  }

  void skip() {
    ++off;
    if (off>=L) seek(pos+off);
  }

  mmap_swriter_t<T, L> clone() const {
    return mmap_swriter_t<T, L>(*this);
  }

};


// L **must** be a multiple of sysconf(_SC_PAGESIZE) (usually 4096)

template <typename T, fuint_t L=PAGE_SIZE*1024/sizeof(T)>
class mmap_reader_t {
 public:
  typedef T t;

 private:
  int fno;
  luint_t pos;
  luint_t off;
  t* base;
  const bool is_clone;

  mmap_reader_t(const int fno_, const luint_t initial_pos)
    : fno(fno_), pos(0), off(0), base(nullptr), is_clone(true)
  {
    seek(initial_pos);
  }

 public:
  mmap_reader_t(const fname_t& fname)
    : fno(open(fname.c_str(), O_RDONLY, 0)), pos(0), off(0), base(nullptr), is_clone(false)
  {
    assert(fno > 0);
  }

  mmap_reader_t(mmap_reader_t<T,L>&& other)
    : fno(other.fno), pos(other.pos), off(other.off), base(other.base), is_clone(other.is_clone)
  {
    other.fno = -1;
    other.base = nullptr;
  }

  mmap_reader_t(const mmap_reader_t<T,L>&) = delete;

  mmap_reader_t<T,L>& operator=(const mmap_reader_t<T,L>&) = delete;

  ~mmap_reader_t() {
    if (base != nullptr) munmap(base, L);
    if ((fno >= 0) & !is_clone) close(fno);
  }

  T read() {
    assert(base != nullptr);
    assert(off < L);
    const T res = base[off];
    ++off;
    if (off>=L) seek(pos+off);
    return res;
  }

  void skip() {
    ++off;
    if (off>=L) seek(pos+off);
  }

  void seek(const luint_t new_pos) {
    if (base == nullptr || new_pos < pos || new_pos >= pos+L) {
      if (base != nullptr) {
        munmap(base, L);
      }
      pos = new_pos - (new_pos % L);
      base = static_cast<t*>(
        mmap(nullptr, sizeof(t)*L, PROT_READ, MAP_PRIVATE, fno, static_cast<__off_t>(pos))
      );
      assert(base != MAP_FAILED);
    }
    off = new_pos - pos;
  }

  mmap_reader_t<T, L> clone() const {
    return mmap_reader_t<T, L>(fno, pos+off);
  }

  mmap_reader_t<T, L> cloneAndSeek(const luint_t new_pos) const {
    assert(fno > 0);
    return mmap_reader_t<T, L>(fno, new_pos);
  }
};

template <typename IO_NUCL_T, unsigned int NEL>
class optifasta_reader_t {
 public:
  typedef enc_nucl_t t;

 private:
  typedef optifasta_reader_t<IO_NUCL_T, NEL> my_t;
  typedef IO_NUCL_T io_nucl_t;

  static constexpr io_nucl_t SHIFT = (sizeof(io_nucl_t) * 8 / NEL) % (sizeof(io_nucl_t) * 8);
  static constexpr io_nucl_t MASK = ((1 << (sizeof(io_nucl_t) * 8 / NEL)) - 1) << ((NEL - 1) * SHIFT);
  static constexpr enc_nucl_t EOFEL = (1 << (sizeof(io_nucl_t) * 8 / NEL)) - 1;

  seq_reader_t<enc_nucl_t> reader;
  unsigned int nel;
  io_nucl_t io_el;
  enc_nucl_t result;
  bool iseof;

  void advance() {
    if (nel == 0 && !iseof) {
      if (reader.has_next()) {
        io_el = reader.read();
        nel = NEL;
      } else {
        iseof = true;
      }
    }
    if (!iseof) {
      result = (io_el & MASK) >> ((NEL - 1) * SHIFT);
      io_el <<= SHIFT;
      --nel;
      if (result == EOFEL) {
        nel = 0;
        iseof = true;
      }
    }
  }

 public:
  optifasta_reader_t(const fname_t& fname, const char* mode = "rb")
      : reader(fname, mode), nel(0),
        io_el(static_cast<io_nucl_t>(0)),
        result(static_cast<enc_nucl_t>(0)), iseof(false)
  {
    iseof = !reader.has_next();
    advance();
  }

  optifasta_reader_t(my_t&& other)
      : reader(std::move(other.reader)), nel(other.nel),
        io_el(other.io_el), result(other.result), iseof(other.iseof)
  {}

  my_t& operator=(my_t&& other) {
    reader = std::move(other.reader);
    nel = other.nel;
    io_el = other.io_el;
    result = other.result;
    iseof = other.iseof;
    return *this;
  }

  enc_nucl_t read() {
    enc_nucl_t res = result;
    advance();
    return res;
  }

  bool has_next() const {
    return !iseof;
  }

};

template <typename IO_NUCL_T, unsigned int NEL>
class optifasta_rnd_reader_t {
 public:
  typedef enc_nucl_t t;

 private:
  typedef IO_NUCL_T io_nucl_t;
  static constexpr io_nucl_t SHIFT = (sizeof(io_nucl_t) * 8 / NEL) % (sizeof(io_nucl_t) * 8);
  static constexpr io_nucl_t MASK = (1 << (sizeof(io_nucl_t) * 8 / NEL)) - 1;

  int fileno;
  IO_NUCL_T* fpl;
  const luint_t real_len;

 public:
  optifasta_rnd_reader_t(const fname_t& fname, const luint_t len)
      : fileno(open(fname.c_str(), O_RDONLY)),
        real_len(len % NEL == 0 ? len / NEL : (len+NEL-(len%NEL)) / NEL)
  {
    fpl = static_cast<IO_NUCL_T*>(mmap(nullptr, real_len, PROT_READ, MAP_PRIVATE, fileno, 0));
  }

  ~optifasta_rnd_reader_t() {
    munmap(fpl, real_len);
    close(fileno);
  }

  enc_nucl_t read(const luint_t pos) const {
    const luint_t index = pos / NEL;
    const luint_t offset = pos % 2;
    // printf("r %d %d %d %x\n", pos, index, offset, fpl[index]);
    // printf("c %d %x %d\n", SHIFT, MASK, SHIFT*(NEL-1-offset));
    return (fpl[index] >> (SHIFT*(NEL-1-offset))) & MASK;
  }

};

template <typename IO_NUCL_T, unsigned int NEL>
class optifasta_writer_t {
 public:
  typedef enc_nucl_t t;

 private:
  typedef IO_NUCL_T io_nucl_t;
  typedef optifasta_writer_t<IO_NUCL_T, NEL> my_t;

  static constexpr io_nucl_t SHIFT = (sizeof(io_nucl_t) * 8 / NEL) % (sizeof(io_nucl_t) * 8);
  static constexpr io_nucl_t MASK = (1 << (sizeof(io_nucl_t) * 8 / NEL)) - 1;

  seq_writer_t<enc_nucl_t> writer;
  unsigned int nel;
  io_nucl_t to_write;

 public:
  optifasta_writer_t(const fname_t& fname, const char* mode = "wb")
      : writer(fname, mode), nel(0), to_write(static_cast<io_nucl_t>(0))
  {}

  optifasta_writer_t(my_t&& other)
      : writer(std::move(other.writer)), nel(other.nel), to_write(other.to_write)
  {}

  ~optifasta_writer_t() {
    if (nel > 0) {
      while (nel < NEL) {
        to_write = static_cast<io_nucl_t>((to_write << SHIFT) | MASK);
        ++nel;
      }
      writer.write(to_write);
    }
  }

  void write(const enc_nucl_t& el) {
    to_write = static_cast<io_nucl_t>((to_write << SHIFT) | (el & MASK));
    ++nel;
    if (nel == NEL) {
      writer.write(to_write);
      nel = 0;
      to_write = static_cast<io_nucl_t>(0);
    }
  }

};


fnames_t
template_generator(const fname_t& tmpl, const suint_t imin, const suint_t imax) {
  fnames_t names;
  char *filepath;
  for (size_t i = imin; i < imax; ++i) {
    asprintf(&filepath, tmpl.c_str(), i);
    names.emplace_back(filepath);
    free(filepath);
  }
  return names;
}

template <typename base_writer_t>
class multi_seq_writer_t {
 private:
  const fnames_t fnames;
  const char* mode;

  std::vector<base_writer_t> writers;
 public:
  multi_seq_writer_t(const fnames_t& fnames_, const char* mode_ = "wb")
      : fnames(fnames_), mode(mode_)
  {
    open();
  }


  void write(const suint_t pos, typename base_writer_t::t el) {
    writers[pos].write(el);
  }

  void open() {
    writers.clear();
    for (const fname_t& fname: fnames) {
      writers.emplace_back(fname, mode);
    }
  }

  void close() {
    writers.clear();
  }
};

template <typename base_reader_t>
class multi_seq_reader_t {
 private:
  typedef multi_seq_reader_t<base_reader_t> my_t;

  const fnames_t fnames;
  const char* mode;
  fnames_t::const_iterator curr_name;

  base_reader_t reader;

  void advance() {
    while (!reader.has_next() && curr_name != fnames.end()) {
      reader = base_reader_t(*curr_name, mode);
      ++curr_name;
    }
  }

 public:
  multi_seq_reader_t(const fnames_t& fnames_, const char* mode_ = "rb")
      : fnames(fnames_), mode(mode_),
        curr_name(fnames.begin()), reader(*curr_name, mode_)
  {
    if (curr_name != fnames.end()) ++curr_name;
  }

  multi_seq_reader_t(my_t&& other)
      : fnames(other.fnames), mode(other.mode),
        curr_name(fnames.begin() + (other.curr_name - other.fnames.begin())),
        reader(std::move(other.reader))
  {
  }

  typename base_reader_t::t read() {
    advance();
    return reader.read();
  }

  bool has_next() {
    advance();
    return reader.has_next();
  }

};


void copy_file(const fname_t& from, const fname_t& to) {
  const std::size_t buf_sz = 32768;
  char* buf = new char[buf_sz];

  int infile = ::open(from.c_str(), O_RDONLY, 0);
  int outfile = ::open(to.c_str(), O_CREAT | O_WRONLY | O_TRUNC, 0666);

  ssize_t sz, sz_read=1, sz_write;
  while ((sz_read > 0) &&
         (sz_read = ::read(infile, buf, buf_sz)) > 0) {
    sz_write = 0;
    do {
      if ((sz = ::write(outfile, buf + sz_write,
                        static_cast<size_t>(sz_read - sz_write))) < 0) {
        sz_read = sz; // cause read loop termination
        break;        //  and error reported after closes
      }
      sz_write += sz;
    } while (sz_write < sz_read);
  }

  delete [] buf;
  ::close(infile);
  ::close(outfile);
}

void remove_file(const fname_t& fname) {
  ::remove(fname.c_str());
}

void remove_files(const fnames_t& fnames) {
  for (const fname_t& fname: fnames) {
    remove_file(fname);
  }
}

void rename_file(const fname_t& from, const fname_t& to) {
  ::rename(from.c_str(), to.c_str());
}

template <typename T>
void create_empty_file(const fname_t& fname, const luint_t n) {
  int fno = creat(fname.c_str(), 0666);
#ifdef FALLOC_FL_ZERO_RANGE
  fallocate(fno, FALLOC_FL_ZERO_RANGE, 0, static_cast<__off_t>(n*sizeof(T)));
#else
  fallocate(fno, 0, 0, static_cast<__off_t>(n*sizeof(T)));
#endif
  close(fno);
}

#endif
