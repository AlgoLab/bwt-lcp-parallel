#ifndef COMMON_TYPES_HPP
#define COMMON_TYPES_HPP

#include "base_types.hpp"
#include "io.hpp"

#include <array>
#include <vector>
#include <tuple>

typedef optifasta_writer_t<unsigned char, 2> enc_fasta_writer_t;
typedef optifasta_reader_t<unsigned char, 2> enc_fasta_reader_t;
typedef optifasta_rnd_reader_t<unsigned char, 2> enc_fasta_rnd_reader_t;

typedef seq_writer_t<enc_nucl_t> plain_fasta_writer_t;
typedef seq_reader_t<enc_nucl_t> plain_fasta_reader_t;

typedef multi_seq_writer_t<enc_fasta_writer_t> multi_enc_fasta_writer_t;

typedef seq_writer_t<strlen_t> strlen_writer_t;
typedef seq_writer_t<strnum_t> strnum_writer_t;
typedef seq_reader_t<strlen_t> strlen_reader_t;
typedef seq_reader_t<strnum_t> strnum_reader_t;

typedef multi_seq_writer_t<strnum_writer_t> multi_strnum_writer_t;
typedef multi_seq_reader_t<strnum_reader_t> multi_strnum_reader_t;

typedef std::array<luint_t, ALPHSIZE> occ_count_t;
typedef std::vector<occ_count_t> occ_countv_t;

typedef std::tuple<
  luint_t,   // b
  luint_t,   // e
  strnum_t *,  // Qppos
  strlen_t   // p
> segment_t;

typedef mmap_swriter_t<strlen_t> lcp_writer_t;
typedef mmap_reader_t<enc_nucl_t> q_reader_t;
typedef std::vector<q_reader_t> qp_readers_t;
typedef std::vector<qp_readers_t> qpl_readers_t;


#endif
