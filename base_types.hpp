#ifndef BASE_TYPES_HPP
#define BASE_TYPES_HPP

#include <stdint.h>
#include <cstring>
#include <string>
#include <vector>

#define ALPHSIZE 6

#ifndef PAGE_SIZE
#define PAGE_SIZE 4096
#endif

typedef uint8_t strlen_t;
typedef uint32_t strnum_t;
typedef uint32_t suint_t;
typedef uint64_t luint_t;
typedef uint_fast32_t fuint_t;

typedef char nucl_t;
typedef uint8_t enc_nucl_t;

typedef std::string fname_t;
typedef std::vector<fname_t> fnames_t;

#endif
