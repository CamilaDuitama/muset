#pragma once

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>
#include <type_traits>
#include <unistd.h>

#include <fmt/format.h>

#include "config.h"

namespace fs = std::filesystem;


namespace kmat {

constexpr int KL[MUSET_KMER_N] = {MUSET_KMER_LIST};

static const int isnuc[256] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   1,   0,   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   1,   0,
    0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   1,   0,   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   1,   0,
    0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
};


static const unsigned char rctable[256] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
   16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
   32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
   48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
   64, 'T',  66, 'G',  68,  69,  70, 'C',  72,  73,  74,  75,  76,  77, 'N',  79,
   80,  81,  82,  83, 'A',  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,
   96, 't',  98, 'g', 100, 101, 102, 'c', 104, 105, 106, 107, 108, 109, 'n', 111,
  112, 113, 114, 115, 'a', 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
  128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
  144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
  160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
  176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
  192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
  208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
  224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
  240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};


static const int n2kt[256] = {
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 0, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 0, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};


static const unsigned char kt2n[4] = { 'A', 'C', 'T', 'G' };


static inline int actg_compare(const char *a, const char *b, size_t n) {
  for (; n != 0; ++a, ++b, --n) {
    int r = n2kt[*a] - n2kt[*b];
    if (r != 0) return r;
  }
  return 0;
}


static inline int actg_compare(const std::string &a, const std::string &b) {
  const auto a_size = a.size();
  const auto b_size = b.size();
  const int cmp = actg_compare(a.data(), b.data(), std::min(a_size,b_size));
  return cmp != 0 ? cmp : (a_size == b_size ? 0 : a_size < b_size ? -1 : 1 );
}


static inline bool is_valid_kmer(const std::string& kmer) {
  return std::all_of(kmer.begin(), kmer.end(), [](char c){return kmat::isnuc[c];});
}


static inline void reverse_complement_inplace(std::string &seq) {
    std::reverse(seq.begin(), seq.end());
    std::for_each(seq.begin(), seq.end(), [](char& c){ c = rctable[c]; });
}


static inline bool is_canonical(const std::string &seq, bool actg_order=false) {

  auto left = seq.cbegin();
  for (auto right = seq.crbegin(); right != seq.crend(); ++right) {
      auto fc = (unsigned char) *left++;
      auto rc = (unsigned char) rctable[*right++];
      if (fc != rc) {
        return actg_order ? (n2kt[fc] < n2kt[rc]) : (fc < rc);
      }
  }
  return true;
}


static size_t get_nb_samples(const std::string_view line, bool skip_first = false) {
  size_t nb_samples{0};
  size_t idx = line.find_first_not_of(" \t");

  if (skip_first && idx != std::string_view::npos) { 
    idx = line.find_first_of(" \t", idx);
    idx = line.find_first_not_of(" \t", idx);
  }

  while (idx != std::string_view::npos) {
    nb_samples++;
    idx = line.find_first_of(" \t", idx);
    if (idx != std::string_view::npos) {
      idx = line.find_first_not_of(" \t", idx);
    }
  }

  return nb_samples;
}


static inline bool is_kmtricks_dir(const fs::path& dir) {
  return fs::is_directory(dir)
    && fs::is_regular_file(dir/"kmtricks.fof")
    && fs::is_regular_file(dir/"run_infos.txt");
};


static inline bool remove_file(const fs::path& path) {
  if(!path.empty() && fs::is_regular_file(path)) {
    fs::remove(path);
  }
  return false;
}


// https://locklessinc.com/articles/sat_arithmetic/
template<typename T, class = typename std::enable_if<std::is_unsigned_v<T>>::type>
inline T add_sat(T a, T b) noexcept {
  T res = a + b;
	res |= -(res < a);
	return res;
}

}; // namespace kmat
