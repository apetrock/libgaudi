#ifndef GAUDI_WINDYCHIEN_BRAID_HPP
#define GAUDI_WINDYCHIEN_BRAID_HPP

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace gaudi {
namespace windychien {

/// KnotTheory-style braid representative BR[k, l]:
///   strands = k
///   word[j] = ±i  → crossing between strands i and i+1 (1-based);
///              + = right-handed, - = left-handed.
/// Optional layer_counts: partition word into simultaneous crossing frames
/// (one planar/sphere column per layer). Empty ⇒ one generator per column.
struct braid {
  std::string name;
  int strands = 0;
  std::vector<int> word;
  std::vector<int> layer_counts;
};

inline int braid_n_columns(const braid &b) {
  if (b.layer_counts.empty()) {
    return static_cast<int>(b.word.size());
  }
  return static_cast<int>(b.layer_counts.size());
}

inline void validate_braid(const braid &b) {
  if (b.strands < 2) {
    throw std::runtime_error("windychien::braid: need strands >= 2 (got " +
                             std::to_string(b.strands) + ")");
  }
  for (int g : b.word) {
    const int i = std::abs(g);
    if (i < 1 || i >= b.strands) {
      throw std::runtime_error(
          "windychien::braid '" + b.name + "': generator " + std::to_string(g) +
          " invalid for " + std::to_string(b.strands) + " strands");
    }
  }
  if (!b.layer_counts.empty()) {
    int sum = 0;
    for (int c : b.layer_counts) {
      if (c < 0) {
        throw std::runtime_error("windychien::braid '" + b.name +
                                 "': negative layer_counts entry");
      }
      sum += c;
    }
    if (sum != static_cast<int>(b.word.size())) {
      throw std::runtime_error(
          "windychien::braid '" + b.name +
          "': layer_counts sum != word size");
    }
  }
}

/// Plain / tabby weave: n_frames axial steps; even frames cross pairs
/// (0,1),(2,3),… ; odd frames cross (1,2),(3,4),… . Signs checkerboard.
inline braid make_plain_weave(int strands, int n_frames,
                              bool checkerboard_signs = true) {
  if (strands < 2) {
    throw std::runtime_error("make_plain_weave: strands >= 2");
  }
  if (n_frames < 1) {
    throw std::runtime_error("make_plain_weave: n_frames >= 1");
  }
  braid b;
  b.name = "plain_weave_" + std::to_string(strands) + "x" +
           std::to_string(n_frames);
  b.strands = strands;
  b.word.reserve(static_cast<size_t>(n_frames) *
                 static_cast<size_t>((strands + 1) / 2));
  b.layer_counts.reserve(static_cast<size_t>(n_frames));

  for (int f = 0; f < n_frames; ++f) {
    // Even frame: generators 1,3,5,… (pairs 0-1, 2-3, …)
    // Odd frame:  generators 2,4,6,… (pairs 1-2, 3-4, …)
    const int g0 = (f % 2 == 0) ? 1 : 2;
    int k = 0;
    const size_t before = b.word.size();
    for (int g = g0; g < strands; g += 2) {
      int sign = 1;
      if (checkerboard_signs) {
        sign = ((f + k) % 2 == 0) ? 1 : -1;
      }
      b.word.push_back(sign * g);
      ++k;
    }
    b.layer_counts.push_back(static_cast<int>(b.word.size() - before));
  }
  validate_braid(b);
  return b;
}

/// After the word, pos[s] = row occupied by starting strand s (0-based).
inline std::vector<int> braid_end_rows(const braid &b) {
  validate_braid(b);
  std::vector<int> pos(static_cast<size_t>(b.strands));
  for (int s = 0; s < b.strands; ++s) {
    pos[static_cast<size_t>(s)] = s;
  }
  for (int g : b.word) {
    const int a = std::abs(g) - 1;
    const int b_row = a + 1;
    for (int s = 0; s < b.strands; ++s) {
      int &r = pos[static_cast<size_t>(s)];
      if (r == a) {
        r = b_row;
      } else if (r == b_row) {
        r = a;
      }
    }
  }
  return pos;
}

} // namespace windychien
} // namespace gaudi

#endif
