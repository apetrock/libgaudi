#ifndef GAUDI_WINDYCHIEN_BRAID_HPP
#define GAUDI_WINDYCHIEN_BRAID_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace gaudi {
namespace windychien {

/// One knot source. `name` selects it.
/// `"plain_weave"` or `"plane_weave"` builds a tabby cloth from `strands` and
/// `frames`. `"random"` builds an irregular braid with `strands` strings and
/// `frames` crossings, repeated from `seed`. Any other name is loaded from the
/// braid catalog (`"K11a359"`, `"3_1"`, `"T8_7"`, …), and `strands` / `frames`
/// are unused.
struct knot_spec {
  std::string name = "K11a359";
  int strands = 50;
  int frames = 50;
  std::uint32_t seed = 1;
};

inline bool knot_is_plain_weave(const std::string &name) {
  return name == "plain_weave" || name == "plane_weave";
}

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

/// Link components of the closure. Signs do not matter: a generator only swaps
/// the two strands in its rows.
inline int braid_component_count(const braid &b) {
  const std::vector<int> row = braid_end_rows(b);
  std::vector<char> seen(static_cast<size_t>(b.strands), 0);
  int cycles = 0;
  for (int s = 0; s < b.strands; ++s) {
    if (seen[static_cast<size_t>(s)]) {
      continue;
    }
    ++cycles;
    int j = s;
    do {
      seen[static_cast<size_t>(j)] = 1;
      j = row[static_cast<size_t>(j)];
    } while (j != s);
  }
  return cycles;
}

/// Irregular braid. `n_crossings` is the requested number of crossings.
/// Each column is one slide: the number of rows it crosses is chosen
/// uniformly from a single pair up to the whole stack, and the slide can
/// sit anywhere in the stack. Up-and-over is never paired with down-and-over
/// (same for under): one chirality ties direction to side, and a slide that
/// crosses at least two rows weaves the other way at one or two crossings.
/// Every adjacent pair is used, and the closure is a single component. A
/// single n-cycle needs the same parity as `strands - 1`, so one crossing
/// is added when the request has the other parity. Same `seed` repeats the
/// word.
inline braid make_random_braid(int strands, int n_crossings,
                               std::uint32_t seed = 1) {
  if (strands < 2) {
    throw std::runtime_error("make_random_braid: strands >= 2");
  }
  if (n_crossings < strands - 1) {
    throw std::runtime_error(
        "make_random_braid: need at least strands-1 crossings so every pair "
        "is used");
  }
  // Each generator is a transposition, so the closure permutation has sign
  // (-1)^{n_crossings}. An n-cycle has sign (-1)^{strands-1}.
  if ((n_crossings % 2) != ((strands - 1) % 2)) {
    ++n_crossings;
  }

  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> coin(0, 1);

  // lo, hi are inclusive 0-based rows. `up` sends the strand at lo to hi.
  // `weave[i]` flips over/under at that crossing along the slide.
  const auto push_slide = [](std::vector<int> &col, int lo, int hi, bool up,
                             bool over, const std::vector<char> &weave) {
    int which = 0;
    const auto emit = [&](int g) {
      const bool woven = weave[static_cast<size_t>(which)] != 0;
      const bool side_over = woven ? !over : over;
      const int sign = (up == side_over) ? 1 : -1;
      col.push_back(sign * g);
      ++which;
    };
    if (up) {
      for (int g = lo + 1; g <= hi; ++g) {
        emit(g);
      }
    } else {
      for (int g = hi; g >= lo + 1; --g) {
        emit(g);
      }
    }
  };

  constexpr int k_max_attempts = 4096;
  for (int attempt = 0; attempt < k_max_attempts; ++attempt) {
    std::vector<std::vector<int>> columns;
    std::vector<std::vector<int>> built_cols;
    bool undone = false;
    int left = n_crossings;
    // Up-over cancels down-over, and up-under cancels down-under.
    // Tie the side to the direction so a return trip is not the inverse.
    const bool up_is_over = coin(rng) != 0;

    const auto column_is_inverse = [&](const std::vector<int> &col) {
      for (const std::vector<int> &prev : built_cols) {
        if (prev.size() != col.size()) {
          continue;
        }
        bool inv = true;
        for (size_t i = 0; i < col.size(); ++i) {
          if (col[i] != -prev[prev.size() - 1 - i]) {
            inv = false;
            break;
          }
        }
        if (inv) {
          return true;
        }
      }
      return false;
    };

    const auto take_slide = [&](std::vector<int> &col, int lo, int hi) {
      const int steps = hi - lo;
      int n_weave = 0;
      if (steps >= 2) {
        n_weave = 1;
        if (steps >= 3 && coin(rng) != 0) {
          n_weave = 2;
        }
        if (n_weave > steps - 1) {
          n_weave = steps - 1;
        }
      }
      std::vector<int> slots(static_cast<size_t>(steps));
      std::iota(slots.begin(), slots.end(), 0);
      bool kept = false;
      for (int try_i = 0; try_i < 8; ++try_i) {
        const bool up = coin(rng) != 0;
        const bool over = up ? up_is_over : !up_is_over;
        std::shuffle(slots.begin(), slots.end(), rng);
        std::vector<char> weave(static_cast<size_t>(steps), 0);
        for (int w = 0; w < n_weave; ++w) {
          weave[static_cast<size_t>(slots[static_cast<size_t>(w)])] = 1;
        }
        col.clear();
        push_slide(col, lo, hi, up, over, weave);
        if (!column_is_inverse(col)) {
          kept = true;
          break;
        }
      }
      if (!kept) {
        undone = true;
      }
      built_cols.push_back(col);
      left -= steps;
    };

    while (left > 0) {
      const int max_steps = std::min(strands - 1, left);
      std::uniform_int_distribution<int> step_dist(1, max_steps);
      const int steps = step_dist(rng);
      const int slack = strands - (steps + 1);
      std::uniform_int_distribution<int> off_dist(0, std::max(slack, 0));
      const int lo = off_dist(rng);
      std::vector<int> col;
      take_slide(col, lo, lo + steps);
      columns.push_back(std::move(col));
    }

    braid b;
    b.name = "random_" + std::to_string(strands) + "x" +
             std::to_string(n_crossings);
    b.strands = strands;
    b.layer_counts.reserve(columns.size());
    for (const std::vector<int> &col : columns) {
      b.word.insert(b.word.end(), col.begin(), col.end());
      b.layer_counts.push_back(static_cast<int>(col.size()));
    }
    if (undone || static_cast<int>(b.word.size()) != n_crossings) {
      continue;
    }
    std::vector<char> covered(static_cast<size_t>(strands - 1), 0);
    for (int g : b.word) {
      covered[static_cast<size_t>(std::abs(g) - 1)] = 1;
    }
    if (std::find(covered.begin(), covered.end(), 0) != covered.end()) {
      continue;
    }
    validate_braid(b);
    if (braid_component_count(b) == 1) {
      return b;
    }
  }
  throw std::runtime_error(
      "make_random_braid: no single-component braid in " +
      std::to_string(k_max_attempts) + " draws (strands " +
      std::to_string(strands) + ", crossings " + std::to_string(n_crossings) +
      ", seed " + std::to_string(seed) + ")");
}

} // namespace windychien
} // namespace gaudi

#endif
