
#ifndef __ASAWA_FACE_LIST_LOADER__
#define __ASAWA_FACE_LIST_LOADER__

#include <cassert>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"

namespace gaudi {
namespace asawa {

bool patch_hole(index_t i, std::vector<index_t> &face_start,
                std::vector<index_t> &vert_start,
                std::vector<index_t> &corners_next,
                std::vector<index_t> &corners_vert,
                std::vector<index_t> &corners_face) {

  auto other = [](int id) { return 2 * (id / 2) + (id + 1) % 2; };
  auto next = [&corners_next](int id) { return corners_next[id]; };
  auto swing = [&corners_next, other, next](int id) { return next(other(id)); };

  index_t c0 = i;
  index_t c1 = other(i);
  assert(next(c1) == -1);

  std::vector<index_t> patch;

  index_t j1 = c0;
  index_t j0 = next(j1);

  bool it = true;
  int k = 0;
  int maxk = 100;
  while (it && k < maxk) {
    it = j0 != j1;
    bool its = true;
    if (swing(j0) < 0) {
      patch.push_back(j0);
      j0 = next(j0);
    } else {
      j0 = swing(j0);
    }
    k++;
  }
#if 0
  for (int i = 0; i < patch.size(); i++) {
    std::cout << other(patch[i]) << " ";
  }
  std::cout << std::endl;
  for (int i = 0; i < patch.size(); i++) {
    std::cout << next(other(patch[i])) << " ";
  }
  std::cout << std::endl;
#endif
  if (k < maxk) {
    int fs = face_start.size();
    for (int i = 0; i < patch.size(); i++) {
      int ip = (i + 1) % patch.size();

      index_t cn = other(patch[i]);
      index_t cp = other(patch[ip]);
      index_t v = corners_vert[swing(cp)];
      corners_next[cp] = cn;
      corners_vert[cp] = v;
      corners_face[cp] = fs;
    }
    face_start.push_back(other(patch[0]));
  }
  return true;
}

void fill_holes(std::vector<index_t> &face_start,
                std::vector<index_t> &vert_start,
                std::vector<index_t> &corners_next,
                std::vector<index_t> &corners_vert,
                std::vector<index_t> &corners_face) {
  int kk = 0;
  int max = 40;

  for (int i = 0; i < corners_next.size(); i += 2) {
    index_t other = i + 1;
    if (corners_next[other] < 0) {
      patch_hole(i, face_start, vert_start, corners_next, corners_vert,
                 corners_face);
    }
  }
}

void assemble_table(const std::vector<vec3> &vertices,
                    const std::vector<std::vector<index_t>> &faces,
                    std::vector<index_t> &corners_next,
                    std::vector<index_t> &corners_vert,
                    std::vector<index_t> &corners_face) {
  std::vector<std::vector<index_t>> face_corners;
  std::vector<index_t> flat_corners;
  std::vector<std::vector<index_t>> corners_on_vertices(vertices.size());

  index_t idx = 0;
  for (index_t i = 0; i < faces.size(); i++) {
    const std::vector<index_t> &face = faces[i];
    std::vector<index_t> corners(face);

    for (index_t j = 0; j < face.size(); j++) {
      index_t corner_idx = idx;

      corners[j] = corner_idx;
      flat_corners.push_back(corner_idx);
      corners_on_vertices[face[j]].push_back(corner_idx);

      corners_face.push_back(i);
      corners_vert.push_back(face[j]);

      idx++;
    }
    face_corners.push_back(corners);
  }

  // build corners next
  corners_next = std::vector<index_t>(flat_corners.size(), -1);
  std::vector<index_t> corners_prev =
      std::vector<index_t>(flat_corners.size(), -1);

  for (index_t i = 0; i < face_corners.size(); i++) {
    const std::vector<index_t> &face = face_corners[i];
    for (index_t j = 1; j < face.size() + 1; j++) {
      index_t jm = (j - 1);
      index_t j0 = j % face.size();
      index_t jp = (j0 + 1) % face.size();

      index_t jjm = face[jm];
      index_t jj0 = face[j0];
      index_t jjp = face[jp];

      corners_next[jj0] = jjp;
      corners_prev[jj0] = jjm;
    }
  }

#if 0
  for (int i = 0; i < corners_next.size(); i++) {
    std::cout << " cn " << corners_next[i] << " c0 " << flat_corners[i]
              << " cp " << corners_prev[i] << std::endl;
  }
#endif

  // build corners next
  std::vector<index_t> interaction_list;
  for (index_t idx_vert = 0; idx_vert < corners_on_vertices.size();
       idx_vert++) {
    const std::vector<index_t> &faces_on_vert = corners_on_vertices[idx_vert];
    for (index_t i = 0; i < faces_on_vert.size(); i++) {
      for (index_t j = i + 1; j < faces_on_vert.size(); j++) {
        interaction_list.push_back(faces_on_vert[i]);
        interaction_list.push_back(faces_on_vert[j]);
      }
    }
  }

  std::vector<index_t> corner_index(flat_corners);
  std::vector<index_t> corners(flat_corners.size(), -1);
  std::vector<index_t> corners_twin =
      std::vector<index_t>(flat_corners.size(), -1);

  for (index_t k = 0; k < interaction_list.size(); k += 2) {
    index_t ci = interaction_list[k + 0];
    index_t cin = corners_next[ci];
    index_t cip = corners_prev[ci];

    index_t cj = interaction_list[k + 1];
    index_t cjn = corners_next[cj];
    index_t cjp = corners_prev[cj];

    index_t vi = corners_vert[ci];
    index_t vin = corners_vert[cin];
    index_t vip = corners_vert[cip];

    index_t vj = corners_vert[cj];
    index_t vjn = corners_vert[cjn];
    index_t vjp = corners_vert[cjp];
    std::cout << ci << " " << cj << std::endl;
    std::cout << vi << " " << vj << std::endl;
    std::cout << " " << vi << " " << vin << " " << vip << std::endl;
    std::cout << " " << vj << " " << vjn << " " << vjp << std::endl;

    if (vip == vjn) {
      corners_twin[cj] = cip;
      corners_twin[cip] = cj;

      std::cout << "   ci:" << cj << " " << cip << std::endl;
      std::cout << "   vi:" << vj << " " << vip << std::endl;
      // std::cout << "    v: " << vi << " " << vjn << " " << vin <<
      // std::endl;
    }
    if (vjp == vin) {
      corners_twin[ci] = cjp;
      corners_twin[cjp] = ci;

      std::cout << "   cj:" << ci << " " << cjp << std::endl;
      std::cout << "   vj:" << vi << " " << vjp << std::endl;

      // std::cout << "    v: " << vi << " " << vjn << " " << vin <<
      // std::endl;
    }
  }

  for (int i = 0; i < corners.size(); i++) {
    std::cout << " ci " << corner_index[i] << " ct " << corners_twin[i]
              << " cn " << corners_next[i] << " cf " << corners_face[i]
              << " cv " << corners_vert[i] << std::endl;
  }
  std::cout << std::endl;

  idx = 0;
  std::vector<bool> corner_allocated(flat_corners.size(), false);
  std::vector<index_t> new_corners;
  for (int i = 0; i < flat_corners.size(); i++) {
    if (corner_allocated[i])
      continue;
    index_t ct = corners_twin[i];
    new_corners.push_back(i);
    new_corners.push_back(ct);
    corner_allocated[i] = true;
    if (ct > 0)
      corner_allocated[ct] = true;
  }

  std::vector<index_t> corner_map(corner_index.size());
  std::vector<index_t> new_corners_next(new_corners.size(), -1);
  std::vector<index_t> new_corners_vert(new_corners.size(), -1);
  std::vector<index_t> new_corners_face(new_corners.size(), -1);
  std::cout << new_corners.size() << std::endl;
  for (int i = 0; i < new_corners.size(); i++) {
    int ii = new_corners[i];
    if (ii > -1) {
      new_corners_next[i] = corners_next[ii];
      new_corners_vert[i] = corners_vert[ii];
      new_corners_face[i] = corners_face[ii];
      corner_map[ii] = i;
    }
  }

  for (int i = 0; i < new_corners_next.size(); i++) {
    if (new_corners_next[i] > 0) {
      std::cout << " " << new_corners_next[i] << " "
                << corner_map[new_corners_next[i]] << std::endl;
      new_corners_next[i] = corner_map[new_corners_next[i]];
    }
  }

  corners_next = new_corners_next;
  corners_face = new_corners_face;
  corners_vert = new_corners_vert;

  for (int i = 0; i < new_corners.size(); i++) {
    std::cout << i << " ci " << new_corners[i] << " cn " << new_corners_next[i]
              << " cf " << new_corners_face[i] << " cv " << new_corners_vert[i]
              << std::endl;
  }

  std::vector<index_t> face_start(faces.size(), -1);
  std::vector<index_t> vert_start(vertices.size(), -1);
  for (int i = 0; i < corners_next.size(); i++) {
    if (corners_face[i] > 0)
      face_start[corners_face[i]] = i;
    if (corners_vert[i] > 0)
      vert_start[corners_vert[i]] = i;
  }
  fill_holes(face_start, vert_start, corners_next, corners_vert, corners_face);

  return;
#if 0
  for (int i = 0; i < face_start.size(); i++) {
    int j_end = corners[face_start[i]];
    int j0 = corners_next[j_end];
    bool it = true;
    while (it) {
      it = j0 != j_end;
      corners_prev[corners_next[j0]] = j0;
      j0 = corners_next[j0];
    }
  }
#endif
}
} // namespace asawa
} // namespace gaudi

#endif