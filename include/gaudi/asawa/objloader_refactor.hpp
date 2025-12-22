/*
 *  obj_loader.h
 *  Manifold
 *
 *  Created by John Delaney on 3/23/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <sstream>
#include <vector>

#include <cmath>

#include "gaudi/vec_addendum.h"

#ifndef __ASAWA_OBJ_LOADER__
#define __ASAWA_OBJ_LOADER__
namespace gaudi {
namespace asawa {

typedef double real;
typedef int index_t;
typedef Eigen::Matrix<real, 3, 1> vec3;

// Core parsing from stream (shared logic)
void loadObjStream(std::istream &stream, std::vector<vec3> &vertices,
                 std::vector<std::vector<int>> &faces) {
  std::string line;

  while (std::getline(stream, line)) {
      std::istringstream ss(line);
      std::vector<std::string> vals;
      while (ss) {
        std::string s;
      if (!std::getline(ss, s, ' '))
          break;
        vals.push_back(s);
      }
      if (vals.size() == 0)
        continue;

      if (vals[0] == "v") {
        std::string fvs = line.substr(2);
        float vx, vy, vz;
        sscanf(fvs.c_str(), "%f %f %f", &vx, &vy, &vz);
        vec3 v(vx, vy, vz);
        vertices.push_back(v);
      } else if (vals[0] == "f") {
        std::vector<int> curFace;
      for (size_t i = 1; i < vals.size(); i++) {
          int index;
          std::istringstream iss(vals[i]);
          if (iss >> index)
            curFace.push_back(index - 1);
        }
        faces.push_back(curFace);
    }
    // Skip normals (vn) and other directives for now
  }
}

// Load OBJ from file path (original interface)
void loadObjfile(const std::string &s, std::vector<vec3> &vertices,
                 std::vector<std::vector<int>> &faces) {
  std::string fname = s;
  std::ifstream myfile(fname.c_str());

  std::cout << "loading " << s << " from " << fname << " . . . " << std::endl;
  if (myfile.is_open()) {
    std::cout << ".obj file open" << std::endl;
    loadObjStream(myfile, vertices, faces);
  }
  myfile.close();
}

// Load OBJ from string content (for WASM/JS)
void loadObjFromString(const std::string &content, std::vector<vec3> &vertices,
                       std::vector<std::vector<int>> &faces) {
  std::istringstream stream(content);
  loadObjStream(stream, vertices, faces);
}
} // namespace asawa
} // namespace gaudi

#endif
