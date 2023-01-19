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

#include "manifold/vec_addendum.h"

#ifndef __ASAWA_OBJ_LOADER__
#define __ASAWA_OBJ_LOADER__

namespace asawa {

typedef double real;
typedef int index_t;
typedef Eigen::Matrix<real, 3, 1> vec3;

void loadObjfile(const std::string &s, std::vector<vec3> &vertices,
                 std::vector<std::vector<int>> &faces) {

  std::string line;
  std::string fname = s;
  std::ifstream myfile(fname.c_str());

  bool bNormalsIncluded = 0;
  bool bFacesIncluded = 0;
  std::cout << "loading " << s << " from " << fname << " . . . " << endl;
  // Parse File and push_back vectors for FACE, VERTS, MAP, and NORMALS
  if (myfile.is_open()) {
    std::cout << ".obj file open" << endl;
    while (!myfile.eof()) {

      getline(myfile, line);

      std::istringstream ss(line);
      std::vector<std::string> vals;
      while (ss) {
        std::string s;
        if (!getline(ss, s, ' '))
          break;
        vals.push_back(s);
      }
      if (vals.size() == 0)
        continue;
      int in = line.find_first_of(" ");
      std::string a = line.substr(0, in);

      int vIdx, fIdx = 0;

      if (vals[0] == "v") {

        std::string fvs = line.substr(2);
        float vx, vy, vz;
        sscanf(fvs.c_str(), "%f %f %f", &vx, &vy, &vz);
        vec3 v(vx, vy, vz);

        vertices.push_back(v);

        vIdx++;
        // cout << mPoints[mPoints.size()-1] << endl;

      } else if (vals[0] == "f") {
        bFacesIncluded = 1;
        int faceSize = vals.size() - 1;
        std::vector<int> curFace;

        for (int i = 1; i < vals.size(); i++) {

          int index;
          std::istringstream iss(vals[i]);
          if (iss >> index)
            curFace.push_back(index - 1);
        }

        faces.push_back(curFace);
        fIdx++;
      } else if (vals[0] == "vn") {
        bNormalsIncluded = 1;

        std::string fvs = line.substr(3);
        int fi = fvs.find_first_of(" ");
        std::string fv1 = fvs.substr(0, fi);
        std::string fvs2 = fvs.substr(fi + 1);
        fi = fvs2.find_first_of(" ");
        std::string fv2 = fvs2.substr(0, fi);
        std::string fv3 = fvs2.substr(fi + 1, 10);

        vec3 normal;

        std::stringstream(fv1) >> normal[0];
        std::stringstream(fv2) >> normal[1];
        std::stringstream(fv3) >> normal[2];

        // mNormal.push_back(normal);
      }
    }
  }
  myfile.close();
}

} // namespace asawa

#endif
