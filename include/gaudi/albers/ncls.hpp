#ifndef NORMAL_CONSTRAINED_LEAST_SQUARES_H
#define NORMAL_CONSTRAINED_LEAST_SQUARES_H
#include <Eigen/Dense>
#include "gaudi/common.h"

namespace gaudi
{
    namespace albers
    {
        vec4 mk_N(const vec3 &N) { return vec4(0.0, N[0], N[1], N[2]); }
    }
}

#endif
