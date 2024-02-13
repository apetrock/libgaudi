#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <vector>
#include <zlib.h>
#include <random>

#ifndef __DE_COLORES__
#define __DE_COLORES__

namespace gaudi
{
    namespace duchamp
    {

        vec2 random_vec2_with_angle()
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<double> dis(0.0f, 2.0f * M_PI);
            double angle = dis(gen);
            return vec2(std::cos(angle), std::sin(angle));
        }

        double random_normal(double mean, double std_dev)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<double> dis(mean, std_dev);
            return dis(gen);
        }

        vec2 rotate_on_disk(const vec2 &v, double angle)
        {
            double radius = v.norm();
            double theta = std::atan2(v.y(), v.x());
            theta += angle;
            double x = radius * std::cos(theta);
            double y = radius * std::sin(theta);
            return vec2(x, y);
        }

        vec3 compose_color(const vec2 &v)
        {
            double angle = std::atan2(v.y(), v.x());
            double red = std::cos(angle);
            double green = std::cos(angle + M_PI * 2.0f / 3.0f);
            double blue = std::cos(angle + M_PI * 4.0f / 3.0f);
            return vec3(red, green, blue) * 0.5f + vec3(0.5f, 0.5f, 0.5f);
        }

        std::array<vec3, 3> get_rand_colors()
        {
            double gd = (3.0 - sqrt(5.0)) * M_PI;
            double phi = (1.0 + sqrt(5.0)) / 2.0;
            vec2 c0 = random_vec2_with_angle();
            c0.normalize();
            gaudi::vec2 c1 = rotate_on_disk(c0, random_normal(1.0, 0.75) * M_PI);
            c1.normalize();
            gaudi::vec2 c2 = rotate_on_disk(c0, random_normal(1.0, 0.75) * M_PI);
            c2.normalize();

            vec3 c03 = compose_color(c0);
            vec3 c13 = compose_color(c1);
            vec3 c23 = compose_color(c2);

            double D = random_normal(0.1, 0.05);
            if (random_normal(0.0, 0.5) > 0)
            {
                c03 *= 1.0 + D;
                c13 *= 1.0 - D;
            }
            else
            {
                c03 *= 1.0 - D;
                c13 *= 1.0 + D;
            }
            return {c03, c13, c23};
        }

    } // namespace duchamp
} // namespace gaudi

#endif