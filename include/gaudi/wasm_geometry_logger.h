#ifndef __GAUDI_WASM_GEOMETRY_LOGGER_H__
#define __GAUDI_WASM_GEOMETRY_LOGGER_H__

#include "common.h"
#include <vector>
#include <Eigen/Dense>

namespace gaudi {

class wasm_geometry_logger {
public:
    static wasm_geometry_logger& get_instance() {
        static wasm_geometry_logger instance;
        return instance;
    }

    static void clear() {
        get_instance()._clear();
    }

    static void point(const vec3& p0, const vec4& color) {
        get_instance()._point(p0, color);
    }

    static void line(const vec3& p0, const vec3& p1, const vec4& color) {
        get_instance()._line(p0, p1, color);
    }

    static void line4(const vec4& p0, const vec4& p1, const vec4& color) {
        get_instance()._line(p0.head<3>(), p1.head<3>(), color);
    }

    static void lines(const std::vector<vec3>& p0, const std::vector<vec3>& p1, const std::vector<vec4>& colors) {
        get_instance()._lines(p0, p1, colors);
    }

    static void box(const vec3 &cen, const vec3 &h, const vec4 &col) {
        get_instance()._box(cen, h, col);
    }

    static void ext(const vec3 &mn, const vec3 &mx, const vec4 &col) {
        get_instance()._ext(mn, mx, col);
    }

    static void field(const std::vector<vec3> &p, const std::vector<vec3> &dirs,
                    double D = 0.1) {
        get_instance()._field(p, dirs, D);
    }

    static void frame(const mat3 &M, const vec3 &c, double C) {
        get_instance()._frame(M, c, C);
    }

    const std::vector<vec3>& get_points() const { return _points; }
    const std::vector<vec4>& get_point_colors() const { return _point_colors; }
    const std::vector<vec3>& get_lines() const { return _lines; }
    const std::vector<vec4>& get_line_colors() const { return _line_colors; }

private:
    wasm_geometry_logger() = default;
    wasm_geometry_logger(const wasm_geometry_logger&) = delete;
    wasm_geometry_logger& operator=(const wasm_geometry_logger&) = delete;

    void _clear() {
        _points.clear();
        _point_colors.clear();
        _lines.clear();
        _line_colors.clear();
    }

    void _point(const vec3& p0, const vec4& color) {
        _points.push_back(p0);
        _point_colors.push_back(color);
    }

    void _line(const vec3& p0, const vec3& p1, const vec4& color) {
        _lines.push_back(p0);
        _lines.push_back(p1);
        _line_colors.push_back(color);
        _line_colors.push_back(color);
    }

    void _lines(const std::vector<vec3>& p0, const std::vector<vec3>& p1, const std::vector<vec4>& colors) {
      for (size_t i = 0; i < p0.size(); ++i) {
        _lines.push_back(p0[i]);
        _lines.push_back(p1[i]);
        _line_colors.push_back(colors[i]);
        _line_colors.push_back(colors[i]);
      }
    }

    void _box(const vec3 &cen, const vec3 &h, const vec4 &col) {
        vec3 mn = cen - h;
        vec3 mx = cen + h;
        _ext(mn, mx, col);
    }

    void _ext(const vec3 &mn, const vec3 &mx, const vec4 &col) {
        _line(vec3(mn[0], mn[1], mn[2]), vec3(mx[0], mn[1], mn[2]), col);
        _line(vec3(mn[0], mn[1], mn[2]), vec3(mn[0], mx[1], mn[2]), col);
        _line(vec3(mn[0], mn[1], mn[2]), vec3(mn[0], mn[1], mx[2]), col);

        _line(vec3(mx[0], mx[1], mx[2]), vec3(mn[0], mx[1], mx[2]), col);
        _line(vec3(mx[0], mx[1], mx[2]), vec3(mx[0], mn[1], mx[2]), col);
        _line(vec3(mx[0], mx[1], mx[2]), vec3(mx[0], mx[1], mn[2]), col);

        _line(vec3(mn[0], mx[1], mn[2]), vec3(mx[0], mx[1], mn[2]), col);
        _line(vec3(mn[0], mn[1], mx[2]), vec3(mx[0], mn[1], mx[2]), col);

        _line(vec3(mx[0], mn[1], mn[2]), vec3(mx[0], mn[1], mx[2]), col);
        _line(vec3(mn[0], mx[1], mn[2]), vec3(mn[0], mx[1], mx[2]), col);

        _line(vec3(mn[0], mn[1], mx[2]), vec3(mn[0], mx[1], mx[2]), col);
        _line(vec3(mx[0], mn[1], mn[2]), vec3(mx[0], mx[1], mn[2]), col);
    }

    void _field(const std::vector<vec3> &p, const std::vector<vec3> &dirs, double D) {
        for (size_t i = 0; i < p.size(); i++) {
            _line(p[i], p[i] + D * dirs[i], vec4(0.5, 0.5, 0.5, 1.0));
        }
    }

    void _frame(const mat3 &M, const vec3 &c, double C) {
        _line(c, c + C * M.col(0), vec4(1, 0, 0, 1));
        _line(c, c + C * M.col(1), vec4(0, 1, 0, 1));
        _line(c, c + C * M.col(2), vec4(0, 0, 1, 1));
    }

    std::vector<vec3> _points;
    std::vector<vec4> _point_colors;
    std::vector<vec3> _lines;
    std::vector<vec4> _line_colors;
};

} 

#endif
