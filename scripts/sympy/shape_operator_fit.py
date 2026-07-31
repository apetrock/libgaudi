"""Shape-operator tip-in linearizers: vech(P_n H(Q) P_n) = M Q."""

from __future__ import annotations

import sympy as sp
from sympy import Matrix, symbols

from albers_printer import AlbersPrinter, write_header
from darboux import C as darboux_C
from darboux import H_expr as darboux_H
from implicit import x, y, z
from quadric import H_expr as quadric_H
from quadric import Q as quadric_Q


def _vech3(M: Matrix) -> Matrix:
    return Matrix([M[0, 0], M[1, 1], M[2, 2], M[0, 1], M[0, 2], M[1, 2]])


def _projector(n: Matrix) -> Matrix:
    return sp.eye(3) - n * n.T


def mk_quad_W_rows_symbolic(n: Matrix) -> Matrix:
    """6x10: vech(P H(Q) P) linear in Q (H independent of x for quadrics)."""
    H = quadric_H()
    P = _projector(n)
    Php = P * H * P
    v = _vech3(Php)
    rows = []
    for i in range(6):
        row = [sp.simplify(v[i].diff(quadric_Q[j])) for j in range(10)]
        rows.append(row)
    return Matrix(rows)


def mk_darboux_W_rows_symbolic(n: Matrix, xv: Matrix) -> Matrix:
    """6x14 at sample x: vech(P H(C,x) P) linear in C."""
    H = darboux_H().subs({x: xv[0], y: xv[1], z: xv[2]})
    P = _projector(n)
    Php = P * H * P
    v = _vech3(Php)
    rows = []
    for i in range(6):
        row = [sp.simplify(v[i].diff(darboux_C[j])) for j in range(14)]
        rows.append(row)
    return Matrix(rows)


def _emit_mat_assign(printer: AlbersPrinter, M: Matrix, var: str) -> str:
    lines = [f"  {var}.setZero();"]
    for i in range(M.rows):
        for j in range(M.cols):
            e = sp.simplify(M[i, j])
            if e == 0:
                continue
            lines.append(f"  {var}({i},{j}) = {printer.doprint(e)};")
    return "\n".join(lines)


def emit_shape_operator_fit() -> str:
    """Body only — write_header adds banner/footer."""
    nx, ny, nz = symbols("nx ny nz", real=True)
    n = Matrix([nx, ny, nz])
    xx, yy, zz = symbols("xx yy zz", real=True)
    xv = Matrix([xx, yy, zz])

    printer = AlbersPrinter(
        symbol_map={
            nx: "n[0]",
            ny: "n[1]",
            nz: "n[2]",
            xx: "x[0]",
            yy: "x[1]",
            zz: "x[2]",
        }
    )

    Mq = mk_quad_W_rows_symbolic(n)
    Md = mk_darboux_W_rows_symbolic(n, xv)

    body = []
    body.append(
        "inline Eigen::Matrix<real, 6, 10> mk_quad_W_rows(const vec3 &n) {\n"
        "  Eigen::Matrix<real, 6, 10> M;\n"
        + _emit_mat_assign(printer, Mq, "M")
        + "\n  return M;\n}\n"
    )
    body.append(
        "inline Eigen::Matrix<real, 6, 14> mk_darboux_W_rows(const vec3 &n, const vec3 &x) {\n"
        "  Eigen::Matrix<real, 6, 14> M;\n"
        + _emit_mat_assign(printer, Md, "M")
        + "\n  return M;\n}\n"
    )
    return "\n".join(body)


if __name__ == "__main__":
    from pathlib import Path

    out = Path(__file__).resolve().parent / "generated"
    write_header(out / "shape_operator_fit_generated.hpp", emit_shape_operator_fit())
    print(f"Wrote {out / 'shape_operator_fit_generated.hpp'}")
