"""Matrix-aware sympy → Eigen C++ printer (mat3, vec3, no flat scalar CSE blobs)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import sympy as sp
from sympy import Matrix, Symbol
from sympy.printing.c import C99CodePrinter


@dataclass(frozen=True)
class EigenCodegenConfig:
    mat3_type: str = "mat3"
    vec3_type: str = "vec3"
    real_type: str = "real"
    identity: str = "mat3::Identity()"
    zero_mat3: str = "mat3::Zero()"


class EigenPrinter(C99CodePrinter):
    """Emit Eigen-style mat3/vec3 expressions from sympy matrix algebra."""

    def __init__(
        self,
        symbol_map: Optional[Dict[Symbol, str]] = None,
        config: Optional[EigenCodegenConfig] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self._symbol_map = symbol_map or {}
        self._cfg = config or EigenCodegenConfig()

    def _print_Symbol(self, expr: Symbol) -> str:
        if expr in self._symbol_map:
            return self._symbol_map[expr]
        return super()._print_Symbol(expr)

    def _print_Rational(self, expr) -> str:
        if expr.q == 1:
            return f"{float(expr.p)}.0"
        return f"{float(expr.p)}.0/{float(expr.q)}.0"

    def _print_Integer(self, expr) -> str:
        return f"{int(expr)}.0"

    def _print_Pow(self, expr) -> str:
        base, exp = expr.as_base_exp()
        if exp == sp.S.Half:
            return f"std::sqrt({self._print(base)})"
        if exp == sp.Rational(-1, 2):
            return f"(1.0/std::sqrt({self._print(base)}))"
        if exp == -1:
            return f"(1.0/({self._print(base)}))"
        if exp.is_Integer:
            return f"std::pow({self._print(base)}, {int(exp)}.0)"
        if exp.is_Rational and exp.q != 1:
            return f"std::pow({self._print(base)}, {float(exp.p)}/{float(exp.q)})"
        return super()._print_Pow(expr)

    def _print_MatMul(self, expr) -> str:
        parts = [self._print(a) for a in expr.args]
        if len(parts) == 1:
            return parts[0]
        return "(" + " * ".join(parts) + ")"

    def _print_Transpose(self, expr) -> str:
        return f"({self._print(expr.args[0])}).transpose()"

    def _print_Trace(self, expr) -> str:
        return f"({self._print(expr.args[0])}).trace()"

    def _print_Identity(self, expr) -> str:
        del expr
        return self._cfg.identity

    def _print_ZeroMatrix(self, expr) -> str:
        rows, cols = expr.shape
        if rows == 3 and cols == 3:
            return self._cfg.zero_mat3
        return super()._print_ZeroMatrix(expr)

    def _print_MatrixSymbol(self, expr) -> str:
        if expr in self._symbol_map:
            return self._symbol_map[expr]
        return str(expr)

    def _print_MatrixElement(self, expr) -> str:
        parent = expr.parent
        if parent in self._symbol_map:
            base = self._symbol_map[parent]
            return f"{base}({expr.i},{expr.j})"
        return super()._print_MatrixElement(expr)


def emit_stmt(printer: EigenPrinter, lhs: str, expr: sp.Basic, indent: str = "  ") -> str:
    return f"{indent}{lhs} = {printer.doprint(expr)};"


def emit_mat3_out_fill(
    printer: EigenPrinter,
    out_var: str,
    matrix: Matrix,
    *,
    indent: str = "  ",
) -> str:
    lines: List[str] = []
    for i in range(3):
        for j in range(3):
            lines.append(
                f"{indent}{out_var}({i},{j}) = {printer.doprint(matrix[i, j])};"
            )
    return "\n".join(lines)


def emit_vec3_function(
    name: str,
    args: str,
    components: Sequence[sp.Expr],
    printer: EigenPrinter,
    preamble: str = "",
) -> str:
    body = "\n".join(
        f"  {printer._cfg.real_type} {c} = {printer.doprint(expr)};"
        for c, expr in zip(("dx", "dy", "dz"), components)
    )
    return f"""
inline {printer._cfg.vec3_type} {name}({args}) {{
{preamble}{body}
  return {printer._cfg.vec3_type}(dx, dy, dz);
}}
"""
