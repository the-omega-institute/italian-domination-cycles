#!/usr/bin/env python3
"""Exact rational-generating-function validation of the deficit formulas."""
from __future__ import annotations
import json
from math import comb
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent
n, y, z = sp.symbols("n y z")
Q = sp.Matrix([
    [0, 0, 0, 0, 1],
    [0, 0, 0, 1, 1],
    [1, 0, 0, 1, 1],
    [0, 1, 0, 1, 1],
    [0, 0, 1, 1, 1],
])
M = sp.diag(y**2, y**2, y**2, y, 1) * Q
S = sp.factor(sp.trace(z * M * (sp.eye(5) - z * M).inv()))
rows = {int(k): v for k, v in json.loads((ROOT / "coefficients.json").read_text()).items()}
thresholds = {0: 3, 1: 3, 2: 3, 3: 3, 4: 3, 5: 4, 6: 4, 7: 4, 8: 5, 9: 6, 10: 7}
out = {}
for j in range(11):
    G = sp.factor(sp.series(S, y, 0, j + 1).removeO().coeff(y, j))
    N = sp.Poly(sp.cancel(G * (1 - z) ** (j + 1)), z)
    coeffs = {int(exp[0]): int(value) for exp, value in N.terms()}
    def gf_value(nn):
        return sum(value * (comb(nn - r + j, j) if nn >= r else 0)
                   for r, value in coeffs.items())
    for nn in range(3, 101):
        k = 2 * nn - j
        actual = rows[nn][k] if 0 <= k < len(rows[nn]) else 0
        assert gf_value(nn) == actual, (j, nn)
    n0 = thresholds[j]
    poly = sp.factor(sp.interpolate(
        [(nn, rows[nn][2 * nn - j] if 0 <= 2 * nn - j < len(rows[nn]) else 0) for nn in range(n0, n0 + j + 1)], n))
    for nn in range(n0, 101):
        k = 2 * nn - j
        actual = rows[nn][k] if 0 <= k < len(rows[nn]) else 0
        assert int(poly.subs(n, nn)) == actual, (j, nn)
    if n0 > 3:
        k = 2 * (n0 - 1) - j
        actual = rows[n0 - 1][k] if 0 <= k < len(rows[n0 - 1]) else 0
        assert int(poly.subs(n, n0 - 1)) != actual
    out[str(j)] = {
        "n0": n0,
        "polynomial": str(sp.expand(poly)),
        "polynomial_factored": str(poly),
        "N_coefficients_low_to_high": [coeffs.get(r, 0) for r in range(max(coeffs) + 1)],
    }
(ROOT / "deficit_formulas.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("validated G_j coefficient extraction and thresholded polynomials for j=0..10, n=3..100")
