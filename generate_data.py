#!/usr/bin/env python3
"""Generate the exact Italian-domination coefficient table for C_n, 3<=n<=100.

The five-state weighted transfer matrix gives the primary table.  The same
polynomial powers are checked against the nine-state labelled transfer matrix
and against exhaustive labelled assignments for n<=12.
"""
from __future__ import annotations
import csv, hashlib, json
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# Polynomial is represented by a low-to-high list of Python integers.
def matmul(A, B):
    nr, mid, nc = len(A), len(B), len(B[0])
    C = [[[0] for _ in range(nc)] for _ in range(nr)]
    for i in range(nr):
        for k in range(mid):
            a = A[i][k]
            if a == [0]:
                continue
            for j in range(nc):
                b = B[k][j]
                if b == [0]:
                    continue
                out = C[i][j]
                need = len(a) + len(b) - 1
                if len(out) < need:
                    out.extend([0] * (need - len(out)))
                for u, au in enumerate(a):
                    if au:
                        for v, bv in enumerate(b):
                            if bv:
                                out[u + v] += au * bv
    return C

def trace(A):
    out = [0]
    for i in range(len(A)):
        a = A[i][i]
        if len(out) < len(a):
            out.extend([0] * (len(a) - len(out)))
        for k, ak in enumerate(a):
            out[k] += ak
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out

def identity(n):
    return [[[int(i == j)] for j in range(n)] for i in range(n)]

def weighted_matrix(Q, row_weights):
    A = [[[0] for _ in row] for row in Q]
    for i, row in enumerate(Q):
        for j, value in enumerate(row):
            if value:
                A[i][j] = [0] * row_weights[i] + [value]
    return A

def five_state_table():
    Q = [
        [0, 0, 0, 0, 1],
        [0, 0, 0, 1, 1],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 1, 1],
        [0, 0, 1, 1, 1],
    ]
    A = identity(5)
    weighted = weighted_matrix(Q, [0, 0, 0, 1, 2])
    rows = {}
    for n in range(1, 101):
        A = matmul(A, weighted)
        if n >= 3:
            rows[n] = trace(A)
    return rows

def nine_state_table_through_12():
    states = list(product(range(3), repeat=2))
    T = [[0] * 9 for _ in range(9)]
    for i, (a, b) in enumerate(states):
        for j, (bp, c) in enumerate(states):
            if bp == b and (b != 0 or a + c >= 2):
                T[i][j] = 1
    A = identity(9)
    weighted = weighted_matrix(T, [b for a, b in states])
    out = {}
    for n in range(1, 13):
        A = matmul(A, weighted)
        if n >= 3:
            out[n] = trace(A)
    return out

def brute_row(n):
    out = [0] * (2 * n + 1)
    for values in product(range(3), repeat=n):
        if all(values[i] != 0 or values[(i - 1) % n] + values[(i + 1) % n] >= 2
                   for i in range(n)):
            out[sum(values)] += 1
    return out

def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    rows = five_state_table()
    direct = nine_state_table_through_12()
    direct_checks = []
    for n in range(3, 13):
        assert rows[n] == direct[n]
        assert rows[n] == brute_row(n)
        direct_checks.append(n)

    csv_path = ROOT / "coefficients.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["n", "k", "d_I(C_n,k)"])
        for n in range(3, 101):
            for k, value in enumerate(rows[n]):
                writer.writerow([n, k, value])

    # Arrays retain zero coefficients and make the row length explicit.
    json_path = ROOT / "coefficients.json"
    json_path.write_text(json.dumps({str(n): rows[n] for n in range(3, 101)},
                                    sort_keys=True, separators=(",", ":")) + "\n")

    cert = {
        "object": "d_I(C_n,k)",
        "range": {"n_min": 3, "n_max": 100, "k_min": 0, "k_max": "2*n"},
        "primary_method": "exact integer polynomial powering of D(x)Q, D(x)=diag(1,1,1,x,x^2)",
        "state_matrix": "Q=[[0,0,0,0,1],[0,0,0,1,1],[1,0,0,1,1],[0,1,0,1,1],[0,0,1,1,1]]",
        "independent_method": "exact integer polynomial powering of the labelled 9-state transfer matrix",
        "bruteforce_method": "exhaustive values in {0,1,2}^n for n=3..12",
        "independent_checks_n": direct_checks,
        "sha256": {
            "coefficients.csv": sha256(csv_path),
            "coefficients.json": sha256(json_path),
        },
        "row_lengths": "2*n+1; CSV includes zero coefficients; JSON stores full arrays",
    }
    (ROOT / "certificate.json").write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
    print(json.dumps(cert, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
