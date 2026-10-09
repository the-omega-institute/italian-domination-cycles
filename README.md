# Italian domination polynomials of cycles

Data and code for the Italian domination polynomial D_I(C_n, x) of the cycle C_n (vertices labelled by Z/nZ, rotations distinguished), accompanying correspondence on

> P. Shao and C. Zhao, *Counting Weight-k Italian Dominating Sets on Trees and Cycles*, arXiv:2610.00108v1.

An Italian dominating function is f : Z/nZ -> {0,1,2} with f(i) = 0 implying f(i-1) + f(i+1) >= 2; d_I(C_n, k) counts such f of total weight k.

## Files

| File | Content |
|---|---|
| `coefficients.csv` | d_I(C_n, k) for 3 <= n <= 100 and 0 <= k <= 2n, one row per (n, k), zeros included |
| `coefficients.json` | the same table, one array of length 2n+1 per n |
| `certificate.json` | records of the independent checks (nine-state transfer matrix and exhaustive enumeration, 3 <= n <= 12) |
| `deficit_formulas.json` | the polynomials d_I(C_n, 2n-j), j = 0..10, with their exact thresholds |
| `generate_data.py` | computes the table from D_I(C_n, x) = tr(Q(x)^n) with exact integer polynomial arithmetic and runs the checks |
| `validate_deficits.py` | extracts the deficit generating functions and verifies the formulas and thresholds |

## Method

With the five-state matrix Q (row classes {(0,0)}, {(1,0)}, {(2,0)}, {(a,1)}, {(a,2)}) and Q(x) = diag(1,1,1,x,x^2) Q, one has D_I(C_n, x) = tr(Q(x)^n). The total counts a(n) = D_I(C_n, 1) = 23, 60, 167, 467, 1297, ... satisfy a(n+5) = 2a(n+4) + 2a(n+3) + a(n+2) - a(n+1) - a(n). All arithmetic is exact; no floating point is used.

Run `python3 generate_data.py` and `python3 validate_deficits.py` (Python 3.10+, standard library).

Authors: Haobo Ma (ChronoAI Pte Ltd; The Omega Institute) and Wenlin Zhang (National University of Singapore; The Omega Institute).
