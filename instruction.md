I am solving a square linear system over a prime field, and I need a solver that can also prove when the system is not uniquely solvable.

Write `/app/gf_solve.py` to solve square systems \(A x = b\) over \(GF(p)\), and to produce a deterministic certificate when the system is not uniquely solvable.

Input is `/app/system.json` with keys `prime`, `rows`, `rhs`: \(p\) is prime (`2 <= p <= 1_000_000_007`), `rows` is an `n x n` integer matrix \(A\), and `rhs` is a length-`n` integer list \(b\) (`n >= 1`). Reduce all numbers modulo \(p\) before doing any arithmetic.

Run it as `python3 /app/gf_solve.py [--input PATH] [--output PATH]`. Defaults are `/app/system.json` and `/app/result.json`. Write JSON only to the output file, and keep keys exactly as specified below.

If there is exactly one solution, write `{"status":"ok","solution":[...]}` with residues in `[0,p)`.

Otherwise write `{"status":"singular","kind":K,"certificate":C}`.

If `K="inconsistent"`, then `C={"y":[...]}` where `y` has length `n`, residues in `[0,p)`, and satisfies \(y^T A = 0\) but \(y^T b \ne 0\) mod \(p\). If multiple `y` work, output the canonical `y` from this exact Gauss–Jordan run on `[A|b]`: for each column, take the leftmost available pivot by scanning rows top-to-bottom; track each current row as a linear combination of the original rows; and when elimination first produces a contradiction row `[0...0|c]` with `c!=0`, output that row’s combination coefficients as `y`.

If `K="non-unique"`, then `C={"x0":[...],"z":[...]}` where `x0` is a solution (residues in `[0,p)`), and `z` is a nonzero nullspace vector (residues in `[0,p)`) with \(A z = 0\) mod \(p\). If multiple choices exist, output the canonical pair from the same elimination run: set all free variables to 0 to get `x0`, then pick the smallest-index free column `f`, set `z[f]=1` (all other free vars 0), and solve for pivot variables from the reduced system.

Use only the Python standard library. Do not call `eval`, `exec`, `compile`, or `__import__`. Do not import `fractions`, `decimal`, `numpy`, `sympy`, `gmpy2`, `importlib`, `inspect`, `ctypes`, `subprocess`, `multiprocessing`, `pickle`, `marshal`, `os`, `builtins`, `types`, `code`, `sqlite3`, `zlib`, `base64`, `ssl`, `socket`. For modular inverse of nonzero values mod \(p\), `pow(a, p - 2, p)` is fine. Keep everything exact (no floats).