You’re solving a square linear system over a prime field, and I need a solver that can also prove when the system is *not* uniquely solvable.

The input is `/app/system.json`. It is a single JSON object with only three top-level fields: `prime`, `rows`, and `rhs`. The value `prime` is a prime \(p\) with `2 <= p <= 1_000_000_007`. The field `rows` is an `n x n` integer matrix and `rhs` is a length-`n` integer list (with `n >= 1`). Together they mean \(A x = b\) modulo \(p\). Numbers may be negative or bigger than \(p\), so reduce everything modulo \(p\) before doing any arithmetic.

Write `/app/gf_solve.py` and run it with `python3`. It must accept `--input PATH` and `--output PATH`, and if you omit them it must default to `/app/system.json` and `/app/result.json`, so `python3 /app/gf_solve.py` works from `/app`.

Write JSON only to the output file, and keep the top-level keys exactly as specified.

If there is exactly one solution, write `{\"status\": \"ok\", \"solution\": [...]}` where every entry is an integer residue in `[0, p)`.

If the system does not have exactly one solution, write `{\"status\": \"singular\", \"kind\": ..., \"certificate\": ...}` and provide a certificate that matches the rules below.

If the system is inconsistent (no solutions), set `kind` to `inconsistent` and write `{\"certificate\": {\"y\": [...]}}`. The vector `y` must have length `n`, must use residues in `[0, p)`, and must satisfy \(y^T A = 0\) but \(y^T b \\ne 0\) modulo \(p\). If there are many valid `y`, you must output the canonical one produced by a deterministic elimination run: do Gauss–Jordan elimination on the augmented matrix `[A|b]` using the leftmost-available pivot in each column (scan rows top-to-bottom), track each current row as a linear combination of the original input rows, and when the elimination first produces a contradiction row `[0 ... 0 | c]` with `c != 0`, output that row’s combination coefficients as `y`.

If the system is non-unique (rank-deficient but consistent), set `kind` to `non-unique` and write `{\"certificate\": {\"x0\": [...], \"z\": [...]}}`. The vector `x0` must be one valid solution in residues `[0, p)`. The vector `z` must be a nonzero nullspace vector in residues `[0, p)`, meaning \(A z = 0\) modulo \(p\), and `z` must not be the all-zero vector. If there are many valid choices, you must output the canonical pair from the same elimination run: set all free variables to 0 to get `x0`, then pick the smallest-index free column `f`, set `z[f]=1`, set all other free variables to 0, and solve for the pivot variables from the reduced system.

Use only the Python standard library. Do not call `eval`, `exec`, `compile`, or `__import__`. Do not import: `fractions`, `decimal`, `numpy`, `sympy`, `gmpy2`, `importlib`, `inspect`, `ctypes`, `subprocess`, `multiprocessing`, `pickle`, `marshal`, `os`, `builtins`, `types`, `code`, `sqlite3`, `zlib`, `base64`, `ssl`, `socket`. You may use `math.gcd` and `pow`; for modular inverse when a value is nonzero mod \(p\), `pow(a, p - 2, p)` is fine.

The included `/app/system.json` is only a smoke test. Hidden tests use larger systems, large primes, negative coefficients, and lots of singular cases; the certificate rules are what make this task tricky, so keep everything exact and modular (no floats).