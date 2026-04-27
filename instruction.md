We need a small solver for linear equations over a prime modulus, but we also need it to explain *why* a system fails to have a unique solution.

The input file is `/app/system.json`. 
It is a single JSON object with only three top-level fields: `prime`, `rows`, and `rhs`. 
Here `prime` is a prime number \(p\) with `2 <= p <= 1_000_000_007`. 
The field `rows` is an `n x n` integer matrix and `rhs` is a length-`n` integer list (with `n >= 1`). 
Treat this as the modular system \(A x = b \pmod p\). 
Coefficients in the JSON may be negative or larger than \(p\), so reduce everything modulo \(p\) before doing the math.

Please write `/app/gf_solve.py` (run it with `python3`). 
It should accept `--input PATH` and `--output PATH`. If you omit those flags, it must default to `/app/system.json` and `/app/result.json`, so `python3 /app/gf_solve.py` works from `/app`.

Write JSON only to the output file, and keep the top-level keys exactly as specified below.

If there is exactly one solution vector, write `{\"status\": \"ok\", \"solution\": [...]}` where every entry is an integer in `[0, p)`.

If the system does *not* have exactly one solution, you must still explain which case you are in and provide a certificate.

If the system is inconsistent (no solutions), write `{\"status\": \"singular\", \"kind\": \"inconsistent\", \"certificate\": {\"y\": [...]}}`. 
The vector `y` must have length `n`, must use residues in `[0, p)`, and must satisfy \(y^T A = 0\) and \(y^T b \\ne 0\) modulo \(p\).

If the system is non-unique (rank-deficient but consistent), write `{\"status\": \"singular\", \"kind\": \"non-unique\", \"certificate\": {\"x0\": [...], \"z\": [...]}}`. The vector `x0` must be a valid solution with entries in `[0, p)`. The vector `z` must be a nonzero nullspace vector (entries in `[0, p)`), meaning \(A z = 0\) modulo \(p\), and it must not be the all-zero vector.

Use only the Python standard library. 
Do not call `eval`, `exec`, `compile`, or `__import__`. Do not import: `fractions`, `decimal`, `numpy`, `sympy`, `gmpy2`, `importlib`, `inspect`, `ctypes`, `subprocess`, `multiprocessing`, `pickle`, `marshal`, `os`, `builtins`, `types`, `code`, `sqlite3`, `zlib`, `base64`, `ssl`, `socket`. You may use `math.gcd` and `pow`; for modular inverse when a value is nonzero mod \(p\), `pow(a, p - 2, p)` is fine.

The included `/app/system.json` is only a small smoke test. 
Hidden tests include larger systems, large primes, negative coefficients, and many singular cases; the certificate rules above are what make this task tricky. Make sure the solver uses exact modular arithmetic, not floats.