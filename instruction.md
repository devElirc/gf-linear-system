We need a small solver for linear equations over a prime modulus.

The input file is `/app/system.json`. It should contain only these top-level fields: `prime`, `rows`, and `rhs`.

`prime` is a real prime number between `2` and `1_000_000_007`.

`rows` is an `n x n` matrix of integers, and `rhs` is a list with the same length. `n` is at least `1`.

Treat this as solving:

`A x = b (mod p)`

The numbers in the input may be negative or larger than `p`. Reduce everything modulo `p` before doing the math.

Please write `/app/gf_solve.py`. It should run with `python3`.

The script should accept:

`--input PATH`  
`--output PATH`

If those flags are not provided, use:

`/app/system.json`  
`/app/result.json`

So this should work from `/app`:

`python3 /app/gf_solve.py`

Write JSON only to the output file.

If the system has exactly one solution, write:

`{"status": "ok", "solution": [...]}`

Each value in `solution` must be an integer from `0` to `p - 1`.

If the system has no solution, or if it has more than one solution, write exactly:

`{"status": "singular"}`

Do not include a `solution` key or any extra fields in that case.

Use only the Python standard library.

Do not call `eval`, `exec`, `compile`, or `__import__`.

Do not import these modules: `fractions`, `decimal`, `numpy`, `sympy`, `gmpy2`, `importlib`, `inspect`, `ctypes`, `subprocess`, `multiprocessing`, `pickle`, `marshal`, `os`, `builtins`, `types`, `code`, `sqlite3`, `zlib`, `base64`, `ssl`, or `socket`.

You can use `math.gcd` and `pow`. For modular inverse with a nonzero value, `pow(a, p - 2, p)` is fine.

The included `/app/system.json` is only a small smoke test. The hidden tests will include larger systems, large primes, negative coefficients, values bigger than `p`, singular systems, and inconsistent systems. Make sure the solver uses exact modular arithmetic, not floats.