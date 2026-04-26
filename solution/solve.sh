#!/bin/bash
set -euo pipefail

cat > /app/gf_solve.py << 'PY'
#!/usr/bin/env python3
"""Reference solver: Ax = b over GF(p) via Gauss–Jordan elimination."""

from __future__ import annotations

import argparse
import json
from typing import Any


def _solve_mod(A: list[list[int]], b: list[int], p: int) -> list[int] | None:
    n = len(A)
    if n == 0 or any(len(row) != n for row in A) or len(b) != n:
        return None
    M: list[list[int]] = []
    for i in range(n):
        M.append([(A[i][j] % p + p) % p for j in range(n)] + [(b[i] % p + p) % p])
    for col in range(n):
        pivot: int | None = None
        for r in range(col, n):
            if M[r][col] % p != 0:
                pivot = r
                break
        if pivot is None:
            return None
        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]
        inv = pow(M[col][col], p - 2, p)
        for j in range(col, n + 1):
            M[col][j] = (M[col][j] * inv) % p
        for r in range(n):
            if r == col:
                continue
            f = M[r][col] % p
            if f:
                for j in range(col, n + 1):
                    M[r][j] = (M[r][j] - f * M[col][j]) % p
    x = [M[i][n] % p for i in range(n)]
    for i in range(n):
        s = sum((A[i][j] % p + p) % p * x[j] for j in range(n)) % p
        if s != (b[i] % p + p) % p:
            return None
    return x


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="/app/system.json")
    ap.add_argument("--output", default="/app/result.json")
    args = ap.parse_args()
    with open(args.input, encoding="utf-8") as f:
        spec: dict[str, Any] = json.load(f)
    p = int(spec["prime"])
    rows = spec["rows"]
    rhs = spec["rhs"]
    sol = _solve_mod(rows, rhs, p)
    if sol is None:
        out = {"status": "singular"}
    else:
        out = {"status": "ok", "solution": sol}
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(out, f)
        f.write("\n")


if __name__ == "__main__":
    main()
PY

chmod +x /app/gf_solve.py
python3 /app/gf_solve.py
