#!/bin/bash
set -euo pipefail

cat > /app/gf_solve.py << 'PY'
#!/usr/bin/env python3
"""Reference solver: Ax = b over GF(p) via Gauss–Jordan elimination."""

from __future__ import annotations

import argparse
import json
from typing import Any


def _mod(x: int, p: int) -> int:
    return (x % p + p) % p


def _classify_with_certificate(
    A: list[list[int]], b: list[int], p: int
) -> tuple[str, dict[str, object]]:
    n = len(A)
    if n == 0 or any(len(row) != n for row in A) or len(b) != n:
        return "inconsistent", {"certificate": {"y": [0] * max(n, 1)}}

    M = [[_mod(A[i][j], p) for j in range(n)] + [_mod(b[i], p)] for i in range(n)]
    T = [[0] * n for _ in range(n)]
    for i in range(n):
        T[i][i] = 1

    pivot_cols: list[int] = []
    row = 0
    for col in range(n):
        pivot: int | None = None
        for r in range(row, n):
            if M[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        if pivot != row:
            M[row], M[pivot] = M[pivot], M[row]
            T[row], T[pivot] = T[pivot], T[row]
        inv = pow(M[row][col], p - 2, p)
        for j in range(col, n + 1):
            M[row][j] = (M[row][j] * inv) % p
        for j in range(n):
            T[row][j] = (T[row][j] * inv) % p
        for r in range(n):
            if r == row:
                continue
            f = M[r][col]
            if f:
                for j in range(col, n + 1):
                    M[r][j] = (M[r][j] - f * M[row][j]) % p
                for j in range(n):
                    T[r][j] = (T[r][j] - f * T[row][j]) % p
        pivot_cols.append(col)
        row += 1
        if row == n:
            break

    for r in range(n):
        if all(M[r][c] == 0 for c in range(n)) and M[r][n] != 0:
            y = [int(T[r][i]) for i in range(n)]
            return "inconsistent", {"certificate": {"y": y}}

    rank = len(pivot_cols)
    if rank == n:
        x = [0] * n
        for r, col in enumerate(pivot_cols):
            x[col] = int(M[r][n] % p)
        return "ok", {"solution": x}

    x0 = [0] * n
    for r, col in enumerate(pivot_cols):
        x0[col] = int(M[r][n] % p)

    free_cols = [c for c in range(n) if c not in pivot_cols]
    fc = free_cols[0]
    z = [0] * n
    z[fc] = 1
    for r, col in enumerate(pivot_cols):
        z[col] = int((-M[r][fc]) % p)
    return "non-unique", {"certificate": {"x0": x0, "z": z}}


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
    kind, payload = _classify_with_certificate(rows, rhs, p)
    if kind == "ok":
        out = {"status": "ok", "solution": payload["solution"]}
    elif kind == "inconsistent":
        out = {"status": "singular", "kind": "inconsistent", "certificate": payload["certificate"]}
    else:
        out = {"status": "singular", "kind": "non-unique", "certificate": payload["certificate"]}
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(out, f)
        f.write("\n")


if __name__ == "__main__":
    main()
PY

chmod +x /app/gf_solve.py
python3 /app/gf_solve.py
