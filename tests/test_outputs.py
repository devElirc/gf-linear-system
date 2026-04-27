"""Verifier for gf-linear-system — GF(p) linear solve with strict JSON and AST rules."""

import ast
import json
import random
import subprocess
import tempfile
from pathlib import Path


GF_SOLVE = Path("/app/gf_solve.py")
DEFAULT_SYSTEM = Path("/app/system.json")
DEFAULT_RESULT = Path("/app/result.json")

FORBIDDEN_IMPORTS = {
    "fractions",
    "decimal",
    "numpy",
    "sympy",
    "gmpy2",
    "importlib",
    "inspect",
    "ctypes",
    "subprocess",
    "multiprocessing",
    "pickle",
    "marshal",
    "os",
    "builtins",
    "types",
    "code",
    "sqlite3",
    "zlib",
    "base64",
    "ssl",
    "socket",
}
FORBIDDEN_CALLS = {"eval", "exec", "compile", "__import__"}


def _mod(x: int, p: int) -> int:
    return (x % p + p) % p


def _dot(u: list[int], v: list[int], p: int) -> int:
    s = 0
    for a, b in zip(u, v):
        s = (s + _mod(a, p) * _mod(b, p)) % p
    return s


def _matmul_vec_left(y: list[int], A: list[list[int]], p: int) -> list[int]:
    """Return y^T A as a length-n vector."""
    n = len(A)
    out = [0] * n
    for j in range(n):
        s = 0
        for i in range(n):
            s = (s + _mod(y[i], p) * _mod(A[i][j], p)) % p
        out[j] = s
    return out


def _oracle_classify(
    A: list[list[int]], b: list[int], p: int
) -> tuple[str, dict[str, object]]:
    """
    Return (kind, certificate_payload) where kind is one of:
    - "ok" with {"solution": x}
    - "inconsistent" with {"certificate": {"y": y}}
    - "non-unique" with {"certificate": {"x0": x0, "z": z}}

    Oracle constructs certificates by tracking row operations.
    """
    n = len(A)
    if n == 0 or any(len(row) != n for row in A) or len(b) != n:
        return "inconsistent", {"certificate": {"y": [0] * max(n, 1)}}

    # Augmented matrix [A | b] plus a row-tracking matrix T so each row is a combination of original rows:
    # current_rows = T * original_rows
    M = [[_mod(A[i][j], p) for j in range(n)] + [_mod(b[i], p)] for i in range(n)]
    T = [[0] * n for _ in range(n)]
    for i in range(n):
        T[i][i] = 1

    pivot_cols: list[int] = []
    row = 0
    for col in range(n):
        pivot = None
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

    # Detect inconsistency: a zero row in A but nonzero in b.
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

    # Non-unique but consistent. Build one solution x0 (free vars = 0), and a nullspace vector z.
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


def _matvec(A: list[list[int]], x: list[int], p: int) -> list[int]:
    n = len(A)
    out: list[int] = []
    for i in range(n):
        s = 0
        for j in range(n):
            s = (s + (A[i][j] % p + p) % p * (x[j] % p + p) % p) % p
        out.append(s)
    return out


def _assert_ok_payload(data: object) -> None:
    assert isinstance(data, dict), "result must be a JSON object"
    assert set(data.keys()) == {"status", "solution"}, (
        "ok result must contain only status and solution, got "
        f"{sorted(data.keys())!r}"
    )
    assert data["status"] == "ok"
    sol = data["solution"]
    assert isinstance(sol, list), "solution must be a JSON array"
    assert all(isinstance(v, int) for v in sol), "solution entries must be integers"


def _assert_singular_payload(data: object) -> None:
    assert isinstance(data, dict), "result must be a JSON object"
    assert data.get("status") == "singular", f"expected status=singular, got {data!r}"
    assert set(data.keys()) == {"status", "kind", "certificate"}, (
        "singular result must contain only status/kind/certificate keys, got "
        f"{sorted(data.keys())!r}"
    )
    assert data["kind"] in ("inconsistent", "non-unique")
    cert = data["certificate"]
    assert isinstance(cert, dict), "certificate must be an object"
    if data["kind"] == "inconsistent":
        assert set(cert.keys()) == {"y"}
        y = cert["y"]
        assert isinstance(y, list) and all(isinstance(v, int) for v in y)
    else:
        assert set(cert.keys()) == {"x0", "z"}
        x0 = cert["x0"]
        z = cert["z"]
        assert isinstance(x0, list) and all(isinstance(v, int) for v in x0)
        assert isinstance(z, list) and all(isinstance(v, int) for v in z)


def _check_ast_constraints() -> None:
    src = GF_SOLVE.read_text(encoding="utf-8", errors="strict")
    tree = ast.parse(src, filename=str(GF_SOLVE))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                top = alias.name.split(".")[0]
                assert top not in FORBIDDEN_IMPORTS, f"forbidden import: {top}"
        elif isinstance(node, ast.ImportFrom) and node.module:
            top = node.module.split(".")[0]
            assert top not in FORBIDDEN_IMPORTS, f"forbidden import from: {top}"
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            assert node.func.id not in FORBIDDEN_CALLS, (
                f"forbidden call: {node.func.id}()"
            )


def _run_gf_solve(
    inp: Path,
    out: Path,
    *,
    timeout_sec: float = 120.0,
) -> dict:
    proc = subprocess.run(
        [
            "python3",
            str(GF_SOLVE),
            "--input",
            str(inp),
            "--output",
            str(out),
        ],
        capture_output=True,
        text=True,
        timeout=timeout_sec,
    )
    assert proc.returncode == 0, f"gf_solve.py failed:\n{proc.stdout}\n{proc.stderr}"
    assert out.is_file(), "output not written"
    data = json.loads(out.read_text(encoding="utf-8"))
    assert isinstance(data, dict), "output must be a JSON object"
    if data.get("status") == "ok":
        _assert_ok_payload(data)
    elif data.get("status") == "singular":
        _assert_singular_payload(data)
    else:
        raise AssertionError(f"invalid status in output: {data!r}")
    return data


def test_gf_solve_script_exists():
    """Verify /app/gf_solve.py is present before running any test."""
    assert GF_SOLVE.is_file(), "Missing /app/gf_solve.py"


def test_ast_forbidden_imports_and_calls():
    """AST-check that gf_solve.py uses no forbidden modules/builtins."""
    assert GF_SOLVE.is_file()
    _check_ast_constraints()


def test_bundled_default_paths():
    """Shipped system.json is identity mod 17; default CLI must write canonical solution."""
    assert GF_SOLVE.is_file()
    if DEFAULT_RESULT.exists():
        DEFAULT_RESULT.unlink()
    proc = subprocess.run(
        ["python3", str(GF_SOLVE)],
        capture_output=True,
        text=True,
        timeout=120,
        cwd="/app",
    )
    assert proc.returncode == 0, proc.stderr
    data = json.loads(DEFAULT_RESULT.read_text(encoding="utf-8"))
    _assert_ok_payload(data)
    assert data["solution"] == [3, 5, 11]


def test_cli_override_paths():
    """--input / --output flags redirect I/O away from defaults."""
    spec = {"prime": 13, "rows": [[2, 0], [0, 3]], "rhs": [1, 1]}
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        inp, outp = td / "in.json", td / "out.json"
        inp.write_text(json.dumps(spec), encoding="utf-8")
        data = _run_gf_solve(inp, outp)
    x = data["solution"]
    assert _matvec(spec["rows"], x, 13) == [(spec["rhs"][i] % 13) for i in range(2)]


def test_singular_duplicate_equations():
    """Rank-deficient consistent system has many solutions → singular."""
    spec = {
        "prime": 11,
        "rows": [[1, 1], [1, 1]],
        "rhs": [2, 2],
    }
    assert _oracle_solve(spec["rows"], spec["rhs"], spec["prime"]) is None
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        inp, outp = td / "s.json", td / "r.json"
        inp.write_text(json.dumps(spec), encoding="utf-8")
        data = _run_gf_solve(inp, outp)
    _assert_singular_payload(data)


def test_singular_inconsistent():
    """Parallel rows disagree → no solution."""
    spec = {
        "prime": 7,
        "rows": [[1, 1], [1, 1]],
        "rhs": [1, 2],
    }
    assert _oracle_solve(spec["rows"], spec["rhs"], 7) is None
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        inp, outp = td / "s.json", td / "r.json"
        inp.write_text(json.dumps(spec), encoding="utf-8")
        data = _run_gf_solve(inp, outp)
    _assert_singular_payload(data)


def test_one_by_one_zero_rhs():
    """1×1 system with trivial solution x=0 mod 2."""
    spec = {"prime": 2, "rows": [[1]], "rhs": [0]}
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        inp, outp = td / "s.json", td / "r.json"
        inp.write_text(json.dumps(spec), encoding="utf-8")
        data = _run_gf_solve(inp, outp)
    _assert_ok_payload(data)
    assert data["solution"] == [0]
    assert _matvec(spec["rows"], data["solution"], 2) == [0]


def test_negative_entries_canonicalized():
    """Negative coefficients and rhs values are reduced mod p first."""
    spec = {
        "prime": 19,
        "rows": [[-1, 0], [0, 1]],
        "rhs": [5, -40],
    }
    exp = _oracle_solve(spec["rows"], spec["rhs"], 19)
    assert exp is not None
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        inp, outp = td / "s.json", td / "r.json"
        inp.write_text(json.dumps(spec), encoding="utf-8")
        data = _run_gf_solve(inp, outp)
    assert data["solution"] == exp
    for v in data["solution"]:
        assert 0 <= v < 19


def test_large_prime_unimodular_triangular():
    """Strictly upper perturbation of I has determinant 1 mod p."""
    p = 1_000_000_007
    n = 14
    rng = random.Random(77_001)
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = 1
    for i in range(n):
        for j in range(i + 1, n):
            A[i][j] = rng.randrange(p)
    x = [rng.randrange(p) for _ in range(n)]
    b = _matvec(A, x, p)
    spec = {"prime": p, "rows": A, "rhs": b}
    exp = _oracle_solve(A, b, p)
    assert exp is not None
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        inp, outp = td / "s.json", td / "r.json"
        inp.write_text(json.dumps(spec), encoding="utf-8")
        data = _run_gf_solve(inp, outp, timeout_sec=180.0)
    assert data["solution"] == exp == x


def _random_invertible_lu(rng: random.Random, n: int, p: int) -> list[list[int]]:
    L = [[0] * n for _ in range(n)]
    U = [[0] * n for _ in range(n)]
    for i in range(n):
        L[i][i] = 1
        for j in range(i):
            L[i][j] = rng.randrange(p)
        while True:
            U[i][i] = rng.randrange(1, p)
            if U[i][i] % p != 0:
                break
        for j in range(i + 1, n):
            U[i][j] = rng.randrange(p)
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = 0
            for k in range(n):
                s = (s + L[i][k] * U[k][j]) % p
            A[i][j] = s
    return A


def test_random_invertible_systems():
    """Random LU-factored invertible systems across several primes."""
    primes = [97, 251, 1009, 10_007, 65_537, 1_000_003, 1_000_000_007]
    assert GF_SOLVE.is_file()
    for trial in range(110):
        rng = random.Random(900_000 + trial)
        p = rng.choice(primes)
        n = rng.randint(5, 26)
        A = _random_invertible_lu(rng, n, p)
        x = [rng.randrange(p) for _ in range(n)]
        b = _matvec(A, x, p)
        kind, payload = _oracle_classify(A, b, p)
        assert kind == "ok", f"trial {trial}: expected ok, got {kind}"
        exp = payload["solution"]
        assert exp == x, f"trial {trial}: internal oracle mismatch"
        spec = {"prime": p, "rows": A, "rhs": b}
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            inp, outp = td / "s.json", td / "r.json"
            inp.write_text(json.dumps(spec), encoding="utf-8")
            data = _run_gf_solve(inp, outp, timeout_sec=90.0)
        assert data["status"] == "ok"
        got = data["solution"]
        assert got == exp, f"trial {trial}: mismatch p={p} n={n}"
        assert _matvec(A, got, p) == [(b[i] % p + p) % p for i in range(n)]


def test_random_singular_systems_with_certificates():
    """Random singular systems must provide a correct inconsistency/nullspace certificate."""
    assert GF_SOLVE.is_file()
    primes = [97, 1009, 65_537, 1_000_003]
    for trial in range(160):
        rng = random.Random(2_100_000 + trial)
        p = rng.choice(primes)
        n = rng.randint(8, 34)
        # Build A with a forced dependency: last row = sum of a few previous rows.
        A = [[rng.randrange(p) for _ in range(n)] for _ in range(n)]
        picks = rng.sample(range(n - 1), k=rng.randint(2, min(6, n - 1)))
        for j in range(n):
            A[n - 1][j] = sum(A[i][j] for i in picks) % p

        if rng.random() < 0.5:
            # Consistent non-unique: choose x, set b = A x.
            x = [rng.randrange(p) for _ in range(n)]
            b = _matvec(A, x, p)
            expected_kind = "non-unique"
        else:
            # Inconsistent: pick b with the same dependency violated.
            b = [rng.randrange(p) for _ in range(n)]
            lhs = sum(b[i] for i in picks) % p
            b[n - 1] = (lhs + rng.randrange(1, p)) % p
            expected_kind = "inconsistent"

        spec = {"prime": p, "rows": A, "rhs": b}
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            inp, outp = td / "s.json", td / "r.json"
            inp.write_text(json.dumps(spec), encoding="utf-8")
            data = _run_gf_solve(inp, outp, timeout_sec=120.0)

        assert data["status"] == "singular"
        assert data["kind"] == expected_kind
        cert = data["certificate"]
        if expected_kind == "inconsistent":
            y = cert["y"]
            assert len(y) == n
            ya = _matmul_vec_left(y, A, p)
            assert all(v % p == 0 for v in ya), "y^T A must be 0"
            assert _dot(y, b, p) % p != 0, "y^T b must be nonzero"
        else:
            x0 = cert["x0"]
            z = cert["z"]
            assert len(x0) == n and len(z) == n
            assert _matvec(A, x0, p) == [v % p for v in b], "A x0 must equal b"
            assert _matvec(A, z, p) == [0] * n, "A z must be 0"
            assert any(v % p != 0 for v in z), "z must be nonzero"


def test_singular_near_full_rank_mod2():
    """Four identical rows force deficiency mod 2."""
    spec = {
        "prime": 2,
        "rows": [
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [1, 1, 0, 0],
        ],
        "rhs": [0, 0, 0, 0],
    }
    assert _oracle_solve(spec["rows"], spec["rhs"], 2) is None
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        inp, outp = td / "s.json", td / "r.json"
        inp.write_text(json.dumps(spec), encoding="utf-8")
        data = _run_gf_solve(inp, outp)
    _assert_singular_payload(data)
