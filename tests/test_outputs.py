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


def _oracle_solve(A: list[list[int]], b: list[int], p: int) -> list[int] | None:
    """Must match reference implementation in solution/solve.sh (Gauss–Jordan)."""
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
    assert data == {"status": "singular"}, (
        "singular result must be exactly {\"status\": \"singular\"}, got "
        f"{data!r}"
    )


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
    assert GF_SOLVE.is_file(), "Missing /app/gf_solve.py"


def test_ast_forbidden_imports_and_calls():
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
    primes = [97, 251, 1009, 10_007, 65_537, 1_000_003, 1_000_000_007]
    assert GF_SOLVE.is_file()
    for trial in range(110):
        rng = random.Random(900_000 + trial)
        p = rng.choice(primes)
        n = rng.randint(5, 26)
        A = _random_invertible_lu(rng, n, p)
        x = [rng.randrange(p) for _ in range(n)]
        b = _matvec(A, x, p)
        exp = _oracle_solve(A, b, p)
        assert exp is not None, f"trial {trial}: oracle failed unexpectedly"
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
