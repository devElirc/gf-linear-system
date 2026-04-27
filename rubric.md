## Rubric 1

Agent reads the problem statement and identifies all three distinct output cases (unique solution, inconsistent, non-unique) before writing any code, +2
Agent creates /app/gf_solve.py with --input and --output CLI arguments that correctly default to /app/system.json and /app/result.json, +3
Agent implements Gauss–Jordan elimination using only Python standard library modules (e.g., argparse, json, pow for modular inverse) and no forbidden imports, +3
Agent tracks row-combination coefficients through every elimination step to derive the canonical y certificate when the first contradiction row appears in the inconsistent case, +3
Agent constructs the canonical (x0, z) certificate for non-unique systems by setting all free variables to zero for x0 and choosing the smallest-index free column to build z, +3
Agent runs python3 /app/gf_solve.py against the smoke-test input and inspects /app/result.json to confirm correct output before declaring the task complete, +2
Agent imports a forbidden library (numpy, sympy, fractions, decimal, or any other module listed in the constraints) inside gf_solve.py, -5
Agent uses floating-point arithmetic anywhere in the solver instead of exact integer modular arithmetic, -3
Agent outputs a singular certificate that is not canonical (even if it is mathematically valid), -4

