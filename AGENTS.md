# AGENTS.md

Instructions for coding agents working on this repository. The public API and
its conventions are in `llms.txt`; read it before changing behaviour.

## Layout

| Path | Role |
|---|---|
| `src/cubatures.f90` | **Source of truth**: the whole library, module `cubatures` |
| `test/test_cubatures.f90` | The tests: one program, `check(cond, name)` per assertion |
| `example/*.f90` | Demo programs; each is copied verbatim into `llms.txt` |
| `app/main.f90` | Error and point count against degree for exp(x + y + z) on every element, the fpm app |
| `llms.txt` | API reference for agents, and the online user manual |
| `man/cubatures.3`, `README.md`, `ford.md` | Human documentation |
| `tools/check_docs.py` | CI check: every public name and example is in `llms.txt` |

## Commands

```sh
pip install ford                     # once
fpm test --profile debug             # fast check, gfortran
fpm test --profile release
fpm run --example --all
fpm run                              # the app
python3 tools/check_docs.py
./ci.sh                              # gfortran, ifx, nvfortran, flang if installed
make test FC=ifx                     # make build
make doc                             # FORD HTML in build/doc
```

Done means: the tests pass for every compiler you have in both profiles, the
examples and the app run, and `tools/check_docs.py` passes.

## Decided conventions (do not change without the maintainer)

- One type, `cubature`, with public components `elm`, `dim`, `degree(3)`, `npoints`, `abscissae(dim, npoints)`, `weights(npoints)`. Users read the components directly.
- Element constants `CUB_LIN` to `CUB_WED` are 1 to 6 so they can index arrays of rules.
- Reference elements: [-1, 1]^d for line, quadrilateral and hexahedron; the unit simplex for triangle and tetrahedron; triangle x [-1, 1] for the wedge.
- Rules are selected by degree of exactness, never by number of points. Every rule has positive weights and interior points; a new rule that breaks this is not accepted.
- Point order: first coordinate fastest; 1D Gauss points from +1 down to -1; wedge triangle point fastest.
- Error policy: invalid arguments to the constructor or `set` stop with `error stop` and a message that starts with `cubature%set:`. The Gauss-Jacobi Newton iteration also stops, with the same prefix, if it hits its cap of 100 iterations (13 suffice up to 1000 points). The policy is also stated in `llms.txt` and the man page NOTES.
- Kinds are not exported (`rk` was removed in 0.2.0). Internally `rk` aliases `real64`.
- No module state. Everything is pure except the output procedures `summary` and `show`.

## Adding or changing a procedure

1. Implement it in `src/cubatures.f90` with a `!!` line; rule builders are private `pure` subroutines.
2. Bind it to the type (with a `!!` line) or add it to the `public` statement of its group.
3. Add checks to `test/test_cubatures.f90`; a new rule goes through `check_rule` over a sweep of degrees.
4. Document it in `llms.txt` (the CI check requires the name), the man page and `CHANGELOG.md`.
5. Run the tests in both profiles, the examples and `tools/check_docs.py`.

## Test rules

- Exactness is checked against closed-form monomial integrals with a relative tolerance of `64*epsilon`.
- Error stops cannot be tested.
- The tests must run under floating-point trapping and the sanitizers of the debug profile.

## Compiler quirks that shaped the code

None so far: the code needs no compiler-specific workaround.

## Commits

Commit as the maintainer, Will Klausler, with the GitHub no-reply address
`226649647+willklausler@users.noreply.github.com` (set in the repository's git config).
Do not add `Co-authored-by` trailers: the maintainer is the only contributor.
