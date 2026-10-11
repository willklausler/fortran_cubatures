# Changelog

All notable changes are listed here. The project follows [Semantic Versioning](https://semver.org).

## [0.2.0] - 2026-10-11

- Fixed: `call q%set(q%elm, p)` stopped with an invalid index (debug builds) or read out of bounds, because `set` reset `q%elm` before using it.
- The Gauss-Jacobi Newton iteration stops with a message if it does not converge, instead of returning an inaccurate rule without any signal.
- The kind parameter `rk` is no longer public. Take `real64` from `iso_fortran_env` instead.
- New layout following the maintainer's package standard: `llms.txt` API reference, man page, `AGENTS.md`, `CITATION.cff`, `CONTRIBUTING.md`, three example programs and an app.
- fpm profiles with per-compiler flags; new make build that installs the library, module file and man page.
- CI adds ifx on Windows and LLVM flang; documentation is deployed by a separate workflow.

## [0.1.0] - 2025-09-04

First public release.

- Type `cubature` with constructor and `set` for lines, triangles, quadrilaterals, tetrahedra, hexahedra and wedges, selected by degree of exactness, with per-direction degrees on tensor-product elements.
- Tensor-product Gauss-Legendre rules, tabulated symmetric triangle and tetrahedron rules, and collapsed Gauss-Jacobi rules for any degree.
- `is_valid`, `summary`, `show` and `destroy`.
