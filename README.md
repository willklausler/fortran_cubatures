# fortran_cubatures

[![CI](https://github.com/willklausler/fortran_cubatures/actions/workflows/ci.yml/badge.svg)](https://github.com/willklausler/fortran_cubatures/actions/workflows/ci.yml)
[![Docs](https://github.com/willklausler/fortran_cubatures/actions/workflows/docs.yml/badge.svg)](https://willklausler.github.io/fortran_cubatures/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Gaussian integration rules on the reference elements of the finite element
method: one derived type, `cubature`, that holds the points and weights of a
rule chosen by its degree of exactness.

```fortran
use cubatures
type(cubature) :: q

q = cubature(CUB_TRI, 4)                     ! exact for total degree <= 4
print *, q%npoints, sum(q%weights)           ! 6  0.5: the area of the triangle
print *, sum(q%weights*q%abscissae(1,:)**4)  ! 3.3333333333333340E-002 = 1/30

q = cubature(CUB_HEX, [3, 3, 1])             ! degree 3, 3 and 1 per direction
print *, q%npoints, q%dim                    ! 4  3
```

> Coding agents: the complete API reference is [llms.txt](llms.txt), and
> contributor instructions are in [AGENTS.md](AGENTS.md).

## Features

- Lines, triangles, quadrilaterals, tetrahedra, hexahedra and wedges (prisms).
- Rules are selected by degree of exactness, for any degree.
- Separate degrees per direction on quadrilaterals, hexahedra and wedges.
- All weights positive, all points strictly inside the element.
- Tabulated symmetric rules for low-degree triangles and tetrahedra, collapsed Gauss-Jacobi rules above.
- Pure procedures, no module state, no dependencies.

Conventions: [-1, 1]^d for tensor-product elements, the unit simplex for triangles and tetrahedra, first coordinate fastest, and invalid input stops the program.

## Installation

With [fpm](https://fpm.fortran-lang.org) 0.13 or later, add to your `fpm.toml`:

```toml
[dependencies]
fortran_cubatures = { git = "https://github.com/willklausler/fortran_cubatures" }
```

With make, which builds `build/make/libfortran_cubatures.a` and the module file:

```sh
make FC=gfortran && make test
make install PREFIX=$HOME/.local    # library, .mod file and man page
```

## Documentation

- [User manual and API documentation](https://willklausler.github.io/fortran_cubatures/), generated with FORD (`make doc`)
- [llms.txt](llms.txt): the API reference as one plain-text file, also the manual's source
- `man 3 cubatures` after `make install` ([man/cubatures.3](man/cubatures.3))
- Example programs in [example/](example/) (`fpm run --example --all`) and an app in [app/](app/) (`fpm run`) that tabulates the error and point count against degree for exp(x + y + z) on every element

## Compilers

CI tests every push on:

| Compiler | Linux | macOS | Windows |
| --- | --- | --- | --- |
| gfortran 13, 14 | yes | yes | yes |
| ifx 2026 | yes | | yes |
| LLVM flang 21 | yes | | |
| nvfortran 26 | yes | | |

There are no known compiler limitations; see [llms.txt](llms.txt#compiler-notes).

## Why another cubature package?

No other Fortran package gives finite-element cubature rules on the standard
reference elements.
[stdlib](https://github.com/fortran-lang/stdlib) (`stdlib_quadrature`),
[quadpack](https://github.com/jacobwilliams/quadpack) and
[GaussJacobiQuad](https://github.com/HaoZeke/GaussJacobiQuad) cover one
dimension only.
[quadrature-fortran](https://github.com/jacobwilliams/quadrature-fortran)
integrates a user function adaptively over boxes. It does not return points
and weights, and it has no simplices.
[triasymq](https://github.com/zgimbutas/triasymq) has high-order symmetric
rules for triangles, squares, cubes and tetrahedra. It is fixed-form
FORTRAN 77 with a separate makefile per element, and it has no line or wedge
rules and no fpm support. This package covers all six FE elements, including
wedges, through one derived type. It takes any degree of exactness, with
anisotropic degrees on tensor elements. Every rule has positive weights and
interior points, is verified by monomial-exactness tests, and is packaged for
fpm and tested on several compilers.

## References

1. M. Abramowitz, I. A. Stegun, *Handbook of Mathematical Functions*, §22 and §25.4, National Bureau of Standards, 1964.
2. A. H. Stroud, *Approximate Calculation of Multiple Integrals*, Prentice-Hall, 1971.
3. J. Radon, Zur mechanischen Kubatur, *Monatshefte für Mathematik* 52, 286–300, 1948.
4. D. A. Dunavant, High degree efficient symmetrical Gaussian quadrature rules for the triangle, *International Journal for Numerical Methods in Engineering* 21(6), 1129–1148, 1985.
5. G. E. Karniadakis, S. J. Sherwin, *Spectral/hp Element Methods for Computational Fluid Dynamics*, 2nd ed., Oxford University Press, 2005.
6. O. C. Zienkiewicz, R. L. Taylor, J. Z. Zhu, *The Finite Element Method: Its Basis and Fundamentals*, 6th ed., Butterworth-Heinemann, 2005.

## Status and support

Version 0.2.0. The API may still change before 1.0; changes are recorded in
[CHANGELOG.md](CHANGELOG.md). Bug reports, questions and pull requests are
welcome on [GitHub](https://github.com/willklausler/fortran_cubatures/issues);
see [CONTRIBUTING.md](CONTRIBUTING.md). To cite this package, use
[CITATION.cff](CITATION.cff).

## License

[MIT](LICENSE)
