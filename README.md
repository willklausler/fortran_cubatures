# fortran_cubatures

[![CI](https://github.com/willklausler/fortran_cubatures/actions/workflows/ci.yml/badge.svg)](https://github.com/willklausler/fortran_cubatures/actions/workflows/ci.yml)

Numerical integration rules (quadrature / cubature) on the reference elements of
the finite element method, packaged as one derived type, `cubature`.

- Lines, triangles, quadrilaterals, tetrahedra, hexahedra and wedges (prisms)
- Rules are selected by **degree of exactness**, for any degree
- Anisotropic degrees for quadrilaterals, hexahedra and wedges
- All weights positive, all points strictly inside the element
- Pure procedures, no dependencies, Fortran 2008

## Installation

Add the package to your `fpm.toml`:

```toml
[dependencies]
fortran_cubatures = { git = "https://github.com/willklausler/fortran_cubatures" }
```

Without fpm, compile the single source file `src/cubatures.f90` with your project.

## Quick start

```fortran
use cubatures, only: cubature, rk, CUB_HEX

type(cubature) :: q
real(rk) :: integral
integer :: k

q = cubature(CUB_HEX, 3)          ! exact for degree <= 3 in each coordinate
integral = 0
do k = 1, q%npoints
  integral = integral + f(q%abscissae(:,k))*q%weights(k)
end do
```

A typical element loop evaluates shape functions at each point and scales by the
Jacobian determinant. For example, the consistent mass matrix of a bilinear
quadrilateral:

```fortran
q = cubature(CUB_QUA, 3)          ! N_i N_j det(J) has degree 3 per direction
m = 0
do k = 1, q%npoints
  xi = q%abscissae(:,k)
  call shape_functions(xi, n, dn) ! N(4), dN/dxi(4,2)
  jac = matmul(xnodes, dn)        ! dx/dxi(2,2)
  m = m + spread(n,2,4)*spread(n,1,4)*det2(jac)*q%weights(k)
end do
```

Codes with mixed meshes can keep one rule per element type, indexed by the
element constants:

```fortran
type(cubature) :: rules(6)

rules(CUB_TET) = cubature(CUB_TET, 4)
rules(CUB_WED) = cubature(CUB_WED, [4, 2])   ! degree 4 in the triangle, 2 along the axis
...
associate (q => rules(element_type(e)))
  ...
end associate
```

The runnable version of these snippets is in
[example/cubatures_example.f90](example/cubatures_example.f90):

```sh
fpm run --example
```

## API

| Entity | Description |
| --- | --- |
| `cubature(elm, degree)` | Constructor. `degree` is an integer or an array with one degree per direction |
| `call q%set(elm, degree)` | Build in place, reusing the object |
| `q%elm` | Element type, one of the `CUB_*` constants |
| `q%dim` | Spatial dimension |
| `q%degree(3)` | Degrees of exactness, padded with 0 |
| `q%npoints` | Number of points |
| `q%abscissae(dim, npoints)` | Point coordinates in the reference element |
| `q%weights(npoints)` | Weights |
| `q%is_valid()` | `.true.` if the rule is set and consistent |
| `call q%summary([unit])` | Print element, dimension, degrees and number of points |
| `call q%show([unit])` | Print the summary, every point and weight, and the weight sum |
| `call q%destroy()` | Deallocate and reset |
| `rk` | Real kind of the abscissae and weights (`real64`) |

Invalid input (unknown element, negative degree, wrong number of degrees)
stops with `error stop` and a message.

## Elements and rules

| Constant | Element | Reference domain | Measure | `degree` | Rule | Points |
| --- | --- | --- | --- | --- | --- | --- |
| `CUB_LIN` | line | [-1, 1] | 2 | `p` | Gauss–Legendre | n(p) |
| `CUB_QUA` | quadrilateral | [-1, 1]² | 4 | `p` or `[px, py]` | tensor Gauss–Legendre | n(px)·n(py) |
| `CUB_HEX` | hexahedron | [-1, 1]³ | 8 | `p` or `[px, py, pz]` | tensor Gauss–Legendre | n(px)·n(py)·n(pz) |
| `CUB_TRI` | triangle | x, y ≥ 0, x + y ≤ 1 | 1/2 | `p` | see below | see below |
| `CUB_TET` | tetrahedron | x, y, z ≥ 0, x + y + z ≤ 1 | 1/6 | `p` | see below | see below |
| `CUB_WED` | wedge / prism | triangle × [-1, 1] | 1 | `p` or `[ptri, pz]` | triangle ⊗ Gauss–Legendre | T(ptri)·n(pz) |

Here n(p) = ⌊p/2⌋ + 1 is the number of Gauss points that integrate degree p exactly.

What "degree p" means depends on the element:

- **Triangle, tetrahedron:** all polynomials of total degree ≤ p.
- **Line, quadrilateral, hexahedron:** all polynomials of degree ≤ pᵢ in each
  coordinate i. This includes total degree ≤ min pᵢ.
- **Wedge:** total degree ≤ ptri in (x, y) times degree ≤ pz in z.

Simplex rules use the most efficient positive rule available:

| Degree p | Triangle points T(p) | Rule | Tetrahedron points | Rule |
| --- | --- | --- | --- | --- |
| 0, 1 | 1 | centroid | 1 | centroid |
| 2 | 3 | symmetric, interior [6] | 4 | symmetric [6] |
| 3 | 4 | collapsed Gauss–Jacobi | 8 | collapsed Gauss–Jacobi |
| 4 | 6 | symmetric [4] | 27 | collapsed Gauss–Jacobi |
| 5 | 7 | symmetric [3] | 27 | collapsed Gauss–Jacobi |
| ≥ 6 | n(p)² | collapsed Gauss–Jacobi [2, 5] | n(p)³ | collapsed Gauss–Jacobi [2, 5] |

The collapsed (conical product) rules map the unit square or cube onto the
simplex and use Gauss–Jacobi points to absorb the Jacobian of that map. Gauss
points are computed at run time by Newton iteration on Legendre and Jacobi
polynomials [1], so there is no upper limit on the degree.

### Conventions

- Triangle and tetrahedron abscissae are Cartesian coordinates on the unit
  simplex. These equal the leading barycentric (area or volume) coordinates, and
  the last barycentric coordinate is 1 − Σ ξᵢ.
- Tensor-product rules vary the first coordinate fastest. Wedge rules vary the
  triangle point fastest and the axial point slowest.
- For an element of polynomial order k with an affine map, the mass matrix needs
  degree 2k and the stiffness matrix needs 2k − 2. Curved or distorted elements
  need more.

## Building, testing and documentation

```sh
fpm build
fpm test --profile debug     # bounds checks, sanitizers, FP traps
fpm test --profile release
fpm run --example
ford ford.md                 # API documentation in docs/
```

The unit tests check every rule over a sweep of degrees. Each rule must
integrate every monomial in its polynomial space exactly (compared against
closed-form integrals), have positive weights, and have points inside the
element. The tests also cover point counts, point ordering, the constructor and
`set` interfaces, reuse and output.

CI covers gfortran 13–15, Intel ifx and NVIDIA nvfortran on Linux, and gfortran
on macOS and Windows. It runs both profiles on every push.

## Pedigree and support

The tensor-product and low-degree simplex rules are the standard rules of
finite element textbooks [6]. The degree 4 and 5 triangle rules are those of
Dunavant [4] and Radon [3]. Every rule is verified by the unit tests described
above. The code has not been independently reviewed. It is maintained by the
author on a best-effort basis.

Please report bugs and request features through
[GitHub issues](https://github.com/willklausler/fortran_cubatures/issues).
Pull requests are welcome; please include a test.

## References

1. Abramowitz, M., Stegun, I. A. (1964). *Handbook of Mathematical Functions*, §22, §25.4. National Bureau of Standards.
2. Stroud, A. H. (1971). *Approximate Calculation of Multiple Integrals*. Prentice-Hall.
3. Radon, J. (1948). Zur mechanischen Kubatur. *Monatshefte für Mathematik*, 52, 286–300.
4. Dunavant, D. A. (1985). High degree efficient symmetrical Gaussian quadrature rules for the triangle. *International Journal for Numerical Methods in Engineering*, 21(6), 1129–1148.
5. Karniadakis, G. E., Sherwin, S. J. (2005). *Spectral/hp Element Methods for Computational Fluid Dynamics*, 2nd ed. Oxford University Press.
6. Zienkiewicz, O. C., Taylor, R. L., Zhu, J. Z. (2005). *The Finite Element Method: Its Basis and Fundamentals*, 6th ed. Butterworth-Heinemann.

## License

[MIT](LICENSE) © 2025 Will Klausler
