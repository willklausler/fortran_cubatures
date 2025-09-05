# Fortran cubatures

## About

Numerical integration derived type for finite elements in

- 1D - lineature, line elements,
- 2D - quadrature, triangle and quadrilateral elements, and
- 3D - cubature, tetrahedron, hexahedron, and prism/wedge elements.

Unit tests verify sum of weights, abscissae means, and polynomial integration.

## Usage

```fortran
! Declare type
type(cubature) :: scheme
character(3) :: geo = "HEX"
integer :: order = 2

! Set cubature with geometry and order
! geometry: LINe, TRIangle, QUAdrilateral, TETrahedron, HEXahedron, WEJ (prism)
! order: Integration order (e.g. 1 = single point)
call scheme%set(geo, order)

! Access element type
print *, scheme%elmtype

! Access number of points
print *, scheme%points

! Access dimension
print *, scheme%dime

! Access weights, shape = [points]
print *, scheme%weights

! Access abscissae aka coordinates, shape = [dime, points]
print *, scheme%abscissae

! Display scheme
call scheme%show()
```

## To do

- Unit testing for prisms

- Polynomial integration testing with non-separable functions

- Anisotropic integration for quadrilaterals, hexahedrons, and prisms
