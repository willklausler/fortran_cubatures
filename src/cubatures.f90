module cubatures
!! Numerical integration rules (cubatures) on the reference elements of
!! the finite element method.
!!
!! A [[cubature]] holds the abscissae \(\xi_q\) and weights \(w_q\) of a
!! rule such that
!! \[ \int_{\hat\Omega} f(\xi)\,d\xi \approx \sum_{q=1}^{n} w_q\, f(\xi_q) \]
!! on the reference element \(\hat\Omega\).
!!
!! ## Element types and reference domains
!!
!! | Constant  | Element       | Reference domain                          | Measure |
!! |-----------|---------------|-------------------------------------------|---------|
!! | `CUB_LIN` | line          | \([-1,1]\)                                | 2       |
!! | `CUB_TRI` | triangle      | \(x,y \ge 0,\ x+y \le 1\)                 | 1/2     |
!! | `CUB_QUA` | quadrilateral | \([-1,1]^2\)                              | 4       |
!! | `CUB_TET` | tetrahedron   | \(x,y,z \ge 0,\ x+y+z \le 1\)             | 1/6     |
!! | `CUB_HEX` | hexahedron    | \([-1,1]^3\)                              | 8       |
!! | `CUB_WED` | wedge / prism | triangle \(\times\ [-1,1]\)               | 1       |
!!
!! Triangle and tetrahedron abscissae are the Cartesian coordinates on the
!! unit simplex, which equal the leading barycentric (area / volume)
!! coordinates; the last barycentric coordinate is \(1 - \sum_i \xi_i\).
!!
!! ## Degree of exactness
!!
!! Rules are selected by their *degree of exactness* \(p\):
!!
!! - `CUB_TRI`, `CUB_TET`: all polynomials of total degree \(\le p\).
!! - `CUB_LIN`, `CUB_QUA`, `CUB_HEX`: all polynomials of degree \(\le p_i\)
!!   in each coordinate \(i\) (tensor-product Gauss–Legendre with
!!   \(p_i/2 + 1\) points per direction). One degree or one per direction.
!! - `CUB_WED`: total degree \(\le p_1\) in the triangle coordinates and
!!   degree \(\le p_2\) in the axial coordinate. One or two degrees.
!!
!! All rules have positive weights and interior abscissae. Triangle rules of
!! degree 2, 4 and 5 and the tetrahedron rule of degree 2 are the tabulated
!! symmetric rules; all other simplex rules are collapsed (conical product)
!! Gauss–Jacobi rules, available for any degree.
!!
!! ## Point ordering
!!
!! Tensor-product rules vary the first coordinate fastest. Wedge rules
!! vary the triangle point fastest and the axial point slowest.
!!
!! ## References
!!
!! 1. Abramowitz, M., Stegun, I. A. (1964). *Handbook of Mathematical
!!    Functions*, §22 and §25.4. National Bureau of Standards.
!! 2. Stroud, A. H. (1971). *Approximate Calculation of Multiple Integrals*.
!!    Prentice-Hall.
!! 3. Radon, J. (1948). Zur mechanischen Kubatur. *Monatshefte für
!!    Mathematik*, 52, 286–300.
!! 4. Dunavant, D. A. (1985). High degree efficient symmetrical Gaussian
!!    quadrature rules for the triangle. *International Journal for
!!    Numerical Methods in Engineering*, 21(6), 1129–1148.
!! 5. Karniadakis, G. E., Sherwin, S. J. (2005). *Spectral/hp Element
!!    Methods for Computational Fluid Dynamics*, 2nd ed.
!!    Oxford University Press.
!! 6. Zienkiewicz, O. C., Taylor, R. L., Zhu, J. Z. (2005). *The Finite
!!    Element Method: Its Basis and Fundamentals*, 6th ed.
!!    Butterworth-Heinemann.

  use iso_fortran_env, only: rk => real64, output_unit

  implicit none

  private

  public :: rk
  public :: cubature
  public :: CUB_LIN, CUB_TRI, CUB_QUA, CUB_TET, CUB_HEX, CUB_WED

  integer, parameter :: CUB_LIN = 1 !! Line element
  integer, parameter :: CUB_TRI = 2 !! Triangle element
  integer, parameter :: CUB_QUA = 3 !! Quadrilateral element
  integer, parameter :: CUB_TET = 4 !! Tetrahedron element
  integer, parameter :: CUB_HEX = 5 !! Hexahedron element
  integer, parameter :: CUB_WED = 6 !! Wedge (prism) element

  character(3), parameter :: names(6) = ["LIN", "TRI", "QUA", "TET", "HEX", "WED"]
    !! Element names, indexed by element type
  integer, parameter :: dims(6) = [1, 2, 2, 3, 3, 3]
    !! Spatial dimension, indexed by element type
  integer, parameter :: ndegrees(6) = [1, 1, 2, 1, 3, 2]
    !! Number of independent degrees, indexed by element type

  real(rk), parameter :: pi = acos(-1.0_rk)

  type :: cubature
    !! Abscissae and weights of an integration rule on a reference element.
    !!
    !! Create with the constructor, `q = cubature(CUB_HEX, 3)`, or in place
    !! with `call q%set(CUB_HEX, 3)`. Integrate with
    !! `sum(q%weights * f(q%abscissae))` or a loop over `q%npoints`.

    integer :: elm = 0                        !! Element type, `CUB_*`
    integer :: dim = 0                        !! Spatial dimension
    integer :: degree(3) = 0                  !! Degree of exactness, see module documentation
    integer :: npoints = 0                    !! Number of points
    real(rk), allocatable :: abscissae(:,:)   !! Abscissae, shape `[dim, npoints]`
    real(rk), allocatable :: weights(:)       !! Weights, shape `[npoints]`

  contains

    generic :: set => set_iso, set_aniso      !! Build the rule in place
    procedure, private :: set_iso
    procedure, private :: set_aniso
    procedure :: is_valid
    procedure :: summary
    procedure :: show
    procedure :: destroy

  end type cubature

  interface cubature
    !! Construct a rule, e.g. `cubature(CUB_TRI, 4)` or `cubature(CUB_HEX, [3, 3, 1])`
    module procedure new_iso
    module procedure new_aniso
  end interface cubature

contains

!***********************************************************************

pure function new_iso(elm, degree) result(self)
!! Construct a rule of degree `degree` in every direction

  integer, intent(in) :: elm       !! Element type, `CUB_*`
  integer, intent(in) :: degree    !! Degree of exactness
  type(cubature) :: self

  call self%set(elm, [degree])

end function new_iso

!***********************************************************************

pure function new_aniso(elm, degree) result(self)
!! Construct a rule with per-direction degrees

  integer, intent(in) :: elm         !! Element type, `CUB_*`
  integer, intent(in) :: degree(:)   !! Degrees of exactness, size 1 or number of directions
  type(cubature) :: self

  call self%set(elm, degree)

end function new_aniso

!***********************************************************************

pure subroutine set_iso(self, elm, degree)
!! Build a rule of degree `degree` in every direction

  class(cubature), intent(inout) :: self
  integer, intent(in) :: elm       !! Element type, `CUB_*`
  integer, intent(in) :: degree    !! Degree of exactness

  call self%set(elm, [degree])

end subroutine set_iso

!***********************************************************************

pure subroutine set_aniso(self, elm, degree)
!! Build a rule with per-direction degrees.
!!
!! `degree` has size 1 (same degree in every direction) or `ndegrees(elm)`:
!! 1 for `CUB_LIN`, `CUB_TRI`, `CUB_TET`; 2 for `CUB_QUA`, `CUB_WED`
!! (triangle, axial); 3 for `CUB_HEX`. Stops on invalid input.

  class(cubature), intent(inout) :: self
  integer, intent(in) :: elm         !! Element type, `CUB_*`
  integer, intent(in) :: degree(:)   !! Degrees of exactness

  integer :: p(3)   ! Degrees, padded with 0
  integer :: n      ! Number of independent degrees

  if (elm < 1 .or. elm > size(names)) error stop "cubature%set: invalid element type"
  if (any(degree < 0)) error stop "cubature%set: degree must be non-negative"

  n = ndegrees(elm)
  p = 0
  if (size(degree) == 1) then
    p(1:n) = degree(1)
  else if (size(degree) == n) then
    p(1:n) = degree
  else
    error stop "cubature%set: size(degree) must be 1 or the number of element directions"
  end if

  call self%destroy()
  self%elm = elm
  self%dim = dims(elm)
  self%degree = p

  select case (elm)
  case (CUB_LIN, CUB_QUA, CUB_HEX)
    call tensor_rule(p(1:n), self%abscissae, self%weights)
  case (CUB_TRI)
    call triangle_rule(p(1), self%abscissae, self%weights)
  case (CUB_TET)
    call tetrahedron_rule(p(1), self%abscissae, self%weights)
  case (CUB_WED)
    call wedge_rule(p(1), p(2), self%abscissae, self%weights)
  end select

  self%npoints = size(self%weights)

end subroutine set_aniso

!***********************************************************************

pure logical function is_valid(self)
!! True if the rule is set and its arrays are consistent

  class(cubature), intent(in) :: self

  is_valid = allocated(self%abscissae) .and. allocated(self%weights)
  if (.not. is_valid) return
  is_valid = self%npoints > 0 &
       .and. size(self%weights) == self%npoints &
       .and. all(shape(self%abscissae) == [self%dim, self%npoints])

end function is_valid

!***********************************************************************

subroutine summary(self, unit)
!! Write element type, dimension, degrees and number of points

  class(cubature), intent(in) :: self
  integer, intent(in), optional :: unit   !! Output unit, default `output_unit`

  integer :: u

  u = output_unit
  if (present(unit)) u = unit

  if (self%elm == 0) then
    write(u,"(A)") "cubature: not set"
    return
  end if

  write(u,"(A,A)")              "Element:   ", names(self%elm)
  write(u,"(A,I0)")             "Dimension: ", self%dim
  write(u,"(A,*(I0,:,', '))")   "Degree:    ", self%degree(1:ndegrees(self%elm))
  write(u,"(A,I0)")             "Points:    ", self%npoints

end subroutine summary

!***********************************************************************

subroutine show(self, unit)
!! Write the summary, every abscissa and weight, and the weight sum

  class(cubature), intent(in) :: self
  integer, intent(in), optional :: unit   !! Output unit, default `output_unit`

  integer :: u, q

  u = output_unit
  if (present(unit)) u = unit

  call self%summary(u)
  if (.not. self%is_valid()) return

  do q = 1, self%npoints
    write(u,"(I4,')',*(1X,ES23.15E3))") q, self%abscissae(:,q), self%weights(q)
  end do
  write(u,"(A,ES23.15E3)") "Sum of weights: ", sum(self%weights)

end subroutine show

!***********************************************************************

pure subroutine destroy(self)
!! Deallocate and reset to the unset state

  class(cubature), intent(inout) :: self

  self%elm     = 0
  self%dim     = 0
  self%degree  = 0
  self%npoints = 0
  if (allocated(self%abscissae)) deallocate(self%abscissae)
  if (allocated(self%weights))   deallocate(self%weights)

end subroutine destroy

!***********************************************************************

pure subroutine tensor_rule(p, x, w)
!! Tensor product of Gauss–Legendre rules, first direction fastest

  integer, intent(in) :: p(:)                        !! Degree per direction
  real(rk), allocatable, intent(out) :: x(:,:)       !! Abscissae
  real(rk), allocatable, intent(out) :: w(:)         !! Weights

  integer :: n(size(p))                              ! Points per direction
  real(rk) :: x1(maxval(p)/2+1, size(p))             ! 1D abscissae per direction
  real(rk) :: w1(maxval(p)/2+1, size(p))             ! 1D weights per direction
  integer :: i, j, q, r

  n = p/2 + 1
  do i = 1, size(p)
    call gauss_jacobi(n(i), 0, x1(:n(i),i), w1(:n(i),i))
  end do

  allocate(x(size(p), product(n)), w(product(n)))
  do q = 1, size(w)
    r = q - 1
    w(q) = 1
    do i = 1, size(p)
      j = mod(r, n(i)) + 1
      r = r/n(i)
      x(i,q) = x1(j,i)
      w(q) = w(q)*w1(j,i)
    end do
  end do

end subroutine tensor_rule

!***********************************************************************

pure subroutine triangle_rule(p, x, w)
!! Rule of total degree `p` on the unit triangle

  integer, intent(in) :: p                           !! Degree of exactness
  real(rk), allocatable, intent(out) :: x(:,:)       !! Abscissae, shape `[2, n]`
  real(rk), allocatable, intent(out) :: w(:)         !! Weights

  select case (p)

  ! 3 points, interior (Zienkiewicz et al. 2005)
  case (2)
    allocate(x(2,3), w(3))
    x = reshape([4, 1, 1, 4, 1, 1], [2, 3])/6.0_rk
    w = 1.0_rk/6

  ! 6 points (Dunavant 1985)
  case (4)
    block
      real(rk), parameter :: a = 0.44594849091596488632_rk, wa = 0.22338158967801146570_rk/2
      real(rk), parameter :: b = 0.09157621350977074346_rk, wb = 0.10995174365532186764_rk/2
      allocate(x(2,6), w(6))
      x(:,1:3) = orbit3(a)
      x(:,4:6) = orbit3(b)
      w = [wa, wa, wa, wb, wb, wb]
    end block

  ! 7 points (Radon 1948)
  case (5)
    block
      real(rk), parameter :: s15 = sqrt(15.0_rk)
      real(rk), parameter :: a = (6 - s15)/21, wa = (155 - s15)/2400
      real(rk), parameter :: b = (6 + s15)/21, wb = (155 + s15)/2400
      allocate(x(2,7), w(7))
      x(:,1) = 1.0_rk/3
      x(:,2:4) = orbit3(a)
      x(:,5:7) = orbit3(b)
      w = [9.0_rk/80, wa, wa, wa, wb, wb, wb]
    end block

  case default
    call collapsed_triangle(p, x, w)

  end select

end subroutine triangle_rule

!***********************************************************************

pure function orbit3(a) result(x)
!! The 3 triangle points with barycentric coordinates \((a, a, 1-2a)\) permuted

  real(rk), intent(in) :: a
  real(rk) :: x(2,3)

  x = reshape([a, a, 1-2*a, a, a, 1-2*a], [2, 3])

end function orbit3

!***********************************************************************

pure subroutine collapsed_triangle(p, x, w)
!! Collapsed Gauss–Jacobi rule of total degree `p` on the unit triangle.
!!
!! Maps \((u, v) \in [0,1]^2\) to \((x, y) = (u, (1-u)v)\), Jacobian
!! \(1-u\), absorbed by Gauss–Jacobi \(\alpha = 1\) in \(u\)
!! (Stroud 1971; Karniadakis & Sherwin 2005).

  integer, intent(in) :: p                           !! Degree of exactness
  real(rk), allocatable, intent(out) :: x(:,:)       !! Abscissae, shape `[2, n**2]`
  real(rk), allocatable, intent(out) :: w(:)         !! Weights

  real(rk) :: u(p/2+1), wu(p/2+1), v(p/2+1), wv(p/2+1)
  integer :: i, j, n, q

  n = p/2 + 1
  call gauss_jacobi(n, 1, u, wu)
  call gauss_jacobi(n, 0, v, wv)
  u = (1 + u)/2;  wu = wu/4
  v = (1 + v)/2;  wv = wv/2

  allocate(x(2,n*n), w(n*n))
  q = 0
  do i = 1, n
    do j = 1, n
      q = q + 1
      x(:,q) = [u(i), (1 - u(i))*v(j)]
      w(q) = wu(i)*wv(j)
    end do
  end do

end subroutine collapsed_triangle

!***********************************************************************

pure subroutine tetrahedron_rule(p, x, w)
!! Rule of total degree `p` on the unit tetrahedron

  integer, intent(in) :: p                           !! Degree of exactness
  real(rk), allocatable, intent(out) :: x(:,:)       !! Abscissae, shape `[3, n]`
  real(rk), allocatable, intent(out) :: w(:)         !! Weights

  integer :: i

  select case (p)

  ! 4 points, interior (Zienkiewicz et al. 2005)
  case (2)
    block
      real(rk), parameter :: a = (5 + 3*sqrt(5.0_rk))/20
      real(rk), parameter :: b = (5 -   sqrt(5.0_rk))/20
      allocate(x(3,4), w(4))
      x = b
      do i = 1, 3
        x(i,i) = a
      end do
      w = 1.0_rk/24
    end block

  case default
    call collapsed_tetrahedron(p, x, w)

  end select

end subroutine tetrahedron_rule

!***********************************************************************

pure subroutine collapsed_tetrahedron(p, x, w)
!! Collapsed Gauss–Jacobi rule of total degree `p` on the unit tetrahedron.
!!
!! Maps \((u, v, t) \in [0,1]^3\) to
!! \((x, y, z) = (u, (1-u)v, (1-u)(1-v)t)\), Jacobian \((1-u)^2(1-v)\),
!! absorbed by Gauss–Jacobi \(\alpha = 2\) in \(u\) and \(\alpha = 1\) in
!! \(v\) (Stroud 1971; Karniadakis & Sherwin 2005).

  integer, intent(in) :: p                           !! Degree of exactness
  real(rk), allocatable, intent(out) :: x(:,:)       !! Abscissae, shape `[3, n**3]`
  real(rk), allocatable, intent(out) :: w(:)         !! Weights

  real(rk), dimension(p/2+1) :: u, wu, v, wv, t, wt
  integer :: i, j, k, n, q

  n = p/2 + 1
  call gauss_jacobi(n, 2, u, wu)
  call gauss_jacobi(n, 1, v, wv)
  call gauss_jacobi(n, 0, t, wt)
  u = (1 + u)/2;  wu = wu/8
  v = (1 + v)/2;  wv = wv/4
  t = (1 + t)/2;  wt = wt/2

  allocate(x(3,n**3), w(n**3))
  q = 0
  do i = 1, n
    do j = 1, n
      do k = 1, n
        q = q + 1
        x(:,q) = [u(i), (1 - u(i))*v(j), (1 - u(i))*(1 - v(j))*t(k)]
        w(q) = wu(i)*wv(j)*wt(k)
      end do
    end do
  end do

end subroutine collapsed_tetrahedron

!***********************************************************************

pure subroutine wedge_rule(ptri, plin, x, w)
!! Product of a triangle rule and a Gauss–Legendre rule, triangle fastest

  integer, intent(in) :: ptri                        !! Total degree in the triangle
  integer, intent(in) :: plin                        !! Degree along the axis
  real(rk), allocatable, intent(out) :: x(:,:)       !! Abscissae, shape `[3, n]`
  real(rk), allocatable, intent(out) :: w(:)         !! Weights

  real(rk), allocatable :: xt(:,:), wt(:)
  real(rk) :: z(plin/2+1), wz(plin/2+1)
  integer :: k, nt

  call triangle_rule(ptri, xt, wt)
  call gauss_jacobi(size(z), 0, z, wz)

  nt = size(wt)
  allocate(x(3, nt*size(z)), w(nt*size(z)))
  do k = 1, size(z)
    x(1:2, (k-1)*nt+1:k*nt) = xt
    x(3,   (k-1)*nt+1:k*nt) = z(k)
    w(     (k-1)*nt+1:k*nt) = wt*wz(k)
  end do

end subroutine wedge_rule

!***********************************************************************

pure subroutine gauss_jacobi(n, alpha, x, w)
!! `n`-point Gauss–Jacobi rule on \([-1,1]\) for the weight
!! \((1-x)^\alpha\), exact to degree \(2n-1\). `alpha = 0` is Gauss–Legendre.
!!
!! Roots of \(P_n^{(\alpha,0)}\) by Newton iteration with deflation;
!! weights \(w_i = 2^{\alpha+1} / \big((1-x_i^2)\,P_n'(x_i)^2\big)\)
!! (Abramowitz & Stegun 1964, 25.4.33 with \(\beta = 0\)).

  integer, intent(in) :: n                 !! Number of points, \(\ge 1\)
  integer, intent(in) :: alpha             !! Jacobi exponent, \(\ge 0\)
  real(rk), intent(out) :: x(n)            !! Abscissae
  real(rk), intent(out) :: w(n)            !! Weights

  real(rk) :: z, dz, pn, dpn
  integer :: i, it

  do i = 1, n
    z = cos(pi*(i - 0.25_rk)/(n + 0.5_rk))   ! Legendre root estimate
    do it = 1, 100
      call jacobi(n, alpha, z, pn, dpn)
      dz = pn/dpn
      dz = dz/(1 - dz*sum(1/(z - x(1:i-1))))   ! Deflate roots already found
      z = z - dz
      if (abs(dz) <= 2*epsilon(z)) exit
    end do
    call jacobi(n, alpha, z, pn, dpn)
    x(i) = z
    w(i) = 2.0_rk**(alpha+1)/((1 - z**2)*dpn**2)
  end do

end subroutine gauss_jacobi

!***********************************************************************

pure subroutine jacobi(n, alpha, x, pn, dpn)
!! Jacobi polynomial \(P_n^{(\alpha,0)}(x)\), \(n \ge 1\), and its
!! derivative, by three-term recurrence (Abramowitz & Stegun 1964, 22.7.1
!! and 22.8.1)

  integer, intent(in) :: n                 !! Degree
  integer, intent(in) :: alpha             !! Jacobi exponent
  real(rk), intent(in) :: x                !! Argument, \(|x| < 1\)
  real(rk), intent(out) :: pn              !! \(P_n(x)\)
  real(rk), intent(out) :: dpn             !! \(P_n'(x)\)

  real(rk) :: a, c, p0, p1
  integer :: k

  a = alpha
  p0 = 1
  p1 = ((a + 2)*x + a)/2
  do k = 2, n
    c = 2*k + a
    pn = ((c - 1)*(c*(c - 2)*x + a*a)*p1 - 2*(k + a - 1)*(k - 1)*c*p0) &
       / (2*k*(k + a)*(c - 2))
    p0 = p1
    p1 = pn
  end do
  pn = p1
  c = 2*n + a
  dpn = (n*(a - c*x)*pn + 2*(n + a)*n*p0)/(c*(1 - x**2))

end subroutine jacobi

!***********************************************************************

end module cubatures
