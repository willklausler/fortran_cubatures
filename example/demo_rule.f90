program demo_rule
  !! Build rules, read their points and weights, and keep one rule per element type
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures
  implicit none
  integer, parameter :: mesh(*) = [CUB_TET, CUB_HEX, CUB_WED, CUB_TET]
  type(cubature) :: q, rules(6)
  real(dp) :: measure
  integer :: e

  q = cubature(CUB_TRI, 2)                   ! exact for total degree <= 2
  call q%show()                              ! 3 points and weights, sum 0.5
  print *, q%dim, q%npoints, q%is_valid()    ! 2 3 T

  call q%set(CUB_HEX, [1, 3, 5])             ! reuse q; one degree per direction
  print *, q%degree, q%npoints               ! 1 3 5 6: 1 x 2 x 3 points
  print *, q%abscissae(:, 2)                 ! 0 -0.577 0.775: first coordinate fastest

  ! One rule per element type, indexed by the CUB_* constants
  rules(CUB_TET) = cubature(CUB_TET, 4)      ! 27 points
  rules(CUB_HEX) = cubature(CUB_HEX, 4)      ! 27 points
  rules(CUB_WED) = cubature(CUB_WED, [4, 2]) ! 6 x 2 = 12 points
  do e = 1, size(mesh)
    associate (r => rules(mesh(e)))
      measure = sum(r%weights)               ! 1/6, 8, 1, 1/6
      print '(a, i0, a, i3, a, f9.6)', 'element ', e, ':', r%npoints, ' points, measure', measure
    end associate
  end do
end program demo_rule
