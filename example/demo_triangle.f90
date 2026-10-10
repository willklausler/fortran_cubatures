program demo_triangle
  !! Integrate over a physical triangle through the affine map from the reference triangle
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures
  implicit none
  real(dp), parameter :: xv(2,3) = reshape([1._dp, 0._dp, &
                                            4._dp, 1._dp, &
                                            2._dp, 3._dp], [2, 3])  ! vertices
  type(cubature) :: q
  real(dp) :: jac(2,2), detj, x(2), total, exact
  integer :: k

  q = cubature(CUB_TRI, 2)                   ! x*y has degree 2 after an affine map
  jac(:,1) = xv(:,2) - xv(:,1)               ! x = xv(:,1) + jac*xi
  jac(:,2) = xv(:,3) - xv(:,1)               ! det(jac) = 8 = twice the area
  detj = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)

  total = 0
  do k = 1, q%npoints
    x = xv(:,1) + matmul(jac, q%abscissae(:,k))
    total = total + x(1)*x(2)*detj*q%weights(k)
  end do

  ! Closed form: area/12*(sum(x_k*y_k) + sum(x_k)*sum(y_k))
  exact = detj/24*(sum(xv(1,:)*xv(2,:)) + sum(xv(1,:))*sum(xv(2,:)))
  print *, total, exact                      ! 12.666666666666664 12.666666666666666
end program demo_triangle
