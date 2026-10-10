program demo_mass_matrix
  !! Consistent mass matrix of a distorted bilinear quadrilateral (isoparametric map)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures
  implicit none
  real(dp), parameter :: xn(2,4) = reshape([0._dp, 0._dp, &
                                            2._dp, 0._dp, &
                                            2.5_dp, 1.5_dp, &
                                            0._dp, 1._dp], [2, 4])  ! nodes
  real(dp), parameter :: sn(2,4) = reshape([-1, -1, 1, -1, 1, 1, -1, 1], [2, 4])  ! reference nodes
  type(cubature) :: q
  real(dp) :: n(4), dn(4,2), jac(2,2), detj, m(4,4)
  integer :: k

  q = cubature(CUB_QUA, 3)                   ! N_i*N_j*det(J) has degree 3 per direction
  m = 0
  do k = 1, q%npoints                        ! 2 x 2 Gauss points
    associate (xi => q%abscissae(:,k))
      n = (1 + sn(1,:)*xi(1))*(1 + sn(2,:)*xi(2))/4
      dn(:,1) = sn(1,:)*(1 + sn(2,:)*xi(2))/4
      dn(:,2) = sn(2,:)*(1 + sn(1,:)*xi(1))/4
    end associate
    jac = matmul(xn, dn)
    detj = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
    m = m + spread(n, 2, 4)*spread(n, 1, 4)*detj*q%weights(k)
  end do

  print '(4f10.6)', m                        ! one column per line; M is symmetric
  print *, sum(m)                            ! 2.75: the area of the element
end program demo_mass_matrix
