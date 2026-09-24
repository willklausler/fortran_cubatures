program cubatures_example
!! Use [[cubatures]] the way a finite element code does:
!!
!! 1. inspect a rule,
!! 2. integrate over a physical triangle through an affine map,
!! 3. assemble a bilinear quadrilateral (Q1) mass matrix through an
!!    isoparametric map,
!! 4. keep one rule per element type in a table indexed by `CUB_*`.

  use cubatures

  implicit none

  call inspect_rule()
  call triangle_moment()
  call q1_mass_matrix()
  call rule_table()

contains

!***********************************************************************

subroutine inspect_rule()
!! Build a rule and print its points and weights

  type(cubature) :: q

  write(*,"(/,A)") "1. Degree-2 triangle rule"
  q = cubature(CUB_TRI, 2)
  call q%show()

end subroutine inspect_rule

!***********************************************************************

subroutine triangle_moment()
!! \(\int_T xy\,dA\) over a physical triangle \(T\), mapped affinely
!! from the reference triangle, compared to the closed form
!! \(\frac{|T|}{12}\left(\sum_k x_k y_k + \sum_k x_k \sum_k y_k\right)\)

  real(rk), parameter :: xv(2,3) = reshape([1.0_rk, 0.0_rk, &
                                            4.0_rk, 1.0_rk, &
                                            2.0_rk, 3.0_rk], [2, 3])  ! Vertices
  type(cubature) :: q
  real(rk) :: jac(2,2), detj, x(2), total, area, exact
  integer :: k

  write(*,"(/,A)") "2. Integral of x*y over a physical triangle"

  ! x = xv1 + J xi, J constant for an affine map; integrand is degree 2
  q = cubature(CUB_TRI, 2)
  jac(:,1) = xv(:,2) - xv(:,1)
  jac(:,2) = xv(:,3) - xv(:,1)
  detj = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)

  total = 0
  do k = 1, q%npoints
    x = xv(:,1) + matmul(jac, q%abscissae(:,k))
    total = total + x(1)*x(2)*detj*q%weights(k)
  end do

  area = detj/2
  exact = area/12*(sum(xv(1,:)*xv(2,:)) + sum(xv(1,:))*sum(xv(2,:)))
  write(*,"(A,ES23.15)") "  cubature: ", total
  write(*,"(A,ES23.15)") "  exact:    ", exact

end subroutine triangle_moment

!***********************************************************************

subroutine q1_mass_matrix()
!! Consistent mass matrix \(M_{ij} = \int_\Omega N_i N_j\,dA\) of a
!! distorted bilinear quadrilateral. With the isoparametric map,
!! \(N_i N_j \det J\) has degree 3 per direction, so a degree-3 rule
!! (2 x 2 Gauss points) is exact. The entries of \(M\) sum to the area.

  real(rk), parameter :: xn(2,4) = reshape([0.0_rk, 0.0_rk, &
                                            2.0_rk, 0.0_rk, &
                                            2.5_rk, 1.5_rk, &
                                            0.0_rk, 1.0_rk], [2, 4])  ! Nodes
  real(rk), parameter :: sn(2,4) = reshape([-1, -1,  1, -1,  1,  1, -1,  1], [2, 4])
    ! Reference node coordinates
  type(cubature) :: q
  real(rk) :: n(4), dn(4,2), jac(2,2), detj, m(4,4), area
  integer :: k, i

  write(*,"(/,A)") "3. Q1 mass matrix of a distorted quadrilateral"

  q = cubature(CUB_QUA, 3)
  m = 0
  do k = 1, q%npoints
    associate (xi => q%abscissae(:,k))
      n = (1 + sn(1,:)*xi(1))*(1 + sn(2,:)*xi(2))/4
      dn(:,1) = sn(1,:)*(1 + sn(2,:)*xi(2))/4
      dn(:,2) = sn(2,:)*(1 + sn(1,:)*xi(1))/4
    end associate
    jac = matmul(xn, dn)
    detj = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
    m = m + spread(n, 2, 4)*spread(n, 1, 4)*detj*q%weights(k)
  end do

  do i = 1, 4
    write(*,"(2X,4F10.6)") m(i,:)
  end do

  area = 0.5_rk*abs(sum(xn(1,:)*cshift(xn(2,:), 1) - cshift(xn(1,:), 1)*xn(2,:)))
  write(*,"(A,F10.6,A,F10.6)") "  sum(M) =", sum(m), ",  area =", area

end subroutine q1_mass_matrix

!***********************************************************************

subroutine rule_table()
!! One rule per element type, looked up by the element's `CUB_*` constant

  type(cubature) :: rules(6)
  integer, parameter :: mesh(*) = [CUB_TET, CUB_HEX, CUB_WED, CUB_TET]
    ! Element types of a small mixed mesh
  integer :: e

  write(*,"(/,A)") "4. Rule table for a mixed mesh, degree 4"

  rules(CUB_TET) = cubature(CUB_TET, 4)
  rules(CUB_HEX) = cubature(CUB_HEX, 4)
  rules(CUB_WED) = cubature(CUB_WED, 4)

  do e = 1, size(mesh)
    associate (q => rules(mesh(e)))
      write(*,"(2X,'element ',I0,': ',I3,' points, reference measure ',F8.6)") &
        e, q%npoints, sum(q%weights)
    end associate
  end do

end subroutine rule_table

!***********************************************************************

end program cubatures_example
