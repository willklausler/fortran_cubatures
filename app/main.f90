program degree_convergence
  !! Integrate exp(x + y + z) over every reference element with rules of
  !! increasing degree and print the number of points and the relative error,
  !! the table a finite element user needs to choose a rule for a
  !! non-polynomial integrand.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures, only: cubature, CUB_LIN, CUB_TRI, CUB_QUA, CUB_TET, CUB_HEX, CUB_WED
  implicit none

  integer, parameter :: elements(6) = [CUB_LIN, CUB_TRI, CUB_QUA, CUB_TET, CUB_HEX, CUB_WED]
  character(*), parameter :: names(6) = ['line ', 'tri  ', 'quad ', 'tet  ', 'hex  ', 'prism']
  integer, parameter :: pmax = 12
  real(dp), parameter :: e = exp(1._dp)
  ! Exact integrals of exp(x + y + z): (e - 1/e)**d on [-1,1]**d, 1 on the
  ! unit triangle, (e - 2)/2 on the unit tetrahedron, 1*(e - 1/e) on the prism
  real(dp), parameter :: exact(6) = [e - 1/e, 1._dp, (e - 1/e)**2, (e - 2)/2, &
                                     (e - 1/e)**3, e - 1/e]
  type(cubature) :: q
  integer :: i, p

  write(*, '(a)') 'Relative error (points) of the integral of exp(x + y + z)'
  write(*, '(a6, 6a16)') 'degree', names
  do p = 0, pmax, 2
    write(*, '(i6)', advance='no') p
    do i = 1, size(elements)
      q = cubature(elements(i), p)
      write(*, '(es11.2, " (", i3, ")")', advance='no') &
        abs(sum(q%weights*exp(sum(q%abscissae, dim=1))) - exact(i))/exact(i), q%npoints
    end do
    write(*, '(a)') ''
  end do

end program degree_convergence
