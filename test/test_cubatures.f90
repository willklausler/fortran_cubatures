! SPDX-License-Identifier: MIT
program test_cubatures
  !! Unit tests for module cubatures. For every element type and a sweep of
  !! degrees, each rule must integrate every monomial of its polynomial space
  !! exactly, have positive weights and interior abscissae, and be valid.
  !! Also checks point counts, point ordering, the constructor and `set`
  !! interfaces, reuse, `destroy`, and output.
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures
  implicit none

  character(3), parameter :: ename(6) = ['LIN', 'TRI', 'QUA', 'TET', 'HEX', 'WED']
  ! Relative tolerance on integrals; -O3 changes the last bits
  real(dp), parameter :: tol = 64*epsilon(1._dp)
  integer :: ncheck = 0, nfail = 0
  integer :: p, px, py, pz
  type(cubature) :: q, r, rules(6)
  character(200) :: line
  integer :: u, n, ios
  logical :: overflow

  ! Exactness, positivity, interior points --------------------------------
  do p = 0, 25
    call check_rule(cubature(CUB_LIN, p))
  end do
  do py = 0, 9
    do px = 0, 9
      call check_rule(cubature(CUB_QUA, [px, py]))
    end do
  end do
  do pz = 0, 5
    do py = 0, 5
      do px = 0, 5
        call check_rule(cubature(CUB_HEX, [px, py, pz]))
      end do
    end do
  end do
  do p = 0, 20
    call check_rule(cubature(CUB_TRI, p))
  end do
  do p = 0, 15
    call check_rule(cubature(CUB_TET, p))
  end do
  do pz = 0, 7
    do p = 0, 7
      call check_rule(cubature(CUB_WED, [p, pz]))
    end do
  end do

  ! Point counts ----------------------------------------------------------
  call check(npts(CUB_LIN, [0]) == 1, 'LIN 0 has 1 point')
  call check(npts(CUB_LIN, [5]) == 3, 'LIN 5 has 3 points')
  call check(npts(CUB_TRI, [1]) == 1, 'TRI 1 has 1 point')
  call check(npts(CUB_TRI, [2]) == 3, 'TRI 2 has 3 points')
  call check(npts(CUB_TRI, [3]) == 4, 'TRI 3 has 4 points')
  call check(npts(CUB_TRI, [4]) == 6, 'TRI 4 has 6 points')
  call check(npts(CUB_TRI, [5]) == 7, 'TRI 5 has 7 points')
  call check(npts(CUB_TET, [1]) == 1, 'TET 1 has 1 point')
  call check(npts(CUB_TET, [2]) == 4, 'TET 2 has 4 points')
  call check(npts(CUB_TET, [3]) == 8, 'TET 3 has 8 points')
  call check(npts(CUB_HEX, [1, 3, 5]) == 6, 'HEX [1,3,5] has 6 points')
  call check(npts(CUB_WED, [2, 3]) == 6, 'WED [2,3] has 6 points')

  ! Constructor, set, ordering, reuse, destroy ----------------------------
  call check(.not. q%is_valid(), 'default cubature is not valid')

  call q%set(CUB_HEX, 3)
  r = cubature(CUB_HEX, [3, 3, 3])
  call check(same(q, r), 'set(HEX, 3) = cubature(HEX, [3,3,3])')
  call check(all(q%degree == [3, 3, 3]), 'scalar degree broadcast to every direction')

  q = cubature(CUB_QUA, [3, 1])
  call check(q%npoints == 2, 'QUA [3,1] has 2 points')
  call check(q%abscissae(1,1) /= q%abscissae(1,2) .and. q%abscissae(2,1) == q%abscissae(2,2), &
    & 'first coordinate varies fastest')

  q = cubature(CUB_WED, [0, 3])
  call check(q%abscissae(3,1) /= q%abscissae(3,2), 'WED axial coordinate varies slowest')
  call check(all(q%degree == [0, 3, 0]), 'WED [0,3] degree = [0,3,0]')

  call q%set(CUB_HEX, 5)
  call q%set(CUB_LIN, 1)
  call check(q%is_valid() .and. q%dim == 1 .and. q%npoints == 1, 'set reuses q for a smaller rule')

  q = cubature(CUB_TRI, 2)
  call q%set(q%elm, 4)
  call check(q%elm == CUB_TRI .and. q%npoints == 6, 'set(q%elm, 4) raises the degree of q')

  call q%destroy()
  call check(.not. q%is_valid() .and. q%elm == 0 .and. q%dim == 0 .and. q%npoints == 0, &
    & 'destroy resets to the unset state')

  ! Rules indexed by element type, as in an FE element table
  rules(CUB_LIN) = cubature(CUB_LIN, 2)
  rules(CUB_TRI) = cubature(CUB_TRI, 2)
  rules(CUB_QUA) = cubature(CUB_QUA, 2)
  rules(CUB_TET) = cubature(CUB_TET, 2)
  rules(CUB_HEX) = cubature(CUB_HEX, 2)
  rules(CUB_WED) = cubature(CUB_WED, 2)
  call check(all([(rules(p)%elm == p, p = 1, 6)]), 'rules(elm)%elm = elm')

  ! Output ----------------------------------------------------------------
  open(newunit=u, status='scratch', action='readwrite')
  call q%summary(u)
  q = cubature(CUB_HEX, [1, 3, 5])
  call q%show(u)
  rewind(u)
  n = 0
  overflow = .false.
  do
    read(u, '(a)', iostat=ios) line
    if (ios /= 0) exit
    n = n + 1
    overflow = overflow .or. index(line, '*') > 0
  end do
  close(u)
  ! 1 unset + 4 summary + 6 points + 1 sum
  call check(n == 12, 'summary of unset rule and show write 12 lines')
  call check(.not. overflow, 'show has no format overflow')

  print '(i0, a, i0, a)', ncheck - nfail, ' of ', ncheck, ' checks passed'
  if (nfail > 0) error stop 1

contains

subroutine check(cond, name)
  logical, intent(in) :: cond
  character(*), intent(in) :: name
  ncheck = ncheck + 1
  if (.not. cond) then
    nfail = nfail + 1
    print '(a)', 'FAIL: ' // name
  end if
end subroutine check

integer function npts(elm, degree) result(n)
  !! Number of points of a rule
  integer, intent(in) :: elm, degree(:)
  type(cubature) :: q
  q = cubature(elm, degree)
  n = q%npoints
end function npts

subroutine check_rule(q)
  !! Check validity, positivity, interior points and exactness of a rule
  type(cubature), intent(in) :: q
  character(40) :: label
  integer :: a, b, c, e(3)
  real(dp) :: err
  real(dp), allocatable :: f(:)

  write(label, '(a, " [", i0, ",", i0, ",", i0, "]")') ename(q%elm), q%degree
  call check(q%is_valid(), trim(label) // ' is valid')
  call check(all(q%weights > 0), trim(label) // ' has positive weights')
  call check(inside(q), trim(label) // ' has interior abscissae')

  err = 0
  do c = 0, maxval(q%degree)
    do b = 0, maxval(q%degree)
      do a = 0, maxval(q%degree)
        if (.not. in_space(q%elm, q%degree, a, b, c)) cycle
        e = [a, b, c]
        f = product(q%abscissae**spread(e(1:q%dim), 2, q%npoints), dim=1)
        associate (exact => monomial_integral(q%elm, a, b, c))
          err = max(err, abs(sum(q%weights*f) - exact)/max(1._dp, abs(exact)))
        end associate
      end do
    end do
  end do
  call check(err <= tol, trim(label) // ' integrates its polynomial space exactly')
end subroutine check_rule

logical function in_space(elm, p, a, b, c) result(res)
  !! True if x**a * y**b * z**c is in the polynomial space of degree p
  integer, intent(in) :: elm, p(3), a, b, c
  select case (elm)
  case (CUB_LIN); res = a <= p(1) .and. b == 0 .and. c == 0
  case (CUB_QUA); res = a <= p(1) .and. b <= p(2) .and. c == 0
  case (CUB_HEX); res = a <= p(1) .and. b <= p(2) .and. c <= p(3)
  case (CUB_TRI); res = a + b <= p(1) .and. c == 0
  case (CUB_TET); res = a + b + c <= p(1)
  case (CUB_WED); res = a + b <= p(1) .and. c <= p(2)
  case default;   res = .false.
  end select
end function in_space

real(dp) function monomial_integral(elm, a, b, c) result(r)
  !! Exact integral of x**a * y**b * z**c over the reference element
  integer, intent(in) :: elm, a, b, c
  select case (elm)
  case (CUB_LIN); r = line_integral(a)
  case (CUB_QUA); r = line_integral(a)*line_integral(b)
  case (CUB_HEX); r = line_integral(a)*line_integral(b)*line_integral(c)
  case (CUB_TRI); r = simplex_integral([a, b])
  case (CUB_TET); r = simplex_integral([a, b, c])
  case (CUB_WED); r = simplex_integral([a, b])*line_integral(c)
  case default;   r = huge(r)
  end select
end function monomial_integral

real(dp) function line_integral(a) result(r)
  !! Integral of x**a over [-1, 1]
  integer, intent(in) :: a
  r = merge(2._dp/(a + 1), 0._dp, mod(a, 2) == 0)
end function line_integral

real(dp) function simplex_integral(e) result(r)
  !! Integral of prod(x_i**e_i) over the unit simplex, prod(e_i!)/(d + sum(e))!
  integer, intent(in) :: e(:)
  r = product(gamma(e + 1._dp))/gamma(size(e) + sum(e) + 1._dp)
end function simplex_integral

logical function inside(q) result(res)
  !! True if every abscissa is strictly inside the reference element
  type(cubature), intent(in) :: q
  associate (x => q%abscissae)
    select case (q%elm)
    case (CUB_LIN, CUB_QUA, CUB_HEX)
      res = all(abs(x) < 1)
    case (CUB_TRI, CUB_TET)
      res = all(x > 0) .and. all(sum(x, dim=1) < 1)
    case (CUB_WED)
      res = all(x(1:2,:) > 0) .and. all(sum(x(1:2,:), dim=1) < 1) .and. all(abs(x(3,:)) < 1)
    case default
      res = .false.
    end select
  end associate
end function inside

logical function same(q, r) result(res)
  !! True if two rules are identical
  type(cubature), intent(in) :: q, r
  res = q%elm == r%elm .and. q%dim == r%dim .and. all(q%degree == r%degree) &
    & .and. q%npoints == r%npoints
  if (res) res = all(q%abscissae == r%abscissae) .and. all(q%weights == r%weights)
end function same

end program test_cubatures
