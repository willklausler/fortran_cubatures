program cubatures_test
!! Unit tests for [[cubatures]].
!!
!! For every element type and a sweep of degrees, checks that each rule
!!
!! - integrates every monomial in its polynomial space exactly,
!! - has positive weights and abscissae strictly inside the element,
!! - is consistent (`is_valid`).
!!
!! Also checks point counts of the tabulated rules, point ordering, the
!! constructor and `set` interfaces, reuse, `destroy`, and output.

  use cubatures

  implicit none

  real(rk), parameter :: tol = 1.0e-13_rk   !! Relative tolerance on integrals
  integer :: nfail = 0                      !! Number of failed checks
  integer :: p, px, py, pz

  call section("Line")
  do p = 0, 25
    call check_rule(cubature(CUB_LIN, p))
  end do

  call section("Quadrilateral")
  do py = 0, 9
    do px = 0, 9
      call check_rule(cubature(CUB_QUA, [px, py]))
    end do
  end do

  call section("Hexahedron")
  do pz = 0, 5
    do py = 0, 5
      do px = 0, 5
        call check_rule(cubature(CUB_HEX, [px, py, pz]))
      end do
    end do
  end do

  call section("Triangle")
  do p = 0, 20
    call check_rule(cubature(CUB_TRI, p))
  end do

  call section("Tetrahedron")
  do p = 0, 15
    call check_rule(cubature(CUB_TET, p))
  end do

  call section("Wedge")
  do pz = 0, 7
    do p = 0, 7
      call check_rule(cubature(CUB_WED, [p, pz]))
    end do
  end do

  call section("Point counts")
  call expect(npts(CUB_LIN, [0]) == 1,       "LIN 0")
  call expect(npts(CUB_LIN, [5]) == 3,       "LIN 5")
  call expect(npts(CUB_TRI, [1]) == 1,       "TRI 1")
  call expect(npts(CUB_TRI, [2]) == 3,       "TRI 2")
  call expect(npts(CUB_TRI, [3]) == 4,       "TRI 3")
  call expect(npts(CUB_TRI, [4]) == 6,       "TRI 4")
  call expect(npts(CUB_TRI, [5]) == 7,       "TRI 5")
  call expect(npts(CUB_TET, [1]) == 1,       "TET 1")
  call expect(npts(CUB_TET, [2]) == 4,       "TET 2")
  call expect(npts(CUB_TET, [3]) == 8,       "TET 3")
  call expect(npts(CUB_HEX, [1, 3, 5]) == 6, "HEX [1,3,5]")
  call expect(npts(CUB_WED, [2, 3]) == 6,    "WED [2,3]")

  call section("Interface")
  call test_interface()

  call section("Output")
  call test_output()

  write(*,*)
  if (nfail > 0) then
    write(*,"(I0,A)") nfail, " check(s) failed"
    error stop 1
  end if
  write(*,"(A)") "All tests passed"

contains

!***********************************************************************

subroutine section(name)
!! Start a group of checks
  character(*), intent(in) :: name
  write(*,"(A)") name
end subroutine section

!***********************************************************************

subroutine expect(cond, msg)
!! Record a failed check
  logical, intent(in) :: cond
  character(*), intent(in) :: msg
  if (cond) return
  nfail = nfail + 1
  write(*,"(2X,'FAIL: ',A)") msg
end subroutine expect

!***********************************************************************

integer function npts(elm, degree)
!! Number of points of a rule
  integer, intent(in) :: elm, degree(:)
  type(cubature) :: q
  q = cubature(elm, degree)
  npts = q%npoints
end function npts

!***********************************************************************

subroutine check_rule(q)
!! Check exactness, positivity, interior points and consistency of a rule

  type(cubature), intent(in) :: q

  character(40) :: label
  integer :: a, b, c, e(3)
  real(rk) :: numerical, exact
  real(rk), allocatable :: f(:)

  write(label,"(I0,' [',I0,',',I0,',',I0,']')") q%elm, q%degree

  call expect(q%is_valid(), trim(label)//" is_valid")
  call expect(all(q%weights > 0), trim(label)//" positive weights")
  call expect(inside(q), trim(label)//" interior abscissae")

  do c = 0, maxval(q%degree)
    do b = 0, maxval(q%degree)
      do a = 0, maxval(q%degree)
        if (.not. in_space(q%elm, q%degree, a, b, c)) cycle
        e = [a, b, c]
        f = product(q%abscissae**spread(e(1:q%dim), 2, q%npoints), dim=1)
        numerical = sum(q%weights*f)
        exact = monomial_integral(q%elm, a, b, c)
        if (abs(numerical - exact) > tol*max(1.0_rk, abs(exact))) then
          call expect(.false., trim(label)//" exactness")
          write(*,"(4X,'x^',I0,' y^',I0,' z^',I0,': ',2ES24.16)") a, b, c, numerical, exact
          return
        end if
      end do
    end do
  end do

end subroutine check_rule

!***********************************************************************

logical function in_space(elm, p, a, b, c)
!! True if \(x^a y^b z^c\) is in the polynomial space of degree `p`

  integer, intent(in) :: elm, p(3), a, b, c

  select case (elm)
  case (CUB_LIN); in_space = a <= p(1) .and. b == 0 .and. c == 0
  case (CUB_QUA); in_space = a <= p(1) .and. b <= p(2) .and. c == 0
  case (CUB_HEX); in_space = a <= p(1) .and. b <= p(2) .and. c <= p(3)
  case (CUB_TRI); in_space = a + b <= p(1) .and. c == 0
  case (CUB_TET); in_space = a + b + c <= p(1)
  case (CUB_WED); in_space = a + b <= p(1) .and. c <= p(2)
  case default;   in_space = .false.
  end select

end function in_space

!***********************************************************************

real(rk) function monomial_integral(elm, a, b, c) result(r)
!! Exact integral of \(x^a y^b z^c\) over the reference element

  integer, intent(in) :: elm, a, b, c

  select case (elm)
  case (CUB_LIN); r = line(a)
  case (CUB_QUA); r = line(a)*line(b)
  case (CUB_HEX); r = line(a)*line(b)*line(c)
  case (CUB_TRI); r = simplex([a, b])
  case (CUB_TET); r = simplex([a, b, c])
  case (CUB_WED); r = simplex([a, b])*line(c)
  case default;   r = huge(r)
  end select

end function monomial_integral

!***********************************************************************

real(rk) function line(a)
!! \(\int_{-1}^{1} x^a\,dx\)
  integer, intent(in) :: a
  line = merge(2.0_rk/(a + 1), 0.0_rk, mod(a, 2) == 0)
end function line

!***********************************************************************

real(rk) function simplex(e)
!! \(\int x_1^{e_1} \cdots x_d^{e_d}\) over the unit simplex,
!! \(= \prod e_i! / (d + \sum e_i)!\)
  integer, intent(in) :: e(:)
  simplex = product(gamma(e + 1.0_rk))/gamma(size(e) + sum(e) + 1.0_rk)
end function simplex

!***********************************************************************

logical function inside(q)
!! True if every abscissa is strictly inside the reference element

  type(cubature), intent(in) :: q

  associate (x => q%abscissae)
    select case (q%elm)
    case (CUB_LIN, CUB_QUA, CUB_HEX)
      inside = all(abs(x) < 1)
    case (CUB_TRI, CUB_TET)
      inside = all(x > 0) .and. all(sum(x, dim=1) < 1)
    case (CUB_WED)
      inside = all(x(1:2,:) > 0) .and. all(sum(x(1:2,:), dim=1) < 1) &
         .and. all(abs(x(3,:)) < 1)
    case default
      inside = .false.
    end select
  end associate

end function inside

!***********************************************************************

subroutine test_interface()
!! Constructor, `set`, scalar vs array degree, ordering, reuse, destroy

  type(cubature) :: q, r
  type(cubature) :: rules(6)

  call expect(.not. q%is_valid(), "default is not valid")

  call q%set(CUB_HEX, 3)
  r = cubature(CUB_HEX, [3, 3, 3])
  call expect(same(q, r), "set(scalar) == cubature(array)")
  call expect(all(q%degree == [3, 3, 3]), "scalar degree broadcast")

  q = cubature(CUB_QUA, [3, 1])
  call expect(q%npoints == 2, "QUA [3,1] points")
  call expect(q%abscissae(1,1) /= q%abscissae(1,2) &
        .and. q%abscissae(2,1) == q%abscissae(2,2), "first coordinate fastest")

  q = cubature(CUB_WED, [0, 3])
  call expect(q%abscissae(3,1) /= q%abscissae(3,2), "WED axial slowest")
  call expect(all(q%degree == [0, 3, 0]), "WED degrees")

  call q%set(CUB_HEX, 5)
  call q%set(CUB_LIN, 1)
  call expect(q%is_valid() .and. q%dim == 1 .and. q%npoints == 1, "reuse with smaller rule")

  call q%destroy()
  call expect(.not. q%is_valid() .and. q%elm == 0 .and. q%dim == 0 &
        .and. q%npoints == 0, "destroy")

  ! Rules indexed by element type, as in an FE element table
  rules(CUB_LIN) = cubature(CUB_LIN, 2)
  rules(CUB_TRI) = cubature(CUB_TRI, 2)
  rules(CUB_QUA) = cubature(CUB_QUA, 2)
  rules(CUB_TET) = cubature(CUB_TET, 2)
  rules(CUB_HEX) = cubature(CUB_HEX, 2)
  rules(CUB_WED) = cubature(CUB_WED, 2)
  call expect(all([(rules(p)%elm == p, p = 1, 6)]), "rule table")

end subroutine test_interface

!***********************************************************************

logical function same(q, r)
!! True if two rules are identical
  type(cubature), intent(in) :: q, r
  same = q%elm == r%elm .and. q%dim == r%dim .and. all(q%degree == r%degree) &
   .and. q%npoints == r%npoints
  if (same) same = all(q%abscissae == r%abscissae) .and. all(q%weights == r%weights)
end function same

!***********************************************************************

subroutine test_output()
!! `summary` and `show` write complete, well-formatted output

  type(cubature) :: q
  character(200) :: line
  integer :: u, n, ios
  logical :: overflow

  open(newunit=u, status="scratch", action="readwrite")
  call q%summary(u)
  q = cubature(CUB_HEX, [1, 3, 5])
  call q%show(u)

  rewind(u)
  n = 0
  overflow = .false.
  do
    read(u,"(A)",iostat=ios) line
    if (ios /= 0) exit
    n = n + 1
    overflow = overflow .or. index(line, "*") > 0
  end do
  close(u)

  ! 1 unset + 4 summary + 6 points + 1 sum
  call expect(n == 12, "show line count")
  call expect(.not. overflow, "no format overflow")

end subroutine test_output

!***********************************************************************

end program cubatures_test
