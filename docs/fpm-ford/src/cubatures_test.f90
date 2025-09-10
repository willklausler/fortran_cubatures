program cubatures_test
!! Test cubatures

  use iso_fortran_env, only: rk => real64
  use cubatures

  implicit none

  integer :: g, h
  integer :: ind
  integer :: ord
  integer :: maxo

  real(rk) :: vol
  real(rk) :: coeff(10,10,10)
  real(rk) :: poly
  real(rk) :: anasol, numsol

  character(1) :: order

  character(*), parameter :: fmt1 = "(A15,': ',$)"

  type(cubature) :: scheme

  write(*,"(A)") "Lineatures"
  write(*,fmt1) "Lines"
  ind = maxval(maxloc(index(elmtypes, "LIN")))
  vol = volumes(ind)
  maxo = maxorders(ind)
  do g = 1,maxo
    write(order,"(I1)") g
    call scheme%set("LIN", [g])

    ! Check weight sum
    if (.not.is_zero(sum(scheme%weights) - vol)) then
      error stop "order "//order//": sum of weights failed"
    end if

    ! Check abscissae balance
    if (.not.is_zero(sum(scheme%abscissae(1,:)))) then
      error stop "order "//order//": sum of abscissae failed"
    end if

    ord = 2*g-1
    coeff = 0
    call random_number(coeff(:,1,1))

    ! Analytical integral
    anasol = ipolyline(coeff(:,1,1),ord)

    ! Numerical integral
    numsol = 0
    do h = 1,scheme%points
      poly = polyfunc(coeff, &
                      scheme%abscissae(1,h), &
                      0.0_rk, &
                      0.0_rk, &
                      ord)
      numsol = numsol + poly*scheme%weights(h)
    end do ! h

    if (.not.is_zero(anasol-numsol)) then
      write(*,*)
      write(*,*) anasol, numsol
      error stop "order "//order//": integration failed"
    end if

  end do ! g
  write(*,"(A)") "passed"

  write(*,"(/,A)") "Quadratures"

  write(*,fmt1) "Quadrilaterals"
  ind = maxval(maxloc(index(elmtypes, "QUA")))
  vol = volumes(ind)
  maxo = maxorders(ind)
  do g = 1,maxo
    write(order,"(I1)") g
    call scheme%set("QUA", [g])

    if (.not.is_zero(sum(scheme%weights) - vol)) then
      error stop "order "//order//": sum of weights failed"
    end if
    do h = 1,2
      if (.not.is_zero(sum(scheme%abscissae(h,:)))) then
        error stop "order "//order//": sum of abscissae failed"
      end if
    end do ! j

    ord = 2*g-1
    coeff = 0
    call random_number(coeff(:,:,1))

    ! Analytical integral
    anasol = ipolysquare(coeff(:,:,1),ord)

    ! Numerical integral
    numsol = 0
    do h = 1,scheme%points
      poly = polyfunc(coeff(:,:,1), &
                      scheme%abscissae(1,h), &
                      scheme%abscissae(2,h), &
                      0.0_rk, &
                      ord)
      numsol = numsol + poly*scheme%weights(h)
    end do ! j

    if (.not.is_zero(anasol-numsol)) then
      write(*,*)
      write(*,*) anasol, numsol
      error stop "order "//order//": integration failed"
    end if

  end do ! g
  write(*,"(A)") "passed"

  write(*,fmt1) "Triangles"
  ind = maxval(maxloc(index(elmtypes, "TRI")))
  vol = volumes(ind)
  maxo = maxorders(ind)
  do g = 1,maxo
    write(order,"(I1)") g
    call scheme%set("TRI", [g])

    if (.not.is_zero(sum(scheme%weights) - vol)) then
      write(*,*)
      write(*,*) sum(scheme%weights)
      error stop "order "//order//": sum of weights failed"
    end if
    do h = 1,2
      if (.not.is_zero(sum(scheme%abscissae(h,:))/scheme%points - 1.0_rk/3)) then
        write(*,*)
        write(*,*) sum(scheme%abscissae(h,:))
        error stop "order "//order//": sum of abscissae failed"
      end if
    end do ! h

    ord = g
    coeff = 0
    call random_number(coeff(:,:,1))

    ! Analytical integral
    anasol = ipolytriangle(coeff(:,:,1),ord)

    ! Numerical integral
    numsol = 0
    do h = 1,scheme%points
      poly = polyfunc(coeff(:,:,1), &
                      scheme%abscissae(1,h), &
                      scheme%abscissae(2,h), &
                      0.0_rk, &
                      ord)
      numsol = numsol + poly*scheme%weights(h)
    end do ! h

    if (.not.is_zero(anasol-numsol)) then
      write(*,*)
      write(*,*) anasol, numsol
      error stop "order "//order//": integration failed"
    end if

  end do ! i
  write(*,"(A)") "passed"

  write(*,"(/,A)") "Cubatures"

  write(*,fmt1) "Hexahedrons"
  ind = maxval(maxloc(index(elmtypes, "HEX")))
  vol = volumes(ind)
  maxo = maxorders(ind)
  do g = 1,maxo
    write(order,"(I1)") g
    call scheme%set("HEX", [g])
    if (.not.is_zero(sum(scheme%weights) - vol)) then
      error stop "order "//order//": sum of weights failed"
    end if
    do h = 1,3
      if (.not.is_zero(sum(scheme%abscissae(h,:)))) then
        error stop "order "//order//": sum of abscissae failed"
      end if
    end do ! h

    ord = 2*g-1
    call random_number(coeff)

    ! Analytical integral
    anasol = ipolycube(coeff,ord)

    ! Numerical integral
    numsol = 0
    do h = 1,scheme%points
      poly = polyfunc(coeff, &
                      scheme%abscissae(1,h), &
                      scheme%abscissae(2,h), &
                      scheme%abscissae(3,h), &
                      ord)
      numsol = numsol + poly*scheme%weights(h)
    end do ! h

    if (.not.is_zero(anasol-numsol)) then
      write(*,*)
      write(*,*) anasol, numsol
      error stop "order "//order//": integration failed"
    end if

  end do ! g
  write(*,"(A)") "passed"

  write(*,fmt1) "Tetrahedrons"
  ind = maxval(maxloc(index(elmtypes, "TET")))
  vol = volumes(ind)
  maxo = maxorders(ind)
  do g = 1,maxo
    write(order,"(I1)") g
    call scheme%set("TET", [g])

    if (.not.is_zero(sum(scheme%weights) - vol)) then
      write(*,*)
      write(*,*) sum(scheme%weights)
      error stop "order "//order//": sum of weights failed"
    end if
    do h = 1,3
      if (.not.is_zero(sum(scheme%abscissae(h,:))/scheme%points - 1.0_rk/4)) then
        write(*,*)
        write(*,*) sum(scheme%abscissae(h,:))
        error stop "order "//order//": sum of abscissae failed"
      end if
    end do ! h

    ord = g
    call random_number(coeff)

    ! Analytical integral
    anasol = ipolytet(coeff,ord)

    ! Numerical integral
    numsol = 0
    do h = 1,scheme%points
      poly = polyfunc(coeff(:,:,1), &
                      scheme%abscissae(1,h), &
                      scheme%abscissae(2,h), &
                      scheme%abscissae(3,h), &
                      ord)
      numsol = numsol + poly*scheme%weights(h)
    end do ! h

    if (.not.is_zero(anasol-numsol)) then
      write(*,*)
      write(*,*) anasol, numsol
      error stop "order "//order//": integration failed"
    end if

  end do ! g
  write(*,"(A)") "passed"

  write(*,fmt1) "Prisms"
  ind = maxval(maxloc(index(elmtypes, "WEJ")))
  vol = volumes(ind)
  maxo = maxorders(ind)
  do g = 1,maxo
    write(order,"(I1)") g
    call scheme%set("WEJ", [g])

    if (.not.is_zero(sum(scheme%weights) - vol)) then
      write(*,*)
      write(*,*) sum(scheme%weights)
      error stop "order "//order//": sum of weights failed"
    end if
    do h = 1,2
      if (.not.is_zero(sum(scheme%abscissae(h,:))/scheme%points - 1.0_rk/3)) then
        write(*,*)
        write(*,*) sum(scheme%abscissae(h,:))
        error stop "order "//order//": sum of abscissae failed"
      end if
    end do ! j
    if (.not.is_zero(sum(scheme%abscissae(3,:)))) then
      write(*,*)
      write(*,*) sum(scheme%abscissae(3,:))
      error stop "order "//order//": sum of abscissae failed"
    end if

    ord = g
    call random_number(coeff)

    ! Analytical integral
    anasol = ipolywedge(coeff,ord)

    ! Numerical integral
    numsol = 0
    do h = 1,scheme%points
      poly = polyfunc(coeff(:,:,1), &
                      scheme%abscissae(1,h), &
                      scheme%abscissae(2,h), &
                      scheme%abscissae(3,h), &
                      ord)
      numsol = numsol + poly*scheme%weights(h)
    end do ! h

    if (.not.is_zero(anasol-numsol)) then
      write(*,*)
      write(*,*) anasol, numsol
      error stop "order "//order//": integration failed"
    end if

  end do ! g
  write(*,"(A)") "passed"

contains

!***********************************************************************

pure elemental function is_zero(r) result(z)
  real(rk), intent(in) :: r
  logical :: z
  z = abs(r) < 10.0_rk**(-12)
end function is_zero

!***********************************************************************

pure function polyfunc(c,x,y,z,o) result(v)
!! Evaluate polynomial at x(3) with coefficients c to order o
!! v = sum_i,j,k=0^{i+j+k=o} c_ijk*x^i*y^j*z^k
  real(rk), intent(in) :: c(10,10,10)
  real(rk), intent(in) :: x, y, z
  integer, intent(in) :: o
  real(rk) :: v

  integer :: i, j, k

  v = 0
  do i = 0,o
    do j = 0,o-i
      do k = 0,o-i-j
        v = v + c(i+1,j+1,k+1)*(x**i)*(y**j)*(z**k)
      end do ! k
    end do ! j
  end do ! i

end function polyfunc

!***********************************************************************

pure function ipolyline(c,o) result(v)
!! Evaluate polynomial integral over domain [xlo, xhi] with coefficients
!! c to order o
!! v = sum_i=0^o c_i*(xhi^i - xlo^i)/i
  real(rk), intent(in) :: c(10)
  integer, intent(in) :: o
  real(rk) :: v
  integer :: i

  v = 0
  do i = 0,o
    v = v + c(i+1)*polyline(i)
  end do ! i
end function ipolyline

!***********************************************************************

pure function ipolysquare(c,o) result(v)
!! Evaluate polynomial integral over square domain
  real(rk), intent(in) :: c(10,10)
  integer, intent(in) :: o
  real(rk) :: v

  integer :: i, j

  v = 0
  do i = 0,o
    do j = 0,o-i
      v = v + c(i+1,j+1)*polyline(i)*polyline(j)
    end do ! j
  end do ! i

end function ipolysquare

!***********************************************************************

pure function ipolytriangle(c,o) result(v)
!! Evaluate polynomial integral over triangular domain
  real(rk), intent(in) :: c(10,10)
  integer, intent(in) :: o
  real(rk) :: v

  integer :: i, j

  v = 0
  do i = 0,o
    do j = 0,o-i
      v = v + c(i+1,j+1)*polytri(i, j)
    end do ! j
  end do ! i

end function ipolytriangle

!***********************************************************************

pure function ipolycube(c,o) result(v)
!! Evaluate polynomial integral over cube domain
  real(rk), intent(in) :: c(10,10,10)
  integer, intent(in) :: o
  real(rk) :: v

  integer :: i, j, k

  v = 0
  do i = 0,o
    do j = 0,o-i
      do k = 0,o-i-j
        v = v + c(i+1,j+1,k+1)*polyline(i)*polyline(j)*polyline(k)
      end do ! k
    end do ! j
  end do ! i

end function ipolycube

!***********************************************************************

pure function ipolytet(c,o) result(v)
!! Evaluate polynomial integral over tetrahedral domain
  real(rk), intent(in) :: c(10,10,10)
  integer, intent(in) :: o
  real(rk) :: v

  integer :: i, j, k

  v = 0
  do i = 0,o
    do j = 0,o-i
      do k = 0,o-i-j
        v = v + c(i+1,j+1,k+1)*polytet(i, j, k)
      end do ! k
    end do ! j
  end do ! i

end function ipolytet

!***********************************************************************

pure function ipolywedge(c,o) result(v)
!! Evaluate polynomial integral over prismatic domain
  real(rk), intent(in) :: c(10,10,10)
  integer, intent(in) :: o
  real(rk) :: v

  integer :: i, j, k

  v = 0
  do i = 0,o
    do j = 0,o-i
      do k = 0,o-i-j
        v = v + c(i+1,j+1,k+1)*polytri(i,j)*polyline(k)
      end do ! k
    end do ! j
  end do ! i

end function ipolywedge

!***********************************************************************

pure elemental integer function factorial(n) result(f)
!! Basic factorial function
  integer, intent(in) :: n
  integer :: i
  f = 1
  do i = 1,n
    f = f*i
  end do ! i
end function factorial

!***********************************************************************

pure elemental real(rk) function polytri(m,n) result(f)
!! Value of polynomial integral over triangular domain
!! int_0^1 int_0^(1-x) x^m*y^n dy dx
  integer, intent(in) :: m, n
  f = 1.0_rk*factorial(m)*factorial(n)/factorial(m+n+2)
end function polytri

!***********************************************************************

pure elemental real(rk) function polytet(m,n,p) result(f)
!! Value of polynomial integral over tetrahedral domain
!! int_0^1 int_0^(1-x) int_0^(1-x-y) x^m*y^n*z^p dz dy dx
  integer, intent(in) :: m, n, p
  f = 1.0_rk*factorial(m)*factorial(n)*factorial(p)/factorial(m+n+p+3)
end function polytet

!***********************************************************************

pure elemental real(rk) function polyline(m) result(f)
!! Value of polynomial integral over linear domain
!! int_-1^1 x^m dx
  integer, intent(in) :: m

  if (mod(m,2) == 0) then
    f = 2.0_rk/(m+1)
  else
    f = 0
  end if

end function polyline

!***********************************************************************

end program cubatures_test