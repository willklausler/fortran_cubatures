program cubatures_test
!! Test cubatures

  use cubatures

  implicit none

  integer :: i, j
  integer :: ord
  integer :: maxo

  real(rk) :: vol
  real(rk) :: coeff(10)
  real(rk) :: poly
  real(rk) :: anasol, numsol

  character(1) :: order

  character(*), parameter :: fmt1 = "(A15,': ',$)"

  type(cubature) :: scheme

  call random_number(coeff)

  write(*,"(A)") "Lineatures"
  write(*,fmt1) "Lines"
  j = maxval(maxloc(index(elmtypes, "LIN")))
  vol = volumes(j)
  maxo = maxorders(j)
  do i = 1,maxo
    write(order,"(I1)") i
    call scheme%set("LIN", [i])

    ! Check weight sum
    if (.not.is_zero(sum(scheme%weights) - vol)) then
      error stop "order "//order//": sum of weights failed"
    end if

    ! Check abscissae balance
    if (.not.is_zero(sum(scheme%abscissae(1,:)))) then
      error stop "order "//order//": sum of abscissae failed"
    end if

    ord = 2*i-1

    ! Analytical integral
    anasol = ipolynomial(coeff,ord)

    ! Numerical integral
    numsol = 0
    do j = 1,scheme%points
      poly = polynomial(coeff,scheme%abscissae(1,j),ord)
      numsol = numsol + poly*scheme%weights(j)
    end do ! j

    if (.not.is_zero(anasol-numsol)) then
      write(*,*)
      write(*,*) anasol, numsol
      error stop "order "//order//": integration failed"
    end if

  end do ! i
  write(*,"(A)") "passed"

  write(*,"(/,A)") "Quadratures"

  write(*,fmt1) "Quadrilaterals"
  j = maxval(maxloc(index(elmtypes, "QUA")))
  vol = volumes(j)
  maxo = maxorders(j)
  do i = 1,maxo
    write(order,"(I1)") i
    call scheme%set("QUA", [i])

    if (.not.is_zero(sum(scheme%weights) - vol)) then
      error stop "order "//order//": sum of weights failed"
    end if
    do j = 1,2
      if (.not.is_zero(sum(scheme%abscissae(j,:)))) then
        error stop "order "//order//": sum of abscissae failed"
      end if
    end do ! j

    ord = 2*i-1

    ! Analytical integral
    anasol = ipolynomial(coeff,ord)**2

    ! Numerical integral
    numsol = 0
    do j = 1,scheme%points
      poly = polynomial(coeff,scheme%abscissae(1,j),ord) &
            *polynomial(coeff,scheme%abscissae(2,j),ord)
      numsol = numsol + poly*scheme%weights(j)
    end do ! j

    if (.not.is_zero(anasol-numsol)) then
      write(*,*)
      write(*,*) anasol, numsol
      error stop "order "//order//": integration failed"
    end if

  end do ! i
  write(*,"(A)") "passed"

  write(*,fmt1) "Triangles"
  j = maxval(maxloc(index(elmtypes, "TRI")))
  vol = volumes(j)
  maxo = maxorders(j)
  do i = 1,maxo
    write(order,"(I1)") i
    call scheme%set("TRI", [i])

    if (.not.is_zero(sum(scheme%weights) - vol)) then
      write(*,*)
      write(*,*) sum(scheme%weights)
      error stop "order "//order//": sum of weights failed"
    end if
    do j = 1,2
      if (.not.is_zero(sum(scheme%abscissae(j,:))/scheme%points - 1.0_rk/3)) then
        write(*,*)
        write(*,*) sum(scheme%abscissae(j,:))
        error stop "order "//order//": sum of abscissae failed"
      end if
    end do ! j

    ! ord = 2*i-1

    ! ! Analytical integral
    ! anasol = ipolynomial(coeff,ord)**2

    ! ! Numerical integral
    ! numsol = 0
    ! do j = 1,scheme%points
    !   poly = polynomial(coeff,scheme%abscissae(1,j),ord) &
    !         *polynomial(coeff,scheme%abscissae(2,j),ord)
    !   numsol = numsol + poly*scheme%weights(j)
    ! end do ! j

    ! if (.not.is_zero(anasol-numsol)) then
    !   write(*,*)
    !   write(*,*) anasol, numsol
    !   error stop "order "//order//": integration failed"
    ! end if

  end do ! i
  write(*,"(A)") "passed"

  write(*,"(/,A)") "Cubatures"

  write(*,fmt1) "Hexahedrons"
  j = maxval(maxloc(index(elmtypes, "HEX")))
  vol = volumes(j)
  maxo = maxorders(j)
  do i = 1,maxo
    write(order,"(I1)") i
    call scheme%set("HEX", [i])
    if (.not.is_zero(sum(scheme%weights) - vol)) then
      error stop "order "//order//": sum of weights failed"
    end if
    do j = 1,3
      if (.not.is_zero(sum(scheme%abscissae(j,:)))) then
        error stop "order "//order//": sum of abscissae failed"
      end if
    end do ! j

    ord = 2*i-1

    ! Analytical integral
    anasol = ipolynomial(coeff,ord)**3

    ! Numerical integral
    numsol = 0
    do j = 1,scheme%points
      poly = polynomial(coeff,scheme%abscissae(1,j),ord) &
            *polynomial(coeff,scheme%abscissae(2,j),ord) &
            *polynomial(coeff,scheme%abscissae(3,j),ord)
      numsol = numsol + poly*scheme%weights(j)
    end do ! j

    if (.not.is_zero(anasol-numsol)) then
      write(*,*)
      write(*,*) anasol, numsol
      error stop "order "//order//": integration failed"
    end if

  end do ! i
  write(*,"(A)") "passed"

  write(*,fmt1) "Tetrahedrons"
  j = maxval(maxloc(index(elmtypes, "TET")))
  vol = volumes(j)
  maxo = maxorders(j)
  do i = 1,maxo
    write(order,"(I1)") i
    call scheme%set("TET", [i])

    if (.not.is_zero(sum(scheme%weights) - vol)) then
      write(*,*)
      write(*,*) sum(scheme%weights)
      error stop "order "//order//": sum of weights failed"
    end if
    do j = 1,3
      if (.not.is_zero(sum(scheme%abscissae(j,:))/scheme%points - 1.0_rk/4)) then
        write(*,*)
        write(*,*) sum(scheme%abscissae(j,:))
        error stop "order "//order//": sum of abscissae failed"
      end if
    end do ! j

    ! ord = 2*i-1

    ! ! Analytical integral
    ! anasol = ipolynomial(coeff,ord)**2

    ! ! Numerical integral
    ! numsol = 0
    ! do j = 1,scheme%points
    !   poly = polynomial(coeff,scheme%abscissae(1,j),ord) &
    !         *polynomial(coeff,scheme%abscissae(2,j),ord)
    !   numsol = numsol + poly*scheme%weights(j)
    ! end do ! j

    ! if (.not.is_zero(anasol-numsol)) then
    !   write(*,*)
    !   write(*,*) anasol, numsol
    !   error stop "order "//order//": integration failed"
    ! end if

  end do ! i
  write(*,"(A)") "passed"

contains

pure elemental function is_zero(r) result(z)
  real(rk), intent(in) :: r
  logical :: z
  z = abs(r) < 10.0_rk**(-12)
end function is_zero

!***********************************************************************

pure function polynomial(c,x,o) result(v)
!! Evaluate polynomial at x with coefficients c to order o
!! v = sum_i=1^o+1 c_i*x^(i-1)
  real(rk), intent(in) :: c(10), x
  integer, intent(in) :: o
  real(rk) :: v
  integer :: i
  v = c(1)
  do i = 2,o+1
    v = v + c(i)*x**(i-1)
  end do ! i
end function polynomial

!***********************************************************************

pure function ipolynomial(c,o) result(v)
!! Evaluate polynomial integral over domain [xlo, xhi] with coefficients
!! c to order o
!! v = sum_i=1^o+1 c_i/i*xhi^i - c_i*xlo^i
  real(rk), intent(in) :: c(10)
  integer, intent(in) :: o
  real(rk) :: v
  real(rk), parameter :: xlo = -1, xhi = 1
  integer :: i
  v = 0
  do i = 1,o+1
    v = v + c(i)*(xhi**i - xlo**i)/i
  end do ! i
end function ipolynomial

end program cubatures_test