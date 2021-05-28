! Subroutines

subroutine p__smoothing_linear(nsmooth, narr, arr)

  implicit none

  integer(4), parameter  :: sp = 4, dp = 8
  integer nsmooth, narr, i, j
  real(dp) arr(narr)
  real(dp) tmparr1(-nsmooth:narr+nsmooth)
  real(dp) tmparr2(-nsmooth:narr+nsmooth)

  ! smoothing
  do i = 1, narr
    tmparr1(i) = arr(i)
  end do

  do i = 0, nsmooth
    tmparr1(-i) = 2_dp*tmparr1(1-i) - tmparr1(2-i)
  end do

  do i = 0, nsmooth
    tmparr1(narr+i+1) = 2_dp*tmparr1(narr+i) - tmparr1(narr+i-1)
  end do

  do i = 1, narr
    tmparr2(i) = 0.0_dp
    do j = -nsmooth, nsmooth
      tmparr2(i) = tmparr2(i) + tmparr1(i+j)
    end do
    arr(i) = tmparr2(i)/dble(2*nsmooth+1)
  end do

end subroutine p__smoothing_linear

subroutine p__smoothing_log(nsmooth, narr, arr)

  implicit none

  integer(4), parameter  :: sp = 4, dp = 8
  integer nsmooth, narr, i, j
  real(dp) arr(narr)
  real(dp) tmp
  real(dp) tmparr1(-nsmooth:narr+nsmooth)
  real(dp) tmparr2(-nsmooth:narr+nsmooth)

  ! smoothing
  do i = 1, narr
    tmparr1(i) = arr(i)
  end do

  do i = 0, nsmooth
    tmparr1(-i) = tmparr1(1-i)**2_dp / tmparr1(2-i)
  end do

  do i = 0, nsmooth
    tmparr1(narr+i+1) = tmparr1(narr+i)**2_dp / tmparr1(narr+i-1)
  end do

  do i = 1, narr
    tmparr2(i) = 0.0_dp
    do j = -nsmooth, nsmooth
      tmparr2(i) = tmparr2(i) + dlog(tmparr1(i+j))
    end do
    arr(i) = dexp(tmparr2(i)/dble(2*nsmooth+1))
  end do

end subroutine p__smoothing_log

