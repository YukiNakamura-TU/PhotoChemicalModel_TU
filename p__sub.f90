! Subroutines

! counting headers
subroutine p__count_header(nh, fname)
  integer nh
  character(len=256) fname
  character(len=256) strm
  !Skip header lines
  open(11, file = fname, status = 'unknown' )
    nh = 0
    do
      read(11,'(A)') strm
      if (strm(1:1) == "#") then
        nh = nh + 1
      else
        exit
      end if
    end do
  close(11)
end subroutine p__count_header

! solving max value
subroutine p__max(arrin, N, imax, max)
  integer N, i, imax
  real(8) arrin(N), tmp, max

  tmp = arrin(1)
  imax = 1
  do i = 2, N
    if (tmp < arrin(i)) then
      tmp = arrin(i)
      imax = i
    end if
  end do
  max = tmp

end subroutine p__max
