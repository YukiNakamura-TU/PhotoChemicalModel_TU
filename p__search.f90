module p__search
  use v__tdec, only : cst_, spl_, var_, grd_, flx_

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: p__search_reactant, p__search_product, sp_index

contains

  function p__search_reactant(spl, ich, species)
    type(spl_),       intent(in) :: spl
    integer,          intent(in) :: ich
    character(len=*), intent(in) :: species
    integer p__search_reactant
    integer i

    p__search_reactant = 0
    do i = 1, spl%reactant_list(ich,1)
      if( spl%species(spl%reactant_list(ich,i+1)) == species ) then
        p__search_reactant = 1
      end if
    end do

  end function p__search_reactant

  function p__search_product(spl, ich, species)
    type(spl_),       intent(in) :: spl
    integer,          intent(in) :: ich
    character(len=*), intent(in) :: species
    integer p__search_product
    integer i

    p__search_product = 0
    do i = 1, spl%product_list(ich,1)
      if( spl%species(spl%product_list(ich,i+1)) == species ) then
        p__search_product = 1
      end if
    end do

  end function p__search_product

  function sp_index(spl, species)
    type(spl_),       intent(in) :: spl
    character(len=*), intent(in) :: species
    integer sp_index
    integer isp

    sp_index = 0
    loop: do isp = 1, spl%nsp
      if( spl%species(isp) == species ) then
        sp_index = isp
        exit loop
      end if
    end do loop

  end function sp_index


end module p__search

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
