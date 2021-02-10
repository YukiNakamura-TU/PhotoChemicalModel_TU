module p__eddy_diffusion
  use v__tdec,        only : spl_, var_, grd_, cst_, set_
  use p__search,      only : p__search_reactant, p__search_product, sp_index
  
  implicit none
  integer(4), parameter :: sp = 4, dp = 8
  
  private
  public :: p__eddy_diffusion__exe
  
contains
  

  subroutine p__eddy_diffusion__exe(spl, cst, grd, & ! in
    &                               var            ) ! inout
    implicit none
    type(cst_),           intent(in)    :: cst
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(var_),           intent(inout) :: var
    integer iz
    real(dp) tmp

    !--------------------------------------------------------------------------------------
    !
    !                                       Mars
    !
    !--------------------------------------------------------------------------------------
    ! var%n_tot(iz) is total neutral density in the unit of [/m^3]
    if (spl%planet == 'Mars') then
      
      ! Chaffin et al., 2017 case
      do iz = 1, grd%nz
          var%K_eddy(iz) = 1.0d2 ![m^2/s]
          tmp = 2.0e9_dp / dsqrt(var%n_tot(iz)*1.0e-6_dp) ![m^2/s]
          if (tmp > 1.0e2_dp) var%K_eddy(iz) = tmp
      end do

    end if

    !--------------------------------------------------------------------------------------
    !
    !                                      Jupiter
    !
    !--------------------------------------------------------------------------------------
    if (spl%planet == 'Jupiter') then

      var%K_eddy = 2.0d2 ![m^2/s]
      
    end if



  end subroutine p__eddy_diffusion__exe


end module p__eddy_diffusion