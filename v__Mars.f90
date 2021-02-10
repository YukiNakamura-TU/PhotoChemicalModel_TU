!------------------------------------------------------------------------------
!
!                                       Mars
!
!------------------------------------------------------------------------------


module v__Mars
  use v__tdec,     only : set_, cst_, spl_, var_, grd_, flx_
  use p__search,   only : p__search_reactant, p__search_product, sp_index

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: v__Mars__ini, v__Mars__exe

contains

  !----------------------------------------------------------------------------------------------------------
  !
  !                                            Calculation settings
  !
  !----------------------------------------------------------------------------------------------------------
  subroutine v__Mars__ini(spl, & ! in
    &                     set  ) ! inout
    type(spl_),   intent(in)      :: spl
    type(set_),   intent(inout)   :: set

    real(dp), parameter :: Mday = 88775.0_dp

    character(len=256) fname

    if ( spl%planet == 'Mars' ) then

      !----------------------------------------------------------------
      !                    SETTING CALCULATION MODE
      !----------------------------------------------------------------
      ! mode selection :
      !      1 : 1D stable solution only (average over LT, latitude is fixed)
      !      2 : 2D stable solution only (average over LT)
      !
      !      3 : 2D rotation (longitude, latitude are fixed)
      !      4 : 3D rotation (longitude is fixed)
      !      5 : 3D global calculation # NOT YET #
      !
      !         mode 1, 3    => Set latitude below !
      !         mode 3, 4, 5 => Select CALCULATE or SKIP stable calculation below !
      !----------------------------------------------------------------
      !        SELECT 'CALCULATE' OR 'SKIP' STABLE CALCULATION
      !----------------------------------------------------------------
      ! calc_stable :
      !      1 : calclate stable solution
      !      0 : skip stable calclation
      !----------------------------------------------------------------
      !               SETTING INPUT STABLE DENSITY PATH
      !----------------------------------------------------------------
      ! if 'calc_stable =  0' (= skip stable calclation)
      !  -> please fill in the path of input density data
      !  -> grids of density data should be the same as the current grids
      !----------------------------------------------------------------
      !                   SETTING PLANETARY DAYS
      !----------------------------------------------------------------
      ! nday  : number of planetary days for rotational calculation
      !         (calculation days = nday - 0.5  [days]
      !----------------------------------------------------------------
      !                   SETTING CALCULATION SCHEME
      !----------------------------------------------------------------
      ! scheme = 1 : 'implicit'      : calculate chemical-flux Jacobian
      !                                to use implicit method
      !          2 : 'semi-implicit' : dt must be small
      !          3 : 'explicit'      : dt must be very small
      !

      ! calculation mode
      !set%mode        =  1
      set%calc_stable =  1
      set%fnamestable = './Mars/input/density/n_stable_Ls90.dat'
      set%test_loc    =  0
      set%sza         =  0.0_dp

      ! Stable calculation settings
      !set%nstep       = 30000 ! for stable solution (you should set large enough to calculate stable mode)
      !set%fin_sec     = 1.0e8_dp*3.0e7_dp ! [sec]: calculation stops when sum of dt reashes this value
      !set%dtime_limit = 1.0e14_dp ! [sec]: dt DO NOT exceed this value
      !set%latitude    = 0.0_dp ! latitude [deg]

      ! Solar longitude
      !set%Ls          = 270

      ! rotational calculation setting
      !set%nday        = 2 ! planetary days

      ! photochemical calculation scheme
      !set%scheme      = 1

      !----------------------------------------------------------------
      !      SETTING HORIZONTAL RESOLUTION FOR 2D, 3D CALCULATION
      !----------------------------------------------------------------
      ! latitude, local time grid number settings
      !  
      !      resx : resolution of local time, longitude [deg]
      !      resy : resolution of latitude [deg]
      !
      !      nx  : numer of local time grids : adapted automatically
      !      ny  : number of latitude grids  : adapted automatically
      !
      !  latitude, local time is calculated in 'p__prm.f90'
      !      local time :  00 - 24 LT (1: 00LT, nx: 24LT)
      !      latitude   : -90 ~ 90    (1: -90,  ny: +90)
      !--------------------------------
      !  nx, ny : MUST BE ODD NUMBER
      !--------------------------------

      set%resx = 3 ! [deg] localtime resolution
      set%resy = 3 ! [deg] latitude resolution

      !                                           END Calculation settings
      !----------------------------------------------------------------------------------------------------------

      !----------------------------------------------------------------
      !                  RESOLUTION ERROR MESSAGE
      !----------------------------------------------------------------
      if ( mod(360,set%resx) /= 0 .or. mod(180,set%resy) /= 0) then
        fname = './progress.dat'
        open(11, file = fname, status = 'unknown' )
          write(11,*) 'change resolution!'
          write(11,*)'!------------------------------------------'
          write(11,*)'!  360/resx, 180/resy : MUST BE INTEGER'
          write(11,*)'!------------------------------------------'
        close(11)
        stop
      else if ( mod(360/set%resx+1,2)==0 .or. mod(360/set%resy+1,2)==0 ) then
        fname = './progress.dat'
        open(11, file = fname, status = 'unknown' )
          write(11,*) 'change resolution!', 360/set%resx+1, 360/set%resy+1
          write(11,*)'!------------------------------------------'
          write(11,*)'!  nx, ny : MUST BE ODD NUMBER'
          write(11,*)'!------------------------------------------'
        close(11)
        stop
      end if

      ! number of latitude and local time grids are automatically adapted with the selected mode
      if ( set%mode == '1D' ) then
        set%nx     =  1
        set%ny     =  1
      else if ( set%mode == '2D Lat' ) then
        set%nx     =  1
        set%ny     =  180/set%resy + 1
      else if( set%mode == '2D Rot' ) then
        set%nx     =  360/set%resx + 1
        set%ny     =  1
      else if( set%mode == '3D Rot' .or. set%mode == '3D Global' ) then
        set%nx     =  360/set%resx + 1
        set%ny     =  180/set%resy + 1
      end if

    end if


  end subroutine v__Mars__ini



  !----------------------------------------------------------------------------------------------------------
  !
  !                                           Special treatment for Mars
  !
  !----------------------------------------------------------------------------------------------------------
  subroutine v__Mars__exe(spl, cst, grd, flx, & ! in
    &                     var                 ) ! inout
    type(spl_),   intent(in)     :: spl
    type(cst_),   intent(in)     :: cst
    type(grd_),   intent(in)     :: grd
    type(flx_),   intent(in)     :: flx
    type(var_),   intent(inout)  :: var
    integer i, j, ix, iy, iz, isp, ich, nspecial
    real(dp) tmp, tmp1, tmp2, tmpzarr1(grd%nz), tmpzarr2(grd%nz)
    character(len=256) fname

    if ( spl%planet == 'Mars' ) then 
      if (var%nspecial == 0) nspecial = 1
      nspecial = var%nspecial
      allocate(var%ich_special(nspecial), var%ki_special(nspecial,grd%nx,grd%ny,grd%nz))
      var%ki_special = 0.0_dp

      do ich = 1, spl%nch

        if ( spl%reaction_type_char(ich) == 'electron impact' ) then

          open(11, file = './Mars/OHtest/input/P_CO2+(B2Sigmau+)_2017_09_13.dat', status = 'unknown' )
            do iz = 1, grd%nz
              read(11,*)
              read(11,*) tmp, tmpzarr1(iz), tmpzarr2(iz)
            end do
          close(11)

          var%ich_special(1) = ich

          do iz = 1, grd%nz
          do iy = 1, grd%ny
          do ix = 1, grd%nx
            var%ki_special(1,ix,iy,iz) = ((tmpzarr1(iz) + tmpzarr2(iz))*1.0e6_dp )*1.0e0_dp
          end do
          end do
          !print *, ich, iz, var%ki_special(1,1,1,iz)
          end do
          !stop
            

        end if

      end do



    end if

  end subroutine v__Mars__exe


end module v__Mars