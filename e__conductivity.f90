!==============================================================================================================
!
!                                         Ionospheric conductivity model
!
!==============================================================================================================

include "v__tdec.f90"
include "c__prm.f90"
include "p__io.f90"
include "p__search.f90"
include "v__Jupiter.f90"
include "v__Mars.f90"
include "p__EUVAC.f90"
include "p__UV.f90"
include "p__eddy_diffusion.f90"
include "p__molecular_diffusion.f90"
include "p__photochem_opticaldepth.f90"
include "p__photochem_rate.f90"
include "p__photochem_transport.f90"
include "p__photochem_scheme.f90"
include "p__airglow.f90" ! not yet

include "./Jupiter/metal_Hill/v__in.f90"

program e__conductivity

  use v__tdec,                   only : set_, grd_, var_, cst_, xct_, spl_, flx_
  use c__prm,                    only : c__prm__ini, c__prm__planet
  use p__search,                 only : p__search_reactant, p__search_product, sp_index
  use v__Jupiter,                only : v__Jupiter__ini, v__Jupiter__exe
  use v__Mars,                   only : v__Mars__ini, v__Mars__exe
  use v__in,                     only : v__in__ini, v__in__exe
  use p__EUVAC,                  only : p__EUVAC_flux, p__EUVAC_cross_section
  use p__UV,                     only : p__UV_flux, p__UV_cross_section_dat, p__UV_cross_section_exe
  use p__photochem_opticaldepth, only : p__photochem_opticaldepth__exe
  use p__photochem_rate,         only : p__photochem_rate__exe
  use p__photochem_transport,    only : p__photochem_transport__exe
  use p__photochem_scheme,       only : p__photochem_scheme__exe
  use p__airglow,                only : p__airglow__exe
  use p__io,                     only : p__io_stable__fin, p__io_rotation__fin, p__io_deallocate, &
    &                                   p__io_progress

  implicit none
  integer(4), parameter  :: sp = 4, dp = 8
  type(set_) :: set ! type of calculation settings
  type(grd_) :: grd ! type of grid
  type(var_) :: var ! type of variables
  type(cst_) :: cst ! type of physical constant
  type(xct_) :: xct ! type of cross section
  type(spl_) :: spl ! type of species list and planet info
  type(flx_) :: flx ! type of solar flux
  integer i, j, k, ix, xs, iy, iz, ich, isp, jsp, is, iday, s, N
  real(dp) tmp
  character(len=256) fname, num

  real(dp) a0, aH2, nu_e_H2, nu_i_H2, omega_e, omega_i, n_H2, m_H2, n_e, n_i, m_i, m_e
  real(dp) tmp1, tmp2, tmp3
  real(dp), allocatable :: B_abs(:,:,:), sinI(:,:,:), cosI(:,:,:)
  real(dp), allocatable :: sigma_H(:,:,:), sigma_P(:,:,:), sigma_0(:,:,:)
  real(dp), allocatable :: hi_sigma_P(:,:), hi_sigma_H(:,:), hi_sigma_0(:,:)
  real(dp), allocatable :: sigma_tt(:,:,:), sigma_tp(:,:,:), sigma_pp(:,:,:)
  real(dp), allocatable :: hi_sigma_tt(:,:), hi_sigma_tp(:,:), hi_sigma_pp(:,:)
  real(dp), allocatable :: sigma_prime_tp(:,:,:), sigma_prime_tz(:,:,:), sigma_prime_zp(:,:,:), sigma_prime_zz(:,:,:)
  real(dp), allocatable :: tmparr(:,:)
  real(dp) B_factor, sigma_min
  

  real(dp) t1, t2

  call cpu_time(t1)

  !-----------------------------------------------------
  ! Magnetic field factor
  !-----------------------------------------------------
  ! i = 0  -> 0.1
  ! i = 2N -> 10
  i = 1

  N = 10
  B_factor = 10.0_dp**(dble(i-N)/dble(N))
  print *, 'B_factor = ', B_factor


  !----------------------------------------------------------------------------------------------------------
  !
  !                                     Initialization : automatically adapted
  !
  !----------------------------------------------------------------------------------------------------------

  !-----------------------------------------------------
  ! read physical constant
  !-----------------------------------------------------
  call c__prm__ini(cst) ! out

  !-----------------------------------------------------
  ! Planet selection: automatically
  !-----------------------------------------------------
  call v__in__ini(spl, set) ! out

  !-----------------------------------------------------
  ! Calculation settings for the selected planet: auto
  !-----------------------------------------------------
  if (      spl%planet == 'Mars' ) then
    call v__Mars__ini(spl, set)
  else if ( spl%planet == 'Jupiter' ) then
    call v__Jupiter__ini(spl, set)
  end if

  !-----------------------------------------------------
  ! initialize variables
  !-----------------------------------------------------
  call v__in__exe(cst,      & ! in
    &             set, spl, & ! inout
    &             var, grd  ) ! out

  set%start_rot = 0
  grd%latitude = set%latitude

  !-----------------------------------------------------
  ! apply physical constant and geometry to the chosen planet
  !-----------------------------------------------------
  call c__prm__planet(spl, set,    & ! in
    &                 cst, grd, flx) ! inout

  !-----------------------------------------------------
  ! Special treatment for selected planet
  !-----------------------------------------------------
  if ( spl%planet == 'Mars' ) then
    call v__Mars__exe(spl, cst, grd, flx, & ! in
      &               var                 ) ! inout
  else if ( spl%planet == 'Jupiter' ) then
    call v__Jupiter__exe(spl, cst, grd, flx, set, & ! in
      &                  var                      ) ! inout
  end if

  !-----------------------------------------------------
  ! EUVAC solar EUV flux Model
  !-----------------------------------------------------
  call p__EUVAC_flux(spl, grd,     & ! in
    &                var, xct, flx ) ! inout

  !-----------------------------------------------------
  ! Woods et al., solar UV flux reference
  !-----------------------------------------------------
  call p__UV_flux(spl, grd, cst,& ! in
    &             var, xct, flx ) ! inout

  !-----------------------------------------------------
  ! UV cross section data 
  !-----------------------------------------------------
  call p__UV_cross_section_dat(xct ) ! inout



  !----------------------------------------------------------------------------------------------------------
  !
  !                                     Ionospheric conductivity calculation
  !
  !----------------------------------------------------------------------------------------------------------

  allocate(B_abs(grd%nx,grd%ny,grd%nz),sinI(grd%nx,grd%ny,grd%nz),cosI(grd%nx,grd%ny,grd%nz))
  allocate(sigma_H(grd%nx,grd%ny,grd%nz))
  allocate(sigma_P(grd%nx,grd%ny,grd%nz))
  allocate(sigma_0(grd%nx,grd%ny,grd%nz))
  allocate(hi_sigma_H(grd%nx,grd%ny))
  allocate(hi_sigma_P(grd%nx,grd%ny))
  allocate(hi_sigma_0(grd%nx,grd%ny))
  allocate(sigma_tt(grd%nx,grd%ny,grd%nz))
  allocate(sigma_tp(grd%nx,grd%ny,grd%nz))
  allocate(sigma_pp(grd%nx,grd%ny,grd%nz))
  allocate(hi_sigma_tt(grd%nx,grd%ny))
  allocate(hi_sigma_tp(grd%nx,grd%ny))
  allocate(hi_sigma_pp(grd%nx,grd%ny))
  allocate(tmparr(grd%nx,grd%ny))
  allocate(sigma_prime_tp(grd%nx,grd%ny,grd%nz))
  allocate(sigma_prime_tz(grd%nx,grd%ny,grd%nz))
  allocate(sigma_prime_zp(grd%nx,grd%ny,grd%nz))
  allocate(sigma_prime_zz(grd%nx,grd%ny,grd%nz))

  print *, 'reading density data...'
  do ix = 1, grd%nx
    write(num,*) ix
    do isp = 1, spl%nsp
      fname = trim(ADJUSTL(set%dir_name))//'/output/density/3Drot/'//trim(ADJUSTL(num))//'/'&
        &   //trim(ADJUSTL(spl%species(isp)))//'.dat'
      open(11, file = fname, status = 'unknown' )
        do iz = 1, grd%nz
          read(11, *) tmp, (var%ni_3d(isp,ix,iy,iz), iy = 1, grd%ny)
        end do
      close(11)
    end do
  end do
  print *, 'finished reading density data!'


  !-----------------------------------------------------
  ! Magnetic field data 
  !-----------------------------------------------------
  if (spl%planet == 'Jupiter') then 

    print *, 'reading B data...'
    open(11, file = './Jupiter/dipole_Jupiter.dat', status = 'unknown')
      do iz = 1, grd%nz
      do iy = 1, grd%ny
      do ix = 1, grd%nx
          !read(11,*) B_abs(ix,iy,iz), sinI(ix,iy,iz), cosI(ix,iy,iz)
          read(11,*) B_abs(ix,iy,iz), sinI(ix,iy,iz), cosI(ix,iy,iz)
          B_abs(ix,iy,iz) = B_abs(ix,iy,iz) * B_factor
      end do
      end do
      end do
    close(11)
    print *, 'finished reading B data!'

  end if


  ! ionospheric conductivity ---------------------------------------------

  a0 = 5.2917d-11 !Bohr radius [m]
  aH2 = 2.7d-10

  sigma_H = 0.0_dp ! Hall conductivity
  sigma_P = 0.0_dp ! Pedersen conductivity
  sigma_0 = 0.0_dp ! parallel conductivity

  print *, 'calculating ionospheric conductivity...'
  do iz = 1, grd%nz
  do iy = 1, grd%ny
  do ix = 1, grd%nx

    jsp  = sp_index(spl, 'H2')
    n_H2 = var%ni_3d(jsp,ix,iy,iz)
    m_H2 = var%m(jsp)

    jsp  = sp_index(spl, 'e-')
    n_e  = var%ni_3d(jsp,ix,iy,iz)
    m_e  = var%m(jsp)

    do isp = 1, spl%nsp

      n_i  = var%ni_3d(isp,ix,iy,iz)
      m_i  = var%m(isp)

      ! electron
      if (spl%species(isp) == 'e-') then 
        !electron and H2
        tmp1 = var%Tn(iz) * cst%eV
        tmp2 = dsqrt( tmp1 )
        tmp3 = tmp1**2.0d0

        ! collision frequency 
        nu_e_H2 = n_H2 * dsqrt( 2.0d0 * cst%k_B * 11600.0d0 / m_e ) &
          &                   * ( - 0.0189d0 * tmp3**2.0d0 * tmp2 &
          &                       + 0.5354d0 * tmp3 * tmp1 * tmp2 &
          &                       - 5.5138d0 * tmp3 * tmp2 &
          &                       + 21.133d0 * tmp1 * tmp2 &
          &                       + 29.475d0 * tmp2 ) * a0**2.0d0
        ! gyro frequency
        omega_e = cst%q_e * B_abs(ix,iy,iz) / m_i

        ! Pedersen conductivity : electron term
        sigma_P(ix,iy,iz) = sigma_P(ix,iy,iz) &
          &               + cst%q_e * n_e / B_abs(ix,iy,iz) &
          &               * omega_e * nu_e_H2 &
          &               / ( nu_e_H2**2.0d0 + omega_e**2.0d0 )

        ! Hall conductivity : electron term
        sigma_H(ix,iy,iz) = sigma_H(ix,iy,iz) &
          &               + cst%q_e * n_e / B_abs(ix,iy,iz) &
          &               * omega_e**2.0d0 &
          &               / ( nu_e_H2**2.0d0 + omega_e**2.0d0 )

        ! Parallel conductivity : electron term
        sigma_0(ix,iy,iz) = sigma_0(ix,iy,iz) &
          &               + cst%q_e * n_e / B_abs(ix,iy,iz) &
          &               * omega_e / nu_e_H2


      ! ion
      else if (spl%species(isp) /= 'e-' .and. var%q(isp) > 0.0_dp) then 

        tmp1 = 2.0d0 * n_H2 * aH2**2.0d0 &
          &          * dsqrt( 2.0d0 * cst%pi * cst%k_B * var%Tn(iz) / m_H2 )

        ! collision frequency
        nu_i_H2 = tmp1 * dsqrt( 1.0d0 + m_H2 / m_i )
        ! gyro frequency
        omega_i = cst%q_e * B_abs(ix,iy,iz) / m_i

        ! Pedersen conductivity : ion term
        sigma_P(ix,iy,iz) = sigma_P(ix,iy,iz) &
          &               + cst%q_e * n_i / B_abs(ix,iy,iz) &
          &               * omega_i * nu_i_H2  &
          &               / ( nu_i_H2**2.0d0 + omega_i**2.0d0 )

        ! Hall conductivity : ion term
        sigma_H(ix,iy,iz) = sigma_H(ix,iy,iz) &
          &               - cst%q_e * n_i / B_abs(ix,iy,iz) &
          &               * omega_i**2.0d0  &
          &               / ( nu_i_H2**2.0d0 + omega_i**2.0d0 )

        ! Parallel conductivity : ion term
        sigma_0(ix,iy,iz) = sigma_0(ix,iy,iz) &
          &               + cst%q_e * n_i / B_abs(ix,iy,iz) &
          &               * omega_i / nu_i_H2 

      end if

    end do

  end do
  end do
  end do

  ! height integrated ionospheric conductivity ---------------------------------------------
  hi_sigma_H = 0.0_dp
  hi_sigma_P = 0.0_dp
  hi_sigma_0 = 0.0_dp
  do iy = 1, grd%ny
  do ix = 1, grd%nx
    do iz = 1, grd%nz
      hi_sigma_H(ix,iy) = hi_sigma_H(ix,iy) + sigma_H(ix,iy,iz)*grd%dalt(iz)
      hi_sigma_P(ix,iy) = hi_sigma_P(ix,iy) + sigma_P(ix,iy,iz)*grd%dalt(iz)
      hi_sigma_0(ix,iy) = hi_sigma_0(ix,iy) + sigma_0(ix,iy,iz)*grd%dalt(iz)
    end do
  end do 
  end do

  ! transforming coordinate --------------------------------------------- 
  do iz = 1, grd%nz
  do iy = 1, grd%ny
  do ix = 1, grd%nx

    tmp1 = sigma_P(ix,iy,iz) * cosI(ix,iy,iz)**2.0d0 &
      &  + sigma_0(ix,iy,iz) * sinI(ix,iy,iz)**2.0d0

    sigma_tt(ix,iy,iz) = sigma_P(ix,iy,iz) * sigma_0(ix,iy,iz) / tmp1

    sigma_tp(ix,iy,iz) = -sigma_H(ix,iy,iz) * sigma_0(ix,iy,iz) * sinI(ix,iy,iz) / tmp1

    sigma_pp(ix,iy,iz) = sigma_P(ix,iy,iz) &
      &                + ( sigma_H(ix,iy,iz) * cosI(ix,iy,iz) )**2.0d0 / tmp1

  end do
  end do
  end do


  ! height integrated ionospheric conductivity ---------------------------------------------
  hi_sigma_tt = 0.0_dp
  hi_sigma_tp = 0.0_dp
  hi_sigma_pp = 0.0_dp
  do iy = 1, grd%ny
  do ix = 1, grd%nx
    do iz = 1, grd%nz
      hi_sigma_tt(ix,iy) = hi_sigma_tt(ix,iy) + sigma_tt(ix,iy,iz)*grd%dalt(iz)
      hi_sigma_tp(ix,iy) = hi_sigma_tp(ix,iy) + sigma_tp(ix,iy,iz)*grd%dalt(iz)
      hi_sigma_pp(ix,iy) = hi_sigma_pp(ix,iy) + sigma_pp(ix,iy,iz)*grd%dalt(iz)
    end do
  end do 
  end do


  !modulation written by Nakamezo et al., 2012
    !theta-theta
  do ix = 1, grd%nx
      hi_sigma_tt(ix,91) = 2.0d0 * hi_sigma_tt(ix,90) - hi_sigma_tt(ix,89)
  end do

      !phi-phi
  do ix = 1, grd%nx
      do iy = 66, 88
          hi_sigma_pp(ix,iy) = ( (dble(iy) - 66.0d0) / 23.0d0 + 1.0d0 ) * hi_sigma_pp(ix,iy)
      end do
      hi_sigma_pp(ix,89) = 2.0d0 * hi_sigma_pp(ix,89)
      hi_sigma_pp(ix,90) = 1.8d0 * hi_sigma_pp(ix,90)
      hi_sigma_pp(ix,91) = 1.2d0 * hi_sigma_pp(ix,91)
      hi_sigma_pp(ix,92) = 1.8d0 * hi_sigma_pp(ix,92)
      hi_sigma_pp(ix,93) = 2.0d0 * hi_sigma_pp(ix,93)
      do j = 94, 116
          hi_sigma_pp(ix,iy) = ( -(dble(iy) - 93.0d0) / 23.0d0 + 2.0d0 ) * hi_sigma_pp(ix,iy)
      end do
  end do

      !theta-phi
  do iz = 1, grd%nz
  do iy = 62, 120
  do ix = 1, grd%nx
      sigma_prime_tp(ix,iy,iz) = -sigma_H(ix,iy,iz) * sinI(ix,iy,iz)

      sigma_prime_tz(ix,iy,iz) = (sigma_0(ix,iy,iz) - sigma_P(ix,iy,iz)) &
  &                             * sinI(ix,iy,iz) * cosI(ix,iy,iz)

      sigma_prime_zp(ix,iy,iz) = sigma_H(ix,iy,iz) * cosI(ix,iy,iz)

      sigma_prime_zz(ix,iy,iz) = sigma_0(ix,iy,iz) * sinI(ix,iy,iz)**2.0d0 &
  &                             + sigma_P(ix,iy,iz) * cosI(ix,iy,iz)**2.0d0
  end do
  end do
  end do

  do ix = 1, grd%nx
      do iy = 62, 91
      do iz = 1, grd%nz
          sigma_tp(ix,iy,iz) &
          = sigma_prime_tp(ix,iy,iz) - (-(dble(iy)-61.0d0) / 30.0d0 + 1.0d0) &
  &           * sigma_prime_tz(ix,iy,iz) * sigma_prime_zp(ix,iy,iz) / sigma_prime_zz(ix,iy,iz)
      end do
      end do

      do iy = 92, 120
      do iz = 1, grd%nz
          sigma_tp(ix,iy,iz) &
          = sigma_prime_tp(ix,iy,iz) - (dble(iy)-91.0d0) / 30.0d0 &
  &           * sigma_prime_tz(ix,iy,iz) * sigma_prime_zp(ix,iy,iz) / sigma_prime_zz(ix,iy,iz)
      end do
      end do
  end do

  do iy = 62, 120
  do ix = 1, grd%nx
      hi_sigma_tp(ix,iy) = 0.0d0
  end do
  end do

  do iy = 62, 120
  do ix = 1, grd%nx
      do iz = 1, (grd%nz-1)
          hi_sigma_tp(ix,iy) = hi_sigma_tp(ix,iy) + 0.5d0 * grd%dalt(iz) &
  &                       * (sigma_tp(ix,iy,iz) + sigma_tp(i,j,iz+1))
      end do
  end do
  end do

  !minimum values
  sigma_min = 1.0d100

  do iy = 1, grd%ny
  do ix = 1, grd%nx
      if( hi_sigma_tt(ix,iy) /= 0.0d0 &
  &           .and. hi_sigma_tt(ix,iy) < sigma_min ) then
          sigma_min = hi_sigma_tt(ix,iy)
      end if
      if( hi_sigma_pp(ix,iy) /= 0.0d0 &
  &           .and. hi_sigma_pp(ix,iy) < sigma_min ) then
          sigma_min = hi_sigma_pp(ix,iy)
      end if
      if( hi_sigma_tp(ix,iy) /= 0.0d0 &
  &           .and. hi_sigma_tp(ix,iy) < sigma_min ) then
          sigma_min = hi_sigma_tp(ix,iy)
      end if
  end do
  end do

  do iy = 1, grd%ny
  do ix = 1, grd%nx
      if( hi_sigma_tt(ix,iy) == 0.0d0 ) then
          hi_sigma_tt(ix,iy) = sigma_min * 0.01d0
      end if
      if( hi_sigma_pp(ix,iy) == 0.0d0 ) then
          hi_sigma_pp(ix,iy) = sigma_min * 0.01d0
      end if
      if( hi_sigma_tp(ix,iy) == 0.0d0 ) then
          hi_sigma_tp(ix,iy) = sigma_min * 0.01d0
      end if
  end do
  end do

  write(*,*) sigma_min, sigma_min, sigma_min
  write(*,*) hi_sigma_tt(271,91), hi_sigma_tp(271,91), hi_sigma_pp(271,91)

  !!!!!!!!!!
  ! smoothing
  do s = 1, 3
      do iy = 2, 180
      do ix = 2, 360
          tmparr(ix,iy) = ( 4.0d0 * hi_sigma_tt(ix,iy) + hi_sigma_tt(ix-1,iy) + hi_sigma_tt(ix+1,iy) &
  &                    + hi_sigma_tt(ix,iy-1) + hi_sigma_tt(ix,iy+1) ) * 0.125d0
      end do
      end do

      do i = 2, 360
          tmparr(i,1) = ( 2.0d0 * hi_sigma_tt(ix,1) + hi_sigma_tt(ix-1,1) &
  &                    + hi_sigma_tt(ix+1,1) ) * 0.25d0
          tmparr(i,181) = ( 2.0d0 * hi_sigma_tt(ix,181) + hi_sigma_tt(ix-1,181) &
  &                      + hi_sigma_tt(ix+1,181) ) * 0.25d0
      end do

      hi_sigma_tt(2:360, 1:181) = tmparr(2:360, 1:181)
      hi_sigma_tt(1, 1:181) = hi_sigma_tt(360, 1:181)
      hi_sigma_tt(361, 1:181) = hi_sigma_tt(1, 1:181)

      do j = 2, 180
      do i = 2, 360
          tmparr(ix,iy) = ( 4.0d0 * hi_sigma_pp(ix,iy) + hi_sigma_pp(ix-1,iy) + hi_sigma_pp(ix+1,iy) &
  &                    + hi_sigma_pp(ix,iy-1) + hi_sigma_pp(ix,iy+1) ) * 0.125d0
      end do
      end do

      do i = 2, 360
          tmparr(i,1) = ( 2.0d0 * hi_sigma_pp(ix,1) + hi_sigma_pp(ix-1,1) &
  &                    + hi_sigma_pp(ix+1,1) ) * 0.25d0
          tmparr(i,181) = ( 2.0d0 * hi_sigma_pp(ix,181) + hi_sigma_pp(ix-1,181) &
  &                      + hi_sigma_pp(ix+1,181) ) * 0.25d0
      end do

      hi_sigma_pp(2:360, 1:181) = tmparr(2:360, 1:181)
      hi_sigma_pp(1, 1:181) = hi_sigma_pp(360, 1:181)
      hi_sigma_pp(361, 1:181) = hi_sigma_pp(1, 1:181)

      do j = 2, 180
      do i = 2, 360
          tmparr(ix,iy) = ( 4.0d0 * hi_sigma_tp(ix,iy) + hi_sigma_tp(ix-1,iy) + hi_sigma_tp(ix+1,iy) &
  &                    + hi_sigma_tp(ix,iy-1) + hi_sigma_tp(ix,iy+1) ) * 0.125d0
      end do
      end do

      do i = 2, 360
          tmparr(i,1) = ( 2.0d0 * hi_sigma_tp(i,1) + hi_sigma_tp(ix-1,1) &
  &                    + hi_sigma_tp(ix+1,1) ) * 0.25d0
          tmparr(i,181) = ( 2.0d0 * hi_sigma_tp(i,181) + hi_sigma_tp(ix-1,181) &
  &                      + hi_sigma_tp(ix+1,181) ) * 0.25d0
      end do

      hi_sigma_tp(2:360, 1:181) = tmparr(2:360, 1:181)
      hi_sigma_tp(1, 1:181) = hi_sigma_tp(360, 1:181)
      hi_sigma_tp(361, 1:181) = hi_sigma_tp(1, 1:181)
  end do
  !!!!!!!!!!

  ! output ----------------------------------

  fname = trim(ADJUSTL(set%dir_name))//'/output/conductivity/sigma_HP0.dat'
  open(11, file = fname, status = 'unknown' )
    do iz = 1, grd%nz 
    do iy = 1, grd%ny
    do ix = 1, grd%nx
      write(11, *) iy-91,ix, sigma_H(ix,iy,iz), sigma_P(ix,iy,iz), sigma_0(ix,iy,iz)
    end do
    end do
    end do
  close(11)

  fname = trim(ADJUSTL(set%dir_name))//'/output/conductivity/hi_sigma_HP0.dat'
  open(11, file = fname, status = 'unknown' )
    do iy = 1, grd%ny
    do ix = 1, grd%nx
      write(11, *) iy-91,ix, hi_sigma_H(ix,iy), hi_sigma_P(ix,iy), hi_sigma_0(ix,iy)
    end do
    end do
  close(11)

  fname = trim(ADJUSTL(set%dir_name))//'/output/conductivity/hi_sigma_HP0_ij.dat'
  open(11, file = fname, status = 'unknown' )
    do iy = 1, grd%ny
    do ix = 1, grd%nx
      write(11, *) iy-91,ix, log10(hi_sigma_H(ix,iy)), log10(hi_sigma_P(ix,iy)), log10(hi_sigma_0(ix,iy))
    end do
      write(11, *)
    end do
  close(11)

  fname = trim(ADJUSTL(set%dir_name))//'/output/conductivity/hi_sigma_t-p.dat'
  open(11, file = fname, status = 'unknown' )
    do iy = 1, grd%ny
    do ix = 1, grd%nx
      write(11, *) iy-91,ix, hi_sigma_tt(ix,iy), hi_sigma_tp(ix,iy), hi_sigma_pp(ix,iy)
    end do
    end do
  close(11)




  deallocate(B_abs,sinI,cosI)
  deallocate(sigma_H)
  deallocate(sigma_P)
  deallocate(sigma_0)
  deallocate(hi_sigma_H)
  deallocate(hi_sigma_P)
  deallocate(hi_sigma_0)
  deallocate(sigma_tt)
  deallocate(sigma_tp)
  deallocate(sigma_pp)
  deallocate(hi_sigma_tt)
  deallocate(hi_sigma_tp)
  deallocate(hi_sigma_pp)
  deallocate(tmparr)
  deallocate(sigma_prime_tp)
  deallocate(sigma_prime_tz)
  deallocate(sigma_prime_zp)
  deallocate(sigma_prime_zz)

  call p__io_deallocate(var, spl, xct, flx) ! inout

  print *, 'finish'


end program e__conductivity
