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

include "./Jupiter/no_metal_Hill/v__in.f90"

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
  integer i, j, k, ix, xs, iy, iz, ich, isp, jsp, is, iday, s, NB, iB
  real(dp) tmp
  character(len=256) fname, num, char, command, outdir, dir

  real(dp) a0, aH2, nu_e_H2, nu_i_H2, omega_e, omega_i, n_H2, m_H2, n_e, n_i, m_i, m_e
  real(dp) tmp1, tmp2, tmp3
  real(dp), allocatable :: B_abs(:,:,:), sinI(:,:,:), cosI(:,:,:)
  real(dp), allocatable :: sigma_H(:), sigma_P(:), sigma_0(:)
  real(dp)  hi_sigma_P, hi_sigma_H, hi_sigma_0
  real(dp), allocatable :: tmparr(:,:)
  real(dp) B_factor, sigma_min, jpara

  integer ij
  

  real(dp) t1, t2

  call cpu_time(t1)


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

  jpara = 0.40

  !dir = 'metal_x1.0'
  dir = 'no_metal'

  write(num,'(f4.2)') jpara
  !print *, num
  !stop

  print *, 'reading density data...'
    do isp = 1, spl%nsp
      fname = trim(ADJUSTL(set%dir_name))//'/output/density/jparatest/'//trim(ADJUSTL(dir))//'/'//trim(ADJUSTL(num))//&
      &'/num/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
      !print *, trim(ADJUSTL(fname))
      open(11, file = fname, status = 'unknown' )
        do iz = 1, grd%nz
          read(11, *) tmp, var%ni(isp,iz)
        end do
      close(11)
    end do
  print *, 'finished reading density data!'

  !-----------------------------------------------------
  ! Magnetic field factor
  !-----------------------------------------------------
  NB = 5

  ! i = 0  ->  0.1
  ! i = N  ->  1
  ! i = 2N -> 10

  do iB = NB, NB

    B_factor = 10.0_dp**(dble(iB-NB)/dble(NB))
    print *, iB, ':   B_factor = ', B_factor

    allocate(B_abs(361,181,141),sinI(361,181,141),cosI(361,181,141))
    allocate(sigma_H(grd%nz))
    allocate(sigma_P(grd%nz))
    allocate(sigma_0(grd%nz))

    !-----------------------------------------------------
    ! Magnetic field data 
    !-----------------------------------------------------
    if (spl%planet == 'Jupiter') then 

      print *, 'reading B data...'
      open(11, file = './Jupiter/data/dipole_Jupiter.dat', status = 'unknown')
        do iz = 1, 141
        do iy = 1, 181
        do ix = 1, 361
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

    iy = 164
    ix = 181

    do iz = 1, grd%nz

      jsp  = sp_index(spl, 'H2')
      n_H2 = var%ni(jsp,iz)
      m_H2 = var%m(jsp)

      jsp  = sp_index(spl, 'e-')
      n_e  = var%ni(jsp,iz)
      m_e  = var%m(jsp)

      do isp = 1, spl%nsp

        n_i  = var%ni(isp,iz)
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
          sigma_P(iz) = sigma_P(iz) &
            &               + cst%q_e * n_e / B_abs(ix,iy,iz) &
            &               * omega_e * nu_e_H2 &
            &               / ( nu_e_H2**2.0d0 + omega_e**2.0d0 )

          ! Hall conductivity : electron term
          sigma_H(iz) = sigma_H(iz) &
            &               + cst%q_e * n_e / B_abs(ix,iy,iz) &
            &               * omega_e**2.0d0 &
            &               / ( nu_e_H2**2.0d0 + omega_e**2.0d0 )

          ! Parallel conductivity : electron term
          sigma_0(iz) = sigma_0(iz) &
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
          sigma_P(iz) = sigma_P(iz) &
            &               + cst%q_e * n_i / B_abs(ix,iy,iz) &
            &               * omega_i * nu_i_H2  &
            &               / ( nu_i_H2**2.0d0 + omega_i**2.0d0 )

          ! Hall conductivity : ion term
          sigma_H(iz) = sigma_H(iz) &
            &               - cst%q_e * n_i / B_abs(ix,iy,iz) &
            &               * omega_i**2.0d0  &
            &               / ( nu_i_H2**2.0d0 + omega_i**2.0d0 )

          ! Parallel conductivity : ion term
          sigma_0(iz) = sigma_0(iz) &
            &               + cst%q_e * n_i / B_abs(ix,iy,iz) &
            &               * omega_i / nu_i_H2 

        end if

      end do

    end do

    ! height integrated ionospheric conductivity ---------------------------------------------
    hi_sigma_H = 0.0_dp
    hi_sigma_P = 0.0_dp
    hi_sigma_0 = 0.0_dp
    do iz = 1, grd%nz
      hi_sigma_H = hi_sigma_H + sigma_H(iz)*grd%dalt(iz)
      hi_sigma_P = hi_sigma_P + sigma_P(iz)*grd%dalt(iz)
      hi_sigma_0 = hi_sigma_0 + sigma_0(iz)*grd%dalt(iz)
    end do


    ! output ----------------------------------

    write(num,'(f4.2)') jpara

    outdir = trim(ADJUSTL(set%dir_name))//'/output/conductivity/'//trim(ADJUSTL(dir))//'/jpara_'//trim(ADJUSTL(num))
    write(command,*) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
    call system(command)

    fname = trim(ADJUSTL(outdir))//'/sigma_HP0.dat'
    open(11, file = fname, status = 'unknown' )
      do iz = 1, grd%nz
        write(11, *) grd%alt(iz), sigma_H(iz), sigma_P(iz), sigma_0(iz)
      end do
    close(11)

    fname = trim(ADJUSTL(outdir))//'/hi_sigma_HP0.dat'
    open(11, file = fname, status = 'unknown' )
        write(11, *) hi_sigma_H, hi_sigma_P, hi_sigma_0
    close(11)

    deallocate(B_abs,sinI,cosI)
    deallocate(sigma_H)
    deallocate(sigma_P)
    deallocate(sigma_0)

  end do

  call p__io_deallocate(var, spl, xct, flx) ! inout

  print *, 'finish'


end program e__conductivity
