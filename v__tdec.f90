module v__tdec
  !---------------------------------------------------------------------
  ! Type declaration statements
  !---------------------------------------------------------------------
  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private
  public :: grd_, var_, cst_, xct_, spl_, flx_, set_

  ! settings
  type set_
    integer(sp)            :: nx, ny, resx, resy            ! number of grid index
    real(dp) latitude, sza, Ls, DOY
    real(dp) fin_sec, dtime_limit
    integer nday, nlat, nstep
    integer calc_stable, start_rot, test_loc, read_stable
    character(len=256) scheme, inversion, mode
    character(len=256) fnamestable, dir_name
  end type set_

  ! grid
  type grd_
    integer(sp)            :: nx, ny, nz               ! number of grid index
    integer(sp)            :: xs, xe, ys, ye, zs, ze   ! grid index range
    integer(sp)            :: ix, iy, iday
    real(dp), allocatable  :: dalt(:), alt(:), lt(:), lat(:), sza(:,:), sza_xact(:,:)   ! altitude, longitude, latitude
    real(dp)               :: tilt, season, latitude ! solar zenith angle, tilt angle, season
  end type grd_

  ! variables
  type var_
    ! 種類ごとに分けて欲しい
    !-------------------------------------------------------------------
    ! nstep     : total number of time steps
    ! istep     : index of time step
    ! species_i : species for searching cross sections
    ! reactants : reactant list
    ! products  : product list
    ! m         : mass of species
    ! m_major   : mass of major species
    ! q         : charge of species
    !
    !-------------------------------------------------------------------
    integer                :: nsteps, istep, nspecial, iter ! number and index of time steps
    real(dp)               :: t1, t2, t3, t4
    character(len=256)     :: species_i, reactants(10), products(10)
    real(dp), allocatable  :: m(:), m_mean(:), q(:)
    real(dp), allocatable  :: ni(:,:), ni_stable(:,:,:), ni_3d(:,:,:,:), ni_new(:,:), n_tot(:)   ! density
    real(dp), allocatable  :: ni_0(:,:)
    real(dp), allocatable  :: clm_ni(:,:) ! column density
    real(dp), allocatable  :: Ti(:), Te(:), Tn(:)
    real(dp), allocatable  :: Ti_3d(:,:,:), Te_3d(:,:,:), Tn_3d(:,:,:)
    real(dp), allocatable  :: tau_EUV(:,:), tau_UV(:,:)
    real(dp), allocatable  :: E_fld(:,:,:), B_fld(:,:,:) ! electric and magnetic field
    real(dp)               :: dtime, sum_time
    real(dp), allocatable  :: ki(:,:), ki_special(:,:,:,:), ich_special(:), Pi(:,:), Li(:,:), Jmtx(:,:), rate(:,:)
    real(dp), allocatable  :: Fluxup(:,:), Fluxdwn(:,:), vFluxup(:,:), vFluxdwn(:,:)
    real(dp), allocatable  :: UpperBC(:,:), LowerBC(:,:) ! upper and lower boundary condition: label, (1:density, 2:flux, 3:velocity)
    real(dp), allocatable  :: K_eddy(:), D_mol(:,:)
    real(dp), allocatable  :: I_EUV(:,:), I_UV(:,:)
    real(dp), allocatable  :: Phip(:,:), Phim(:,:) ! 
    real(dp), allocatable  :: dPhi_dz(:,:)         ! dPhi/dz
    real(dp), allocatable  :: d_dniu_dPhi_dz(:,:)  ! d/dniu * dPhi/dz
    real(dp), allocatable  :: d_dni0_dPhi_dz(:,:)  ! d/dni0 * dPhi/dz
    real(dp), allocatable  :: d_dnil_dPhi_dz(:,:)  ! d/dnil * dPhi/dz
    real(dp), allocatable  :: Amtx(:,:), Lmtx(:,:), Umtx(:,:), barr(:), xarr(:), yarr(:), dxarr(:), rarr(:)
    real(dp)               :: max_dn_n(3)
  end type var_

  ! species list
  type spl_
    !-------------------------------------------------------------------
    ! planet        : selected planet e.g. Mars, Jupiter, ...
    ! nsp           : total number of species
    ! nsp_i         : total number of variable species
    ! nsp_f         : total number of fixed species
    ! nch           : total number of chemical reactions
    ! nch_P         : total number of columns for reading Production list
    ! nch_L         : total number of columns for reading Loss list
    ! nch_J         : total number of columns for reading Jacobian list
    ! n_Jlist       : total number of lines for reading Jacobian list
    !
    ! species       : species list
    ! reactant_list :
    !
    !-------------------------------------------------------------------
    character(len=256)     :: planet
    integer                :: nsp, nsp_i, nch, nch_P, nch_L, nch_J, n_Jlist, nrpn
    character(len=256), allocatable :: species(:), reaction_type_char(:)
    integer,  allocatable  :: reactant_list(:,:), product_list(:,:), Prod_list(:,:), Loss_list(:,:)
    integer,  allocatable  :: reaction_type_list(:), label_fix(:), all_to_var(:), var_to_all(:), major_species(:)
    integer,  allocatable  :: Jmtx_list(:,:)
    integer,  allocatable  :: rate_rpn_label(:,:,:)
    real(dp), allocatable  :: rate_rpn_token(:,:,:), T_range(:,:,:)
    integer,  allocatable  :: rate_cases(:)
  end type spl_

  ! physical constant
  type cst_
    real(dp)               :: pi, c, h
    real(dp)               :: NA  ! Avogadro constant
    real(dp)               :: k_B ! Boltzmann constant
    real(dp)               :: g   ! gravitational acceleration
    real(dp)               :: BigG ! gravitational constant
    real(dp)               :: Mplanet ! gravitational acceleration
    real(dp)               :: R   ! planetary radius
    real(dp)               :: m_u ! atomic mass unit
    real(dp)               :: q_e
    real(dp)               :: eV  ! electron volt
    real(dp)               :: planck 
    real(dp)               :: daysec
  end type cst_

  ! cross sections
  type xct_
    character(len=256)    :: type          ! type of cross section
    real(dp), allocatable :: sigma_a_EUV(:,:)    ! EUV absorption cross section
    real(dp), allocatable :: sigma_i_EUV(:,:)    ! EUV ionization cross section
    real(dp), allocatable :: sigma_a_UV(:,:,:) ! UV absorption cross section
    real(dp), allocatable :: sigma_d_UV(:,:,:) ! UV dissociation cross section
    real(dp), allocatable :: sigma_e_UV(:) ! UV excitation cross section
    real(dp), allocatable :: sigma_a_UV_EUV(:,:)    ! EUV absorption cross section

    real(dp)              :: co2xdata(159,3)
    real(dp)              :: h2oxdata(99,2)
    real(dp)              :: h2o2xdata(70,2), h2o2xdata_l(90,2)
    real(dp)              :: h2xdata(277,2)
    real(dp)              :: ohxdata(601,2), oho1Dxdata(601,2)
    real(dp)              :: o3xdata(215+1601,2)
    real(dp)              :: o2xdata(115,2), o2schr130K(30,4),o2schr190K(30,4), o2schr280K(30,4)
  end type xct_

  ! solar flux
  type flx_
    real(dp)              :: orbit ! orbit radius [AU]
    real(dp)              :: mode_factor

    integer               :: nwl_EUV
    real(dp), allocatable :: lambda_EUV(:), F74113(:), Ai(:), F10_7, solar_EUV(:)
    real(dp)              :: multiplying_factor_EUV

    integer               :: nwl_UV
    real(dp), allocatable :: lambda_UV(:), solar_UV(:)
    real(dp)              :: multiplying_factor_UV

  end type flx_


end module v__tdec
