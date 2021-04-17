module v__in
  use v__tdec,   only : set_, var_, cst_, spl_, grd_
  use p__search, only : p__search_reactant, p__search_product, sp_index

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: v__in__ini, v__in__exe

contains


  subroutine v__in__ini(spl, set) ! out
    type(spl_),           intent(out) :: spl
    type(set_),           intent(out) :: set

    ! Planet type
    spl%planet = 'Jupiter'

    ! Calculation settings
    set%mode = '2D Lat'
    set%nstep = 30000
    set%fin_sec = 35729.685e3_dp
    set%dtime_limit = 1.0e5_dp
    set%latitude = 0.0_dp
    set%Ls = 0.0_dp
    set%nday = 3.0_dp
    set%scheme = 'implicit'
    set%inversion = 'Catling'
    ! directory setting
    set%dir_name = './Jupiter/metal_dev10'
    set%fnamestable = './Jupiter/metal_dev10/output/density/n_stable.dat'

  end subroutine v__in__ini


  subroutine v__in__exe(cst,      & ! in
    &                   set, spl, & ! inout
    &                   var, grd  ) ! out
    type(cst_),           intent(in)    :: cst
    type(set_),           intent(inout) :: set
    type(spl_),           intent(inout) :: spl
    type(var_),           intent(out)   :: var
    type(grd_),           intent(out)   :: grd
    integer i, j, ip, isp, jsp, ich, jch, iz, nh
    real(dp) tmp
    character(len = 256) strm, fname


    ! grid setting
    grd%nx    = set%nx
    grd%ny    = set%ny
    grd%nz    = 141
    allocate(grd%alt(grd%nz),grd%dalt(grd%nz))
    grd%dalt(1:141) = 20.0e3_dp ! [m]
    grd%alt(1)      = 200.0e3_dp ! [m]
    do iz = 2, grd%nz
      grd%alt(iz) = grd%alt(iz-1) + grd%dalt(iz)
    end do

    ! reactions, chemical species
    spl%nsp     = 88
    spl%nsp_i   = 69
    spl%nch     = 408
    spl%nch_P   = 203
    spl%nch_L   = 159
    spl%n_Jlist = 746
    spl%nch_J   = 161
    spl%nrpn    = 8

    ! allocate
    allocate(var%ni(0:spl%nsp,grd%nz))
    allocate(var%ni_0(0:spl%nsp,grd%nz))
    allocate(var%ni_new(spl%nsp,grd%nz))
    allocate(var%ni_stable(spl%nsp,grd%ny,grd%nz))
    allocate(var%ni_3d(spl%nsp,grd%nx,grd%ny,grd%nz))
    allocate(var%clm_ni(spl%nsp,grd%nz))
    allocate(var%Ti(grd%nz),var%Te(grd%nz),var%Tn(grd%nz))
    allocate(var%Ti_3d(grd%nx,grd%ny,grd%nz))
    allocate(var%Te_3d(grd%nx,grd%ny,grd%nz))
    allocate(var%Tn_3d(grd%nx,grd%ny,grd%nz))
    allocate(var%m(spl%nsp), var%q(spl%nsp))
    allocate(spl%reactant_list(spl%nch, 10))
    allocate(spl%product_list(spl%nch, 10))
    allocate(spl%species(spl%nsp))
    allocate(spl%label_fix(spl%nsp))
    allocate(spl%all_to_var(spl%nsp))
    allocate(spl%var_to_all(spl%nsp_i))
    allocate(spl%Prod_list(spl%nsp_i, spl%nch_P))
    allocate(spl%Loss_list(spl%nsp_i, spl%nch_L))
    allocate(spl%Jmtx_list(spl%n_Jlist, spl%nch_J))
    allocate(spl%reaction_type_list(spl%nch))
    allocate(spl%reaction_type_char(spl%nch))
    allocate(var%ki(spl%nch,grd%nz))
    allocate(var%rate(spl%nch,grd%nz))
    allocate(var%Pi(spl%nsp_i,grd%nz))
    allocate(var%Li(spl%nsp_i,grd%nz))
    allocate(var%K_eddy(grd%nz),var%D_mol(spl%nsp,grd%nz))
    allocate(var%Fluxup(spl%nsp_i,0:grd%nz),var%Fluxdwn(spl%nsp_i,0:grd%nz))
    allocate(var%vFluxup(spl%nsp_i,grd%nz),var%vFluxdwn(spl%nsp_i,grd%nz))
    allocate(var%UpperBC(spl%nsp,3),var%LowerBC(spl%nsp,3))
    allocate(var%Jmtx(spl%nsp_i, spl%nsp_i))
    allocate(spl%rate_rpn_token(spl%nch,3,spl%nrpn))
    allocate(spl%rate_rpn_label(spl%nch,3,spl%nrpn))
    allocate(spl%rate_cases(spl%nch))
    allocate(spl%T_range(spl%nch,3,3))
    allocate(spl%major_species(grd%nz))
    allocate(var%n_tot(grd%nz),var%m_mean(grd%nz))
    allocate(var%Phip(spl%nsp_i,grd%nz))
    allocate(var%Phim(spl%nsp_i,grd%nz))
    allocate(var%dPhi_dz(spl%nsp_i,grd%nz))
    allocate(var%d_dniu_dPhi_dz(spl%nsp_i,grd%nz))
    allocate(var%d_dni0_dPhi_dz(spl%nsp_i,grd%nz))
    allocate(var%d_dnil_dPhi_dz(spl%nsp_i,grd%nz))
    allocate(var%barr(spl%nsp_i*grd%nz), var%xarr(spl%nsp_i*grd%nz))
    allocate(var%yarr(spl%nsp_i*grd%nz), var%dxarr(spl%nsp_i*grd%nz))
    allocate(var%Amtx(spl%nsp_i*grd%nz,2*spl%nsp_i+1))
    allocate(var%Lmtx(spl%nsp_i*grd%nz,2*spl%nsp_i+1))
    allocate(var%Umtx(spl%nsp_i*grd%nz,2*spl%nsp_i+1))

    ! species
    spl%species(1) = 'H'
    spl%species(2) = 'H+'
    spl%species(3) = 'e-'
    spl%species(4) = 'H2'
    spl%species(5) = 'H2+'
    spl%species(6) = 'He'
    spl%species(7) = 'He+'
    spl%species(8) = 'CH4'
    spl%species(9) = 'CH4+'
    spl%species(10) = 'CH3+'
    spl%species(11) = 'C2H2'
    spl%species(12) = 'C2H2+'
    spl%species(13) = 'C2H4'
    spl%species(14) = 'C2H4+'
    spl%species(15) = 'C2H3+'
    spl%species(16) = 'C2H+'
    spl%species(17) = 'C2H6'
    spl%species(18) = 'C2H6+'
    spl%species(19) = 'C2H5+'
    spl%species(20) = 'HeH+'
    spl%species(21) = 'H3+'
    spl%species(22) = 'C+'
    spl%species(23) = 'C'
    spl%species(24) = 'CH+'
    spl%species(25) = 'CH2+'
    spl%species(26) = 'CH'
    spl%species(27) = 'CH2'
    spl%species(28) = 'CH3'
    spl%species(29) = 'CH5+'
    spl%species(30) = 'C2+'
    spl%species(31) = 'C2'
    spl%species(32) = 'C2H'
    spl%species(33) = 'C2H3'
    spl%species(34) = 'C2H5'
    spl%species(35) = 'C2H7+'
    spl%species(36) = 'C3H+'
    spl%species(37) = 'C3H2+'
    spl%species(38) = 'C3H3+'
    spl%species(39) = 'C3H4+'
    spl%species(40) = 'C3H5+'
    spl%species(41) = 'C3H6+'
    spl%species(42) = 'C3H7+'
    spl%species(43) = 'C3H8+'
    spl%species(44) = 'C3H9+'
    spl%species(45) = 'C4H+'
    spl%species(46) = 'C4H2+'
    spl%species(47) = 'C4H3+'
    spl%species(48) = 'C4H5+'
    spl%species(49) = 'C4H7+'
    spl%species(50) = 'C4H9+'
    spl%species(51) = 'H2(v>=2)'
    spl%species(52) = 'H2(v>=4)'
    spl%species(53) = 'Na'
    spl%species(54) = 'Na+'
    spl%species(55) = 'Fe'
    spl%species(56) = 'Fe+'
    spl%species(57) = 'Mg'
    spl%species(58) = 'Mg+'
    spl%species(59) = 'Si'
    spl%species(60) = 'Si+'
    spl%species(61) = 'NaH'
    spl%species(62) = 'NaCH3'
    spl%species(63) = 'NaH2+'
    spl%species(64) = 'NaCH4+'
    spl%species(65) = 'NaC2H2+'
    spl%species(66) = 'NaC2H4+'
    spl%species(67) = 'H3O+'
    spl%species(68) = 'H2O'
    spl%species(69) = 'FeH'
    spl%species(70) = 'FeH+'
    spl%species(71) = 'FeH2+'
    spl%species(72) = 'FeCH4+'
    spl%species(73) = 'FeC2H2+'
    spl%species(74) = 'FeC2H4+'
    spl%species(75) = 'MgH'
    spl%species(76) = 'MgH+'
    spl%species(77) = 'MgH2+'
    spl%species(78) = 'MgCH4+'
    spl%species(79) = 'MgC2H2+'
    spl%species(80) = 'MgC2H4+'
    spl%species(81) = 'SiH'
    spl%species(82) = 'SiH+'
    spl%species(83) = 'SiH2+'
    spl%species(84) = 'SiCnHm+'
    spl%species(85) = 'SiOH+'
    spl%species(86) = 'H3O'
    spl%species(87) = 'SiCH'
    spl%species(88) = 'SiCH2'

    ! label_fix
    spl%label_fix(1) = 0 ! H: variable
    spl%label_fix(2) = 0 ! H+: variable
    spl%label_fix(3) = 0 ! e-: variable
    spl%label_fix(4) = 1 ! H2: fixed
    spl%label_fix(5) = 0 ! H2+: variable
    spl%label_fix(6) = 1 ! He: fixed
    spl%label_fix(7) = 0 ! He+: variable
    spl%label_fix(8) = 1 ! CH4: fixed
    spl%label_fix(9) = 0 ! CH4+: variable
    spl%label_fix(10) = 0 ! CH3+: variable
    spl%label_fix(11) = 1 ! C2H2: fixed
    spl%label_fix(12) = 0 ! C2H2+: variable
    spl%label_fix(13) = 1 ! C2H4: fixed
    spl%label_fix(14) = 0 ! C2H4+: variable
    spl%label_fix(15) = 0 ! C2H3+: variable
    spl%label_fix(16) = 0 ! C2H+: variable
    spl%label_fix(17) = 1 ! C2H6: fixed
    spl%label_fix(18) = 0 ! C2H6+: variable
    spl%label_fix(19) = 0 ! C2H5+: variable
    spl%label_fix(20) = 0 ! HeH+: variable
    spl%label_fix(21) = 0 ! H3+: variable
    spl%label_fix(22) = 0 ! C+: variable
    spl%label_fix(23) = 1 ! C: fixed
    spl%label_fix(24) = 0 ! CH+: variable
    spl%label_fix(25) = 0 ! CH2+: variable
    spl%label_fix(26) = 1 ! CH: fixed
    spl%label_fix(27) = 0 ! CH2: variable
    spl%label_fix(28) = 0 ! CH3: variable
    spl%label_fix(29) = 0 ! CH5+: variable
    spl%label_fix(30) = 0 ! C2+: variable
    spl%label_fix(31) = 1 ! C2: fixed
    spl%label_fix(32) = 1 ! C2H: fixed
    spl%label_fix(33) = 1 ! C2H3: fixed
    spl%label_fix(34) = 1 ! C2H5: fixed
    spl%label_fix(35) = 0 ! C2H7+: variable
    spl%label_fix(36) = 0 ! C3H+: variable
    spl%label_fix(37) = 0 ! C3H2+: variable
    spl%label_fix(38) = 0 ! C3H3+: variable
    spl%label_fix(39) = 0 ! C3H4+: variable
    spl%label_fix(40) = 0 ! C3H5+: variable
    spl%label_fix(41) = 0 ! C3H6+: variable
    spl%label_fix(42) = 0 ! C3H7+: variable
    spl%label_fix(43) = 0 ! C3H8+: variable
    spl%label_fix(44) = 0 ! C3H9+: variable
    spl%label_fix(45) = 0 ! C4H+: variable
    spl%label_fix(46) = 0 ! C4H2+: variable
    spl%label_fix(47) = 0 ! C4H3+: variable
    spl%label_fix(48) = 0 ! C4H5+: variable
    spl%label_fix(49) = 0 ! C4H7+: variable
    spl%label_fix(50) = 0 ! C4H9+: variable
    spl%label_fix(51) = 1 ! H2(v>=2): fixed
    spl%label_fix(52) = 1 ! H2(v>=4): fixed
    spl%label_fix(53) = 0 ! Na: variable
    spl%label_fix(54) = 0 ! Na+: variable
    spl%label_fix(55) = 0 ! Fe: variable
    spl%label_fix(56) = 0 ! Fe+: variable
    spl%label_fix(57) = 0 ! Mg: variable
    spl%label_fix(58) = 0 ! Mg+: variable
    spl%label_fix(59) = 0 ! Si: variable
    spl%label_fix(60) = 0 ! Si+: variable
    spl%label_fix(61) = 0 ! NaH: variable
    spl%label_fix(62) = 1 ! NaCH3: fixed
    spl%label_fix(63) = 0 ! NaH2+: variable
    spl%label_fix(64) = 0 ! NaCH4+: variable
    spl%label_fix(65) = 0 ! NaC2H2+: variable
    spl%label_fix(66) = 0 ! NaC2H4+: variable
    spl%label_fix(67) = 0 ! H3O+: variable
    spl%label_fix(68) = 0 ! H2O: variable
    spl%label_fix(69) = 0 ! FeH: variable
    spl%label_fix(70) = 0 ! FeH+: variable
    spl%label_fix(71) = 0 ! FeH2+: variable
    spl%label_fix(72) = 0 ! FeCH4+: variable
    spl%label_fix(73) = 0 ! FeC2H2+: variable
    spl%label_fix(74) = 0 ! FeC2H4+: variable
    spl%label_fix(75) = 0 ! MgH: variable
    spl%label_fix(76) = 0 ! MgH+: variable
    spl%label_fix(77) = 0 ! MgH2+: variable
    spl%label_fix(78) = 0 ! MgCH4+: variable
    spl%label_fix(79) = 0 ! MgC2H2+: variable
    spl%label_fix(80) = 0 ! MgC2H4+: variable
    spl%label_fix(81) = 0 ! SiH: variable
    spl%label_fix(82) = 0 ! SiH+: variable
    spl%label_fix(83) = 0 ! SiH2+: variable
    spl%label_fix(84) = 0 ! SiCnHm+: variable
    spl%label_fix(85) = 1 ! SiOH+: fixed
    spl%label_fix(86) = 1 ! H3O: fixed
    spl%label_fix(87) = 1 ! SiCH: fixed
    spl%label_fix(88) = 1 ! SiCH2: fixed

    ! all_to_var
    spl%all_to_var = 0
    spl%all_to_var(1) = 1 ! H: variable
    spl%all_to_var(2) = 2 ! H+: variable
    spl%all_to_var(3) = 3 ! e-: variable
    spl%all_to_var(5) = 4 ! H2+: variable
    spl%all_to_var(7) = 5 ! He+: variable
    spl%all_to_var(9) = 6 ! CH4+: variable
    spl%all_to_var(10) = 7 ! CH3+: variable
    spl%all_to_var(12) = 8 ! C2H2+: variable
    spl%all_to_var(14) = 9 ! C2H4+: variable
    spl%all_to_var(15) = 10 ! C2H3+: variable
    spl%all_to_var(16) = 11 ! C2H+: variable
    spl%all_to_var(18) = 12 ! C2H6+: variable
    spl%all_to_var(19) = 13 ! C2H5+: variable
    spl%all_to_var(20) = 14 ! HeH+: variable
    spl%all_to_var(21) = 15 ! H3+: variable
    spl%all_to_var(22) = 16 ! C+: variable
    spl%all_to_var(24) = 17 ! CH+: variable
    spl%all_to_var(25) = 18 ! CH2+: variable
    spl%all_to_var(27) = 19 ! CH2: variable
    spl%all_to_var(28) = 20 ! CH3: variable
    spl%all_to_var(29) = 21 ! CH5+: variable
    spl%all_to_var(30) = 22 ! C2+: variable
    spl%all_to_var(35) = 23 ! C2H7+: variable
    spl%all_to_var(36) = 24 ! C3H+: variable
    spl%all_to_var(37) = 25 ! C3H2+: variable
    spl%all_to_var(38) = 26 ! C3H3+: variable
    spl%all_to_var(39) = 27 ! C3H4+: variable
    spl%all_to_var(40) = 28 ! C3H5+: variable
    spl%all_to_var(41) = 29 ! C3H6+: variable
    spl%all_to_var(42) = 30 ! C3H7+: variable
    spl%all_to_var(43) = 31 ! C3H8+: variable
    spl%all_to_var(44) = 32 ! C3H9+: variable
    spl%all_to_var(45) = 33 ! C4H+: variable
    spl%all_to_var(46) = 34 ! C4H2+: variable
    spl%all_to_var(47) = 35 ! C4H3+: variable
    spl%all_to_var(48) = 36 ! C4H5+: variable
    spl%all_to_var(49) = 37 ! C4H7+: variable
    spl%all_to_var(50) = 38 ! C4H9+: variable
    spl%all_to_var(53) = 39 ! Na: variable
    spl%all_to_var(54) = 40 ! Na+: variable
    spl%all_to_var(55) = 41 ! Fe: variable
    spl%all_to_var(56) = 42 ! Fe+: variable
    spl%all_to_var(57) = 43 ! Mg: variable
    spl%all_to_var(58) = 44 ! Mg+: variable
    spl%all_to_var(59) = 45 ! Si: variable
    spl%all_to_var(60) = 46 ! Si+: variable
    spl%all_to_var(61) = 47 ! NaH: variable
    spl%all_to_var(63) = 48 ! NaH2+: variable
    spl%all_to_var(64) = 49 ! NaCH4+: variable
    spl%all_to_var(65) = 50 ! NaC2H2+: variable
    spl%all_to_var(66) = 51 ! NaC2H4+: variable
    spl%all_to_var(67) = 52 ! H3O+: variable
    spl%all_to_var(68) = 53 ! H2O: variable
    spl%all_to_var(69) = 54 ! FeH: variable
    spl%all_to_var(70) = 55 ! FeH+: variable
    spl%all_to_var(71) = 56 ! FeH2+: variable
    spl%all_to_var(72) = 57 ! FeCH4+: variable
    spl%all_to_var(73) = 58 ! FeC2H2+: variable
    spl%all_to_var(74) = 59 ! FeC2H4+: variable
    spl%all_to_var(75) = 60 ! MgH: variable
    spl%all_to_var(76) = 61 ! MgH+: variable
    spl%all_to_var(77) = 62 ! MgH2+: variable
    spl%all_to_var(78) = 63 ! MgCH4+: variable
    spl%all_to_var(79) = 64 ! MgC2H2+: variable
    spl%all_to_var(80) = 65 ! MgC2H4+: variable
    spl%all_to_var(81) = 66 ! SiH: variable
    spl%all_to_var(82) = 67 ! SiH+: variable
    spl%all_to_var(83) = 68 ! SiH2+: variable
    spl%all_to_var(84) = 69 ! SiCnHm+: variable

    ! var_to_all
    spl%var_to_all(1) = 1 ! H: variable
    spl%var_to_all(2) = 2 ! H+: variable
    spl%var_to_all(3) = 3 ! e-: variable
    spl%var_to_all(4) = 5 ! H2+: variable
    spl%var_to_all(5) = 7 ! He+: variable
    spl%var_to_all(6) = 9 ! CH4+: variable
    spl%var_to_all(7) = 10 ! CH3+: variable
    spl%var_to_all(8) = 12 ! C2H2+: variable
    spl%var_to_all(9) = 14 ! C2H4+: variable
    spl%var_to_all(10) = 15 ! C2H3+: variable
    spl%var_to_all(11) = 16 ! C2H+: variable
    spl%var_to_all(12) = 18 ! C2H6+: variable
    spl%var_to_all(13) = 19 ! C2H5+: variable
    spl%var_to_all(14) = 20 ! HeH+: variable
    spl%var_to_all(15) = 21 ! H3+: variable
    spl%var_to_all(16) = 22 ! C+: variable
    spl%var_to_all(17) = 24 ! CH+: variable
    spl%var_to_all(18) = 25 ! CH2+: variable
    spl%var_to_all(19) = 27 ! CH2: variable
    spl%var_to_all(20) = 28 ! CH3: variable
    spl%var_to_all(21) = 29 ! CH5+: variable
    spl%var_to_all(22) = 30 ! C2+: variable
    spl%var_to_all(23) = 35 ! C2H7+: variable
    spl%var_to_all(24) = 36 ! C3H+: variable
    spl%var_to_all(25) = 37 ! C3H2+: variable
    spl%var_to_all(26) = 38 ! C3H3+: variable
    spl%var_to_all(27) = 39 ! C3H4+: variable
    spl%var_to_all(28) = 40 ! C3H5+: variable
    spl%var_to_all(29) = 41 ! C3H6+: variable
    spl%var_to_all(30) = 42 ! C3H7+: variable
    spl%var_to_all(31) = 43 ! C3H8+: variable
    spl%var_to_all(32) = 44 ! C3H9+: variable
    spl%var_to_all(33) = 45 ! C4H+: variable
    spl%var_to_all(34) = 46 ! C4H2+: variable
    spl%var_to_all(35) = 47 ! C4H3+: variable
    spl%var_to_all(36) = 48 ! C4H5+: variable
    spl%var_to_all(37) = 49 ! C4H7+: variable
    spl%var_to_all(38) = 50 ! C4H9+: variable
    spl%var_to_all(39) = 53 ! Na: variable
    spl%var_to_all(40) = 54 ! Na+: variable
    spl%var_to_all(41) = 55 ! Fe: variable
    spl%var_to_all(42) = 56 ! Fe+: variable
    spl%var_to_all(43) = 57 ! Mg: variable
    spl%var_to_all(44) = 58 ! Mg+: variable
    spl%var_to_all(45) = 59 ! Si: variable
    spl%var_to_all(46) = 60 ! Si+: variable
    spl%var_to_all(47) = 61 ! NaH: variable
    spl%var_to_all(48) = 63 ! NaH2+: variable
    spl%var_to_all(49) = 64 ! NaCH4+: variable
    spl%var_to_all(50) = 65 ! NaC2H2+: variable
    spl%var_to_all(51) = 66 ! NaC2H4+: variable
    spl%var_to_all(52) = 67 ! H3O+: variable
    spl%var_to_all(53) = 68 ! H2O: variable
    spl%var_to_all(54) = 69 ! FeH: variable
    spl%var_to_all(55) = 70 ! FeH+: variable
    spl%var_to_all(56) = 71 ! FeH2+: variable
    spl%var_to_all(57) = 72 ! FeCH4+: variable
    spl%var_to_all(58) = 73 ! FeC2H2+: variable
    spl%var_to_all(59) = 74 ! FeC2H4+: variable
    spl%var_to_all(60) = 75 ! MgH: variable
    spl%var_to_all(61) = 76 ! MgH+: variable
    spl%var_to_all(62) = 77 ! MgH2+: variable
    spl%var_to_all(63) = 78 ! MgCH4+: variable
    spl%var_to_all(64) = 79 ! MgC2H2+: variable
    spl%var_to_all(65) = 80 ! MgC2H4+: variable
    spl%var_to_all(66) = 81 ! SiH: variable
    spl%var_to_all(67) = 82 ! SiH+: variable
    spl%var_to_all(68) = 83 ! SiH2+: variable
    spl%var_to_all(69) = 84 ! SiCnHm+: variable

    ! mass
    var%m(1) = 1.00794000_dp * cst%m_u !H
    var%m(2) = 1.00739142_dp * cst%m_u !H+
    var%m(3) = 0.00054858_dp * cst%m_u !e-
    var%m(4) = 2.01588000_dp * cst%m_u !H2
    var%m(5) = 2.01533142_dp * cst%m_u !H2+
    var%m(6) = 4.00260200_dp * cst%m_u !He
    var%m(7) = 4.00205342_dp * cst%m_u !He+
    var%m(8) = 16.04246000_dp * cst%m_u !CH4
    var%m(9) = 16.04191142_dp * cst%m_u !CH4+
    var%m(10) = 15.03397142_dp * cst%m_u !CH3+
    var%m(11) = 26.03728000_dp * cst%m_u !C2H2
    var%m(12) = 26.03673142_dp * cst%m_u !C2H2+
    var%m(13) = 28.05316000_dp * cst%m_u !C2H4
    var%m(14) = 28.05261142_dp * cst%m_u !C2H4+
    var%m(15) = 27.04467142_dp * cst%m_u !C2H3+
    var%m(16) = 25.02879142_dp * cst%m_u !C2H+
    var%m(17) = 30.06904000_dp * cst%m_u !C2H6
    var%m(18) = 30.06849142_dp * cst%m_u !C2H6+
    var%m(19) = 29.06055142_dp * cst%m_u !C2H5+
    var%m(20) = 5.00999342_dp * cst%m_u !HeH+
    var%m(21) = 3.02327142_dp * cst%m_u !H3+
    var%m(22) = 12.01015142_dp * cst%m_u !C+
    var%m(23) = 12.01070000_dp * cst%m_u !C
    var%m(24) = 13.01809142_dp * cst%m_u !CH+
    var%m(25) = 14.02603142_dp * cst%m_u !CH2+
    var%m(26) = 13.01864000_dp * cst%m_u !CH
    var%m(27) = 14.02658000_dp * cst%m_u !CH2
    var%m(28) = 15.03452000_dp * cst%m_u !CH3
    var%m(29) = 17.04985142_dp * cst%m_u !CH5+
    var%m(30) = 24.02085142_dp * cst%m_u !C2+
    var%m(31) = 24.02140000_dp * cst%m_u !C2
    var%m(32) = 25.02934000_dp * cst%m_u !C2H
    var%m(33) = 27.04522000_dp * cst%m_u !C2H3
    var%m(34) = 29.06110000_dp * cst%m_u !C2H5
    var%m(35) = 31.07643142_dp * cst%m_u !C2H7+
    var%m(36) = 37.03949142_dp * cst%m_u !C3H+
    var%m(37) = 38.04743142_dp * cst%m_u !C3H2+
    var%m(38) = 39.05537142_dp * cst%m_u !C3H3+
    var%m(39) = 40.06331142_dp * cst%m_u !C3H4+
    var%m(40) = 41.07125142_dp * cst%m_u !C3H5+
    var%m(41) = 42.07919142_dp * cst%m_u !C3H6+
    var%m(42) = 43.08713142_dp * cst%m_u !C3H7+
    var%m(43) = 44.09507142_dp * cst%m_u !C3H8+
    var%m(44) = 45.10301142_dp * cst%m_u !C3H9+
    var%m(45) = 49.05019142_dp * cst%m_u !C4H+
    var%m(46) = 50.05813142_dp * cst%m_u !C4H2+
    var%m(47) = 51.06607142_dp * cst%m_u !C4H3+
    var%m(48) = 53.08195142_dp * cst%m_u !C4H5+
    var%m(49) = 55.09783142_dp * cst%m_u !C4H7+
    var%m(50) = 57.11371142_dp * cst%m_u !C4H9+
    var%m(51) = 2.01588000_dp * cst%m_u !H2(v>=2)
    var%m(52) = 2.01588000_dp * cst%m_u !H2(v>=4)
    var%m(53) = 22.98900000_dp * cst%m_u !Na
    var%m(54) = 22.98845142_dp * cst%m_u !Na+
    var%m(55) = 55.84500000_dp * cst%m_u !Fe
    var%m(56) = 55.84445142_dp * cst%m_u !Fe+
    var%m(57) = 24.30500000_dp * cst%m_u !Mg
    var%m(58) = 24.30445142_dp * cst%m_u !Mg+
    var%m(59) = 28.08550000_dp * cst%m_u !Si
    var%m(60) = 28.08495142_dp * cst%m_u !Si+
    var%m(61) = 23.99694000_dp * cst%m_u !NaH
    var%m(62) = 38.02352000_dp * cst%m_u !NaCH3
    var%m(63) = 25.00433142_dp * cst%m_u !NaH2+
    var%m(64) = 39.03091142_dp * cst%m_u !NaCH4+
    var%m(65) = 49.02573142_dp * cst%m_u !NaC2H2+
    var%m(66) = 51.04161142_dp * cst%m_u !NaC2H4+
    var%m(67) = 19.02267142_dp * cst%m_u !H3O+
    var%m(68) = 18.01528000_dp * cst%m_u !H2O
    var%m(69) = 56.85294000_dp * cst%m_u !FeH
    var%m(70) = 56.85239142_dp * cst%m_u !FeH+
    var%m(71) = 57.86033142_dp * cst%m_u !FeH2+
    var%m(72) = 71.88691142_dp * cst%m_u !FeCH4+
    var%m(73) = 81.88173142_dp * cst%m_u !FeC2H2+
    var%m(74) = 83.89761142_dp * cst%m_u !FeC2H4+
    var%m(75) = 25.31294000_dp * cst%m_u !MgH
    var%m(76) = 25.31239142_dp * cst%m_u !MgH+
    var%m(77) = 26.32033142_dp * cst%m_u !MgH2+
    var%m(78) = 40.34691142_dp * cst%m_u !MgCH4+
    var%m(79) = 50.34173142_dp * cst%m_u !MgC2H2+
    var%m(80) = 52.35761142_dp * cst%m_u !MgC2H4+
    var%m(81) = 29.09344000_dp * cst%m_u !SiH
    var%m(82) = 29.09289142_dp * cst%m_u !SiH+
    var%m(83) = 30.10083142_dp * cst%m_u !SiH2+
    var%m(84) = 28.08495142_dp * cst%m_u !SiCnHm+
    var%m(85) = 45.09229142_dp * cst%m_u !SiOH+
    var%m(86) = 19.02322000_dp * cst%m_u !H3O
    var%m(87) = 41.10414000_dp * cst%m_u !SiCH
    var%m(88) = 42.11208000_dp * cst%m_u !SiCH2

    ! mass zero error
    do isp = 1, spl%nsp
      if ( var%m(isp) == 0_dp ) then 
        write(*,*) 'mass zero error'
        write(*,*) 'please check species list.'
        stop
      end  if
    end do

    ! charge
    var%q(1) = 0.0_dp * cst%q_e !H
    var%q(2) = 1.0_dp * cst%q_e !H+
    var%q(3) = -1.0_dp * cst%q_e !e-
    var%q(4) = 0.0_dp * cst%q_e !H2
    var%q(5) = 1.0_dp * cst%q_e !H2+
    var%q(6) = 0.0_dp * cst%q_e !He
    var%q(7) = 1.0_dp * cst%q_e !He+
    var%q(8) = 0.0_dp * cst%q_e !CH4
    var%q(9) = 1.0_dp * cst%q_e !CH4+
    var%q(10) = 1.0_dp * cst%q_e !CH3+
    var%q(11) = 0.0_dp * cst%q_e !C2H2
    var%q(12) = 1.0_dp * cst%q_e !C2H2+
    var%q(13) = 0.0_dp * cst%q_e !C2H4
    var%q(14) = 1.0_dp * cst%q_e !C2H4+
    var%q(15) = 1.0_dp * cst%q_e !C2H3+
    var%q(16) = 1.0_dp * cst%q_e !C2H+
    var%q(17) = 0.0_dp * cst%q_e !C2H6
    var%q(18) = 1.0_dp * cst%q_e !C2H6+
    var%q(19) = 1.0_dp * cst%q_e !C2H5+
    var%q(20) = 1.0_dp * cst%q_e !HeH+
    var%q(21) = 1.0_dp * cst%q_e !H3+
    var%q(22) = 1.0_dp * cst%q_e !C+
    var%q(23) = 0.0_dp * cst%q_e !C
    var%q(24) = 1.0_dp * cst%q_e !CH+
    var%q(25) = 1.0_dp * cst%q_e !CH2+
    var%q(26) = 0.0_dp * cst%q_e !CH
    var%q(27) = 0.0_dp * cst%q_e !CH2
    var%q(28) = 0.0_dp * cst%q_e !CH3
    var%q(29) = 1.0_dp * cst%q_e !CH5+
    var%q(30) = 1.0_dp * cst%q_e !C2+
    var%q(31) = 0.0_dp * cst%q_e !C2
    var%q(32) = 0.0_dp * cst%q_e !C2H
    var%q(33) = 0.0_dp * cst%q_e !C2H3
    var%q(34) = 0.0_dp * cst%q_e !C2H5
    var%q(35) = 1.0_dp * cst%q_e !C2H7+
    var%q(36) = 1.0_dp * cst%q_e !C3H+
    var%q(37) = 1.0_dp * cst%q_e !C3H2+
    var%q(38) = 1.0_dp * cst%q_e !C3H3+
    var%q(39) = 1.0_dp * cst%q_e !C3H4+
    var%q(40) = 1.0_dp * cst%q_e !C3H5+
    var%q(41) = 1.0_dp * cst%q_e !C3H6+
    var%q(42) = 1.0_dp * cst%q_e !C3H7+
    var%q(43) = 1.0_dp * cst%q_e !C3H8+
    var%q(44) = 1.0_dp * cst%q_e !C3H9+
    var%q(45) = 1.0_dp * cst%q_e !C4H+
    var%q(46) = 1.0_dp * cst%q_e !C4H2+
    var%q(47) = 1.0_dp * cst%q_e !C4H3+
    var%q(48) = 1.0_dp * cst%q_e !C4H5+
    var%q(49) = 1.0_dp * cst%q_e !C4H7+
    var%q(50) = 1.0_dp * cst%q_e !C4H9+
    var%q(51) = 0.0_dp * cst%q_e !H2(v>=2)
    var%q(52) = 0.0_dp * cst%q_e !H2(v>=4)
    var%q(53) = 0.0_dp * cst%q_e !Na
    var%q(54) = 1.0_dp * cst%q_e !Na+
    var%q(55) = 0.0_dp * cst%q_e !Fe
    var%q(56) = 1.0_dp * cst%q_e !Fe+
    var%q(57) = 0.0_dp * cst%q_e !Mg
    var%q(58) = 1.0_dp * cst%q_e !Mg+
    var%q(59) = 0.0_dp * cst%q_e !Si
    var%q(60) = 1.0_dp * cst%q_e !Si+
    var%q(61) = 0.0_dp * cst%q_e !NaH
    var%q(62) = 0.0_dp * cst%q_e !NaCH3
    var%q(63) = 1.0_dp * cst%q_e !NaH2+
    var%q(64) = 1.0_dp * cst%q_e !NaCH4+
    var%q(65) = 1.0_dp * cst%q_e !NaC2H2+
    var%q(66) = 1.0_dp * cst%q_e !NaC2H4+
    var%q(67) = 1.0_dp * cst%q_e !H3O+
    var%q(68) = 0.0_dp * cst%q_e !H2O
    var%q(69) = 0.0_dp * cst%q_e !FeH
    var%q(70) = 1.0_dp * cst%q_e !FeH+
    var%q(71) = 1.0_dp * cst%q_e !FeH2+
    var%q(72) = 1.0_dp * cst%q_e !FeCH4+
    var%q(73) = 1.0_dp * cst%q_e !FeC2H2+
    var%q(74) = 1.0_dp * cst%q_e !FeC2H4+
    var%q(75) = 0.0_dp * cst%q_e !MgH
    var%q(76) = 1.0_dp * cst%q_e !MgH+
    var%q(77) = 1.0_dp * cst%q_e !MgH2+
    var%q(78) = 1.0_dp * cst%q_e !MgCH4+
    var%q(79) = 1.0_dp * cst%q_e !MgC2H2+
    var%q(80) = 1.0_dp * cst%q_e !MgC2H4+
    var%q(81) = 0.0_dp * cst%q_e !SiH
    var%q(82) = 1.0_dp * cst%q_e !SiH+
    var%q(83) = 1.0_dp * cst%q_e !SiH2+
    var%q(84) = 1.0_dp * cst%q_e !SiCnHm+
    var%q(85) = 1.0_dp * cst%q_e !SiOH+
    var%q(86) = 0.0_dp * cst%q_e !H3O
    var%q(87) = 0.0_dp * cst%q_e !SiCH
    var%q(88) = 0.0_dp * cst%q_e !SiCH2

    ! read P, L, J list
    open(11, file = './Jupiter/metal_dev10/input/PLJ_list/Production_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Prod_list(isp,ich), ich = 1, spl%nch_P)
      end do
    close(11)

    open(11, file = './Jupiter/metal_dev10/input/PLJ_list/Loss_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Loss_list(isp,ich), ich = 1, spl%nch_L)
      end do
    close(11)

    open(11, file = './Jupiter/metal_dev10/input/PLJ_list/Jacobian_list.dat', status = 'unknown' )
      do i = 1, spl%n_Jlist
        read(11,*) (spl%Jmtx_list(i,j), j = 1, spl%nch_J)
      end do
    close(11)

    open(11, file = './Jupiter/metal_dev10/input/PLJ_list/reactant_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%reactant_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Jupiter/metal_dev10/input/PLJ_list/product_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%product_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Jupiter/metal_dev10/input/PLJ_list/reaction_type_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) spl%reaction_type_list(ich)
      end do
    close(11)

    open(11, file = './Jupiter/metal_dev10/input/PLJ_list/rate_rpn_token.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_token(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Jupiter/metal_dev10/input/PLJ_list/rate_rpn_label.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_label(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Jupiter/metal_dev10/input/PLJ_list/rate_cases.dat', status = 'unknown' )
      read(11,*) (spl%rate_cases(ich), ich = 1, spl%nch)
    close(11)

    open(11, file = './Jupiter/metal_dev10/input/PLJ_list/T_range.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%T_range(ich,i,j), j = 1, 3)
        end do
      end do
    close(11)

    ! reaction type characters
    var%nspecial = 0
    do ich = 1, spl%nch
      if (      spl%reaction_type_list(ich) == 1 ) then
        spl%reaction_type_char(ich) = 'photoionization'
      else if ( spl%reaction_type_list(ich) == 2 ) then
        spl%reaction_type_char(ich) = 'photodissociation'
      else if ( spl%reaction_type_list(ich) == 11 ) then
        spl%reaction_type_char(ich) = 'electron impact'
        var%nspecial = var%nspecial + 1
      else if ( spl%reaction_type_list(ich) == 12 ) then
        spl%reaction_type_char(ich) = 'proton impact'
        var%nspecial = var%nspecial + 1
      else if ( spl%reaction_type_list(ich) == 13 ) then
        spl%reaction_type_char(ich) = 'H impact'
        var%nspecial = var%nspecial + 1
      else if ( spl%reaction_type_list(ich) == 30 ) then
        spl%reaction_type_char(ich) = 'pressure_dependent_3body'
      else if ( spl%reaction_type_list(ich) == 31 ) then
        spl%reaction_type_char(ich) = 'pressure_dependent_3bodyM'
      else if ( spl%reaction_type_list(ich) == 41 ) then
        spl%reaction_type_char(ich) = 'Meteoroid ablation'
        var%nspecial = var%nspecial + 1
      end if
    end do

    ! input Temperature profiles

    open(11, file = './Jupiter/metal_dev10/input/Temperature/T_e.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Te(iz)
      end do
    close(11)

    open(11, file = './Jupiter/metal_dev10/input/Temperature/T_i.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Ti(iz)
      end do
    close(11)

    open(11, file = './Jupiter/metal_dev10/input/Temperature/T_n.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Tn(iz)
      end do
    close(11)

    ! input density profiles
    var%ni   = 1.0e-20_dp

    isp = sp_index(spl, 'H2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/metal/input/density/H2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'He')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/metal/input/density/He.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'CH4')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/metal/input/density/CH4.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'C2H2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/metal/input/density/C2H2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'C2H4')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/metal/input/density/C2H4.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'C2H6')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/metal/input/density/C2H6.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'H2(v>=2)')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/metal/input/density/H2_v2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'H2(v>=4)')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/metal/input/density/H2_v4.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    var%ni_0 = var%ni

    ! Lower boundary condition
    var%LowerBC = 0.0_dp


    ! Upper boundary condition
    var%UpperBC = 0.0_dp


  end subroutine v__in__exe

end module v__in
