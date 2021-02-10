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
    spl%planet = 'Mars'

    ! Calculation settings
    set%mode = '1D'
    set%nstep = 10000
    set%fin_sec = 3.0e5_dp
    set%dtime_limit = 1.0e5_dp
    set%latitude = 0.0_dp
    set%Ls = 0.0_dp
    set%nday = 3
    set%scheme = 'implicit'
    set%inversion = 'Catling'

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

    ! directory setting
    set%dir_name = './Mars/ionosphere'

    ! grid setting
    grd%nx    = set%nx
    grd%ny    = set%ny
    grd%nz    = 101
    allocate(grd%alt(grd%nz),grd%dalt(grd%nz))
    grd%dalt(1:101) = 2.0e3_dp ! [m]
    grd%alt(1)      = 0.0e3_dp ! [m]
    do iz = 2, grd%nz
      grd%alt(iz) = grd%alt(iz-1) + grd%dalt(iz)
    end do

    ! reactions, chemical species
    spl%nsp     = 55
    spl%nsp_i   = 39
    spl%nch     = 285
    spl%nch_P   = 129
    spl%nch_L   = 85
    spl%n_Jlist = 561
    spl%nch_J   = 97
    spl%nrpn    = 13

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
    spl%species(1) = 'CO2'
    spl%species(2) = 'CO2+'
    spl%species(3) = 'e-'
    spl%species(4) = 'CO+'
    spl%species(5) = 'O'
    spl%species(6) = 'O+(4S)'
    spl%species(7) = 'CO'
    spl%species(8) = 'C+'
    spl%species(9) = 'O2'
    spl%species(10) = 'O2+'
    spl%species(11) = 'O+(2D)'
    spl%species(12) = 'O+(2P)'
    spl%species(13) = 'C'
    spl%species(14) = 'N2'
    spl%species(15) = 'N2+'
    spl%species(16) = 'N+'
    spl%species(17) = 'N'
    spl%species(18) = 'H'
    spl%species(19) = 'H+'
    spl%species(20) = 'H2'
    spl%species(21) = 'H2+'
    spl%species(22) = 'He'
    spl%species(23) = 'He+'
    spl%species(24) = 'Ar'
    spl%species(25) = 'Ar+'
    spl%species(26) = 'O(1D)'
    spl%species(27) = 'H2O'
    spl%species(28) = 'OH'
    spl%species(29) = 'O3'
    spl%species(30) = 'HO2'
    spl%species(31) = 'H2O2'
    spl%species(32) = 'H3+'
    spl%species(33) = 'NO'
    spl%species(34) = 'NO+'
    spl%species(35) = 'N(2D)'
    spl%species(36) = 'HCO2+'
    spl%species(37) = 'HCO+'
    spl%species(38) = 'N2H+'
    spl%species(39) = 'N(4S)'
    spl%species(40) = 'O(3P)'
    spl%species(41) = 'OH+'
    spl%species(42) = 'CH+'
    spl%species(43) = 'ArH+'
    spl%species(44) = 'NH'
    spl%species(45) = 'N(2P)'
    spl%species(46) = 'O(1S)'
    spl%species(47) = 'CN'
    spl%species(48) = 'C(1D)'
    spl%species(49) = 'HNO+'
    spl%species(50) = 'HNO'
    spl%species(51) = 'H2O+'
    spl%species(52) = 'H3O+'
    spl%species(53) = 'HO2+'
    spl%species(54) = 'HOC+'
    spl%species(55) = 'M'

    ! label_fix
    spl%label_fix(1) = 1 ! CO2: fixed
    spl%label_fix(2) = 0 ! CO2+: variable
    spl%label_fix(3) = 0 ! e-: variable
    spl%label_fix(4) = 0 ! CO+: variable
    spl%label_fix(5) = 0 ! O: variable
    spl%label_fix(6) = 0 ! O+(4S): variable
    spl%label_fix(7) = 0 ! CO: variable
    spl%label_fix(8) = 0 ! C+: variable
    spl%label_fix(9) = 0 ! O2: variable
    spl%label_fix(10) = 0 ! O2+: variable
    spl%label_fix(11) = 0 ! O+(2D): variable
    spl%label_fix(12) = 0 ! O+(2P): variable
    spl%label_fix(13) = 0 ! C: variable
    spl%label_fix(14) = 1 ! N2: fixed
    spl%label_fix(15) = 0 ! N2+: variable
    spl%label_fix(16) = 0 ! N+: variable
    spl%label_fix(17) = 0 ! N: variable
    spl%label_fix(18) = 0 ! H: variable
    spl%label_fix(19) = 0 ! H+: variable
    spl%label_fix(20) = 0 ! H2: variable
    spl%label_fix(21) = 0 ! H2+: variable
    spl%label_fix(22) = 0 ! He: variable
    spl%label_fix(23) = 0 ! He+: variable
    spl%label_fix(24) = 0 ! Ar: variable
    spl%label_fix(25) = 0 ! Ar+: variable
    spl%label_fix(26) = 0 ! O(1D): variable
    spl%label_fix(27) = 1 ! H2O: fixed
    spl%label_fix(28) = 0 ! OH: variable
    spl%label_fix(29) = 0 ! O3: variable
    spl%label_fix(30) = 0 ! HO2: variable
    spl%label_fix(31) = 0 ! H2O2: variable
    spl%label_fix(32) = 0 ! H3+: variable
    spl%label_fix(33) = 0 ! NO: variable
    spl%label_fix(34) = 0 ! NO+: variable
    spl%label_fix(35) = 0 ! N(2D): variable
    spl%label_fix(36) = 0 ! HCO2+: variable
    spl%label_fix(37) = 0 ! HCO+: variable
    spl%label_fix(38) = 0 ! N2H+: variable
    spl%label_fix(39) = 1 ! N(4S): fixed
    spl%label_fix(40) = 1 ! O(3P): fixed
    spl%label_fix(41) = 0 ! OH+: variable
    spl%label_fix(42) = 1 ! CH+: fixed
    spl%label_fix(43) = 1 ! ArH+: fixed
    spl%label_fix(44) = 1 ! NH: fixed
    spl%label_fix(45) = 0 ! N(2P): variable
    spl%label_fix(46) = 0 ! O(1S): variable
    spl%label_fix(47) = 1 ! CN: fixed
    spl%label_fix(48) = 1 ! C(1D): fixed
    spl%label_fix(49) = 0 ! HNO+: variable
    spl%label_fix(50) = 1 ! HNO: fixed
    spl%label_fix(51) = 1 ! H2O+: fixed
    spl%label_fix(52) = 1 ! H3O+: fixed
    spl%label_fix(53) = 1 ! HO2+: fixed
    spl%label_fix(54) = 1 ! HOC+: fixed
    spl%label_fix(55) = 1 ! M: fixed

    ! all_to_var
    spl%all_to_var = 0
    spl%all_to_var(2) = 1 ! CO2+: variable
    spl%all_to_var(3) = 2 ! e-: variable
    spl%all_to_var(4) = 3 ! CO+: variable
    spl%all_to_var(5) = 4 ! O: variable
    spl%all_to_var(6) = 5 ! O+(4S): variable
    spl%all_to_var(7) = 6 ! CO: variable
    spl%all_to_var(8) = 7 ! C+: variable
    spl%all_to_var(9) = 8 ! O2: variable
    spl%all_to_var(10) = 9 ! O2+: variable
    spl%all_to_var(11) = 10 ! O+(2D): variable
    spl%all_to_var(12) = 11 ! O+(2P): variable
    spl%all_to_var(13) = 12 ! C: variable
    spl%all_to_var(15) = 13 ! N2+: variable
    spl%all_to_var(16) = 14 ! N+: variable
    spl%all_to_var(17) = 15 ! N: variable
    spl%all_to_var(18) = 16 ! H: variable
    spl%all_to_var(19) = 17 ! H+: variable
    spl%all_to_var(20) = 18 ! H2: variable
    spl%all_to_var(21) = 19 ! H2+: variable
    spl%all_to_var(22) = 20 ! He: variable
    spl%all_to_var(23) = 21 ! He+: variable
    spl%all_to_var(24) = 22 ! Ar: variable
    spl%all_to_var(25) = 23 ! Ar+: variable
    spl%all_to_var(26) = 24 ! O(1D): variable
    spl%all_to_var(28) = 25 ! OH: variable
    spl%all_to_var(29) = 26 ! O3: variable
    spl%all_to_var(30) = 27 ! HO2: variable
    spl%all_to_var(31) = 28 ! H2O2: variable
    spl%all_to_var(32) = 29 ! H3+: variable
    spl%all_to_var(33) = 30 ! NO: variable
    spl%all_to_var(34) = 31 ! NO+: variable
    spl%all_to_var(35) = 32 ! N(2D): variable
    spl%all_to_var(36) = 33 ! HCO2+: variable
    spl%all_to_var(37) = 34 ! HCO+: variable
    spl%all_to_var(38) = 35 ! N2H+: variable
    spl%all_to_var(41) = 36 ! OH+: variable
    spl%all_to_var(45) = 37 ! N(2P): variable
    spl%all_to_var(46) = 38 ! O(1S): variable
    spl%all_to_var(49) = 39 ! HNO+: variable

    ! var_to_all
    spl%var_to_all(1) = 2 ! CO2+: variable
    spl%var_to_all(2) = 3 ! e-: variable
    spl%var_to_all(3) = 4 ! CO+: variable
    spl%var_to_all(4) = 5 ! O: variable
    spl%var_to_all(5) = 6 ! O+(4S): variable
    spl%var_to_all(6) = 7 ! CO: variable
    spl%var_to_all(7) = 8 ! C+: variable
    spl%var_to_all(8) = 9 ! O2: variable
    spl%var_to_all(9) = 10 ! O2+: variable
    spl%var_to_all(10) = 11 ! O+(2D): variable
    spl%var_to_all(11) = 12 ! O+(2P): variable
    spl%var_to_all(12) = 13 ! C: variable
    spl%var_to_all(13) = 15 ! N2+: variable
    spl%var_to_all(14) = 16 ! N+: variable
    spl%var_to_all(15) = 17 ! N: variable
    spl%var_to_all(16) = 18 ! H: variable
    spl%var_to_all(17) = 19 ! H+: variable
    spl%var_to_all(18) = 20 ! H2: variable
    spl%var_to_all(19) = 21 ! H2+: variable
    spl%var_to_all(20) = 22 ! He: variable
    spl%var_to_all(21) = 23 ! He+: variable
    spl%var_to_all(22) = 24 ! Ar: variable
    spl%var_to_all(23) = 25 ! Ar+: variable
    spl%var_to_all(24) = 26 ! O(1D): variable
    spl%var_to_all(25) = 28 ! OH: variable
    spl%var_to_all(26) = 29 ! O3: variable
    spl%var_to_all(27) = 30 ! HO2: variable
    spl%var_to_all(28) = 31 ! H2O2: variable
    spl%var_to_all(29) = 32 ! H3+: variable
    spl%var_to_all(30) = 33 ! NO: variable
    spl%var_to_all(31) = 34 ! NO+: variable
    spl%var_to_all(32) = 35 ! N(2D): variable
    spl%var_to_all(33) = 36 ! HCO2+: variable
    spl%var_to_all(34) = 37 ! HCO+: variable
    spl%var_to_all(35) = 38 ! N2H+: variable
    spl%var_to_all(36) = 41 ! OH+: variable
    spl%var_to_all(37) = 45 ! N(2P): variable
    spl%var_to_all(38) = 46 ! O(1S): variable
    spl%var_to_all(39) = 49 ! HNO+: variable

    ! mass
    var%m(1) = 44.00950000_dp * cst%m_u !CO2
    var%m(2) = 44.00895142_dp * cst%m_u !CO2+
    var%m(3) = 0.00054858_dp * cst%m_u !e-
    var%m(4) = 28.00955142_dp * cst%m_u !CO+
    var%m(5) = 15.99940000_dp * cst%m_u !O
    var%m(6) = 15.99885142_dp * cst%m_u !O+(4S)
    var%m(7) = 28.01010000_dp * cst%m_u !CO
    var%m(8) = 12.01015142_dp * cst%m_u !C+
    var%m(9) = 31.99880000_dp * cst%m_u !O2
    var%m(10) = 31.99825142_dp * cst%m_u !O2+
    var%m(11) = 15.99885142_dp * cst%m_u !O+(2D)
    var%m(12) = 15.99885142_dp * cst%m_u !O+(2P)
    var%m(13) = 12.01070000_dp * cst%m_u !C
    var%m(14) = 28.01340000_dp * cst%m_u !N2
    var%m(15) = 28.01285142_dp * cst%m_u !N2+
    var%m(16) = 14.00615142_dp * cst%m_u !N+
    var%m(17) = 14.00670000_dp * cst%m_u !N
    var%m(18) = 1.00794000_dp * cst%m_u !H
    var%m(19) = 1.00739142_dp * cst%m_u !H+
    var%m(20) = 2.01588000_dp * cst%m_u !H2
    var%m(21) = 2.01533142_dp * cst%m_u !H2+
    var%m(22) = 4.00260200_dp * cst%m_u !He
    var%m(23) = 4.00205342_dp * cst%m_u !He+
    var%m(24) = 39.94800000_dp * cst%m_u !Ar
    var%m(25) = 39.94745142_dp * cst%m_u !Ar+
    var%m(26) = 15.99940000_dp * cst%m_u !O(1D)
    var%m(27) = 18.01528000_dp * cst%m_u !H2O
    var%m(28) = 17.00734000_dp * cst%m_u !OH
    var%m(29) = 47.99820000_dp * cst%m_u !O3
    var%m(30) = 33.00674000_dp * cst%m_u !HO2
    var%m(31) = 34.01468000_dp * cst%m_u !H2O2
    var%m(32) = 3.02327142_dp * cst%m_u !H3+
    var%m(33) = 30.00610000_dp * cst%m_u !NO
    var%m(34) = 30.00555142_dp * cst%m_u !NO+
    var%m(35) = 14.00670000_dp * cst%m_u !N(2D)
    var%m(36) = 45.01689142_dp * cst%m_u !HCO2+
    var%m(37) = 29.01749142_dp * cst%m_u !HCO+
    var%m(38) = 29.02079142_dp * cst%m_u !N2H+
    var%m(39) = 14.00670000_dp * cst%m_u !N(4S)
    var%m(40) = 15.99940000_dp * cst%m_u !O(3P)
    var%m(41) = 17.00679142_dp * cst%m_u !OH+
    var%m(42) = 13.01809142_dp * cst%m_u !CH+
    var%m(43) = 40.95539142_dp * cst%m_u !ArH+
    var%m(44) = 15.01464000_dp * cst%m_u !NH
    var%m(45) = 14.00670000_dp * cst%m_u !N(2P)
    var%m(46) = 15.99940000_dp * cst%m_u !O(1S)
    var%m(47) = 26.01740000_dp * cst%m_u !CN
    var%m(48) = 12.01070000_dp * cst%m_u !C(1D)
    var%m(49) = 31.01349142_dp * cst%m_u !HNO+
    var%m(50) = 31.01404000_dp * cst%m_u !HNO
    var%m(51) = 18.01473142_dp * cst%m_u !H2O+
    var%m(52) = 19.02267142_dp * cst%m_u !H3O+
    var%m(53) = 33.00619142_dp * cst%m_u !HO2+
    var%m(54) = 29.01749142_dp * cst%m_u !HOC+
    var%m(55) = 10.00000000_dp * cst%m_u !M

    ! mass zero error
    do isp = 1, spl%nsp
      if ( var%m(isp) == 0_dp ) then 
        write(*,*) 'mass zero error'
        write(*,*) 'please check species list.'
        stop
      end  if
    end do

    ! charge
    var%q(1) = 0.0_dp * cst%q_e !CO2
    var%q(2) = 1.0_dp * cst%q_e !CO2+
    var%q(3) = -1.0_dp * cst%q_e !e-
    var%q(4) = 1.0_dp * cst%q_e !CO+
    var%q(5) = 0.0_dp * cst%q_e !O
    var%q(6) = 1.0_dp * cst%q_e !O+(4S)
    var%q(7) = 0.0_dp * cst%q_e !CO
    var%q(8) = 1.0_dp * cst%q_e !C+
    var%q(9) = 0.0_dp * cst%q_e !O2
    var%q(10) = 1.0_dp * cst%q_e !O2+
    var%q(11) = 1.0_dp * cst%q_e !O+(2D)
    var%q(12) = 1.0_dp * cst%q_e !O+(2P)
    var%q(13) = 0.0_dp * cst%q_e !C
    var%q(14) = 0.0_dp * cst%q_e !N2
    var%q(15) = 1.0_dp * cst%q_e !N2+
    var%q(16) = 1.0_dp * cst%q_e !N+
    var%q(17) = 0.0_dp * cst%q_e !N
    var%q(18) = 0.0_dp * cst%q_e !H
    var%q(19) = 1.0_dp * cst%q_e !H+
    var%q(20) = 0.0_dp * cst%q_e !H2
    var%q(21) = 1.0_dp * cst%q_e !H2+
    var%q(22) = 0.0_dp * cst%q_e !He
    var%q(23) = 1.0_dp * cst%q_e !He+
    var%q(24) = 0.0_dp * cst%q_e !Ar
    var%q(25) = 1.0_dp * cst%q_e !Ar+
    var%q(26) = 0.0_dp * cst%q_e !O(1D)
    var%q(27) = 0.0_dp * cst%q_e !H2O
    var%q(28) = 0.0_dp * cst%q_e !OH
    var%q(29) = 0.0_dp * cst%q_e !O3
    var%q(30) = 0.0_dp * cst%q_e !HO2
    var%q(31) = 0.0_dp * cst%q_e !H2O2
    var%q(32) = 1.0_dp * cst%q_e !H3+
    var%q(33) = 0.0_dp * cst%q_e !NO
    var%q(34) = 1.0_dp * cst%q_e !NO+
    var%q(35) = 0.0_dp * cst%q_e !N(2D)
    var%q(36) = 1.0_dp * cst%q_e !HCO2+
    var%q(37) = 1.0_dp * cst%q_e !HCO+
    var%q(38) = 1.0_dp * cst%q_e !N2H+
    var%q(39) = 0.0_dp * cst%q_e !N(4S)
    var%q(40) = 0.0_dp * cst%q_e !O(3P)
    var%q(41) = 1.0_dp * cst%q_e !OH+
    var%q(42) = 1.0_dp * cst%q_e !CH+
    var%q(43) = 1.0_dp * cst%q_e !ArH+
    var%q(44) = 0.0_dp * cst%q_e !NH
    var%q(45) = 0.0_dp * cst%q_e !N(2P)
    var%q(46) = 0.0_dp * cst%q_e !O(1S)
    var%q(47) = 0.0_dp * cst%q_e !CN
    var%q(48) = 0.0_dp * cst%q_e !C(1D)
    var%q(49) = 1.0_dp * cst%q_e !HNO+
    var%q(50) = 0.0_dp * cst%q_e !HNO
    var%q(51) = 1.0_dp * cst%q_e !H2O+
    var%q(52) = 1.0_dp * cst%q_e !H3O+
    var%q(53) = 1.0_dp * cst%q_e !HO2+
    var%q(54) = 1.0_dp * cst%q_e !HOC+
    var%q(55) = 0.0_dp * cst%q_e !M

    ! read P, L, J list
    open(11, file = './Mars/ionosphere/input/PLJ_list/Production_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Prod_list(isp,ich), ich = 1, spl%nch_P)
      end do
    close(11)

    open(11, file = './Mars/ionosphere/input/PLJ_list/Loss_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Loss_list(isp,ich), ich = 1, spl%nch_L)
      end do
    close(11)

    open(11, file = './Mars/ionosphere/input/PLJ_list/Jacobian_list.dat', status = 'unknown' )
      do i = 1, spl%n_Jlist
        read(11,*) (spl%Jmtx_list(i,j), j = 1, spl%nch_J)
      end do
    close(11)

    open(11, file = './Mars/ionosphere/input/PLJ_list/reactant_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%reactant_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Mars/ionosphere/input/PLJ_list/product_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%product_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Mars/ionosphere/input/PLJ_list/reaction_type_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) spl%reaction_type_list(ich)
      end do
    close(11)

    open(11, file = './Mars/ionosphere/input/PLJ_list/rate_rpn_token.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_token(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Mars/ionosphere/input/PLJ_list/rate_rpn_label.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_label(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Mars/ionosphere/input/PLJ_list/rate_cases.dat', status = 'unknown' )
      read(11,*) (spl%rate_cases(ich), ich = 1, spl%nch)
    close(11)

    open(11, file = './Mars/ionosphere/input/PLJ_list/T_range.dat', status = 'unknown' )
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

    open(11, file = './Mars/ionosphere/input/Temperature/T_e.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Te(iz)
      end do
    close(11)

    open(11, file = './Mars/ionosphere/input/Temperature/T_i.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Ti(iz)
      end do
    close(11)

    open(11, file = './Mars/ionosphere/input/Temperature/T_n.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Tn(iz)
      end do
    close(11)

    ! input density profiles
    var%ni   = 1.0e-20_dp

    isp = sp_index(spl, 'CO2')
    open(11, file = './Mars/ionosphere/input/density/CO2.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    isp = sp_index(spl, 'O')
    open(11, file = './Mars/ionosphere/input/density/O.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    isp = sp_index(spl, 'CO')
    open(11, file = './Mars/ionosphere/input/density/CO.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    isp = sp_index(spl, 'O2')
    open(11, file = './Mars/ionosphere/input/density/O2.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    isp = sp_index(spl, 'N2')
    open(11, file = './Mars/ionosphere/input/density/N2.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    isp = sp_index(spl, 'H')
    open(11, file = './Mars/ionosphere/input/density/H.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    isp = sp_index(spl, 'H2')
    open(11, file = './Mars/ionosphere/input/density/H2.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    isp = sp_index(spl, 'O(1D)')
    open(11, file = './Mars/ionosphere/input/density/O(1D).dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    isp = sp_index(spl, 'H2O')
    open(11, file = './Mars/ionosphere/input/density/H2O.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    isp = sp_index(spl, 'OH')
    open(11, file = './Mars/ionosphere/input/density/OH.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    isp = sp_index(spl, 'O3')
    open(11, file = './Mars/ionosphere/input/density/O3.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    isp = sp_index(spl, 'HO2')
    open(11, file = './Mars/ionosphere/input/density/HO2.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    isp = sp_index(spl, 'H2O2')
    open(11, file = './Mars/ionosphere/input/density/H2O2.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp, var%ni(isp,iz)
      end do
    close(11)

    var%ni_0 = var%ni

    ! Lower boundary condition
    var%LowerBC = 0.0_dp


    ! Upper boundary condition
    var%UpperBC = 0.0_dp

    isp = sp_index(spl, 'O')
    var%UpperBC(isp,1) = 2.0_dp
    var%UpperBC(isp,2) = 1.2e12_dp

    isp = sp_index(spl, 'H')
    var%UpperBC(isp,1) = 10.0_dp ! Jeans escape

    isp = sp_index(spl, 'H2')
    var%UpperBC(isp,1) = 10.0_dp ! Jeans escape


  end subroutine v__in__exe

end module v__in
