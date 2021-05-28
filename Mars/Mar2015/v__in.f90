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
    set%fin_sec = 8.0e4_dp
    set%dtime_limit = 1.0e2_dp
    set%latitude = 0.0_dp
    set%Ls = 0.0_dp
    set%nday = 3_dp
    set%scheme = 'implicit'
    set%inversion = 'Catling'
    ! directory setting
    set%dir_name = './Mars/Mar2015'
    set%fnamestable = './Mars/Mar2015/output/density/n_stable.dat'

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
    grd%nz    = 101
    allocate(grd%alt(grd%nz),grd%dalt(grd%nz))
    grd%dalt(1:101) = 2.0e3_dp ! [m]
    grd%alt(1)      = 0.0e3_dp ! [m]
    do iz = 2, grd%nz
      grd%alt(iz) = grd%alt(iz-1) + grd%dalt(iz)
    end do

    ! reactions, chemical species
    spl%nsp     = 45
    spl%nsp_i   = 31
    spl%nch     = 117
    spl%nch_P   = 49
    spl%nch_L   = 31
    spl%n_Jlist = 212
    spl%nch_J   = 33
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
    allocate(var%d_dneu_dPhi_dz_add(spl%nsp_i,grd%nz))
    allocate(var%d_dne0_dPhi_dz_add(spl%nsp_i,grd%nz))
    allocate(var%d_dnel_dPhi_dz_add(spl%nsp_i,grd%nz))
    allocate(var%barr(spl%nsp_i*grd%nz), var%xarr(spl%nsp_i*grd%nz))
    allocate(var%yarr(spl%nsp_i*grd%nz), var%dxarr(spl%nsp_i*grd%nz))
    allocate(var%Amtx(spl%nsp_i*grd%nz,2*spl%nsp_i+1))
    allocate(var%Lmtx(spl%nsp_i*grd%nz,2*spl%nsp_i+1))
    allocate(var%Umtx(spl%nsp_i*grd%nz,2*spl%nsp_i+1))

    ! species
    spl%species(1) = 'CO2+'
    spl%species(2) = 'e-'
    spl%species(3) = 'CO2'
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
    spl%species(14) = 'H'
    spl%species(15) = 'H+'
    spl%species(16) = 'H2'
    spl%species(17) = 'H2+'
    spl%species(18) = 'H2O'
    spl%species(19) = 'H2O+'
    spl%species(20) = 'OH+'
    spl%species(21) = 'OH'
    spl%species(22) = 'O+'
    spl%species(23) = 'O(1D)'
    spl%species(24) = 'O3'
    spl%species(25) = 'HO2'
    spl%species(26) = 'H2O2'
    spl%species(27) = 'H3O+'
    spl%species(28) = 'N2'
    spl%species(29) = 'NO+'
    spl%species(30) = 'N'
    spl%species(31) = 'HCO+'
    spl%species(32) = 'M'
    spl%species(33) = 'CO2+CO2'
    spl%species(34) = 'O2+CO2'
    spl%species(35) = 'O2+H2O'
    spl%species(36) = 'O2+O2'
    spl%species(37) = 'H3O+OH'
    spl%species(38) = 'H3O+(H2O)'
    spl%species(39) = 'H3O+(H2O)2'
    spl%species(40) = 'H3O+(H2O)3'
    spl%species(41) = 'H3O+(H2O)4'
    spl%species(42) = 'H3O+(H2O)5'
    spl%species(43) = 'H3O+(H2O)6'
    spl%species(44) = 'H3O+(H2O)7'
    spl%species(45) = 'HOCO'

    ! label_fix
    spl%label_fix(1) = 0 ! CO2+: variable
    spl%label_fix(2) = 0 ! e-: variable
    spl%label_fix(3) = 1 ! CO2: fixed
    spl%label_fix(4) = 1 ! CO+: fixed
    spl%label_fix(5) = 0 ! O: variable
    spl%label_fix(6) = 0 ! O+(4S): variable
    spl%label_fix(7) = 0 ! CO: variable
    spl%label_fix(8) = 1 ! C+: fixed
    spl%label_fix(9) = 0 ! O2: variable
    spl%label_fix(10) = 0 ! O2+: variable
    spl%label_fix(11) = 1 ! O+(2D): fixed
    spl%label_fix(12) = 1 ! O+(2P): fixed
    spl%label_fix(13) = 1 ! C: fixed
    spl%label_fix(14) = 0 ! H: variable
    spl%label_fix(15) = 0 ! H+: variable
    spl%label_fix(16) = 0 ! H2: variable
    spl%label_fix(17) = 1 ! H2+: fixed
    spl%label_fix(18) = 1 ! H2O: fixed
    spl%label_fix(19) = 0 ! H2O+: variable
    spl%label_fix(20) = 1 ! OH+: fixed
    spl%label_fix(21) = 0 ! OH: variable
    spl%label_fix(22) = 1 ! O+: fixed
    spl%label_fix(23) = 0 ! O(1D): variable
    spl%label_fix(24) = 0 ! O3: variable
    spl%label_fix(25) = 0 ! HO2: variable
    spl%label_fix(26) = 0 ! H2O2: variable
    spl%label_fix(27) = 0 ! H3O+: variable
    spl%label_fix(28) = 0 ! N2: variable
    spl%label_fix(29) = 1 ! NO+: fixed
    spl%label_fix(30) = 1 ! N: fixed
    spl%label_fix(31) = 1 ! HCO+: fixed
    spl%label_fix(32) = 1 ! M: fixed
    spl%label_fix(33) = 0 ! CO2+CO2: variable
    spl%label_fix(34) = 0 ! O2+CO2: variable
    spl%label_fix(35) = 0 ! O2+H2O: variable
    spl%label_fix(36) = 0 ! O2+O2: variable
    spl%label_fix(37) = 0 ! H3O+OH: variable
    spl%label_fix(38) = 0 ! H3O+(H2O): variable
    spl%label_fix(39) = 0 ! H3O+(H2O)2: variable
    spl%label_fix(40) = 0 ! H3O+(H2O)3: variable
    spl%label_fix(41) = 0 ! H3O+(H2O)4: variable
    spl%label_fix(42) = 0 ! H3O+(H2O)5: variable
    spl%label_fix(43) = 0 ! H3O+(H2O)6: variable
    spl%label_fix(44) = 0 ! H3O+(H2O)7: variable
    spl%label_fix(45) = 0 ! HOCO: variable

    ! all_to_var
    spl%all_to_var = 0
    spl%all_to_var(1) = 1 ! CO2+: variable
    spl%all_to_var(2) = 2 ! e-: variable
    spl%all_to_var(5) = 3 ! O: variable
    spl%all_to_var(6) = 4 ! O+(4S): variable
    spl%all_to_var(7) = 5 ! CO: variable
    spl%all_to_var(9) = 6 ! O2: variable
    spl%all_to_var(10) = 7 ! O2+: variable
    spl%all_to_var(14) = 8 ! H: variable
    spl%all_to_var(15) = 9 ! H+: variable
    spl%all_to_var(16) = 10 ! H2: variable
    spl%all_to_var(19) = 11 ! H2O+: variable
    spl%all_to_var(21) = 12 ! OH: variable
    spl%all_to_var(23) = 13 ! O(1D): variable
    spl%all_to_var(24) = 14 ! O3: variable
    spl%all_to_var(25) = 15 ! HO2: variable
    spl%all_to_var(26) = 16 ! H2O2: variable
    spl%all_to_var(27) = 17 ! H3O+: variable
    spl%all_to_var(28) = 18 ! N2: variable
    spl%all_to_var(33) = 19 ! CO2+CO2: variable
    spl%all_to_var(34) = 20 ! O2+CO2: variable
    spl%all_to_var(35) = 21 ! O2+H2O: variable
    spl%all_to_var(36) = 22 ! O2+O2: variable
    spl%all_to_var(37) = 23 ! H3O+OH: variable
    spl%all_to_var(38) = 24 ! H3O+(H2O): variable
    spl%all_to_var(39) = 25 ! H3O+(H2O)2: variable
    spl%all_to_var(40) = 26 ! H3O+(H2O)3: variable
    spl%all_to_var(41) = 27 ! H3O+(H2O)4: variable
    spl%all_to_var(42) = 28 ! H3O+(H2O)5: variable
    spl%all_to_var(43) = 29 ! H3O+(H2O)6: variable
    spl%all_to_var(44) = 30 ! H3O+(H2O)7: variable
    spl%all_to_var(45) = 31 ! HOCO: variable

    ! var_to_all
    spl%var_to_all(1) = 1 ! CO2+: variable
    spl%var_to_all(2) = 2 ! e-: variable
    spl%var_to_all(3) = 5 ! O: variable
    spl%var_to_all(4) = 6 ! O+(4S): variable
    spl%var_to_all(5) = 7 ! CO: variable
    spl%var_to_all(6) = 9 ! O2: variable
    spl%var_to_all(7) = 10 ! O2+: variable
    spl%var_to_all(8) = 14 ! H: variable
    spl%var_to_all(9) = 15 ! H+: variable
    spl%var_to_all(10) = 16 ! H2: variable
    spl%var_to_all(11) = 19 ! H2O+: variable
    spl%var_to_all(12) = 21 ! OH: variable
    spl%var_to_all(13) = 23 ! O(1D): variable
    spl%var_to_all(14) = 24 ! O3: variable
    spl%var_to_all(15) = 25 ! HO2: variable
    spl%var_to_all(16) = 26 ! H2O2: variable
    spl%var_to_all(17) = 27 ! H3O+: variable
    spl%var_to_all(18) = 28 ! N2: variable
    spl%var_to_all(19) = 33 ! CO2+CO2: variable
    spl%var_to_all(20) = 34 ! O2+CO2: variable
    spl%var_to_all(21) = 35 ! O2+H2O: variable
    spl%var_to_all(22) = 36 ! O2+O2: variable
    spl%var_to_all(23) = 37 ! H3O+OH: variable
    spl%var_to_all(24) = 38 ! H3O+(H2O): variable
    spl%var_to_all(25) = 39 ! H3O+(H2O)2: variable
    spl%var_to_all(26) = 40 ! H3O+(H2O)3: variable
    spl%var_to_all(27) = 41 ! H3O+(H2O)4: variable
    spl%var_to_all(28) = 42 ! H3O+(H2O)5: variable
    spl%var_to_all(29) = 43 ! H3O+(H2O)6: variable
    spl%var_to_all(30) = 44 ! H3O+(H2O)7: variable
    spl%var_to_all(31) = 45 ! HOCO: variable

    ! mass
    var%m(1) = 44.00895142_dp * cst%m_u !CO2+
    var%m(2) = 0.00054858_dp * cst%m_u !e-
    var%m(3) = 44.00950000_dp * cst%m_u !CO2
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
    var%m(14) = 1.00794000_dp * cst%m_u !H
    var%m(15) = 1.00739142_dp * cst%m_u !H+
    var%m(16) = 2.01588000_dp * cst%m_u !H2
    var%m(17) = 2.01533142_dp * cst%m_u !H2+
    var%m(18) = 18.01528000_dp * cst%m_u !H2O
    var%m(19) = 18.01473142_dp * cst%m_u !H2O+
    var%m(20) = 17.00679142_dp * cst%m_u !OH+
    var%m(21) = 17.00734000_dp * cst%m_u !OH
    var%m(22) = 15.99885142_dp * cst%m_u !O+
    var%m(23) = 15.99940000_dp * cst%m_u !O(1D)
    var%m(24) = 47.99820000_dp * cst%m_u !O3
    var%m(25) = 33.00674000_dp * cst%m_u !HO2
    var%m(26) = 34.01468000_dp * cst%m_u !H2O2
    var%m(27) = 19.02267142_dp * cst%m_u !H3O+
    var%m(28) = 28.01340000_dp * cst%m_u !N2
    var%m(29) = 30.00555142_dp * cst%m_u !NO+
    var%m(30) = 14.00670000_dp * cst%m_u !N
    var%m(31) = 29.01749142_dp * cst%m_u !HCO+
    var%m(32) = 10.00000000_dp * cst%m_u !M
    var%m(33) = 88.01845142_dp * cst%m_u !CO2+CO2
    var%m(34) = 76.00775142_dp * cst%m_u !O2+CO2
    var%m(35) = 50.01353142_dp * cst%m_u !O2+H2O
    var%m(36) = 63.99705142_dp * cst%m_u !O2+O2
    var%m(37) = 36.03001142_dp * cst%m_u !H3O+OH
    var%m(38) = 5.03915142_dp * cst%m_u !H3O+(H2O)
    var%m(39) = 55.05323142_dp * cst%m_u !H3O+(H2O)2
    var%m(40) = 73.06851142_dp * cst%m_u !H3O+(H2O)3
    var%m(41) = 91.08379142_dp * cst%m_u !H3O+(H2O)4
    var%m(42) = 109.09907142_dp * cst%m_u !H3O+(H2O)5
    var%m(43) = 127.11435142_dp * cst%m_u !H3O+(H2O)6
    var%m(44) = 145.12963142_dp * cst%m_u !H3O+(H2O)7
    var%m(45) = 45.01744000_dp * cst%m_u !HOCO

    ! mass zero error
    do isp = 1, spl%nsp
      if ( var%m(isp) == 0.0_dp ) then 
        write(*,*) 'mass zero error'
        write(*,*) 'please check species list.'
        stop
      end  if
    end do

    ! charge
    var%q(1) = 1.0_dp * cst%q_e !CO2+
    var%q(2) = -1.0_dp * cst%q_e !e-
    var%q(3) = 0.0_dp * cst%q_e !CO2
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
    var%q(14) = 0.0_dp * cst%q_e !H
    var%q(15) = 1.0_dp * cst%q_e !H+
    var%q(16) = 0.0_dp * cst%q_e !H2
    var%q(17) = 1.0_dp * cst%q_e !H2+
    var%q(18) = 0.0_dp * cst%q_e !H2O
    var%q(19) = 1.0_dp * cst%q_e !H2O+
    var%q(20) = 1.0_dp * cst%q_e !OH+
    var%q(21) = 0.0_dp * cst%q_e !OH
    var%q(22) = 1.0_dp * cst%q_e !O+
    var%q(23) = 0.0_dp * cst%q_e !O(1D)
    var%q(24) = 0.0_dp * cst%q_e !O3
    var%q(25) = 0.0_dp * cst%q_e !HO2
    var%q(26) = 0.0_dp * cst%q_e !H2O2
    var%q(27) = 1.0_dp * cst%q_e !H3O+
    var%q(28) = 0.0_dp * cst%q_e !N2
    var%q(29) = 1.0_dp * cst%q_e !NO+
    var%q(30) = 0.0_dp * cst%q_e !N
    var%q(31) = 1.0_dp * cst%q_e !HCO+
    var%q(32) = 0.0_dp * cst%q_e !M
    var%q(33) = 1.0_dp * cst%q_e !CO2+CO2
    var%q(34) = 1.0_dp * cst%q_e !O2+CO2
    var%q(35) = 1.0_dp * cst%q_e !O2+H2O
    var%q(36) = 1.0_dp * cst%q_e !O2+O2
    var%q(37) = 1.0_dp * cst%q_e !H3O+OH
    var%q(38) = 1.0_dp * cst%q_e !H3O+(H2O)
    var%q(39) = 1.0_dp * cst%q_e !H3O+(H2O)2
    var%q(40) = 1.0_dp * cst%q_e !H3O+(H2O)3
    var%q(41) = 1.0_dp * cst%q_e !H3O+(H2O)4
    var%q(42) = 1.0_dp * cst%q_e !H3O+(H2O)5
    var%q(43) = 1.0_dp * cst%q_e !H3O+(H2O)6
    var%q(44) = 1.0_dp * cst%q_e !H3O+(H2O)7
    var%q(45) = 0.0_dp * cst%q_e !HOCO

    ! read P, L, J list
    open(11, file = './Mars/Mar2015/input/PLJ_list/Production_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Prod_list(isp,ich), ich = 1, spl%nch_P)
      end do
    close(11)

    open(11, file = './Mars/Mar2015/input/PLJ_list/Loss_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Loss_list(isp,ich), ich = 1, spl%nch_L)
      end do
    close(11)

    open(11, file = './Mars/Mar2015/input/PLJ_list/Jacobian_list.dat', status = 'unknown' )
      do i = 1, spl%n_Jlist
        read(11,*) (spl%Jmtx_list(i,j), j = 1, spl%nch_J)
      end do
    close(11)

    open(11, file = './Mars/Mar2015/input/PLJ_list/reactant_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%reactant_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Mars/Mar2015/input/PLJ_list/product_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%product_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Mars/Mar2015/input/PLJ_list/reaction_type_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) spl%reaction_type_list(ich)
      end do
    close(11)

    open(11, file = './Mars/Mar2015/input/PLJ_list/rate_rpn_token.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_token(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Mars/Mar2015/input/PLJ_list/rate_rpn_label.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_label(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Mars/Mar2015/input/PLJ_list/rate_cases.dat', status = 'unknown' )
      read(11,*) (spl%rate_cases(ich), ich = 1, spl%nch)
    close(11)

    open(11, file = './Mars/Mar2015/input/PLJ_list/T_range.dat', status = 'unknown' )
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

    open(11, file = './Mars/Mar2015/input/Temperature/T_e.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Te(iz)
      end do
    close(11)

    open(11, file = './Mars/Mar2015/input/Temperature/T_i.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Ti(iz)
      end do
    close(11)

    open(11, file = './Mars/Mar2015/input/Temperature/T_n.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Tn(iz)
      end do
    close(11)

    ! input density profiles
    var%ni   = 1.0e-20_dp

    isp = sp_index(spl, 'CO2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/CO2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'O')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/O.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'CO')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/CO.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'O2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/O2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'H')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/H.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'H2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/H2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'H2O')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/H2O.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'OH')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/OH.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'O(1D)')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/O(1D).dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'O3')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/O3.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'HO2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/HO2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'H2O2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/H2O2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'N2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/OHtest/input/density_3day/N2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    var%ni_0 = var%ni

    ! Lower boundary condition
    var%LowerBC = 0.0_dp

    isp = sp_index(spl, 'CO2')
    var%LowerBC(isp,1) = 1.0_dp
    var%LowerBC(isp,2) = 2.1e23_dp


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
