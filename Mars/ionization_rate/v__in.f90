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
    set%nstep = 1
    set%fin_sec = 1.0e6_dp
    set%dtime_limit = 1.0e3_dp
    set%latitude = 0.0_dp
    set%Ls = 60.0_dp
    set%nday = 3_dp
    set%scheme = 'implicit'
    set%inversion = 'Catling'
    ! directory setting
    set%dir_name = './Mars/ionization_rate'
    set%fnamestable = './Mars/ionization_rate/output/density/n_stable.dat'

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
    grd%nz    = 16
    allocate(grd%alt(grd%nz),grd%dalt(grd%nz))
    grd%dalt(1:16) = 5.0e3_dp ! [m]
    grd%alt(1)      = 122.5e3_dp ! [m]
    do iz = 2, grd%nz
      grd%alt(iz) = grd%alt(iz-1) + grd%dalt(iz)
    end do

    ! reactions, chemical species
    spl%nsp     = 15
    spl%nsp_i   = 7
    spl%nch     = 14
    spl%nch_P   = 11
    spl%nch_L   = 7
    spl%n_Jlist = 22
    spl%nch_J   = 9
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
    allocate(var%d_dneu_dPhi_dz_add(spl%nsp_i,grd%nz))
    allocate(var%d_dne0_dPhi_dz_add(spl%nsp_i,grd%nz))
    allocate(var%d_dnel_dPhi_dz_add(spl%nsp_i,grd%nz))
    allocate(var%barr(spl%nsp_i*grd%nz), var%xarr(spl%nsp_i*grd%nz))
    allocate(var%yarr(spl%nsp_i*grd%nz), var%dxarr(spl%nsp_i*grd%nz))
    allocate(var%Amtx(spl%nsp_i*grd%nz,2*spl%nsp_i+1))
    allocate(var%Lmtx(spl%nsp_i*grd%nz,2*spl%nsp_i+1))
    allocate(var%Umtx(spl%nsp_i*grd%nz,2*spl%nsp_i+1))

    ! species
    spl%species(1) = 'CO2'
    spl%species(2) = 'CO2+'
    spl%species(3) = 'e-'
    spl%species(4) = 'O+(4S)'
    spl%species(5) = 'CO'
    spl%species(6) = 'O'
    spl%species(7) = 'O+(2D)'
    spl%species(8) = 'O+(2P)'
    spl%species(9) = 'O2+'
    spl%species(10) = 'N2'
    spl%species(11) = 'NO+'
    spl%species(12) = 'N'
    spl%species(13) = 'H'
    spl%species(14) = 'H+'
    spl%species(15) = 'HCO+'

    ! label_fix
    spl%label_fix(1) = 1 ! CO2: fixed
    spl%label_fix(2) = 0 ! CO2+: variable
    spl%label_fix(3) = 0 ! e-: variable
    spl%label_fix(4) = 0 ! O+(4S): variable
    spl%label_fix(5) = 1 ! CO: fixed
    spl%label_fix(6) = 1 ! O: fixed
    spl%label_fix(7) = 1 ! O+(2D): fixed
    spl%label_fix(8) = 1 ! O+(2P): fixed
    spl%label_fix(9) = 0 ! O2+: variable
    spl%label_fix(10) = 0 ! N2: variable
    spl%label_fix(11) = 1 ! NO+: fixed
    spl%label_fix(12) = 1 ! N: fixed
    spl%label_fix(13) = 0 ! H: variable
    spl%label_fix(14) = 0 ! H+: variable
    spl%label_fix(15) = 1 ! HCO+: fixed

    ! all_to_var
    spl%all_to_var = 0
    spl%all_to_var(2) = 1 ! CO2+: variable
    spl%all_to_var(3) = 2 ! e-: variable
    spl%all_to_var(4) = 3 ! O+(4S): variable
    spl%all_to_var(9) = 4 ! O2+: variable
    spl%all_to_var(10) = 5 ! N2: variable
    spl%all_to_var(13) = 6 ! H: variable
    spl%all_to_var(14) = 7 ! H+: variable

    ! var_to_all
    spl%var_to_all(1) = 2 ! CO2+: variable
    spl%var_to_all(2) = 3 ! e-: variable
    spl%var_to_all(3) = 4 ! O+(4S): variable
    spl%var_to_all(4) = 9 ! O2+: variable
    spl%var_to_all(5) = 10 ! N2: variable
    spl%var_to_all(6) = 13 ! H: variable
    spl%var_to_all(7) = 14 ! H+: variable

    ! mass
    var%m(1) = 44.00950000_dp * cst%m_u !CO2
    var%m(2) = 44.00895142_dp * cst%m_u !CO2+
    var%m(3) = 0.00054858_dp * cst%m_u !e-
    var%m(4) = 15.99885142_dp * cst%m_u !O+(4S)
    var%m(5) = 28.01010000_dp * cst%m_u !CO
    var%m(6) = 15.99940000_dp * cst%m_u !O
    var%m(7) = 15.99885142_dp * cst%m_u !O+(2D)
    var%m(8) = 15.99885142_dp * cst%m_u !O+(2P)
    var%m(9) = 31.99825142_dp * cst%m_u !O2+
    var%m(10) = 28.01340000_dp * cst%m_u !N2
    var%m(11) = 30.00555142_dp * cst%m_u !NO+
    var%m(12) = 14.00670000_dp * cst%m_u !N
    var%m(13) = 1.00794000_dp * cst%m_u !H
    var%m(14) = 1.00739142_dp * cst%m_u !H+
    var%m(15) = 29.01749142_dp * cst%m_u !HCO+

    ! mass zero error
    do isp = 1, spl%nsp
      if ( var%m(isp) == 0.0_dp ) then 
        write(*,*) 'mass zero error'
        write(*,*) 'please check species list.'
        stop
      end  if
    end do

    ! charge
    var%q(1) = 0.0_dp * cst%q_e !CO2
    var%q(2) = 1.0_dp * cst%q_e !CO2+
    var%q(3) = -1.0_dp * cst%q_e !e-
    var%q(4) = 1.0_dp * cst%q_e !O+(4S)
    var%q(5) = 0.0_dp * cst%q_e !CO
    var%q(6) = 0.0_dp * cst%q_e !O
    var%q(7) = 1.0_dp * cst%q_e !O+(2D)
    var%q(8) = 1.0_dp * cst%q_e !O+(2P)
    var%q(9) = 1.0_dp * cst%q_e !O2+
    var%q(10) = 0.0_dp * cst%q_e !N2
    var%q(11) = 1.0_dp * cst%q_e !NO+
    var%q(12) = 0.0_dp * cst%q_e !N
    var%q(13) = 0.0_dp * cst%q_e !H
    var%q(14) = 1.0_dp * cst%q_e !H+
    var%q(15) = 1.0_dp * cst%q_e !HCO+

    ! read P, L, J list
    open(11, file = './Mars/ionization_rate/input/PLJ_list/Production_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Prod_list(isp,ich), ich = 1, spl%nch_P)
      end do
    close(11)

    open(11, file = './Mars/ionization_rate/input/PLJ_list/Loss_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Loss_list(isp,ich), ich = 1, spl%nch_L)
      end do
    close(11)

    open(11, file = './Mars/ionization_rate/input/PLJ_list/Jacobian_list.dat', status = 'unknown' )
      do i = 1, spl%n_Jlist
        read(11,*) (spl%Jmtx_list(i,j), j = 1, spl%nch_J)
      end do
    close(11)

    open(11, file = './Mars/ionization_rate/input/PLJ_list/reactant_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%reactant_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Mars/ionization_rate/input/PLJ_list/product_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%product_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Mars/ionization_rate/input/PLJ_list/reaction_type_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) spl%reaction_type_list(ich)
      end do
    close(11)

    open(11, file = './Mars/ionization_rate/input/PLJ_list/rate_rpn_token.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_token(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Mars/ionization_rate/input/PLJ_list/rate_rpn_label.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_label(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Mars/ionization_rate/input/PLJ_list/rate_cases.dat', status = 'unknown' )
      read(11,*) (spl%rate_cases(ich), ich = 1, spl%nch)
    close(11)

    open(11, file = './Mars/ionization_rate/input/PLJ_list/T_range.dat', status = 'unknown' )
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

    open(11, file = './Mars/ionization_rate/input/Temperature/T_e.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Te(iz)
      end do
    close(11)

    open(11, file = './Mars/ionization_rate/input/Temperature/T_i.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Ti(iz)
      end do
    close(11)

    open(11, file = './Mars/ionization_rate/input/Temperature/T_n.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Tn(iz)
      end do
    close(11)

    ! input density profiles
    var%ni   = 1.0e-20_dp

    isp = sp_index(spl, 'CO2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/ionization_rate/input/density/Ls60/CO2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'O')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/ionization_rate/input/density/Ls60/O.dat', status = 'unknown' )
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

    isp = sp_index(spl, 'O')
    var%UpperBC(isp,1) = 2.0_dp
    var%UpperBC(isp,2) = 1.2e12_dp

    isp = sp_index(spl, 'H')
    var%UpperBC(isp,1) = 10.0_dp ! Jeans escape

    isp = sp_index(spl, 'H2')
    var%UpperBC(isp,1) = 10.0_dp ! Jeans escape


  end subroutine v__in__exe

end module v__in
