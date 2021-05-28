! p__photochem_transport.f90
!
!
!
!

module p__photochem_transport
  use v__tdec,                only : spl_, var_, grd_, cst_, set_
  use p__search,              only : p__search_reactant, p__search_product, sp_index
  use p__eddy_diffusion,      only : p__eddy_diffusion__exe
  use p__molecular_diffusion, only : p__molecular_diffusion__exe

  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private
  public :: p__photochem_transport__exe

contains


  subroutine p__photochem_transport__exe(spl, cst, grd, set, & ! in
    &                                    var                 ) ! inout
    implicit none
    type(cst_),           intent(in)    :: cst
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(set_),           intent(in)    :: set
    type(var_),           intent(inout) :: var

    integer isp, jsp, iz
    real(dp) tmp, tmp1, tmp2
    real(dp) Texo, lambda, vth

    ! for diffusion vertical transport
    real(dp) ne(0:grd%nz+1)
    real(dp) Ti(0:grd%nz+1), Te(0:grd%nz+1), Tn(0:grd%nz+1)
    real(dp) z0, zu, zl, zp, zm
    real(dp) m0, mu, ml, mp, mm
    real(dp) dz0, dzl, mi, niu, ni0, nil, nip, nim, neu, ne0, nel, nep, nem
    real(dp) Diu, Di0, Dil, Dip, Dim, Kiu, Ki0, Kil, Kip, Kim
    real(dp) Tiu, Ti0, Til, Teu, Te0, Tel, Tep, Tem, Tip, Tim
    real(dp) Te_Tiu, Te_Ti0, Te_Til, Te_Tip, Te_Tim
    real(dp) Peu, Pe0, Pel, Pep, Pem
    real(dp) dni_dzp, dni_dzm, dne_dzp, dne_dzm
    real(dp) dTi_dzp, dTi_dzm, dTe_dzp, dTe_dzm, dPe_dzp, dPe_dzm
    real(dp) Hiu, Hi0, Hil, Hip, Him, Hau, Ha0, Hal, Hap, Ham
    real(dp) gradPep, gradPem
    real(dp) Thermp, Thermm, alpha(spl%nsp)
    real(dp) zetap, zetam

    !-------------------------------------------------------------------------------------------------
    !
    !                                    Photochemical calculation
    !
    !-------------------------------------------------------------------------------------------------

    ! transport ---------------------------------------------------------------------------

      !--------------------------------------
      ! Eddy diffusion coefficient
      !--------------------------------------
      var%K_eddy = 0.0_dp
      call p__eddy_diffusion__exe(spl, cst, grd, & ! in
        &                         var            ) ! inout

      !--------------------------------------
      ! Molecular diffusion coefficient
      !--------------------------------------
      var%D_mol = 0.0_dp
      call p__molecular_diffusion__exe(spl, cst, grd, & ! in
        &                              var            ) ! inout

      !--------------------------------------
      ! Jeans escape
      !--------------------------------------
      do isp = 1, spl%nsp
        if (nint(var%UpperBC(isp,1)) == 10) then 
          Texo = var%Tn(grd%nz)
          lambda = cst%BigG * cst%Mplanet * var%m(isp) / ( cst%k_B * Texo * (cst%R+grd%alt(grd%nz)) )
          vth = dsqrt(2.0_dp*cst%k_B*Texo/var%m(isp))
          var%UpperBC(isp,2) = dexp(-lambda)*vth*(lambda+1.0_dp)/(2.0_dp*cst%pi**0.5_dp) 
        end if
      end do

      ! Thermal diffusion factors
      alpha = 0.0_dp

      isp = sp_index(spl, 'H')
      if (isp >= 1 .and. isp <= spl%nsp) then 
        alpha(isp) = -0.25_dp
      end if

      isp = sp_index(spl, 'H2')
      if (isp >= 1 .and. isp <= spl%nsp) then 
        alpha(isp) = -0.25_dp
      end if

      ! temperature & electron density
      do iz = 1, grd%nz
        ne(iz) = 1.0e-20_dp
        do isp = 1, spl%nsp
          if (spl%species(isp) == 'e-') then
            ne(iz) = var%ni(isp,iz)
          end if
        end do
        Tn(iz) = var%Tn(iz)
        Ti(iz) = var%Ti(iz) 
        Te(iz) = var%Te(iz)
      end do

      !--------------------------------------
      ! Vertical Diffusive Transport Equation 
      !--------------------------------------
      !  Phi = flux
      var%dPhi_dz        = 0.0_dp ! dPhi/dz
      var%d_dniu_dPhi_dz = 0.0_dp ! d/dniu * dPhi/dz
      var%d_dni0_dPhi_dz = 0.0_dp ! d/dni0 * dPhi/dz
      var%d_dnil_dPhi_dz = 0.0_dp ! d/dnil * dPhi/dz

      do isp = 1, spl%nsp

        jsp = spl%all_to_var(isp)

        if (jsp >= 1 .and. jsp <= spl%nsp_i) then 

          if (spl%label_fix(isp) == 0 .and. spl%species(isp) /= 'e-' &
            & .and. spl%species(isp) /= 'M' .and. spl%species(isp) /= 'hv') then

            do iz = 2, grd%nz-1

              ! grid definition

              ! iz+1    : p1
              ! iz+1/2  : ph ＞ dz00
              ! iz      : 00
              ! iz-1/2  : mh ＞ dzm1
              ! iz-1    : m1

              zu  = grd%alt(iz+1)
              z0  = grd%alt(iz)
              zl  = grd%alt(iz-1)
              zp  = ( zu + z0 ) / 2.0_dp
              zm  = ( z0 + zl ) / 2.0_dp

              dz0 = grd%dalt(iz)
              dzl = grd%dalt(iz-1)

              mu  = var%m_mean(iz+1)
              m0  = var%m_mean(iz)
              ml  = var%m_mean(iz-1)
              mp  = ( mu + m0 ) / 2.0_dp
              mm  = ( m0 + ml ) / 2.0_dp

              mi  = var%m(isp)

              niu = var%ni(isp,iz+1)
              ni0 = var%ni(isp,iz)
              nil = var%ni(isp,iz-1)
              nip = ( niu + ni0 ) / 2.0_dp
              nim = ( ni0 + nil ) / 2.0_dp
              neu = ne(iz+1)
              ne0 = ne(iz)
              nel = ne(iz-1)
              nep = ( neu + ne0 ) / 2.0_dp
              nem = ( ne0 + nel ) / 2.0_dp

              Diu = var%D_mol(isp,iz+1)
              Di0 = var%D_mol(isp,iz)
              Dil = var%D_mol(isp,iz-1)
              Dip = ( Diu + Di0 ) / 2.0_dp
              Dim = ( Di0 + Dil ) / 2.0_dp
              Kiu = var%K_eddy(iz+1)
              Ki0 = var%K_eddy(iz)
              Kil = var%K_eddy(iz-1)
              Kip = ( Kiu + Ki0 ) / 2.0_dp
              Kim = ( Ki0 + Kil ) / 2.0_dp

              if ( var%q(isp) /= 0.0_dp ) then
                Tiu = Ti(iz+1)
                Ti0 = Ti(iz)
                Til = Ti(iz-1)
              else if ( var%q(isp) == 0.0_dp ) then
                Tiu = Tn(iz+1)
                Ti0 = Tn(iz)
                Til = Tn(iz-1)
              end if
              Tip = ( Tiu + Ti0 ) / 2.0_dp
              Tim = ( Ti0 + Til ) / 2.0_dp

              Teu = Te(iz+1)
              Te0 = Te(iz)
              Tel = Te(iz-1)
              Tep = ( Teu + Te0 ) / 2.0_dp
              Tem = ( Te0 + Tel ) / 2.0_dp

              Te_Tiu = Te(iz+1) / Ti(iz+1)
              Te_Ti0 = Te(iz  ) / Ti(iz  )
              Te_Til = Te(iz-1) / Ti(iz-1)
              Te_Tip = ( Te_Tiu + Te_Ti0 ) / 2.0_dp
              Te_Tim = ( Te_Ti0 + Te_Til ) / 2.0_dp

              Peu = neu * cst%k_B * Teu
              Pe0 = ne0 * cst%k_B * Te0
              Pel = nel * cst%k_B * Tel
              Pep = ( Peu + Pe0 ) / 2.0_dp
              Pem = ( Pe0 + Pel ) / 2.0_dp

              dni_dzp = ( niu - ni0 ) / dz0
              dni_dzm = ( ni0 - nil ) / dzl
              dne_dzp = ( neu - ne0 ) / dz0
              dne_dzm = ( ne0 - nel ) / dzl
              dTi_dzp = ( Tiu - Ti0 ) / dz0
              dTi_dzm = ( Ti0 - Til ) / dzl
              dTe_dzp = ( Teu - Te0 ) / dz0
              dTe_dzm = ( Te0 - Tel ) / dzl
              dPe_dzp = ( Peu - Pe0 ) / dz0
              dPe_dzm = ( Pe0 - Pel ) / dzl

              if ( neu > 0.0_dp .and. ne0 > 0.0_dp .and. nel > 0.0_dp &
                  & .and. var%q(isp) /= 0.0_dp ) then
                gradPep = Tep/Tip / nep * dne_dzp + 1.0_dp / Tip * dTe_dzp
                gradPem = Tem/Tim / nem * dne_dzm + 1.0_dp / Tim * dTe_dzm
              else
                gradPep = 0.0_dp
                gradPem = 0.0_dp
              end if
              
              Thermp = 1.0_dp / Tip * dTi_dzp
              Thermm = 1.0_dp / Tim * dTi_dzm
              Hiu = cst%k_B * Tiu / mi / (cst%BigG*cst%Mplanet/(cst%R+zu)**2.0_dp)
              Hi0 = cst%k_B * Ti0 / mi / (cst%BigG*cst%Mplanet/(cst%R+z0)**2.0_dp)
              Hil = cst%k_B * Til / mi / (cst%BigG*cst%Mplanet/(cst%R+zl)**2.0_dp)
              Hau = cst%k_B * Tiu / mu / (cst%BigG*cst%Mplanet/(cst%R+zu)**2.0_dp)
              Ha0 = cst%k_B * Ti0 / m0 / (cst%BigG*cst%Mplanet/(cst%R+z0)**2.0_dp)
              Hal = cst%k_B * Til / ml / (cst%BigG*cst%Mplanet/(cst%R+zl)**2.0_dp)
              Hip = ( Hiu + Hi0 ) / 2.0_dp
              Him = ( Hi0 + Hil ) / 2.0_dp
              Hap = ( Hau + Ha0 ) / 2.0_dp
              Ham = ( Ha0 + Hal ) / 2.0_dp

              if ( var%q(isp) == 0.0_dp ) then
                zetap = Dip * ( 1.0_dp/Hip + (1.0_dp+alpha(isp))*Thermp) + Kip * (1.0_dp/Hap + Thermp)
                zetam = Dim * ( 1.0_dp/Him + (1.0_dp+alpha(isp))*Thermm) + Kim * (1.0_dp/Ham + Thermm)
              else if ( var%q(isp) > 0.0_dp ) then 
                zetap = Dip * (gradPep + 1.0_dp/Hip + (1.0_dp+alpha(isp))*Thermp) 
                zetam = Dim * (gradPem + 1.0_dp/Him + (1.0_dp+alpha(isp))*Thermm) 
              end if

              var%dPhi_dz(jsp,iz) = niu * ( - (Dip+Kip)/dz0/dz0 - zetap/2.0_dp/dz0 ) & ! iz+1 -> iz
                &                 + ni0 * (   (Dip+Kip)/dz0/dz0 - zetap/2.0_dp/dz0   & ! iz   -> iz+1
                &                           + (Dim+Kim)/dzl/dzl + zetam/2.0_dp/dzl ) & ! iz   -> iz-1
                &                 + nil * ( - (Dim+Kim)/dzl/dzl + zetam/2.0_dp/dzl )   ! iz-1 -> iz
              var%d_dniu_dPhi_dz(jsp,iz) =  - (Dip+Kip)/dz0/dz0 - zetap/2.0_dp/dz0
              var%d_dni0_dPhi_dz(jsp,iz) =    (Dip+Kip)/dz0/dz0 - zetap/2.0_dp/dz0   &
                &                           + (Dim+Kim)/dzl/dzl + zetam/2.0_dp/dzl
              var%d_dnil_dPhi_dz(jsp,iz) =  - (Dim+Kim)/dzl/dzl + zetam/2.0_dp/dzl

              if (var%q(isp) > 0.0_dp) then 
                var%d_dneu_dPhi_dz_add(jsp,iz) = - (niu+ni0)*Dip*Tep/Tip/nep/dz0/2.0_dp/dz0
                var%d_dne0_dPhi_dz_add(jsp,iz) =   (niu+ni0)*Dip*Tep/Tip/nep/dz0/2.0_dp/dz0 &
                  &                              + (ni0+nil)*Dim*Tem/Tim/nem/dzl/2.0_dp/dzl
                var%d_dnel_dPhi_dz_add(jsp,iz) = - (ni0+nil)*Dim*Tem/Tim/nem/dzl/2.0_dp/dzl
              end if

              var%Phip(jsp,iz) = - (Dip+Kip)*dni_dzp - zetap*nip
              var%Phim(jsp,iz) = - (Dim+Kim)*dni_dzm - zetam*nim

            end do

            ! Lower Boundary Condition
            zu  = grd%alt(2)
            z0  = grd%alt(1)
            zp  = ( zu + z0 ) / 2.0_dp

            dz0 = grd%dalt(1)

            mu  = var%m_mean(2)
            m0  = var%m_mean(1)
            mp  = ( mu + m0 ) / 2.0_dp

            mi  = var%m(isp)

            niu = var%ni(isp,2)
            ni0 = var%ni(isp,1)
            nip = ( niu + ni0 ) / 2.0_dp
            neu = ne(2)
            ne0 = ne(1)
            nep = ( neu + ne0 ) / 2.0_dp

            Diu = var%D_mol(isp,2)
            Di0 = var%D_mol(isp,1)
            Dip = ( Diu + Di0 ) / 2.0_dp
            Kiu = var%K_eddy(2)
            Ki0 = var%K_eddy(1)
            Kip = ( Kiu + Ki0 ) / 2.0_dp

            if ( var%q(isp) /= 0.0_dp ) then
              Tiu = Ti(2)
              Ti0 = Ti(1)
            else if ( var%q(isp) == 0.0_dp ) then
              Tiu = Tn(2)
              Ti0 = Tn(1)
            end if
            Tip = ( Tiu + Ti0 ) / 2.0_dp

            Teu = Te(2)
            Te0 = Te(1)
            Tep = ( Teu + Te0 ) / 2.0_dp

            Te_Tiu = Te(2) / Ti(2)
            Te_Ti0 = Te(1) / Ti(1)
            Te_Tip = ( Te_Tiu + Te_Ti0 ) / 2.0_dp

            Peu = neu * cst%k_B * Teu
            Pe0 = ne0 * cst%k_B * Te0
            Pep = nep * cst%k_B * Tep

            dni_dzp = ( niu - ni0 ) / dz0
            dne_dzp = ( neu - ne0 ) / dz0
            dTi_dzp = ( Tiu - Ti0 ) / dz0
            dTe_dzp = ( Teu - Te0 ) / dz0
            dPe_dzp = ( Peu - Pe0 ) / dz0

            if ( neu > 0.0_dp .and. ne0 > 0.0_dp &
                & .and. var%q(isp) /= 0.0_dp ) then
              gradPep = Tep/Tip / nep * dne_dzp + 1.0_dp / Tip * dTe_dzp
            else
              gradPep = 0.0_dp
            end if
            Thermp = 1.0_dp / Tip * dTi_dzp
            Hiu = cst%k_B * Tiu / mi / (cst%BigG*cst%Mplanet/(cst%R+zu)**2.0_dp)
            Hi0 = cst%k_B * Ti0 / mi / (cst%BigG*cst%Mplanet/(cst%R+z0)**2.0_dp)
            Hau = cst%k_B * Tiu / mu / (cst%BigG*cst%Mplanet/(cst%R+zu)**2.0_dp)
            Ha0 = cst%k_B * Ti0 / m0 / (cst%BigG*cst%Mplanet/(cst%R+z0)**2.0_dp)
            Hip = ( Hiu + Hi0 ) / 2.0_dp
            Hap = ( Hau + Ha0 ) / 2.0_dp

            if ( var%q(isp) == 0.0_dp ) then
              zetap = Dip * ( 1.0_dp/Hip + (1.0_dp+alpha(isp))*Thermp) + Kip * (1.0_dp/Hap + Thermp)
            else if ( var%q(isp) > 0.0_dp ) then 
              zetap = Dip * (gradPep + 1.0_dp/Hip + (1.0_dp+alpha(isp))*Thermp) 
            end if

            var%dPhi_dz(jsp,1) = niu * ( - (Dip+Kip)/dz0/dz0 - zetap/2.0_dp/dz0 ) & ! iz+1 -> iz
              &                + ni0 * (   (Dip+Kip)/dz0/dz0 - zetap/2.0_dp/dz0 )   ! iz   -> iz+1
            var%d_dniu_dPhi_dz(jsp,1) =  - (Dip+Kip)/dz0/dz0 - zetap/2.0_dp/dz0
            var%d_dni0_dPhi_dz(jsp,1) =    (Dip+Kip)/dz0/dz0 - zetap/2.0_dp/dz0
            var%d_dnil_dPhi_dz(jsp,1) = 0.0_dp

            if (var%q(isp) > 0.0_dp) then 
              var%d_dneu_dPhi_dz_add(jsp,1) = - (niu+ni0)*Dip*Tep/Tip/nep/dz0/2.0_dp/dz0
              var%d_dne0_dPhi_dz_add(jsp,1) =   (niu+ni0)*Dip*Tep/Tip/nep/dz0/2.0_dp/dz0 
              var%d_dnel_dPhi_dz_add(jsp,1) = 0.0_dp
            end if

            var%Phip(jsp,1) = - (Dip+Kip)*dni_dzp - zetap*nip

            if (nint(var%LowerBC(isp,1)) == 2) then 
              var%dPhi_dz(jsp,1) = var%dPhi_dz(jsp,1) &
                &                - var%LowerBC(isp,2)/dz0 ! Lower boundary: flux 
              var%Phim(jsp,1) = - var%LowerBC(isp,2)
            else if (nint(var%LowerBC(isp,1)) == 3) then 
              var%dPhi_dz(jsp,1) = var%dPhi_dz(jsp,1) &
                &                - ni0 * var%LowerBC(isp,2)/dz0 ! Lower boundary: velocity
              var%d_dni0_dPhi_dz(jsp,1) =  var%d_dni0_dPhi_dz(jsp,1) &
                &                       - var%LowerBC(isp,2)/dz0
              var%Phim(jsp,1) = - ni0 * var%LowerBC(isp,2)
            end if 

            ! Upper Boundary Condition
            z0  = grd%alt(grd%nz)
            zl  = grd%alt(grd%nz-1)
            zm  = ( z0 + zl ) / 2.0_dp

            dzl = grd%dalt(grd%nz)

            m0  = var%m_mean(grd%nz)
            ml  = var%m_mean(grd%nz-1)
            mm  = ( m0 + ml ) / 2.0_dp

            mi  = var%m(isp)

            ni0 = var%ni(isp,grd%nz)
            nil = var%ni(isp,grd%nz-1)
            nim = ( ni0 + nil ) / 2.0_dp
            ne0 = ne(grd%nz)
            nel = ne(grd%nz-1)
            nem = ( ne0 + nel ) / 2.0_dp

            Di0 = var%D_mol(isp,grd%nz)
            Dil = var%D_mol(isp,grd%nz-1)
            Dim = ( Di0 + Dil ) / 2.0_dp
            Ki0 = var%K_eddy(grd%nz)
            Kil = var%K_eddy(grd%nz-1)
            Kim = ( Ki0 + Kil ) / 2.0_dp

            if ( var%q(isp) /= 0.0_dp ) then
              Ti0 = Ti(grd%nz)
              Til = Ti(grd%nz-1)
            else if ( var%q(isp) == 0.0_dp ) then
              Ti0 = Tn(grd%nz)
              Til = Tn(grd%nz-1)
            end if
            Tim = ( Ti0 + Til ) / 2.0_dp

            Te0 = Te(grd%nz)
            Tel = Te(grd%nz-1)
            Tem = ( Te0 + Tel ) / 2.0_dp

            Te_Ti0 = Te(grd%nz) / Ti(grd%nz)
            Te_Til = Te(grd%nz-1) / Ti(grd%nz-1)
            Te_Tim = ( Te_Ti0 + Te_Til ) / 2.0_dp

            Pe0 = ne0 * cst%k_B * Te0
            Pel = nel * cst%k_B * Tel
            Pem = nem * cst%k_B * Tem

            dni_dzm = ( ni0 - nil ) / dzl
            dne_dzm = ( ne0 - nel ) / dzl
            dTi_dzm = ( Ti0 - Til ) / dzl
            dTe_dzm = ( Te0 - Tel ) / dzl
            dPe_dzm = ( Pe0 - Pel ) / dzl

            if ( ne0 > 0.0_dp .and. nel > 0.0_dp &
                & .and. var%q(isp) /= 0.0_dp ) then
              gradPem = Tem/Tim / nem * dne_dzm + 1.0_dp / Tim * dTe_dzm
            else
              gradPem = 0.0_dp
            end if
            Thermm = 1.0_dp / Tim * dTi_dzm
            Hi0 = cst%k_B * Ti0 / mi / (cst%BigG*cst%Mplanet/(cst%R+z0)**2.0_dp)
            Hil = cst%k_B * Til / mi / (cst%BigG*cst%Mplanet/(cst%R+zl)**2.0_dp)
            Ha0 = cst%k_B * Ti0 / m0 / (cst%BigG*cst%Mplanet/(cst%R+z0)**2.0_dp)
            Hal = cst%k_B * Til / ml / (cst%BigG*cst%Mplanet/(cst%R+zl)**2.0_dp)
            Him = ( Hi0 + Hil ) / 2.0_dp
            Ham = ( Ha0 + Hal ) / 2.0_dp

            if ( var%q(isp) == 0.0_dp ) then
              zetam = Dim * ( 1.0_dp/Him + (1.0_dp+alpha(isp))*Thermm) + Kim * (1.0_dp/Ham + Thermm)
            else if ( var%q(isp) > 0.0_dp ) then 
              zetam = Dim * (gradPem + 1.0_dp/Him + (1.0_dp+alpha(isp))*Thermm) 
            end if

            var%dPhi_dz(jsp,grd%nz) = ni0 * (   (Dim+Kim)/dzl/dzl + zetam/2.0_dp/dzl ) & ! iz   -> iz-1
              &                     + nil * ( - (Dim+Kim)/dzl/dzl + zetam/2.0_dp/dzl )   ! iz-1 -> iz
            var%d_dniu_dPhi_dz(jsp,grd%nz) = 0.0_dp
            var%d_dni0_dPhi_dz(jsp,grd%nz) =    (Dim+Kim)/dzl/dzl + zetam/2.0_dp/dzl
            var%d_dnil_dPhi_dz(jsp,grd%nz) =  - (Dim+Kim)/dzl/dzl + zetam/2.0_dp/dzl

            if (var%q(isp) > 0.0_dp) then 
              var%d_dneu_dPhi_dz_add(jsp,grd%nz) = 0.0_dp
              var%d_dne0_dPhi_dz_add(jsp,grd%nz) =   (ni0+nil)*Dim*Tem/Tim/nem/dzl/2.0_dp/dzl
              var%d_dnel_dPhi_dz_add(jsp,grd%nz) = - (ni0+nil)*Dim*Tem/Tim/nem/dzl/2.0_dp/dzl
            end if

            var%Phim(jsp,iz) = - (Dim+Kim)*dni_dzm - zetam*nim

            if (nint(var%UpperBC(isp,1)) == 2) then 
              var%dPhi_dz(jsp,grd%nz) = var%dPhi_dz(jsp,grd%nz) &
                &                     + var%UpperBC(isp,2)/dzl ! Upper boundary: flux 
              var%Phip(jsp,iz) = var%UpperBC(isp,2)
            else if (nint(var%UpperBC(isp,1)) == 3 .or. nint(var%UpperBC(isp,1)) == 10) then 
              var%dPhi_dz(jsp,grd%nz) = var%dPhi_dz(jsp,grd%nz) &
                &                     + ni0 * var%UpperBC(isp,2)/dzl ! Upper boundary: velocity
              var%d_dni0_dPhi_dz(jsp,grd%nz) = var%d_dni0_dPhi_dz(jsp,grd%nz) &
                &                            + var%UpperBC(isp,2)/dzl
              var%Phip(jsp,iz) = ni0 * var%UpperBC(isp,2)
            end if 

          end if

        end if
        
      end do

      ! electron : NO transport
      do isp = 1, spl%nsp

        jsp = spl%all_to_var(isp)

        if (jsp >= 1 .and. jsp <= spl%nsp_i) then 
          if (spl%species(isp) == 'e-') then
            do iz = 1, grd%nz
              var%dPhi_dz(jsp,iz)        = 0.0_dp
              var%d_dniu_dPhi_dz(jsp,iz) = 0.0_dp
              var%d_dni0_dPhi_dz(jsp,iz) = 0.0_dp
              var%d_dnil_dPhi_dz(jsp,iz) = 0.0_dp
            end do
          end if
        end if
      end do

      var%Fluxup  = 0.0_dp
      var%Fluxdwn = 0.0_dp
      do iz = 1, grd%nz
        do isp = 1, spl%nsp

          jsp = spl%all_to_var(isp)

          if (jsp >= 1 .and. jsp <= spl%nsp_i) then 
            if (var%Phip(jsp,iz) >= 0.0_dp) then
              var%Fluxup(jsp,iz) = var%Phip(jsp,iz)
            else if (var%Phip(spl%all_to_var(isp),iz) < 0.0_dp) then
              var%Fluxdwn(jsp,iz) = - var%Phip(jsp,iz)
            end if
          end if

        end do
      end do

      do isp = 1, spl%nsp

        jsp = spl%all_to_var(isp)

        if (jsp >= 1 .and. jsp <= spl%nsp_i) then 
          if (var%Phim(jsp,1) >= 0.0_dp) then
            var%Fluxup(jsp,0) = var%Phip(jsp,1)
          else if (var%Phim(jsp,1) < 0.0_dp) then
            var%Fluxdwn(jsp,0) = - var%Phim(jsp,1)
          end if
        end if

      end do

      ! NO transport case
      !var%dPhi_dz        = !0.0_dp
      !var%d_dniu_dPhi_dz = !0.0_dp
      !var%d_dni0_dPhi_dz = !0.0_dp
      !var%d_dnil_dPhi_dz = !0.0_dp


  end subroutine p__photochem_transport__exe


end module p__photochem_transport
