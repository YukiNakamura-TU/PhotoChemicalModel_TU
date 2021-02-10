! p__photochem_opticaldepth.f90
!
!
!
!

module p__photochem_opticaldepth
  use v__tdec,        only : spl_, var_, grd_, cst_, xct_, flx_, set_
  use p__EUVAC,       only : p__EUVAC_cross_section
  use p__UV,          only : p__UV_cross_section_exe, p__UV_EUV_xct_convergence_exe

  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private
  public :: p__photochem_opticaldepth__exe

contains


  subroutine p__photochem_opticaldepth__exe(spl, cst, flx, grd, set, & ! in
    &                                       xct, var                 ) ! inout
    type(cst_),           intent(in)    :: cst
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(flx_),           intent(in)    :: flx
    type(set_),           intent(in)    :: set
    type(xct_),           intent(inout) :: xct
    type(var_),           intent(inout) :: var

    integer isp, jsp, ksp, ich, jch, iz, iwl
    integer i, j, k
    real(dp) tau_factor
    real(dp) tmp, tmp1, tmp2, tmpzarr(grd%nz)
    character(len=256) fname

    ! for solar zenith angle near and greater than 90deg [Smith et al., 1972]
    real(dp) Hz, yz, Xz, chiz, Chfunc(grd%nz)

    !-------------------------------------------------------------------------------------------------
    !
    !                                    Photochemical calculation
    !
    !-------------------------------------------------------------------------------------------------

    ! column density -------------------------------------------------------------------------------
      var%clm_ni = 0.0_dp
      do isp = 1, spl%nsp
        var%clm_ni(isp,grd%nz) = var%ni(isp,grd%nz)*grd%dalt(grd%nz)
        do iz = grd%nz, 2, -1
          var%clm_ni(isp,iz-1) = var%clm_ni(isp,iz) &
            &                   + (var%ni(isp,iz-1) + var%ni(isp,iz)) &
            &                   / 2.0_dp * grd%dalt(iz-1)
        end do
      end do

    ! 'M' density : sum of all ---------------------------------------------------------
      do iz = 1, grd%nz
        do isp = 1, spl%nsp
          if (spl%species(isp) == 'M') then
            var%ni(isp,iz) = 0.0_dp
            do jsp = 1, spl%nsp
              var%ni(isp,iz) &
                &  = var%ni(isp,iz) + var%ni(jsp,iz)
            end do
          end if
        end do
      end do

    ! total density, mean mass ---------------------------------------------------------
      do iz = 1, grd%nz
        tmp1 = 0.0_dp
        tmp2 = 0.0_dp
        do isp = 1, spl%nsp
          if (spl%species(isp) /= 'M' .and. spl%species(isp) /= 'products' .and. spl%species(isp) /= 'hv') then
            tmp1 = tmp1 + var%ni(isp,iz) * var%m(isp)
            tmp2 = tmp2 + var%ni(isp,iz)
          end if
        end do
        var%m_mean(iz) = tmp1 / tmp2
        var%n_tot(iz) = tmp2
      end do

    ! Ch*(X,chi) function in Smith et al., 1972: treatment of solar zenith angle for terminator ---
    ! currently, sza <= 1.36  >> Chfunc = 1/cos(sza)
      do iz = 1, grd%nz
        chiz = grd%sza(grd%ix, grd%iy)
        Hz   = cst%k_B * var%Tn(iz) / var%m_mean(iz) / cst%g
        Xz   = (cst%R + grd%alt(iz)) / Hz
        yz   = dsqrt(0.5_dp * Xz) * dcos(chiz)
        if (chiz <= 1.36_dp) then
          Chfunc(iz) = 1_dp/dcos(grd%sza(grd%ix,grd%iy))
        else if (chiz > 1.36_dp .and. chiz <= cst%pi) then
          Chfunc(iz) = dsqrt(cst%pi/2.0_dp*Xz) * dexp(yz*yz) * derfc(yz)
        else if (chiz > cst%pi) then
          Chfunc(iz) = 1.0e2_dp
        end if
        !print *, Chfunc(iz), grd%sza_xact(grd%ix,grd%iy)*180.0_dp/cst%pi
      end do

      ! avoid '0' ionization source at night: 1/100 of sza 1.366 [rad] ~ 77 [deg] value
      tmp  = grd%sza_xact(grd%ix,grd%iy)
      tmp1 = 1.36_dp
      tmp2 = cst%pi
      if ( tmp < tmp1 ) then
        tau_factor = 1.0_dp
      else if ( tmp >= tmp1 .and. tmp < tmp2 ) then
        tau_factor  = 0.1d0*dexp(1.7d0*datan(35.0d0*(11.0d0*cst%pi/24.0d0-grd%sza(grd%ix,grd%iy))))
      else if ( tmp > tmp2 ) then
        tau_factor  = 0.01_dp
      end if

    ! solar flux at each altitude ------------------------------------------------------------------
      var%tau_EUV = 0.0_dp
      var%tau_UV  = 0.0_dp
      xct%type       = 'absorption'
      xct%sigma_a_EUV = 0.0_dp
      xct%sigma_a_UV  = 0.0_dp

      call p__EUVAC_cross_section(spl, & ! in
        &                         xct  ) ! inout
      call p__UV_cross_section_exe(var%Tn, grd%nz, spl, flx, & ! in
        &                          xct                       ) ! inout

      do isp = 1, spl%nsp
        
        do iz = 1, grd%nz
          do iwl = 1, flx%nwl_EUV
            var%tau_EUV(iwl,iz) = var%tau_EUV(iwl,iz) &
              &                 + var%clm_ni(isp,iz) * xct%sigma_a_EUV(iwl,isp) &
              &                 * Chfunc(iz)
            if ( var%tau_EUV(iwl,iz) > 100.0_dp ) then
              var%tau_EUV(iwl,iz) = 100.0_dp
            end if
          end do
        end do

        do iz = 1, grd%nz
          do iwl = 1, flx%nwl_UV
            var%tau_UV(iwl,iz) = var%tau_UV(iwl,iz) &
              &                + var%clm_ni(isp,iz) * xct%sigma_a_UV(iwl,iz,isp) &
              &                * Chfunc(iz)
          end do
        end do

      end do ! isp

      ! cross section convergence
      !   EUVAC cross sections are added to UV cross sections
      call p__UV_EUV_xct_convergence_exe(grd%nz, spl, var, flx, & ! in
        &                                xct                    ) ! inout
      do iz = 1, grd%nz
        do iwl = 1, 105
          do isp = 1, spl%nsp
            var%tau_UV(iwl,iz) = var%tau_UV(iwl,iz) &
              &                + var%clm_ni(isp,iz) * xct%sigma_a_UV_EUV(iwl,isp) &
              &                * Chfunc(iz)
          end do
        end do
      end do

      !fname = './'//trim(ADJUSTL(set%dir_name))//'/uvxct.dat'
      !open(11, file = fname, status = 'unknown' )
      !  do iwl = 1, flx%nwl_UV 
      !    write(11,*) flx%lambda_UV(iwl), &
      !      & (xct%sigma_a_UV(iwl,1,isp)+xct%sigma_a_UV_EUV(iwl,isp), isp = 1, spl%nsp)
      !  end do
      !close(11)
      !stop

      do iz  = 1, grd%nz
        do iwl = 1, flx%nwl_EUV
          var%I_EUV(iwl,iz) = flx%solar_EUV(iwl) &
            &               * dexp( -var%tau_EUV(iwl,iz) ) * flx%mode_factor !* tau_factor
        end do
      end do

      do iz  = 1, grd%nz
        do iwl = 1, flx%nwl_UV
          var%I_UV(iwl,iz) = flx%solar_UV(iwl) &
          &                * dexp( -var%tau_UV(iwl,iz) ) * flx%mode_factor !* tau_factor
        end do
      end do
      !stop

  end subroutine p__photochem_opticaldepth__exe


end module p__photochem_opticaldepth
