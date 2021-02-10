module p__UV
  use v__tdec,    only : var_, grd_, cst_, xct_, spl_, flx_
  use c__prm,     only : c__prm__ini


  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: p__UV_flux, p__UV_cross_section_dat, p__UV_cross_section_exe, &
    &       p__UV_EUV_xct_convergence_exe

contains


  !-------------------------------------------------
  ! solar UV model: Woods et al., 2016
  !-------------------------------------------------
  subroutine p__UV_flux(spl, grd, cst, & ! in
    &                   var, xct, flx  ) ! inout
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(cst_),           intent(in)    :: cst
    type(var_),           intent(inout) :: var
    type(xct_),           intent(inout) :: xct
    type(flx_),           intent(inout) :: flx
    integer iwl, i, j
    real(dp) tmp, tmp1, tmp2, tmp3, sum
    real(dp), allocatable :: wavelength(:), tmparr1(:), tmparr2(:)
    character(len=256) fname

    flx%nwl_UV = 2000
    allocate(flx%solar_UV(flx%nwl_UV), flx%lambda_UV(flx%nwl_UV))
    allocate(tmparr1(flx%nwl_UV), tmparr2(flx%nwl_UV))
    allocate(var%tau_UV(flx%nwl_UV,grd%nz),var%I_UV(flx%nwl_UV,grd%nz))
    allocate(xct%sigma_a_UV(flx%nwl_UV,grd%nz,spl%nsp), xct%sigma_d_UV(flx%nwl_UV,grd%nz,spl%nch))
    allocate(xct%sigma_a_UV_EUV(flx%nwl_UV,spl%nsp))

    do iwl = 1, flx%nwl_UV
      flx%lambda_UV(iwl) = dble(iwl)-0.5_dp ![nm]
    end do

    flx%solar_UV = 0.0_dp

    ! for comparison with Chaffin model
    fname = './UV/marssolarphotonflux.dat'
    open(11, file = fname, status = 'unknown' )
      do iwl =  1, flx%nwl_UV
        read(11,*) tmp, flx%solar_UV(iwl)
        flx%solar_UV(iwl) = flx%solar_UV(iwl)  * 1.0e4_dp  ![/cm^2/s] -> [/m^2/s], /2: according to Chaffin
      end do
    close(11)

    fname = './UV/ref_solar_irradiance_whi-2008_ver2.dat'
    open(11, file = fname, status = 'unknown' )
      do i = 1, 142; read(11,*); end do
      do iwl =  1, flx%nwl_UV
        sum = 0.0_dp
        do j = 1, 10
          read(11,*) tmp, tmp1, tmp2, tmp3
          sum = sum + tmp3 / 10.0_dp 
        end do
        sum = sum * flx%lambda_UV(iwl) * 1.0e-9_dp / cst%planck / cst%c
        flx%solar_UV(iwl) = sum / (flx%orbit * flx%orbit) 
      end do
    close(11)


  end subroutine p__UV_flux


  !-------------------------------------------------
  ! cross section data for UV [Chaffin et al., 2017]
  !-------------------------------------------------
  subroutine p__UV_cross_section_dat(xct            ) ! inout
    type(xct_),           intent(inout) :: xct
    integer iwl, jnu, i, j, nh
    real(dp) tmpu, tmpl, tmparr(10000)
    real(dp) a0u, a1u, a2u, axl, ax0, ax1, ax2
    real(dp) o3xdata(215+1601,2), o3chapxdata(0:20001,2), o3xdataadd(1601,2)
    real(dp) o2schr130K_raw(0:16000,6), o2schr190K_raw(0:16000,6), o2schr280K_raw(0:16000,6)
    character(len=256) fname

    ! CO2 data --------------------------------------------------------
    xct%co2xdata = 0.0_dp
    fname = 'UV/uvxsect/CO2.dat'
    call p__count_header(nh, fname)
    open(11, file = fname, status = 'unknown')
      do i = 1, nh
        read(11,*)
      end do
      do iwl = 1, 159
        read(11,*) xct%co2xdata(iwl,1), xct%co2xdata(iwl,2), xct%co2xdata(iwl,3)
        xct%co2xdata(iwl,2) = xct%co2xdata(iwl,2) * 1.0e-4_dp ! [cm2 -> m2]
        xct%co2xdata(iwl,3) = xct%co2xdata(iwl,3) * 1.0e-4_dp ! [cm2 -> m2]
        !print *, xct%co2xdata(iwl,1), xct%co2xdata(iwl,2), xct%co2xdata(iwl,3)
      end do
    close(11)

    ! H2O2 data --------------------------------------------------------
    xct%h2o2xdata = 0.0_dp
    fname = 'UV/uvxsect/H2O2.dat'
    call p__count_header(nh, fname)
    open(11, file = fname, status = 'unknown')
      do i = 1, nh
        read(11,*)
      end do
      do iwl = 1, 70
        read(11,*) xct%h2o2xdata(iwl,1), xct%h2o2xdata(iwl,2) ! [m2]
      end do
    close(11)

    ! H2O data --------------------------------------------------------
    xct%h2oxdata = 0.0_dp
    fname = 'UV/uvxsect/H2Oavgtbl.dat'
    call p__count_header(nh, fname)
    open(11, file = fname, status = 'unknown')
      do i = 1, nh
        read(11,*)
      end do
      do iwl = 1, 99
        read(11,*) xct%h2oxdata(iwl,1), xct%h2oxdata(iwl,2)
        xct%h2oxdata(iwl,2) = xct%h2oxdata(iwl,2) * 1.0e-4_dp ! [cm2 -> m2]
      end do
    close(11)

    ! H2 data --------------------------------------------------------
    xct%h2xdata = 0.0_dp
    fname = 'UV/uvxsect/binnedH2.csv'
    call p__count_header(nh, fname)
    open(11, file = fname, status = 'unknown')
      do i = 1, nh
        read(11,*)
      end do
      do iwl = 1, 277
        read(11,*) xct%h2xdata(iwl,1), xct%h2xdata(iwl,2) ! [m2]
        xct%h2xdata(iwl,2) = xct%h2xdata(iwl,2) * 1.0e-4_dp ! [cm2 -> m2]
      end do
    close(11)

    ! OH data --------------------------------------------------------
    xct%ohxdata = 0.0_dp
    fname = 'UV/uvxsect/binnedOH.csv'
    call p__count_header(nh, fname)
    open(11, file = fname, status = 'unknown')
      do i = 1, nh
        read(11,*)
      end do
      do iwl = 1, 601
        read(11,*) xct%ohxdata(iwl,1), xct%ohxdata(iwl,2) ! [m2]
        xct%ohxdata(iwl,2) = xct%ohxdata(iwl,2) * 1.0e-4_dp ! [cm2 -> m2]
      end do
    close(11)

    ! OH -> O1D data --------------------------------------------------------
    xct%oho1Dxdata = 0.0_dp
    fname = 'UV/uvxsect/binnedOHo1D.csv'
    call p__count_header(nh, fname)
    open(11, file = fname, status = 'unknown')
      do i = 1, nh
        read(11,*)
      end do
      do iwl = 1, 601
        read(11,*) xct%oho1Dxdata(iwl,1), xct%oho1Dxdata(iwl,2) ! [m2]
        xct%oho1Dxdata(iwl,2) = xct%oho1Dxdata(iwl,2) * 1.0e-4_dp ! [cm2 -> m2]
      end do
    close(11)

    ! O3 data --------------------------------------------------------
    xct%o3xdata = 0.0_dp
    o3chapxdata = 0.0_dp
    fname = 'UV/uvxsect/O3.dat'
    call p__count_header(nh, fname)
    open(11, file = fname, status = 'unknown')
      do i = 1, nh
        read(11,*)
      end do
      do iwl = 1, 215
        read(11,*) xct%o3xdata(iwl,1), xct%o3xdata(iwl,2)
        xct%o3xdata(iwl,2) = xct%o3xdata(iwl,2) * 1.0e-4_dp ! [cm2 -> m2]
        !print *, xct%o3xdata(iwl,1), xct%o3xdata(iwl,2)
      end do
    close(11)

    fname = 'UV/uvxsect/O3Chap.dat'
    call p__count_header(nh, fname)
    open(11, file = fname, status = 'unknown')
      do i = 1, nh
        read(11,*)
      end do
      do iwl = 1, 20001
        read(11,*) o3chapxdata(iwl,1), o3chapxdata(iwl,2)
        o3chapxdata(iwl,2) = o3chapxdata(iwl,2) * 1.0e-4_dp ! [cm2 -> m2]
      end do
    close(11)

    do iwl = 1, 20001
      o3chapxdata(iwl,1) = 1.0e7_dp/o3chapxdata(iwl,1)
    end do
    o3chapxdata(0,1) = 2001.0_dp

    do iwl = 400, 2000
      tmpu = 0.0_dp
      tmpl = 0.0_dp
      do jnu = int(1.0e7_dp/dble(iwl+1))-5000+1, int(1.0e7_dp/dble(iwl))-5000
        tmpu = tmpu + (o3chapxdata(jnu,1) - o3chapxdata(jnu-1,1))*o3chapxdata(jnu,2)
        tmpl = tmpl + (o3chapxdata(jnu,1) - o3chapxdata(jnu-1,1))
      end do
      xct%o3xdata(215+iwl-400+1,1) = dble(iwl)+0.5_dp
      xct%o3xdata(215+iwl-400+1,2) = tmpu/tmpl
      !print *, xct%o3xdata(215+iwl-400+1,1), xct%o3xdata(215+iwl-400+1,2)
    end do

    !stop

    ! O2 data --------------------------------------------------------
    xct%o2xdata = 0.0_dp
    xct%o2schr130K = 0.0_dp
    xct%o2schr190K = 0.0_dp
    xct%o2schr280K = 0.0_dp
    fname = 'UV/uvxsect/O2.dat'
    call p__count_header(nh, fname)
    open(11, file = fname, status = 'unknown')
      do i = 1, nh
        read(11,*)
      end do
      do iwl = 1, 115
        read(11,*) xct%o2xdata(iwl,1), xct%o2xdata(iwl,2)
        xct%o2xdata(iwl,2) = xct%o2xdata(iwl,2) * 1.0e-4_dp ! [cm2 -> m2]
      end do
    close(11)

    fname = 'UV/uvxsect/130-190.cf4'
    open(11, file = fname, status = 'unknown')
      do i = 1, 3
        read(11,*)
      end do
      do iwl = 1, 16000
        read(11,*) (o2schr130K_raw(iwl,j), j = 1, 6)
      end do
    close(11)

    fname = 'UV/uvxsect/190-280.cf4'
    open(11, file = fname, status = 'unknown')
      do i = 1, 3
        read(11,*)
      end do
      do iwl = 1, 16000
        read(11,*) (o2schr190K_raw(iwl,j), j = 1, 6)
      end do
    close(11)

    fname = 'UV/uvxsect/280-500.cf4'
    open(11, file = fname, status = 'unknown')
      do i = 1, 3
        read(11,*)
      end do
      do iwl = 1, 16000
        read(11,*) (o2schr280K_raw(iwl,j), j = 1, 6)
      end do
    close(11)

    do iwl = 1, 16000
      o2schr130K_raw(iwl,1) = 1.0e7_dp/o2schr130K_raw(iwl,1)
      o2schr190K_raw(iwl,1) = 1.0e7_dp/o2schr190K_raw(iwl,1)
      o2schr280K_raw(iwl,1) = 1.0e7_dp/o2schr280K_raw(iwl,1)
    end do

    do iwl = 176, 203
      a0u = 0.0_dp
      a1u = 0.0_dp
      a2u = 0.0_dp
      axl = 0.0_dp
      do jnu = 2*(int(1.0e7_dp/dble(iwl+1))-49000)+1, 2*(int(1.0e7_dp/dble(iwl))-49000)
        a0u = a0u + (o2schr130K_raw(jnu,1) - o2schr130K_raw(jnu-1,1))*o2schr130K_raw(jnu,2)
        a1u = a1u + (o2schr130K_raw(jnu,1) - o2schr130K_raw(jnu-1,1))*o2schr130K_raw(jnu,3)
        a2u = a2u + (o2schr130K_raw(jnu,1) - o2schr130K_raw(jnu-1,1))*o2schr130K_raw(jnu,4)
        axl = axl + (o2schr130K_raw(jnu,1) - o2schr130K_raw(jnu-1,1))
      end do
      xct%o2schr130K(iwl-176+1,1) = dble(iwl)+0.5_dp
      xct%o2schr130K(iwl-176+1,2) = a0u/axl
      xct%o2schr130K(iwl-176+1,3) = a1u/axl
      xct%o2schr130K(iwl-176+1,4) = a2u/axl

      a0u = 0.0_dp
      a1u = 0.0_dp
      a2u = 0.0_dp
      axl = 0.0_dp
      do jnu = 2*(int(1.0e7_dp/dble(iwl+1))-49000)+1, 2*(int(1.0e7_dp/dble(iwl))-49000)
        a0u = a0u + (o2schr190K_raw(jnu,1) - o2schr190K_raw(jnu-1,1))*o2schr190K_raw(jnu,2)
        a1u = a1u + (o2schr190K_raw(jnu,1) - o2schr190K_raw(jnu-1,1))*o2schr190K_raw(jnu,3)
        a2u = a2u + (o2schr190K_raw(jnu,1) - o2schr190K_raw(jnu-1,1))*o2schr190K_raw(jnu,4)
        axl = axl + (o2schr190K_raw(jnu,1) - o2schr190K_raw(jnu-1,1))
      end do
      xct%o2schr190K(iwl-176+1,1) = dble(iwl)+0.5_dp
      xct%o2schr190K(iwl-176+1,2) = a0u/axl
      xct%o2schr190K(iwl-176+1,3) = a1u/axl
      xct%o2schr190K(iwl-176+1,4) = a2u/axl

      a0u = 0.0_dp
      a1u = 0.0_dp
      a2u = 0.0_dp
      axl = 0.0_dp
      do jnu = 2*(int(1.0e7_dp/dble(iwl+1))-49000)+1, 2*(int(1.0e7_dp/dble(iwl))-49000)
        a0u = a0u + (o2schr280K_raw(jnu,1) - o2schr280K_raw(jnu-1,1))*o2schr280K_raw(jnu,2)
        a1u = a1u + (o2schr280K_raw(jnu,1) - o2schr280K_raw(jnu-1,1))*o2schr280K_raw(jnu,3)
        a2u = a2u + (o2schr280K_raw(jnu,1) - o2schr280K_raw(jnu-1,1))*o2schr280K_raw(jnu,4)
        axl = axl + (o2schr280K_raw(jnu,1) - o2schr280K_raw(jnu-1,1))
      end do
      xct%o2schr280K(iwl-176+1,1) = dble(iwl)+0.5_dp
      xct%o2schr280K(iwl-176+1,2) = a0u/axl
      xct%o2schr280K(iwl-176+1,3) = a1u/axl
      xct%o2schr280K(iwl-176+1,4) = a2u/axl

      !print *, xct%o2schr130K(iwl-176+1,1), &
      !  &      xct%o2schr130K(iwl-176+1,2), &
      !  &      xct%o2schr130K(iwl-176+1,3), &
      !  &      xct%o2schr130K(iwl-176+1,4)


    end do

    !stop


  end subroutine p__UV_cross_section_dat


  !-------------------------------------------------
  ! cross section data for UV
  !-------------------------------------------------
  subroutine p__UV_cross_section_exe(T, nz, spl, flx, & ! in
    &                                xct              ) ! inout
    integer,    intent(in)    :: nz
    real(dp),   intent(in)    :: T(nz)
    type(spl_), intent(in)    :: spl
    type(flx_), intent(in)    :: flx
    type(xct_), intent(inout) :: xct
    integer iwl, jwl, i, iz, isp, ich
    real(dp) Tclamp, Tfrac, ratio, lambda
    real(dp) X_O3q(3), w_O3q(3), A_O3q(3), v_O3q(2), c_O3q, R_O3q, q_O3q(2), qsum_O3q
    real(dp) O2schr(100,4), del, l, o2xdata(115,2)
    real(dp) A_H2O2(8), B_H2O2(5), lpowA, lpowB, expfac
    real(dp) a, b, sigmamed, vmed, v
    character(len=256) reactants(10),  products(10)

    !----------------------------------------------------------------------------------------------------
    ! UV photoabsorption cross section

    if (xct%type == 'absorption') then

      do isp = 1, spl%nsp

        ! CO2
        if (spl%species(isp) == 'CO2') then

          do iz = 1, nz
            Tclamp = T(iz)
            if (T(iz)<=195.0_dp) Tclamp = 195.0_dp
            if (T(iz)>=295.0_dp) Tclamp = 295.0_dp
            Tfrac = (Tclamp-195.0_dp)/(295.0_dp-195.0_dp)

            !write(num,*) iz
            !fname = 'UV/xct/CO2/'//trim(ADJUSTL(num))//'.dat'
            !open(11, file = fname, status = 'unknown' )

            do iwl = nint(xct%co2xdata(1,1)+0.5_dp), nint(xct%co2xdata(159,1)+0.5_dp)
              jwl = iwl - nint(xct%co2xdata(1,1)+0.5_dp) + 1
              xct%sigma_a_UV(iwl,iz,isp) = (1.0_dp-Tfrac)*xct%co2xdata(jwl,2) + Tfrac*xct%co2xdata(jwl,3)
              !write(11, *) dble(iwl) - 0.5_dp, xct%sigma_a_UV(iwl,iz,isp)
            end do

            !close(11)

          end do

        ! H2O2
        else if (spl%species(isp) == 'H2O2') then
          do iz = 1, nz

            !write(num,*) iz
            !fname = 'UV/xct/H2O2/'//trim(ADJUSTL(num))//'.dat'
            !open(11, file = fname, status = 'unknown' )

            do iwl = nint(xct%h2o2xdata(1,1)+0.5_dp), nint(xct%h2o2xdata(70,1)+0.5_dp)
              jwl = iwl - nint(xct%h2o2xdata(1,1)+0.5_dp) + 1
              xct%sigma_a_UV(iwl,iz,isp) = xct%h2o2xdata(jwl,2)
              !write(11, *) dble(iwl) - 0.5_dp, xct%sigma_a_UV(iwl,iz,isp)
            end do

            Tclamp = T(iz)
            if (T(iz)<=200.0_dp) Tclamp = 200.0_dp
            if (T(iz)>=400.0_dp) Tclamp = 400.0_dp

            A_H2O2(1:8) = (/64761.0_dp, -921.70972_dp, 4.535649_dp, &
              &             -0.0044589016_dp, -0.00004035101_dp, &
              &             1.6878206e-7_dp, -2.652014e-10_dp, 1.5534675e-13_dp/)
            B_H2O2(1:5) = (/6812.3_dp, -51.351_dp, 0.11522_dp, -0.000030493_dp, -1.0924e-7_dp/)

            expfac = 1.0_dp/(1.0_dp+dexp(-1265.0_dp/Tclamp))

            do iwl = 261, 351
              lambda = dble(iwl) - 0.5_dp
              lpowA = 0.0_dp
              do i = 0, 7
                lpowA = lpowA + A_H2O2(i+1) * lambda ** dble(i)
              end do
              lpowB = 0.0_dp
              do i = 0, 4
                lpowB = lpowB + B_H2O2(i+1) * lambda ** dble(i)
              end do
              xct%sigma_a_UV(iwl,iz,isp) = 1.0e-25_dp * (expfac * lpowA + (1.0_dp-expfac) * lpowB)
              !write(11, *) dble(iwl) - 0.5_dp, xct%sigma_a_UV(iwl,iz,isp)
            end do

            !close(11)

          end do

        ! H2O
        else if (spl%species(isp) == 'H2O') then

          do iz = 1, nz

            !write(num,*) iz
            !fname = 'UV/xct/H2O/'//trim(ADJUSTL(num))//'.dat'
            !open(11, file = fname, status = 'unknown' )

            do iwl = nint(xct%h2oxdata(1,1)+0.5_dp), nint(xct%h2oxdata(99,1)+0.5_dp)
              jwl = iwl - nint(xct%h2oxdata(1,1)+0.5_dp) + 1
              xct%sigma_a_UV(iwl,iz,isp) = xct%h2oxdata(jwl,2)
              !write(11, *) dble(iwl) - 0.5_dp, xct%sigma_a_UV(iwl,iz,isp)
            end do

            !close(11)

          end do

        ! HO2
        else if (spl%species(isp) == 'HO2') then

          do iz = 1, nz

            !write(num,*) iz
            !fname = 'UV/xct/HO2/'//trim(ADJUSTL(num))//'.dat'
            !open(11, file = fname, status = 'unknown' )

            do iwl = 191, 250
              l = dble(iwl) - 0.5_dp
              a = 4.91_dp
              b = 30612.0_dp
              sigmamed = 1.64e-18_dp
              vmed = 50260.0_dp
              v = 1.0e7_dp/l
              xct%sigma_a_UV(iwl,iz,isp) = 1.0e-4_dp * &
              &  sigmamed / ( 1.0_dp - b/v ) * dexp( -a * dlog( (v-b)/(vmed-b) )**2.0_dp )
              !write(11, *) dble(iwl) - 0.5_dp, xct%sigma_a_UV(iwl,iz,isp)
            end do

            !close(11)

          end do

        ! H2
        else if (spl%species(isp) == 'H2') then

          do iz = 1, nz

            !write(num,*) iz
            !fname = 'UV/xct/H2/'//trim(ADJUSTL(num))//'.dat'
            !open(11, file = fname, status = 'unknown' )

            do iwl = nint(xct%h2xdata(1,1)+0.5_dp), nint(xct%h2xdata(277,1)+0.5_dp)
              jwl = iwl - nint(xct%h2xdata(1,1)+0.5_dp) + 1
              xct%sigma_a_UV(iwl,iz,isp) = xct%h2xdata(jwl,2)
              !write(11, *) dble(iwl) - 0.5_dp, xct%sigma_a_UV(iwl,iz,isp)
            end do

            !close(11)

          end do

        ! OH
        else if (spl%species(isp) == 'OH') then

          do iz = 1, nz

            !write(num,*) iz
            !fname = 'UV/xct/OH/'//trim(ADJUSTL(num))//'.dat'
            !open(11, file = fname, status = 'unknown' )

            do iwl = nint(xct%ohxdata(1,1)+0.5_dp), nint(xct%ohxdata(601,1)+0.5_dp)
              jwl = iwl - nint(xct%ohxdata(1,1)+0.5_dp) + 1
              xct%sigma_a_UV(iwl,iz,isp) = xct%ohxdata(jwl,2) + xct%oho1Dxdata(jwl,2)
              !write(*,*) dble(iwl) - 0.5_dp, xct%sigma_a_UV(iwl,iz,isp)
            end do

            !close(11)

          end do

        ! O3
        else if (spl%species(isp) == 'O3') then
          do iz = 1, nz

            !write(num,*) iz
            !fname = 'UV/xct/O3/'//trim(ADJUSTL(num))//'.dat'
            !open(11, file = fname, status = 'unknown' )

            do iwl = nint(xct%o3xdata(1,1)+0.5_dp), nint(xct%o3xdata(215+1601,1)+0.5_dp)
              jwl = iwl - nint(xct%O3xdata(1,1)+0.5_dp) + 1
              xct%sigma_a_UV(iwl,iz,isp) = xct%o3xdata(jwl,2)
              !write(11, *) dble(iwl) - 0.5_dp, xct%sigma_a_UV(iwl,iz,isp)
            end do

            !close(11)

          end do

        ! O2
        else if (spl%species(isp) == 'O2') then
          do iz = 1, nz
            Tclamp = T(iz)
            if (T(iz)<=130.0_dp) Tclamp = 130.0_dp
            if (T(iz)>=500.0_dp) Tclamp = 500.0_dp

            O2schr = 0.0_dp
            if (Tclamp >= 130.0_dp .and. Tclamp < 190.0_dp) then
              do iwl = 177, 204
                O2schr(iwl-176,1) = xct%o2schr130K(iwl-176,1)
                O2schr(iwl-176,2) = xct%o2schr130K(iwl-176,2)
                O2schr(iwl-176,3) = xct%o2schr130K(iwl-176,3)
                O2schr(iwl-176,4) = xct%o2schr130K(iwl-176,4)
              end do
            else if (Tclamp >= 190.0_dp .and. Tclamp < 280.0_dp) then
              do iwl = 177, 204
                O2schr(iwl-176,1) = xct%o2schr190K(iwl-176,1)
                O2schr(iwl-176,2) = xct%o2schr190K(iwl-176,2)
                O2schr(iwl-176,3) = xct%o2schr190K(iwl-176,3)
                O2schr(iwl-176,4) = xct%o2schr190K(iwl-176,4)
              end do
            else if (Tclamp >= 280.0_dp) then
              do iwl = 177, 204
                O2schr(iwl-176,1) = xct%o2schr280K(iwl-176,1)
                O2schr(iwl-176,2) = xct%o2schr280K(iwl-176,2)
                O2schr(iwl-176,3) = xct%o2schr280K(iwl-176,3)
                O2schr(iwl-176,4) = xct%o2schr280K(iwl-176,4)
              end do
            end if

            do iwl = nint(xct%o2xdata(1,1)+0.5_dp), nint(xct%o2xdata(115,1)+0.5_dp)
              jwl = iwl - nint(xct%O2xdata(1,1)+0.5_dp) + 1
              o2xdata(jwl,2) = xct%o2xdata(jwl,2)
            end do

            !fill in the schumann-runge bands according to Minschwaner 1992
            del = ((Tclamp-100.0_dp)/10.0_dp)**2.0_dp
            do iwl = 177, 204
              jwl = iwl - nint(xct%O2xdata(1,1)+0.5_dp) + 1
              o2xdata(jwl,2) = o2xdata(jwl,2) &
                & + 1.0e-24_dp * O2schr(iwl-176,2) * del**2.0_dp &
                & + 1.0e-24_dp * O2schr(iwl-176,3) * del &
                & + 1.0e-24_dp * O2schr(iwl-176,4)
            end do

            ! add in the herzberg continuum (though tiny)
            ! measured by yoshino 1992
            do iwl = 193, 245
              jwl = iwl - nint(xct%O2xdata(1,1)+0.5_dp) + 1
              l = dble(iwl) - 0.5_dp
              o2xdata(jwl,2) = o2xdata(jwl,2) &
                & + 1.0e-28_dp*(-2.3837947e4_dp &
                &             +4.1973085e2_dp *l &
                &             -2.7640139e0_dp *l**2.0_dp &
                &             +8.0723193e-3_dp*l**3.0_dp &
                &             -8.8255447e-6_dp*l**4.0_dp)
            end do

            !write(num,*) iz
            !fname = 'UV/xct/O2/'//trim(ADJUSTL(num))//'.dat'
            !open(11, file = fname, status = 'unknown' )

            do iwl = nint(xct%o2xdata(1,1)+0.5_dp), nint(xct%o2xdata(115,1)+0.5_dp)
              jwl = iwl - nint(xct%O2xdata(1,1)+0.5_dp) + 1
              xct%sigma_a_UV(iwl,iz,isp) = o2xdata(jwl,2)
              !write(11, *) dble(iwl) - 0.5_dp, xct%sigma_a_UV(iwl,iz,isp)
            end do

            !close(11)

          end do
          
        end if
      end do!isp

    end if

    !----------------------------------------------------------------------------------------------------
    ! UV photodissociation cross section

    if (xct%type == 'dissociation') then

      do ich = 1, spl%nch
        if (spl%reaction_type_char(ich) == 'photodissociation') then
          do i = 2, spl%reactant_list(ich,1)+1
            reactants(i-1) = trim(spl%species(spl%reactant_list(ich,i)))
          end do
          do i = 2, spl%product_list(ich,1)+1
            products(i-1) = trim(spl%species(spl%product_list(ich,i)))
          end do
          isp = spl%reactant_list(ich,2)

          ! CO2 -------------------------------------------
          ! CO2 + hv -> CO + O
          if ( reactants(1) == 'CO2' &
          & .and. products(1) == 'CO' .and. products(2) == 'O' ) then

            do iz = 1, nz
              do iwl = nint(xct%co2xdata(1,1)+0.5_dp), nint(xct%co2xdata(159,1)+0.5_dp)
                if ( flx%lambda_UV(iwl) > 167.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 1.0_dp
                else if ( flx%lambda_UV(iwl) < 95.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.5_dp
                end if
              end do
              !print *, xct%sigma_d_UV(170,iz), xct%sigma_a_UV(170,iz,isp)
            end do

          ! CO2 + hv -> CO + O(1D)
          else if ( reactants(1) == 'CO2' &
          & .and. products(1) == 'CO' .and. products(2) == 'O(1D)' ) then

            do iz = 1, nz
              do iwl = nint(xct%co2xdata(1,1)+0.5_dp), nint(xct%co2xdata(159,1)+0.5_dp)
                if ( flx%lambda_UV(iwl) > 95.0_dp .and. flx%lambda_UV(iwl) < 167.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 1.0_dp
                else if ( flx%lambda_UV(iwl) < 95.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.5_dp
                end if
              end do
            end do

          ! H2O -------------------------------------------
          ! H2O + hv -> H + OH
          else if ( reactants(1) == 'H2O' &
          & .and. products(1) == 'H' .and. products(2) == 'OH' ) then

            do iz = 1, nz
              do iwl = nint(xct%h2oxdata(1,1)+0.5_dp), nint(xct%h2oxdata(99,1)+0.5_dp)
                if ( flx%lambda_UV(iwl) < 145.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.89_dp
                else if ( flx%lambda_UV(iwl) > 145.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 1.0_dp
                end if
              end do
            end do

          ! H2O + hv -> H2 + O(1D)
          else if ( reactants(1) == 'H2O' &
          & .and. products(1) == 'H2' .and. products(2) == 'O(1D)' ) then

            do iz = 1, nz
              do iwl = nint(xct%h2oxdata(1,1)+0.5_dp), nint(xct%h2oxdata(99,1)+0.5_dp)
                if ( flx%lambda_UV(iwl) < 145.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.11_dp
                else if ( flx%lambda_UV(iwl) > 145.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.0_dp
                end if
              end do
            end do

          ! H2O + hv -> H + H + O
          else if ( reactants(1) == 'H2O' &
          & .and. products(1) == 'H' .and. products(2) == 'H' .and. products(3) == 'O' ) then

            do iz = 1, nz
              do iwl = nint(xct%h2oxdata(1,1)+0.5_dp), nint(xct%h2oxdata(99,1)+0.5_dp)
                xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.0_dp
              end do
            end do

          ! H2O2 -------------------------------------------
          ! H2O2 + hv -> OH + OH
          else if ( reactants(1) == 'H2O2' &
          & .and. products(1) == 'OH' .and. products(2) == 'OH' ) then

            do iz = 1, nz
              do iwl = nint(xct%h2o2xdata(1,1)+0.5_dp), 351
                if ( flx%lambda_UV(iwl) < 230.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.85_dp
                else if ( flx%lambda_UV(iwl) > 230.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 1.0_dp
                end if
              end do
            end do

          ! H2O2 + hv -> HO2 + H
          else if ( reactants(1) == 'H2O2' &
          & .and. products(1) == 'HO2' .and. products(2) == 'H' ) then

            do iz = 1, nz
              do iwl = nint(xct%h2o2xdata(1,1)+0.5_dp), 351
                if ( flx%lambda_UV(iwl) < 230.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.15_dp
                else if ( flx%lambda_UV(iwl) > 230.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.0_dp
                end if
              end do
            end do

          ! H2O2 + hv -> H2O + O(1D)
          else if ( reactants(1) == 'H2O2' &
          & .and. products(1) == 'H2O' .and. products(2) == 'O(1D)' ) then

            do iz = 1, nz
              do iwl = nint(xct%h2o2xdata(1,1)+0.5_dp), 351
                xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.0_dp
              end do
            end do

          ! O3 -------------------------------------------
          ! O3 + hv -> O2 + O
          else if ( reactants(1) == 'O3' &
          & .and. products(1) == 'O2' .and. products(2) == 'O' ) then

            do iz = 1, nz
              do iwl = nint(xct%o3xdata(1,1)+0.5_dp), nint(xct%o3xdata(215+1601,1)+0.5_dp)
                if ( flx%lambda_UV(iwl) < 193.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) &
                    &                    * (1.0_dp - (1.37e-2_dp*193.0_dp-2.16_dp))
                else if ( flx%lambda_UV(iwl) > 193.0_dp .and. flx%lambda_UV(iwl) < 225.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) &
                    &                    * (1.0_dp - (1.37e-2_dp*flx%lambda_UV(iwl)-2.16_dp))
                else if ( flx%lambda_UV(iwl) > 225.0_dp .and. flx%lambda_UV(iwl) < 306.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.1_dp
                else if ( flx%lambda_UV(iwl) > 306.0_dp .and. flx%lambda_UV(iwl) < 328.0_dp ) then
                  lambda = flx%lambda_UV(iwl)
                  X_O3q(1:3) = (/304.225_dp, 314.957_dp, 310.737_dp/)
                  w_O3q(1:3) = (/5.576_dp, 6.601_dp, 2.187_dp/)
                  A_O3q(1:3) = (/0.8036_dp, 8.9061_dp, 0.1192_dp/)
                  v_O3q(1:2) = (/0.0_dp, 825.518_dp/)
                  c_O3q      = 0.0765_dp
                  R_O3q      = 0.695_dp
                  q_O3q(1)   = dexp(-v_O3q(1)/(R_O3q*T(iz)))
                  q_O3q(2)   = dexp(-v_O3q(2)/(R_O3q*T(iz)))
                  qsum_O3q   = q_O3q(1) + q_O3q(2)
                  ratio = q_O3q(1)/qsum_O3q*A_O3q(1)*dexp(-((X_O3q(1)-lambda)/w_O3q(1))**4.0_dp) &
                    &   + q_O3q(2)/qsum_O3q*A_O3q(2)*(T(iz)/300.0_dp)**2.0_dp*dexp(-((X_O3q(2)-lambda)/w_O3q(2))**2.0_dp) &
                    &   + A_O3q(3)*(T(iz)/300.0_dp)**1.5_dp*dexp(-((X_O3q(3)-lambda)/w_O3q(3))**2.0_dp) &
                    &   + c_O3q
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * (1.0_dp - ratio)
                else if ( flx%lambda_UV(iwl) > 328.0_dp .and. flx%lambda_UV(iwl) < 340.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.92_dp
                else if ( flx%lambda_UV(iwl) > 340.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 1.0_dp
                end if
              end do
            end do

          ! O3 + hv -> O2 + O
          else if ( reactants(1) == 'O3' &
          & .and. products(1) == 'O2' .and. products(2) == 'O(1D)' ) then

            do iz = 1, nz
              do iwl = nint(xct%o3xdata(1,1)+0.5_dp), nint(xct%o3xdata(215+1601,1)+0.5_dp)
                if ( flx%lambda_UV(iwl) < 193.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) &
                    &                    * (1.37e-2_dp*193.0_dp-2.16_dp)
                else if ( flx%lambda_UV(iwl) > 193.0_dp .and. flx%lambda_UV(iwl) < 225.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) &
                    &                    * (1.37e-2_dp*flx%lambda_UV(iwl)-2.16_dp)
                else if ( flx%lambda_UV(iwl) > 225.0_dp .and. flx%lambda_UV(iwl) < 306.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.9_dp
                else if ( flx%lambda_UV(iwl) > 306.0_dp .and. flx%lambda_UV(iwl) < 328.0_dp ) then
                  lambda = flx%lambda_UV(iwl)
                  X_O3q(1:3) = (/304.225_dp, 314.957_dp, 310.737_dp/)
                  w_O3q(1:3) = (/5.576_dp, 6.601_dp, 2.187_dp/)
                  A_O3q(1:3) = (/0.8036_dp, 8.9061_dp, 0.1192_dp/)
                  v_O3q(1:2) = (/0.0_dp, 825.518_dp/)
                  c_O3q      = 0.0765_dp
                  R_O3q      = 0.695_dp
                  q_O3q(1)   = dexp(-v_O3q(1)/(R_O3q*T(iz)))
                  q_O3q(2)   = dexp(-v_O3q(2)/(R_O3q*T(iz)))
                  qsum_O3q   = q_O3q(1) + q_O3q(2)
                  ratio = q_O3q(1)/qsum_O3q*A_O3q(1)*dexp(-((X_O3q(1)-lambda)/w_O3q(1))**4.0_dp) &
                    &   + q_O3q(2)/qsum_O3q*A_O3q(2)*(T(iz)/300.0_dp)**2.0_dp*dexp(-((X_O3q(2)-lambda)/w_O3q(2))**2.0_dp) &
                    &   + A_O3q(3)*(T(iz)/300.0_dp)**1.5_dp*dexp(-((X_O3q(3)-lambda)/w_O3q(3))**2.0_dp) &
                    &   + c_O3q
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * ratio
                else if ( flx%lambda_UV(iwl) > 328.0_dp .and. flx%lambda_UV(iwl) < 340.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.08_dp
                else if ( flx%lambda_UV(iwl) > 340.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.0_dp
                end if
              end do
            end do

          ! O2 -------------------------------------------
          ! O2 + hv -> O + O
          else if ( reactants(1) == 'O2' &
          & .and. products(1) == 'O' .and. products(2) == 'O' ) then

            do iz = 1, nz
              do iwl = nint(xct%o2xdata(1,1)+0.5_dp), nint(xct%o2xdata(115,1)+0.5_dp)
                if ( flx%lambda_UV(iwl) < 175.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.0_dp
                else if ( flx%lambda_UV(iwl) > 175.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 1.0_dp
                end if
              end do
            end do

          ! O2 + hv -> O + O
          else if ( reactants(1) == 'O2' &
          & .and. products(1) == 'O' .and. products(2) == 'O(1D)' ) then

            do iz = 1, nz
              do iwl = nint(xct%o2xdata(1,1)+0.5_dp), nint(xct%o2xdata(115,1)+0.5_dp)
                if ( flx%lambda_UV(iwl) < 175.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 1.0_dp
                else if ( flx%lambda_UV(iwl) > 175.0_dp ) then
                  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 0.0_dp
                end if
              end do
            end do

          ! H2 -------------------------------------------
          ! H2 + hv -> H + H
          else if ( reactants(1) == 'H2' &
          & .and. products(1) == 'H' .and. products(2) == 'H' ) then

            do iz = 1, nz
              do iwl = nint(xct%h2xdata(1,1)+0.5_dp), nint(xct%h2xdata(277,1)+0.5_dp)
                xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 1.0_dp
              end do
            end do

          ! OH -------------------------------------------
          ! OH + hv -> O + H
          else if ( reactants(1) == 'OH' &
          & .and. products(1) == 'O' .and. products(2) == 'H' ) then

          do iz = 1, nz
            do iwl = nint(xct%ohxdata(1,1)+0.5_dp), nint(xct%ohxdata(601,1)+0.5_dp)
              jwl = iwl - nint(xct%ohxdata(1,1)+0.5_dp) + 1
              xct%sigma_d_UV(iwl,iz,ich) = xct%ohxdata(jwl,2)
            end do
          end do

          ! OH + hv -> O(1D) + H
          else if ( reactants(1) == 'OH' &
          & .and. products(1) == 'O(1D)' .and. products(2) == 'H' ) then

            do iz = 1, nz
              do iwl = nint(xct%ohxdata(1,1)+0.5_dp), nint(xct%ohxdata(601,1)+0.5_dp)
                jwl = iwl - nint(xct%ohxdata(1,1)+0.5_dp) + 1
                xct%sigma_d_UV(iwl,iz,ich) = xct%oho1Dxdata(jwl,2)
              end do
            end do

          ! HO2 -------------------------------------------
          ! HO2 + hv -> OH + H
          else if ( reactants(1) == 'HO2' &
          & .and. products(1) == 'OH' .and. products(2) == 'O' ) then

            do iz = 1, nz
              do iwl = 191, 250
                xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * 1.0_dp
              end do
            end do



          end if


        end if 
      end do!ich


    end if
    !stop



  end subroutine p__UV_cross_section_exe



  !-------------------------------------------------
  ! cross section data for UV
  !-------------------------------------------------
  subroutine p__UV_EUV_xct_convergence_exe(nz, spl, var, flx, & ! in
    &                                      xct                ) ! inout
    integer,    intent(in)    :: nz
    type(spl_), intent(in)    :: spl
    type(var_), intent(in)    :: var
    type(flx_), intent(in)    :: flx
    type(xct_), intent(inout) :: xct
    integer iwl, jwl, i, iz, isp
    character(len=256) fname, num

    !--------------------------------
    ! UV absorption cross section
    !--------------------------------
    if (xct%type == 'absorption') then
      do isp = 1, spl%nsp!                                               ! UV               : EUVAC
        xct%sigma_a_UV_EUV(  1: 10,isp)  = xct%sigma_a_EUV(1,isp)        !   0.5 -   9.5 nm :   7.5 nm 
        xct%sigma_a_UV_EUV( 11: 15,isp)  = xct%sigma_a_EUV(2,isp)        !  10.5 -  14.5 nm :  12.5 nm 
        xct%sigma_a_UV_EUV( 16: 20,isp)  = xct%sigma_a_EUV(3,isp)        !  15.5 -  19.5 nm :  17.5 nm 
        xct%sigma_a_UV_EUV( 21: 25,isp)  = xct%sigma_a_EUV(4,isp)        !  20.5 -  24.5 nm :  22.5 nm 
        xct%sigma_a_UV_EUV( 26: 30,isp)  = xct%sigma_a_EUV(7,isp)        !  25.5 -  29.5 nm :  27.5 nm 
        xct%sigma_a_UV_EUV( 31: 35,isp)  = xct%sigma_a_EUV(10,isp)       !  30.5 -  34.5 nm :  32.5 nm 
        xct%sigma_a_UV_EUV( 36: 40,isp)  = xct%sigma_a_EUV(12,isp)       !  35.5 -  39.5 nm :  37.5 nm 
        xct%sigma_a_UV_EUV( 41: 45,isp)  = xct%sigma_a_EUV(13,isp)       !  40.5 -  44.5 nm :  42.5 nm 
        xct%sigma_a_UV_EUV( 46: 50,isp)  = xct%sigma_a_EUV(15,isp)       !  45.5 -  49.5 nm :  47.5 nm 
        xct%sigma_a_UV_EUV( 51: 55,isp)  = xct%sigma_a_EUV(16,isp)       !  50.5 -  54.5 nm :  52.5 nm 
        xct%sigma_a_UV_EUV( 56: 60,isp)  = xct%sigma_a_EUV(19,isp)       !  55.5 -  59.5 nm :  57.5 nm 
        xct%sigma_a_UV_EUV( 61: 65,isp)  = xct%sigma_a_EUV(22,isp)       !  60.5 -  65.5 nm :  62.5 nm 
        xct%sigma_a_UV_EUV( 66: 70,isp)  =(xct%sigma_a_EUV(23,isp) &
          &                              + xct%sigma_a_EUV(24,isp))/2.0d0!  65.5 -  69.5 nm :  67.5 nm 
        xct%sigma_a_UV_EUV( 71: 75,isp)  = xct%sigma_a_EUV(25,isp)       !  70.5 -  74.5 nm :  72.5 nm 
        xct%sigma_a_UV_EUV( 76: 80,isp)  = xct%sigma_a_EUV(29,isp)       !  75.5 -  79.5 nm :  77.5 nm 
        xct%sigma_a_UV_EUV( 81: 85,isp)  = xct%sigma_a_EUV(30,isp)       !  80.5 -  84.5 nm :  82.5 nm 
        xct%sigma_a_UV_EUV( 86: 90,isp)  = xct%sigma_a_EUV(31,isp)       !  85.5 -  89.5 nm :  87.5 nm 
        xct%sigma_a_UV_EUV( 91: 95,isp)  = xct%sigma_a_EUV(32,isp)       !  90.5 -  94.5 nm :  92.5 nm 
        xct%sigma_a_UV_EUV( 96:100,isp)  = xct%sigma_a_EUV(34,isp)       !  95.5 -  99.5 nm :  97.5 nm 
        xct%sigma_a_UV_EUV(101:105,isp)  = xct%sigma_a_EUV(37,isp)       ! 100.5 - 104.5 nm : 102.5 nm 

        xct%sigma_a_UV_EUV( 26,isp)  = xct%sigma_a_EUV(5,isp)        !  25.5 nm :  25.632 nm
        xct%sigma_a_UV_EUV( 29,isp)  = xct%sigma_a_EUV(6,isp)        !  28.5 nm :  28.415 nm
        xct%sigma_a_UV_EUV( 31,isp)  =(xct%sigma_a_EUV(8,isp) &  
          &                          + xct%sigma_a_EUV(9,isp))/2.0_dp!  30.5 nm :  30.331 nm & 30.378 nm
        xct%sigma_a_UV_EUV( 37,isp)  = xct%sigma_a_EUV(11,isp)       !  36.5 nm :  36.807 nm
        xct%sigma_a_UV_EUV( 47,isp)  = xct%sigma_a_EUV(14,isp)       !  46.5 nm :  46.522 nm
        xct%sigma_a_UV_EUV( 56,isp)  = xct%sigma_a_EUV(17,isp)       !  55.5 nm :  55.437 nm
        xct%sigma_a_UV_EUV( 59,isp)  = xct%sigma_a_EUV(18,isp)       !  58.5 nm :  58.433 nm
        xct%sigma_a_UV_EUV( 61,isp)  = xct%sigma_a_EUV(20,isp)       !  60.5 nm :  60.976 nm
        xct%sigma_a_UV_EUV( 62,isp)  = xct%sigma_a_EUV(20,isp)       !  61.5 nm :  60.976 nm
        xct%sigma_a_UV_EUV( 63,isp)  = xct%sigma_a_EUV(21,isp)       !  62.5 nm :  62.973 nm
        xct%sigma_a_UV_EUV( 64,isp)  = xct%sigma_a_EUV(21,isp)       !  63.5 nm :  62.973 nm
        xct%sigma_a_UV_EUV( 71,isp)  = xct%sigma_a_EUV(24,isp)       !  70.5 nm :  70.336 nm
        xct%sigma_a_UV_EUV( 77,isp)  = xct%sigma_a_EUV(26,isp)       !  76.5 nm :  76.515 nm
        xct%sigma_a_UV_EUV( 79,isp)  = xct%sigma_a_EUV(28,isp)       !  78.5 nm :  78.936 nm
        xct%sigma_a_UV_EUV( 80,isp)  = xct%sigma_a_EUV(28,isp)       !  79.5 nm :  78.936 nm
        xct%sigma_a_UV_EUV( 98,isp)  = xct%sigma_a_EUV(28,isp)       !  97.5 nm :  97.702 nm
        xct%sigma_a_UV_EUV(103,isp)  = xct%sigma_a_EUV(36,isp)       ! 102.5 nm : 103.191 nm
        xct%sigma_a_UV_EUV(104,isp)  = xct%sigma_a_EUV(36,isp)       ! 103.5 nm : 103.191 nm
      end do


    end if
    !stop



  end subroutine p__UV_EUV_xct_convergence_exe

end module p__UV
