! p__photochem_scheme.f90
!
!
!
!

module p__photochem_scheme
  use v__tdec,        only : spl_, var_, grd_, cst_, xct_, flx_, set_
  use p__search,      only : p__search_reactant, p__search_product, sp_index
  use p__io,          only : p__io_stable__fin, p__io_rotation__fin, p__io_deallocate, &
    &                        p__io_progress

  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private
  public :: p__photochem_scheme__exe

contains


  subroutine p__photochem_scheme__exe(spl, grd, set, & ! in
    &                                 var            ) ! inout
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(set_),           intent(in)    :: set
    type(var_),           intent(inout) :: var

    integer isp, jsp, ich, jch, iz
    integer i, j, k, ijch, ilist, ii, jj
    integer bij
    real(dp) a0, nj, nk, kij, tmp
    real(dp) dn(spl%nsp_i, grd%nz)

    real(dp) lambda ! implicit factor
    real(dp) eps

    real(dp) dni_ni(spl%nsp,grd%nz)
    character(len=256) model

    !-------------------------------------------------------------------------------------------------
    !
    !                                    Photochemical calculation
    !
    !-------------------------------------------------------------------------------------------------

    model = set%inversion

    if (model == 'Catling') eps = 1.0e-1 ! for convergence
    if (model == 'Chaffin') eps = 1.0e-8 ! for convergence

    ! advance time step scheme ------------------------------------------------------------------

      !---------------------------------------------------------
      ! explicit scheme
      !---------------------------------------------------------
        if (set%scheme == 'explicit') then

          do iz = 1, grd%nz
            do isp = 1, spl%nsp
              jsp = spl%all_to_var(isp)
              if (spl%label_fix(isp) == 0) then
                dn(jsp,iz) = ( var%Pi(jsp,iz) - var%Li(jsp,iz) - var%dPhi_dz(jsp,iz) ) * var%dtime
                var%ni_new(isp,iz) = var%ni(isp,iz) + dn(jsp,iz)
              else if (spl%label_fix(isp) == 1) then
                var%ni_new(isp,iz)  = var%ni(isp,iz)
              end if
            end do
          end do

        end if

      !---------------------------------------------------------
      ! semi-implicit scheme
      !---------------------------------------------------------
        if (set%scheme == 'semi-implicit') then

          do iz = 1, grd%nz
            do isp = 1, spl%nsp
              jsp = spl%all_to_var(isp)
              var%ni_new(isp,iz) = 0.0_dp
              if (spl%label_fix(isp) == 0) then
                if (var%ni(isp,iz) /= 0.0_dp) then
                  var%ni_new(isp,iz) = ( var%ni(isp,iz) + (var%Pi(jsp,iz) - var%dPhi_dz(jsp,iz)) * var%dtime ) &
                    &                  / (1.0_dp + (var%Li(jsp,iz) / var%ni(isp,iz)) * var%dtime )
                end if

                if (var%ni(isp,iz) == 0.0_dp) then
                  var%ni_new(isp,iz) = ( var%ni(isp,iz) + (var%Pi(jsp,iz) - var%dPhi_dz(jsp,iz)) * var%dtime ) &
                    &                  / (1.0_dp)
                end if

              else if (spl%label_fix(isp) == 1) then
                var%ni_new(isp,iz)  = var%ni(isp,iz)
              end if

            end do
          end do

        end if

      !---------------------------------------------------------
      ! implicit scheme
      !---------------------------------------------------------
        if (set%scheme == 'implicit') then

          lambda   = 1.0_dp
          var%Amtx = 0.0_dp
          var%Jmtx = 0.0_dp
          var%barr = 0.0_dp
          var%xarr = 0.0_dp
          do iz = 1, grd%nz

            ! generate Chemical Jacobian matrix
            !   Jmtx = d(Pi-Li)/dni
            var%Jmtx = 0.0_dp
            do ilist = 1, spl%n_Jlist
              i = spl%Jmtx_list(ilist,1)
              j = spl%Jmtx_list(ilist,2)
              ii = spl%var_to_all(i)
              jj = spl%var_to_all(j)
              do ich = 1, spl%Jmtx_list(ilist,3)
                ijch = spl%Jmtx_list(ilist,2+2*ich)
                a0 = dble(spl%Jmtx_list(ilist,2+2*ich+1))
                kij = var%ki(ijch,iz)
                tmp = a0 * kij

                bij = 0
                do isp = 2, spl%reactant_list(ijch,1)+1
                  k = spl%reactant_list(ijch,isp)
                  if (jj /= k) then
                    nk = var%ni(k,iz)
                    tmp = tmp * var%ni(k,iz)
                  else if (jj == k) then
                    bij = bij + 1
                  end if
                end do
                nj = var%ni(jj,iz)
                tmp = tmp * dble(bij) * nj**dble(bij-1)

                var%Jmtx(i,j) = var%Jmtx(i,j) + tmp

              end do
            end do

            do isp = 1, spl%nsp_i
              do jsp = 1, spl%nsp_i
                i = (iz-1)*spl%nsp_i+isp
                j = (iz-1)*spl%nsp_i+jsp + spl%nsp_i + 1 - i
                var%Amtx(i,j) = - var%Jmtx(isp,jsp)
              end do
            end do

            ! vector b = P - L - dPhi_dz 
            do isp = 1, spl%nsp_i
              i = (iz-1)*spl%nsp_i+isp
              var%barr(i) = ( var%Pi(isp,iz) - var%Li(isp,iz) ) - var%dPhi_dz(isp,iz)
            end do

          end do ! z

          ! matrix A = ( I/dt - λ*J )
          !   Combining Chemical Jacobian and Transport Jacobian term
          do iz = 1, grd%nz
            do isp = 1, spl%nsp_i
              i = (iz-1)*spl%nsp_i+isp
              j = (iz-1)*spl%nsp_i+isp + spl%nsp_i + 1 - i
              var%Amtx(i,j) = var%Amtx(i,j) + var%d_dni0_dPhi_dz(isp,iz)
            end do
          end do
          do iz = 2, grd%nz
            do isp = 1, spl%nsp_i
              i = (iz-1)*spl%nsp_i+isp
              j = (iz-2)*spl%nsp_i+isp + spl%nsp_i + 1 - i
              var%Amtx(i,j) = var%d_dnil_dPhi_dz(isp,iz)
              i = (iz-2)*spl%nsp_i+isp
              j = (iz-1)*spl%nsp_i+isp + spl%nsp_i + 1 - i
              var%Amtx(i,j) = var%d_dniu_dPhi_dz(isp,iz-1)
            end do
          end do

          if ( model == 'Catling' ) then
            do iz = 1, grd%nz
              do isp = 1, spl%nsp_i
                i = (iz-1)*spl%nsp_i+isp
                j = (iz-1)*spl%nsp_i+isp + spl%nsp_i + 1 - i
                var%Amtx(i,j) = var%Amtx(i,j) + 1.0_dp / var%dtime
              end do
            end do
          end if

          if ( model == 'Chaffin' ) then
            do iz = 1, grd%nz
              do isp = 1, spl%nsp_i
                i = (iz-1)*spl%nsp_i+isp
                j = (iz-1)*spl%nsp_i+isp + spl%nsp_i + 1 - i
                jsp = spl%var_to_all(isp)
                var%Amtx(i,j) = var%Amtx(i,j) + 1.0_dp / var%dtime
                var%barr(i) = var%barr(i) - ( var%ni(jsp,iz) - var%ni_0(jsp,iz) ) / var%dtime
              end do
            end do
          end if

          ! Solve 'A x = b' using LU decomposition
          CALL p__LU_CMP_CHEM_DIFF( spl, grd, var )
          CALL p__LU_SLV_CHEM_DIFF( spl, grd, var )
          !CALL p__LU_IMPROVE(spl, grd, var)

          !do iz = 1, grd%nz
          !  do isp = 1, spl%nsp_i
          !    i = (iz-1)*spl%nsp_i+isp
          !    jsp = spl%var_to_all(isp)
          !    var%ni_new(jsp,iz) = var%ni(jsp,iz) + var%xarr(i)
          !  end do
          !end do

          do iz = 1, grd%nz
            do isp = 1, spl%nsp
              var%ni_new(isp,iz) = 0.0_dp
              if (spl%label_fix(isp) == 0) then
                var%ni_new(isp,iz) = var%ni(isp,iz) &
                  &                + var%xarr((iz-1)*spl%nsp_i+spl%all_to_var(isp))
              else if (spl%label_fix(isp) == 1) then
                var%ni_new(isp,iz) = var%ni(isp,iz)
              end if
            end do
          end do

        end if
        

    ! electron density : charge neutrality ---------------------------------------------------------
      do iz = 1, grd%nz
        do isp = 1, spl%nsp
          if (spl%species(isp) == 'e-') then
            var%ni_new(isp,iz) = 0.0_dp
            do jsp = 1, spl%nsp
              if (var%q(jsp) /= 0.0_dp .and. jsp /= isp) then
                var%ni_new(isp,iz) &
                  &  = var%ni_new(isp,iz) + var%ni_new(jsp,iz)
              end if
            end do
          end if
        end do
      end do

    ! advance time step ------------------------------------------------------------------
      do iz = 1, grd%nz
        do isp = 1, spl%nsp
          if (var%ni_new(isp,iz) < 0.0_dp) then
            var%ni_new(isp,iz) = 0.0_dp
          end if
          if (var%ni_new(isp,iz) < 1.0e-20_dp) then
            var%ni_new(isp,iz) = 1.0e-20_dp
          end if
        end do
      end do
      ! lower boundary condition: density
      do isp = 1, spl%nsp
        if (nint(var%LowerBC(isp,1)) == 1.0_dp) then 
          var%ni_new(isp,1) = var%LowerBC(isp,2)
        end if
      end do

    ! convergence check and iuncrease time step
      if (set%calc_stable == 1 .and. set%start_rot == 0) then
        if ( var%istep >= 2 ) then

          var%max_dn_n = 0.0_dp
          do iz = 1, grd%nz
            do isp = 1, spl%nsp
              ! delta_N / N
              if ( var%ni(isp,iz) /= 0.0d0 ) then
                dni_ni(isp,iz) = dabs(var%ni_new(isp,iz) - var%ni(isp,iz)) / var%ni(isp,iz)
              else if ( var%ni(isp,iz) == 0.0d0 ) then
                dni_ni(isp,iz) = 0.0d0
              end if
              if ( dni_ni(isp,iz) > var%max_dn_n(3) ) then
                var%max_dn_n(1) = isp
                var%max_dn_n(2) = iz
                var%max_dn_n(3) = dni_ni(isp,iz)
                tmp = dble(isp)
              end if
            end do
          end do

          ! if max(dN/N) < 10 %, then goes to larger dt (*10)
          ! if dt becomes 100000 [sec], unstable at HC layer
          if ( var%dtime < set%dtime_limit * 0.99_dp ) then
            if ( var%max_dn_n(3) < eps ) then
              if (model == 'Catling') var%dtime = var%dtime * 1.1e0_dp
              if (model == 'Chaffin') var%dtime = var%dtime * 1.0e1_dp
              var%iter = 0
              do iz = 1, grd%nz
                do isp = 1, spl%nsp
                  var%ni_0(isp,iz) = var%ni_new(isp,iz)
                end do
              end do
            else if ( var%max_dn_n(3) >= eps ) then
              var%iter = var%iter + 1
              if ( var%iter > 5000 ) then
                var%dtime = 1.0e-8_dp
                var%iter = 0
                var%ni_new(isp,iz) = var%ni_0(isp,iz)
              end if
            end if
          else if ( var%dtime  >= set%dtime_limit * 0.99_dp ) then
            var%dtime = set%dtime_limit ! dt DO NOT excess set%dtime_limit <- set at v__'planet'__ini
          end if

        end if
      end if

    ! advance time step ------------------------------------------------------------------
      do iz = 1, grd%nz
        do isp = 1, spl%nsp
          if (var%ni_new(isp,iz) /= var%ni_new(isp,iz)) then
            var%ni(isp,iz) = var%ni_0(isp,iz)
            var%ni_new(isp,iz) = var%ni_0(isp,iz)
            var%dtime = 1.0e-8_dp
            var%iter = 0
            !print *, iz, isp
            !call p__io_stable__fin(set, spl, var, grd) ! in
            !stop
          end if
          !var%ni_0(isp,iz) = var%ni(isp,iz)
          var%ni(isp,iz) = var%ni_new(isp,iz)
        end do
      end do




  end subroutine p__photochem_scheme__exe


  !--------------------------------------------------------------------
  ! LU decomposition : flux implicit
  !--------------------------------------------------------------------
  subroutine p__LU_CMP_CHEM_DIFF( spl, grd, var )
    implicit none
    type(spl_), intent(in)    :: spl
    type(grd_), intent(in)    :: grd
    type(var_), intent(inout) :: var

    integer i, j, k
    integer N, nsp, nz
    double precision sum, max, ipv

    nsp = spl%nsp_i
    nz  = grd%nz
    N   = nsp * nz

    var%Lmtx = 0.0d0
    var%Umtx = 0.0d0

    do i = 1, N
      var%Lmtx(i,i+nsp+1-i) = 1.0d0
    end do

    do j = 1, N
      do i = j-nsp, j
        if (i>=1) then
          sum = 0.0d0
          if (i>=2) then
            do k = i-nsp, i-1
              if ( k>=1 .and. k<=N &
                & .and. k+nsp+1-i>=1 .and. k+nsp+1-i<=2*nsp+1 &
                & .and. j+nsp+1-k>=1 .and. j+nsp+1-k<=2*nsp+1 ) then
                sum = sum + var%Lmtx(i,k+nsp+1-i) * var%Umtx(k,j+nsp+1-k)
              end if
            end do
          end if
          var%Umtx(i,j+nsp+1-i) = var%Amtx(i,j+nsp+1-i) - sum
        end if
      end do

      do i = j+1, j+nsp
        if (i<=N) then
          sum = 0.0d0
          do k = j-nsp, j-1
            if ( k>=1 .and. k<=N &
            & .and. k+nsp+1-i>=1 .and. k+nsp+1-i<=2*nsp+1 &
            & .and. j+nsp+1-k>=1 .and. j+nsp+1-k<=2*nsp+1 ) then
              sum = sum + var%Lmtx(i,k+nsp+1-i) * var%Umtx(k,j+nsp+1-k)
            end if
          end do
          var%Lmtx(i,j+nsp+1-i) = ( var%Amtx(i,j+nsp+1-i) - sum ) / var%Umtx(j,j+nsp+1-j)
        end if
      end do
    end do

  end subroutine p__LU_CMP_CHEM_DIFF


  ! Solve simultaneous N eqations for N unknowns using LU decomposition
  subroutine p__LU_SLV_CHEM_DIFF( spl, grd, var )
    implicit none
    type(spl_), intent(in)    :: spl
    type(grd_), intent(in)    :: grd
    type(var_), intent(inout) :: var
    integer i, j
    integer N, nsp, nz
    double precision sum

    nsp = spl%nsp_i
    nz  = grd%nz
    N   = nsp * nz

    var%yarr(1) = var%barr(1)
    do i = 2, N
      sum = 0.0d0
        do j = i - nsp, i-1
          if (j>=1) then
            sum = sum + var%Lmtx(i,j+nsp+1-i) * var%yarr(j)
          end if
        end do
        var%yarr(i) = var%barr(i) - sum
    end do

    var%xarr(N) = var%yarr(N) / var%Umtx(N,N+nsp+1-N)
    do i = N-1, 1, -1
      sum = 0.0d0
        do j = i+1, i + nsp
          if (j<=N) then
            sum = sum + var%Umtx(i,j+nsp+1-i) * var%xarr(j)
          end if
        end do
        var%xarr(i) = ( var%yarr(i) - sum )/ var%Umtx(i,i+nsp+1-i)
    end do

  end subroutine p__LU_SLV_CHEM_DIFF


  ! Iteration using LU decomposition
  subroutine p__LU_IMPROVE(spl, grd, var)

    type(spl_), intent(in)    :: spl
    type(grd_), intent(in)    :: grd
    type(var_), intent(inout) :: var

    integer i, j, N, nsp, nz
    double precision sum

    nsp = spl%nsp_i
    nz  = grd%nz
    N   = nsp * nz

    do i = 1, N
      sum = 0.0d0
      do j = i - nsp, i + nsp
        if (j >= 1 .and. j <= N &
          & .and. j+nsp+1-i >= 1 .and. j+nsp+1-i <= 2*nsp+1) then 
          sum = sum + var%Amtx(i,j+nsp+1-i) * var%xarr(j)
        end if 
      end do
      var%rarr(i) = sum - var%barr(i)
    end do

    !CALL p__LU_SLV_CHEM_DIFF( spl, grd, var )
    do i = 1, N
      var%xarr(i) = var%xarr(i) - var%dxarr(i)
    end do

  end subroutine p__LU_IMPROVE




end module p__photochem_scheme
