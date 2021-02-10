module p__io
  use v__tdec, only : set_, var_, spl_, xct_, flx_, grd_

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: p__io_stable__fin, p__io_rotation__fin, p__io_deallocate, &
    &       p__io_progress

contains


  subroutine p__io_stable__fin(set, spl, var, grd)
    type(set_), intent(in) :: set
    type(spl_), intent(in) :: spl
    type(var_), intent(in) :: var
    type(grd_), intent(in) :: grd
    character(len=256) fname, num
    integer i, ip, is, iy, iz, isp, ich

    fname = './'//trim(ADJUSTL(set%dir_name))//'/output/density/n_stable.dat'
    open(11, file = fname, status = 'unknown' )
    do iz = 1, grd%nz
      do iy = 1, grd%ny
        do isp = 1, spl%nsp
          write(11,*) var%ni_stable(isp,iy,iz)
        end do
      end do
    end do
    close(11)

    do isp = 1, spl%nsp
      fname = './'//trim(ADJUSTL(set%dir_name))//'/output/density/num/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
      open(11, file = fname, status = 'unknown' )
        do iz = 1, grd%nz
          write(11, *) grd%alt(iz)/1d3, var%ni(isp,iz)
        end do
      close(11)

      fname = './'//trim(ADJUSTL(set%dir_name))//'/output/density/vmr/vmr_'//trim(ADJUSTL(spl%species(isp)))//'.dat'
      open(11, file = fname, status = 'unknown' )
        do iz = 1, grd%nz
          write(11, *) grd%alt(iz)/1d3, var%ni(isp,iz)/var%n_tot(iz)
        end do
      close(11)

      if (spl%label_fix(isp) == 0) then

        fname = './'//trim(ADJUSTL(set%dir_name))//'/output/flux/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
        open(11, file = fname, status = 'unknown' )
          do iz = 1, grd%nz
            write(11, *) grd%alt(iz)/1d3, &
              & var%Fluxup(spl%all_to_var(isp),iz)/1.0e4_dp, &
              & var%Fluxdwn(spl%all_to_var(isp),iz)/1.0e4_dp
          end do
        close(11)

        fname = './'//trim(ADJUSTL(set%dir_name))//'/output/flux/v_'//trim(ADJUSTL(spl%species(isp)))//'.dat'
        open(11, file = fname, status = 'unknown' )
          do iz = 1, grd%nz
            write(11, *) grd%alt(iz)/1d3, &
              & var%Fluxup(spl%all_to_var(isp),iz)/var%ni(isp,iz), &
              & var%Fluxdwn(spl%all_to_var(isp),iz)/var%ni(isp,iz)
          end do
        close(11)

        fname = './'//trim(ADJUSTL(set%dir_name))//'/output/flux/D_'//trim(ADJUSTL(spl%species(isp)))//'.dat'
        open(11, file = fname, status = 'unknown' )
          do iz = 1, grd%nz
            write(11, *) grd%alt(iz)/1d3, var%D_mol(isp,iz)
          end do
        close(11)

      end if
    end do

    do ich = 1, spl%nch
      write(num,'(I3)') ich
      fname = './'//trim(ADJUSTL(set%dir_name))//'/output/rate/k'//trim(ADJUSTL(num))//'.dat'
      open(11, file = fname, status = 'unknown' )
        do iz = 1, grd%nz
          write(11, *) grd%alt(iz)/1.0d3, var%ki(ich,iz)
        end do
      close(11)
    end do

    do ich = 1, spl%nch
      write(num,'(I3)') ich
      fname = './'//trim(ADJUSTL(set%dir_name))//'/output/rate/rate'//trim(ADJUSTL(num))//'.dat'
      open(11, file = fname, status = 'unknown' )
        do iz = 1, grd%nz
          write(11, *) grd%alt(iz)/1.0d3, var%rate(ich,iz)
        end do
      close(11)
    end do

    do isp = 1, spl%nsp_i
      fname = './'//trim(ADJUSTL(set%dir_name))//'/output/rate/P_'//trim(ADJUSTL(spl%species(isp)))//'.dat'
      open(11, file = fname, status = 'unknown' )
        do iz = 1, grd%nz
          write(11, *) grd%alt(iz)/1d3, var%Pi(spl%all_to_var(isp),iz)
        end do
      close(11)

      fname = './'//trim(ADJUSTL(set%dir_name))//'/output/rate/L_'//trim(ADJUSTL(spl%species(isp)))//'.dat'
      open(11, file = fname, status = 'unknown' )
        do iz = 1, grd%nz
          write(11, *) grd%alt(iz)/1d3, var%Li(spl%all_to_var(isp),iz)
        end do
      close(11)
    end do

    fname = './'//trim(ADJUSTL(set%dir_name))//'/output/flux/K_eddy.dat'
    open(11, file = fname, status = 'unknown' )
      do iz = 1, grd%nz
        write(11, *) grd%alt(iz)/1d3, var%K_eddy(iz)
      end do
    close(11)


  end subroutine p__io_stable__fin

  subroutine p__io_rotation__fin(set, spl, var, grd)
    type(set_), intent(in) :: set
    type(spl_), intent(in) :: spl
    type(var_), intent(in) :: var
    type(grd_), intent(in) :: grd
    character(len=256) fname, num
    integer i, ip, is, ix, iy, iz, isp

    do isp = 1, spl%nsp
      fname = './'//trim(ADJUSTL(set%dir_name))//'/output/density/num/rot_'//trim(ADJUSTL(spl%species(isp)))//'.dat'
      open(11, file = fname, status = 'unknown' )
        do ix = 1, grd%nx
          do iy = 1, grd%ny
            do iz = 1, grd%nz
              write(11, *) ix, iy, grd%alt(iz)/1d3, var%ni_3d(isp,ix,iy,iz)
            end do
          write(11,*)
          end do
        end do
      close(11)

      do iz = 1, grd%nz
        write(num,'(I3)') iz
        fname = './'//trim(ADJUSTL(set%dir_name))//'/output/density/global/global_' &
          &      //trim(ADJUSTL(num))//'_'//trim(ADJUSTL(spl%species(isp)))//'.dat'
        open(11, file = fname, status = 'unknown' )
          do iy = 1, grd%ny
            do ix = 1, grd%nx
              write(11, *) 180_dp*dble(iy-((grd%ny+1)/2))/dble(grd%ny-1), &
                &          360_dp*dble(ix-1)/dble(grd%nx-1), var%ni_3d(isp,ix,iy,iz)
            end do
            write(11,*)
          end do
        close(11)
      end do


      if (spl%label_fix(isp) == 0) then

!        fname = './'//trim(ADJUSTL(spl%planet))//'/output/flux/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
!        open(11, file = fname, status = 'unknown' )
!          do iz = 1, grd%nz
!            write(11, *) grd%alt(iz)/1d3, &
!              & var%Fluxup(spl%all_to_var(isp),iz), &
!              & var%Fluxdwn(spl%all_to_var(isp),iz)
!          end do
!        close(11)
!
!        fname = './'//trim(ADJUSTL(spl%planet))//'/output/flux/v_'//trim(ADJUSTL(spl%species(isp)))//'.dat'
!        open(11, file = fname, status = 'unknown' )
!          do iz = 1, grd%nz
!            write(11, *) grd%alt(iz)/1d3, &
!              & var%Fluxup(spl%all_to_var(isp),iz)/var%ni(isp,iz), &
!              & var%Fluxdwn(spl%all_to_var(isp),iz)/var%ni(isp,iz)
!          end do
!        close(11)
!
!        fname = './'//trim(ADJUSTL(spl%planet))//'/output/flux/D_'//trim(ADJUSTL(spl%species(isp)))//'.dat'
!        open(11, file = fname, status = 'unknown' )
!          do iz = 1, grd%nz
!            write(11, *) grd%alt(iz)/1d3, var%D_mol(isp,iz)
!          end do
!        close(11)
!
      end if
    end do

    fname = './'//trim(ADJUSTL(set%dir_name))//'/output/flux/K_eddy.dat'
    open(11, file = fname, status = 'unknown' )
      do iz = 1, grd%nz
        write(11, *) grd%alt(iz)/1d3, var%K_eddy(iz)
      end do
    close(11)


  end subroutine p__io_rotation__fin


  subroutine p__io_progress(spl, var, grd, set)
    type(spl_), intent(in) :: spl
    type(var_), intent(in) :: var
    type(grd_), intent(in) :: grd
    type(set_), intent(in) :: set
    character(len=256) fname

    fname = './progress.dat'
    open(11, file = fname, status = 'unknown' )

      write(11, *) '########################################'
      write(11, *)
      write(11, *) '            ',trim(ADJUSTL(spl%planet))
      write(11, *)
      write(11, *) '########################################'
      write(11, *)
      write(11, *)

      if ( set%mode == '1D' ) then

        write(11, *) '----------------------------------------'
        write(11, *) '         1D stable calculation          '
        write(11, *) '----------------------------------------'
        write(11, *)

        if ( set%calc_stable == 0 ) then
          write(11, *) 'skipped stable solution...'
        else if ( set%calc_stable == 1 ) then
          write(11, *) 'calculate stable solution...'
          if ( grd%iy < grd%ny ) then
            write(11, *) 'current progress = ', int(100 * (grd%iy-1)/grd%ny), ' %'
          else if ( grd%iy == grd%ny ) then
            write(11, *) 'finished stable calculation!'
          end if
          write(11,*) 'dt     = ',var%dtime, ' [s]'
          write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
          write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
        end if

      else if ( set%mode == '2D Lat' ) then

        write(11, *) '----------------------------------------'
        write(11, *) '         2D stable calculation          '
        write(11, *) '----------------------------------------'
        write(11, *)

        if ( set%calc_stable == 0 ) then
          write(11, *) 'skipped stable solution...'
        else if ( set%calc_stable == 1 ) then
          write(11, *) 'calculate stable solution...'
          if ( grd%iy < grd%ny ) then
            write(11, *) 'current progress = ', int(100 * (grd%iy-1)/grd%ny), ' %'
          else if ( grd%iy == grd%ny ) then
            write(11, *) 'finished stable calculation!'
          end if
          write(11,*) 'dt     = ',var%dtime, ' [s]'
          write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
          write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
        end if

      else if (set%mode == '2D Rot') then
        write(11, *) '----------------------------------------'
        write(11, *) '       2D rotational calculation     '
        write(11, *) '----------------------------------------'
        write(11, *)

        if ( set%calc_stable == 0 ) then
          write(11, *) 'skipped stable solution...'
          if ( set%start_rot == 1 ) then
            write(11, *) 'calculating rotational mode...'
            if ( grd%iy < grd%ny ) then
              write(11, *) 'current progress = ', int(100 * (grd%iy-1)/grd%ny), ' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished rotational calculation!'
            end if
            write(11,*) 'cpu_t  = ',var%t4-var%t3, ' [s]'
          end if
        else if ( set%calc_stable == 1 ) then
          if ( set%start_rot == 0 ) then
            write(11, *) 'calculating stable solution...'
            if ( grd%iy < grd%ny ) then
              write(11, *) 'current progress = ', int(100 * (grd%iy-1)/grd%ny), ' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished stable calculation!'
            end if
            write(11,*) 'dt     = ',var%dtime, ' [s]'
            write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
            write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
            write(11,*)
          end if
          if ( set%start_rot == 1 ) then
            write(11, *) 'finished stable calculation!'
            write(11,*) 'dt     = ',var%dtime, ' [s]'
            write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
            write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
            write(11,*)
            if ( grd%iy < grd%ny ) then
              write(11, *) 'calculating rotational mode...'
              write(11, *) 'current progress = ', int(100 * (grd%iy-1)/grd%ny), ' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished rotational calculation!'
            end if
            write(11,*) 'cpu_t  = ',var%t4-var%t3, ' [s]'
          end if
        end if

      else if (set%mode == '3D Rot') then
        write(11, *) '----------------------------------------'
        write(11, *) '       3D rotational calculation     '
        write(11, *) '----------------------------------------'
        write(11, *)

        if ( set%calc_stable == 0 ) then
          write(11, *) 'skipped stable solution...'
          if ( set%start_rot == 1 ) then
            write(11, *) 'calculating rotational mode...'
            if ( grd%iy < grd%ny ) then
              write(11, *) 'current progress = ', int(100 * (grd%iy-1)/grd%ny), ' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished rotational calculation!'
            end if
            write(11,*) 'cpu_t  = ',var%t4-var%t3, ' [s]'
          end if
        else if ( set%calc_stable == 1 ) then
          if ( set%start_rot == 0 ) then
            write(11, *) 'calculating stable solution...'
            if ( grd%iy < grd%ny ) then
              write(11, *) 'current progress = ', int(100 * (grd%iy-1)/grd%ny), ' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished stable calculation!'
            end if
            write(11,*) 'dt     = ',var%dtime, ' [s]'
            write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
            write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
            write(11,*)
          end if
          if ( set%start_rot == 1 ) then
            write(11, *) 'finished stable calculation!'
            write(11,*) 'dt     = ',var%dtime, ' [s]'
            write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
            write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
            write(11,*)
            if ( grd%iy < grd%ny ) then
              write(11, *) 'calculating rotational mode...'
              write(11, *) 'current progress = ', int(100 * (grd%iy-1)/grd%ny), ' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished rotational calculation!'
            end if
            write(11,*) 'cpu_t  = ',var%t4-var%t3, ' [s]'
          end if
        end if

      else if ( set%mode == '3D Global' ) then
        write(11, *) '----------------------------------------'
        write(11, *) '         3D global calculation          '
        write(11, *) '----------------------------------------'
        write(11, *)

        if ( set%calc_stable == 0 ) then
          write(11, *) 'skipped stable solution...'
          if ( set%start_rot == 1 ) then
            write(11, *) 'calculating rotational mode...'
            if ( grd%iy < grd%ny ) then
              write(11, *) 'current progress = ', int(100 * (grd%iy-1)/grd%ny), ' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished rotational calculation!'
            end if
            write(1,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
          end if
        else if ( set%calc_stable == 1 ) then
          write(11, *) 'calculating stable solution...'
          if ( grd%iy < grd%ny ) then
            write(11, *) 'current progress = ', int(100 * (grd%iy-1)/grd%ny), ' %'
          else if ( grd%iy == grd%ny ) then
            write(11, *) 'finished stable calculation!'
          end if
          write(11,*) 'dt     = ',var%dtime, ' [s]'
          write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
          write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
          write(11,*)
          if ( set%start_rot == 1 ) then
            write(11, *) 'calculating rotational mode...'
            if ( grd%iy < grd%ny ) then
              write(11, *) 'current progress = ', int(100 * (grd%iy-1)/grd%ny), ' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished rotational calculation!'
            end if
            write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
          end if
        end if

      end if


    close(11)


  end subroutine p__io_progress


  subroutine p__io_deallocate(var, spl, xct, flx) ! inout
    type(var_),           intent(inout) :: var
    type(spl_),           intent(inout) :: spl
    type(xct_),           intent(inout) :: xct
    type(flx_),           intent(inout) :: flx
    integer i, ip, isp, jsp, ich, jch, iz, nh
    real(dp) tmp
    character(len = 256) strm, fname


    ! deallocate
    deallocate(var%ni)
    deallocate(var%ni_stable)
    deallocate(var%ni_3d)
    deallocate(var%ni_new)
    deallocate(var%clm_ni)
    deallocate(var%Ti,var%Te,var%Tn)
    deallocate(var%Ti_3d,var%Te_3d,var%Tn_3d)
    deallocate(var%m, var%q)
    deallocate(spl%reactant_list)
    deallocate(spl%product_list)
    deallocate(spl%species)
    deallocate(spl%label_fix)
    deallocate(spl%all_to_var)
    deallocate(spl%var_to_all)
    deallocate(spl%Prod_list)
    deallocate(spl%Loss_list)
    deallocate(spl%Jmtx_list)
    deallocate(spl%reaction_type_list)
    deallocate(var%ki)
    deallocate(var%Pi)
    deallocate(var%Li)
    deallocate(var%K_eddy,var%D_mol)
    deallocate(var%Fluxup,var%Fluxdwn)
    deallocate(var%Jmtx)
    deallocate(spl%rate_rpn_token, spl%rate_rpn_label)
    deallocate(spl%rate_cases)
    deallocate(spl%T_range)
    deallocate(spl%major_species)
    deallocate(var%n_tot,var%m_mean)

    deallocate(flx%F74113,flx%Ai,flx%solar_EUV)
    deallocate(var%tau_EUV,var%I_EUV)
    deallocate(xct%sigma_a_EUV, xct%sigma_i_EUV)
    deallocate(flx%solar_UV)
    deallocate(var%tau_UV,var%I_UV)
    deallocate(xct%sigma_a_UV, xct%sigma_d_UV)


  end subroutine p__io_deallocate


end module p__io