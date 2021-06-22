program MCD_density

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  character(len = 256) fname, species(10)
  integer, parameter :: nz = 35
  integer, parameter :: nz_new = 16
  integer, parameter :: nsp = 2
  integer isp, iz, jz, nh, i
  real(dp) mass(nsp), NA, tmp, vmr
  real(dp) alt(nz), den(nz), T_n(nz), Ni(nsp,nz)
  real(dp) dalt_new(nz_new), alt_new(nz_new), den_new(nz_new), T_n_new(nz_new), Ni_new(nsp,nz_new)

  species(1) = 'CO2'; mass(1) = 44_dp
  !species(2) = 'N2' ; mass(2) = 28_dp
  !species(3) = 'H2' ; mass(3) = 2_dp
  !species(4) = 'O2' ; mass(4) = 32_dp
  !species(5) = 'CO' ; mass(5) = 28_dp
  species(2) = 'O'  ; mass(2) = 16_dp

  NA = 6.022e23_dp

  dalt_new(1:nz_new) = 5.0e3_dp ! [km]
  alt_new(1)      = 122.5e3_dp ! [km]
  do iz = 2, nz_new
    alt_new(iz) = alt_new(iz-1) + dalt_new(iz)
  end do


  fname = './density.dat'
  call p__count_header(nh, fname)
  open(11, file = fname, status = 'unknown' )
    do i = 1, nh; read(11,*); enddo
    do iz = 1, nz
      read(11, *) alt(iz), den(iz)
    end do
  close(11)

  do isp = 1, nsp
    fname = './vmr_'//trim(species(isp))//'.dat'
    call p__count_header(nh, fname)
    open(11, file = fname, status = 'unknown' )
      do i = 1, nh; read(11,*); enddo
      do iz = 1, nz
        read(11, *) tmp, vmr
        Ni(isp,iz) = den(iz) * vmr * 1.0e3_dp / mass(isp) * NA
      end do
    close(11)
  end do

  !fname = './T_n_MCD.dat'
  !call p__count_header(nh, fname)
  !open(11, file = fname, status = 'unknown' )
  !  do i = 1, nh; read(11,*); enddo
  !  do iz = 1, nz
  !    read(11, *) tmp, T_n(iz)
  !  end do
  !close(11)

  do iz = 1, nz
    write(*,*) alt(iz)/1000_dp, Ni(1,iz), Ni(2,iz)
  end do
  write(*,*)

  ! atmospheric density at the location of the particles
  do isp = 1, nsp
    do iz = 1, nz-1
      do jz = 1, nz_new
        if (    alt_new(jz) >= alt(iz) &
        & .and. alt_new(jz) <  alt(iz+1) ) then

          tmp = dlog(Ni(isp,iz)) &
            & - ( dlog(Ni(isp,iz)) - dlog(Ni(isp,iz+1)) ) &
            & * ( alt_new(jz) - alt(iz) )/( alt(iz+1) - alt(iz) )
          Ni_new(isp,jz) = dexp(tmp)

          tmp = dlog(T_n(iz)) &
            & - ( dlog(T_n(iz)) - dlog(T_n(iz+1)) ) &
            & * ( alt_new(jz) - alt(iz) )/( alt(iz+1) - alt(iz) )
          T_n_new(jz) = dexp(tmp)

        end if

      end do
    end do
  end do

  do iz = 1, nz_new
    write(*,*) alt_new(iz)/1000_dp, Ni_new(1,iz), Ni_new(2,iz)
  end do

  do isp = 1, nsp
    fname = './'//trim(species(isp))//'.dat'
    open(11, file = fname, status = 'unknown' )
      do iz = 1, nz_new
        write(11, *) alt_new(iz), Ni_new(isp,iz)
      end do
    close(11)
  end do

  !fname = './T_n.dat'
  !open(11, file = fname, status = 'unknown' )
  !  do iz = 1, nz_new
  !    write(11, *) alt_new(iz), T_n_new(iz)
  !  end do
  !close(11)


end program MCD_density

! counting headers
subroutine p__count_header(nh, fname)
  integer nh
  character(len=256) fname
  character(len=256) strm
  !Skip header lines
  open(11, file = fname, status = 'unknown' )
    nh = 0
    do
      read(11,'(A)') strm
      if(strm(1:1) == "#")then
        nh = nh + 1
      else
        exit
      endif
    enddo
  close(11)
end subroutine p__count_header
