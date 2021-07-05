
subroutine g_int(N, M, Mdot, j_para_i, b, db_dL, g_x)
  implicit none
  integer N, M, i
  double precision g_x(0:2*N*M*240)
  double precision j_para_i(0:2*N*240)
  double precision b(0:2*N*M*240)
  double precision db_dL(0:2*N*M*240)
  double precision K
  double precision dx, xi
  double precision Mdot, Rp, Bp, tmp
  double precision g1, g2, g3
  double precision, parameter :: pi = dacos(-1.0d0)

  Rp = 7.1492d+7 ! planetary radius [m]
  Bp = 4.2d-4 ! magnetic field at equator surface [T] = [kg/A/s^2]
  ! Sigma: [A^2 s^3 /kg/m^2]
  ! Rp^2 Bp^2 Sigma / Mdot = m^2 kg^2 A^-2 s^-4 A^2 s^3 kg^-1 m^-2 kg^-1 s = [no unit]

  tmp = 4.0d0 * pi * Rp**2.0d0 * Bp**2.0d0 / Mdot

  g_x = 0.0d0
  do i = 0, 2 * N * M * 240 - 1
    xi = 1.0d0 + dble(i) / dble(2*N*M)
    dx = 1.0d0 / dble(2*N*M)
    g1 = K(N, M, i, j_para_i, b, db_dL, xi)
    g2 = K(N, M, i, j_para_i, b, db_dL, xi + dx/2.0d0)
    g3 = K(N, M, i, j_para_i, b, db_dL, xi + dx)
    g_x(i+1) = g_x(i) + tmp * (g1 + 4.0d0 * g2 * g3) * dx / 6.0d0
  end do

end subroutine g_int

double precision function K(N, M, j, j_para_i, b, db_dL, x)
  implicit none
  integer N, M, j, i
  double precision Sigma_j, jpara, jpara0, jpara1, Bze, Bp
  double precision j_para_i(0:2*N*240)
  double precision b(0:2*N*M*240)
  double precision db_dL(0:2*N*M*240)
  double precision x

  Bp = 4.2d-4
  
  !b(j) = 1.0d0/x
  !db_dL(j) = -1.0d0/x/x

  i = int( dble(j)/dble(2*M) )
  jpara0 = j_para_i(i)
  jpara1 = j_para_i(i+1)
  jpara  = (jpara0 * dble((i+1)*2*M - j) + jpara1 * dble(j - i*2*M)) / dble(2*M)
  if (jpara < 0.0d0) jpara = 0.0d0

  Bze = 3.335d-4 / x**3.0d0 * dexp( -(x/14.501d0)**2.5d0 ) + 5.4d-5 / x**2.71d0

  K =  - Sigma_j(jpara) * b(j) / x**2.0d0 * dsqrt(1.0d0 - b(j)) * db_dL(j)


  ! Nichols and Cowley 2004
  K =  2.0d0 * Sigma_j(jpara) * b(j) * Bze / Bp / x
  
  return
end function K

double precision function Sigma_j(jpara)
    implicit none
    double precision jpara, x
    double precision a, b, c
    ! Metal x1.0
      a = 2.23295d0       
      b = 1.09797d0       
      c = 345.821d0   
    !! Metal x0.1
    !  a = 1.44487d0   
    !  b = 1.05008d0   
    !  c = 242.836d0   
    !! NO Metal
    !  a = 1.28267d0        
    !  b = 1.00612d0        
    !  c = 173.301d0       

    x = jpara*1.0d6
    if (x < 0.0d0) x = 0.0d0

    Sigma_j = a * dlog(b+c*x*x) * 0.5d0

    return
end function Sigma_j


double precision function Sigma(x)
    implicit none
    double precision x
    double precision a, b, c
    ! Metal x1.0
      a = 2.23295d0       
      b = 1.09797d0       
      c = 345.821d0   
    !! Metal x0.1
    !  a = 1.44487d0   
    !  b = 1.05008d0   
    !  c = 242.836d0   
    !! NO Metal
    !  a = 1.28267d0        
    !  b = 1.00612d0        
    !  c = 173.301d0       

    Sigma = a * dlog(b+c*x*x)

    return
end function Sigma


subroutine calc_j_para_i( N, M, Mdot, f, b, db_dL, j_para_i, j_para_i_new)
  implicit none
  integer N, M, i, j, i0, i1, iL0, iL1
  double precision j_para_i(0:2*N*240)
  double precision j_para_i_new(0:2*N*240)
  double precision b(0:2*N*M*240)
  double precision db_dL(0:2*N*M*240)
  double precision Sigma, Sigma_j
  double precision f(2*N*240)
  double precision K
  double precision dx, xi
  double precision Mdot, Rp, Bp, tmp
  double precision L, dL, Bze, OmegaJ, dSpFef_dL, fj1, fj0, f1, f0, jpara0, jpara1
  double precision b0, b1, bj, Sigma0, Sigma1, theta0, theta1, dtheta, LK, Sj, df_dLj
  double precision, parameter :: pi = dacos(-1.0d0)

  Rp = 7.1492d+7 ! planetary radius [m]
  Bp = 4.2d-4 ! magnetic field at equator surface [T] = [kg/A/s^2]
  ! Sigma: [A^2 s^3 /kg/m^2]
  ! Rp^2 Bp^2 Sigma / Mdot = m^2 kg^2 A^-2 s^-4 A^2 s^3 kg^-1 m^-2 kg^-1 s = [no unit]

  OmegaJ = 2.0d0 * pi / 35729.0d0

  j_para_i(0)     = 0.0d0
  j_para_i(100*N-1:110*N) = -0.04d-6

  do i = 1, N * 110 - 1
    L   = 1.0d0 + (dble(i)+0.5d0) / dble(N)
    dL  = 1.0d0 / dble(N)

    Bze = dabs(3.335d-4 / L**3.0d0 * dexp( -(L/14.501d0)**2.5d0 ) + 5.4d-5 / L**2.71d0)



    f0 = f(i) 
    f1 = f(i+1)

    jpara0 = (j_para_i(i-1)+j_para_i(i  ))/2.0d0
    jpara1 = (j_para_i(i  )+j_para_i(i+1))/2.0d0

    dSpFef_dL = ( Sigma_j(jpara1) * b((i+1)*2*M) * Bp * f1 &
      &         - Sigma_j(jpara0) * b((i  )*2*M) * Bp * f0 )  / dL
! 
!    dSpFef_dL = ( b((i+1)*2*M) * Bp * f1 &
!      &         - b((i  )*2*M) * Bp * f0 ) * Sigma_j((jpara0+jpara1)/2.0d0) / dL

    j_para_i_new(i) = 4.0d0 * Bp * OmegaJ / (pi * L * Bze) * dSpFef_dL
    

!    bj  = b(i*2*M+M)
!    Sj  = Sigma_j((jpara0+jpara1)/2.0d0)
!    LK  = 4.0d0 * pi * Rp**2.0d0 * Bp**2.0d0 / Mdot &
!      & * 2.0d0 * Sj * bj * Bze / Bp
!    df_dLj = (f1 - f0) / dL
!
!    j_para_i_new(i) = 4.0d0 * Bp * OmegaJ / (pi * L * Bze) &
!      &             * (   2.0d0 / (2.0d0 + LK) / df_dLj * ( Sigma_j(jpara1) * b((i+1)*2*M) * Bp &
!      &                                                   - Sigma_j(jpara0) * b((i  )*2*M) * Bp ) &
!      &               + Sj * bj * Bp * df_dLj )


!    ! Hill 2000
!
!    f0 = f(i) 
!    f1 = f(i+1)
!
!    jpara0 = (j_para_i(i-1)+j_para_i(i  ))/2.0d0
!    jpara1 = (j_para_i(i  )+j_para_i(i+1))/2.0d0
!
!    b0 = b((i  )*2*M)
!    b1 = b((i+1)*2*M)
!    bj = (b0+b1)/2.0d0
!
!    Sigma0 = Sigma_j(jpara0)
!    Sigma1 = Sigma_j(jpara1)
!
!    theta0 = dasin(dsqrt(b0))
!    theta1 = dasin(dsqrt(b1))
!    dtheta = theta1 - theta0
!
!    j_para_i_new(i) = - 2.0d0 * Bp * OmegaJ / dsqrt(bj) &
!      &             * ( b1 * dsqrt(1.0d0-b1) * Sigma1 * f1 &
!      &               - b0 * dsqrt(1.0d0-b0) * Sigma0 * f0 ) &
!      &             / dtheta

    if (mod(i,N) == 0) print *, L, j_para_i_new(i)*1.0d6, ' uA/m2'
  end do

  j_para_i_new(100*N-1:110*N) = -0.04d-6

end subroutine calc_j_para_i


!begin calculation of Hill current

program plasma_velocity
  implicit none

  !parameter

  double precision, parameter :: pi = dacos(-1.0d0)

  double precision x              !x
  double precision dx             !integration grid size : 1 / M
  double precision L              !distance normalized at RJ
  double precision, parameter :: L0 = 23.0d0 !corotation lag become significant at L0
  double precision f_u            !numerator of f(x)
  double precision farr(10,5000)       !f = dw/wp : w/wp = 1 - f : corotation lag
  double precision f1             !g(x)
  double precision f2             !g(x+dx/2)
  double precision f3             !g(x+dx)
  double precision theta
  double precision Mdot, tmp, LT
  double precision, parameter :: Rp = 7.1492d+7

  integer, parameter :: N = 1   !grid number of L
  integer, parameter :: M = 100000  !grid number of x in 1 L grid
  integer, parameter :: tg = 10   !grid number of theta

  double precision f(2*N*240)       !f = dw/wp : w/wp = 1 - f : corotation lag

  double precision, allocatable :: g_x(:)
  double precision, allocatable :: Sigma(:)
  double precision, allocatable :: b(:)
  double precision, allocatable :: db_dL(:)
  double precision, allocatable :: L_theta(:)
  double precision, allocatable :: Lb(:)
  double precision, allocatable :: j_para_i(:)
  double precision, allocatable :: j_para_i_new(:)
  double precision, allocatable :: arr(:)
  double precision jarr(10,0:5000)

  double precision Sigma_j

  double precision theta_to_L(181)


  integer i   !a grid number of L
  integer j   !a grid number of x
  integer ix, iy, iB, iL, iter
  integer, parameter :: nx = 361
  integer, parameter :: ny = 181

  double precision hi_sigma_P(nx,ny), hi_sigma_H(nx,ny), hi_sigma_0(nx,ny)

  character(len=256) fname, num, char, dir_name

  !real(8) gamain, uicg_p1, uicg, alpha
  !real(8) Fe, be
  !integer ( kind = 4 ) ifault

  allocate(g_x(0:2*N*M*240))
  allocate(Sigma(0:2*N*M*240))
  allocate(b(0:2*N*M*240))
  allocate(db_dL(0:2*N*M*240))
  allocate(L_theta(0:2*N*M*240))
  allocate(Lb(0:2*N*M*240))
  allocate(j_para_i(0:2*N*240))
  allocate(j_para_i_new(0:2*N*240))
  allocate(arr(2*N*M*240+1))



  ! Parameter settings ----------------------------------------
  dir_name = "./metal_Hill"
  Mdot = 1000.0d0
  LT   = 12.0d0



  ! Read mapping function ------------------------------
  print *, 'Reading b function and mapping data...'
  write(char,*) M*N
  fname = './b_function/M'//trim(ADJUSTL(char))//'/b.dat'
  open(11, file = fname, status = 'unknown' )
    do iL = 0, 2 * M * 240 * N
      read(11, *) Lb(iL), b(iL)
    end do
  close(11)

  fname = './b_function/M'//trim(ADJUSTL(char))//'/db_dL.dat'
  open(11, file = fname, status = 'unknown' )
    do iL = 0, 2 * M * 240 * N
      read(11, *) tmp, db_dL(iL)
    end do
  close(11)

  fname = './b_function/M'//trim(ADJUSTL(char))//'/theta.dat'
  open(11, file = fname, status = 'unknown' )
    do iL = 0, 2 * M * 240 * N
      read(11, *) Lb(iL), L_theta(iL)
    end do
  close(11)
  
  !do iL = 0, 2 * M * 240 * N
  !  loop: do i = 1, 90
  !    if (L_theta(iL) < i .and. i <= L_theta(iL+1)) then 
  !      print *, i, L_theta(iL), L_theta(iL+1), Lb(iL), 1.0d0 / ( dcos(dble(i)*pi/180.0d0)**2.0d0 )
  !      exit loop
  !    end if
  !  end do loop
  !end do

  print *, 'Finished reading b function and mapping data.'

  !stop


  ! Integration ------------------------------------------------

  j_para_i = 0.0d0
  j_para_i_new = 0.0d0
  f = 0.0d0
  g_x = 0.0d0

  do iter = 1, 4

    print *, 'Integration started...'

    !if (iter >= 2) then 
    !  call calc_j_para_i(N, M, Mdot, f, b, db_dL, j_para_i, j_para_i_new)
    !  j_para_i = j_para_i_new
    !  !do i = 1, 2*N*M*240+1
    !  !  arr(i) = j_para_i(i-1)
    !  !  if (i > 2*N*M*100) arr(i) = j_para_i(2*N*M*100-1)
    !  !end do
    !  !do j = 1, 10
    !  !  call p__smoothing_linear(19, 2*N*M*240+1, arr)
    !  !end do
    !  !do i = 1, 2*N*M*240+1
    !  !  if (mod(i,2*N*M) == 0) print *, L, j_para_i(i-1)*1.0d6, arr(i)*1.0d6
    !  !  j_para_i(i-1) = arr(i)
    !  !end do
    !  !stop
    !end if

    CALL g_int(N, M, Mdot, j_para_i, b, db_dL, g_x)

    do i = 1, 110 * N
      L = 1.0d0 + ( dble(i) ) / dble(N)
 
      !integration
      f_u = 0.0d0
      do j = 0, M * i - 1
        x = 1.0d0 + dble(j) / dble(N*M)
        dx = 1.0d0 / dble(N*M)
        f1 = 2.0d0 * x              * dexp( g_x(2*j)   - g_x(2*i*M) )
        f2 = 2.0d0 * (x + dx/2.0d0) * dexp( g_x(2*j+1) - g_x(2*i*M) )
        f3 = 2.0d0 * (x + dx)       * dexp( g_x(2*j+2) - g_x(2*i*M) )
        f_u = f_u + ( f1 + 4.0d0 * f2 + f3 ) * dx / 6.0d0
      enddo
      f(i) = f_u / L / L
      if (mod(i,N)==0) print *, 'iter = ', iter, 'L = ',L, 'f = ', f(i)
    enddo

    print *, 'Integration finished.'

    ! calculating field-aligned current
    do i = 1, 1
      call calc_j_para_i(N, M, Mdot, f, b, db_dL, j_para_i, j_para_i_new)
      j_para_i = j_para_i_new
    end do

    !call p__smoothing_linear(5, 2*N*240, f)
    !call p__smoothing_linear(5, 2*N*240, f)
    
    do i = 1, 100 * N
      farr(iter,i) = f(i)
      jarr(iter,i) = j_para_i(i)
    end do

    ! output ----------------------------------------------------
    write(char,*) nint(Mdot)
    write(num,*) nint(LT)
    open(11, file = trim(ADJUSTL(dir_name))//'/output/plasmav/wi_Mdot'//trim(ADJUSTL(char))//'_LT'&
    & //trim(ADJUSTL(num))//'.dat', status = 'unknown' )
      do i = 1, 100 * N
        L = 1.0d0 + ( dble(i) ) / dble(N)
        write(11, *) L, (1.0d0 - farr(j,i), j = 1, 10 )
      enddo
    close(11)

    open(11, file = trim(ADJUSTL(dir_name))//'/output/plasmav/vi_Mdot'//trim(ADJUSTL(char))//'_LT'&
    & //trim(ADJUSTL(num))//'.dat', status = 'unknown' )
      do i = 1, 100 * N
        L =  1.0d0 + ( dble(i) ) / dble(N)
        write(11, *) L, 2.0d0*pi/35729.0d0 * L*Rp/1.0d3, &
        & ((1.0d0 - farr(j,i))*2.0d0*pi/35729.0d0 * L*Rp/1.0d3, j = 1, 10)
      enddo
    close(11)

    open(11, file = trim(ADJUSTL(dir_name))//'/output/plasmav/j_para_i_Mdot'//trim(ADJUSTL(char))//'_LT'&
    & //trim(ADJUSTL(num))//'.dat', status = 'unknown' )
      do i = 1, 100 * N
        L =  1.0d0 + ( dble(i) ) / dble(N)
        write(11, *) L, (jarr(j,i), j = 1, 10)
      enddo
    close(11)

    open(11, file = trim(ADJUSTL(dir_name))//'/output/plasmav/Sigma_P_Mdot'//trim(ADJUSTL(char))//'_LT'&
    & //trim(ADJUSTL(num))//'.dat', status = 'unknown' )
      do i = 1, 100 * N
        L =  1.0d0 + ( dble(i) ) / dble(N)
        write(11, *) L, (Sigma_j(jarr(j,i)), j = 1, 10)
      enddo
    close(11)

  end do ! iter

  ! output ----------------------------------------------------
  write(char,*) nint(Mdot)
  write(num,*) nint(LT)
  open(11, file = trim(ADJUSTL(dir_name))//'/output/plasmav/wi_Mdot'//trim(ADJUSTL(char))//'_LT'&
  & //trim(ADJUSTL(num))//'.dat', status = 'unknown' )
    do i = 1, 100 * N
      L = 1.0d0 + ( dble(i) ) / dble(N)
      write(11, *) L, (1.0d0 - farr(iter,i), iter = 1, 10 )
    enddo
  close(11)

  open(11, file = trim(ADJUSTL(dir_name))//'/output/plasmav/vi_Mdot'//trim(ADJUSTL(char))//'_LT'&
  & //trim(ADJUSTL(num))//'.dat', status = 'unknown' )
    do i = 1, 100 * N
      L =  1.0d0 + ( dble(i) ) / dble(N)
      write(11, *) L, 2.0d0*pi/35729.0d0 * L*Rp/1.0d3, &
      & ((1.0d0 - farr(iter,i))*2.0d0*pi/35729.0d0 * L*Rp/1.0d3, iter = 1, 10)
    enddo
  close(11)

  open(11, file = trim(ADJUSTL(dir_name))//'/output/plasmav/j_para_i_Mdot'//trim(ADJUSTL(char))//'_LT'&
  & //trim(ADJUSTL(num))//'.dat', status = 'unknown' )
    do i = 1, 100 * N
      L =  1.0d0 + ( dble(i) ) / dble(N)
      write(11, *) L, (jarr(iter,i), iter = 1, 10)
    enddo
  close(11)

  open(11, file = trim(ADJUSTL(dir_name))//'/output/plasmav/Sigma_P_Mdot'//trim(ADJUSTL(char))//'_LT'&
  & //trim(ADJUSTL(num))//'.dat', status = 'unknown' )
    do i = 1, 100 * N
      L =  1.0d0 + ( dble(i) ) / dble(N)
      write(11, *) L, (Sigma_j(jarr(iter,i)), iter = 1, 10)
    enddo
  close(11)



end program



subroutine p__smoothing_linear(nsmooth, narr, arr)

  implicit none

  integer(4), parameter  :: sp = 4, dp = 8
  integer nsmooth, narr, i, j
  real(dp) arr(narr)
  real(dp) tmparr1(-nsmooth:narr+nsmooth)
  real(dp) tmparr2(-nsmooth:narr+nsmooth)

  ! smoothing
  do i = 1, narr
    tmparr1(i) = arr(i)
  end do

  do i = 0, nsmooth
    tmparr1(-i) = 2_dp*tmparr1(1-i) - tmparr1(2-i)
  end do

  do i = 0, nsmooth
    tmparr1(narr+i+1) = 2_dp*tmparr1(narr+i) - tmparr1(narr+i-1)
  end do

  do i = 1, narr
    tmparr2(i) = 0.0_dp
    do j = -nsmooth, nsmooth
      tmparr2(i) = tmparr2(i) + tmparr1(i+j)
    end do
    arr(i) = tmparr2(i)/dble(2*nsmooth+1)
  end do

end subroutine p__smoothing_linear