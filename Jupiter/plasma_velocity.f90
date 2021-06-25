!integration in g(x) ; g(x) = exp( h(x) )　これサブルーチンにして、積分計算して呼び出すだけにしたいよね
subroutine g_int( N, M, Mdot, g_x, Sigma, b, db_dL)
  implicit none
  integer N, M, i
  double precision g_x(0:2*N*M*240)
  double precision Sigma(0:2*N*M*240)
  double precision b(0:2*N*M*240)
  double precision db_dL(0:2*N*M*240)
  double precision K
  double precision dx, xi
  double precision Mdot, Rp, Bp, tmp
  double precision g1, g2, g3
  double precision, parameter :: pi = acos(-1.0d0)

  Rp = 7.1492d+7 ! planetary radius [m]
  Bp = 4.2d-4 ! magnetic field at equator surface [T] = [kg/A/s^2]
  ! Sigma: [A^2 s^3 /kg/m^2]
  ! Rp^2 Bp^2 Sigma / Mdot = m^2 kg^2 A^-2 s^-4 A^2 s^3 kg^-1 m^-2 kg^-1 s = [no unit]

  tmp = 4.0d0 * pi * Rp**2.0d0 * Bp**2.0d0 / Mdot

  g_x = 0.0d0
  do i = 0, 2 * N * M * 240 - 1
    xi = 1.0d0 + dble(i) / dble(2*N*M)
    dx = 1.0d0 / dble(2*N*M)
    g1 = K(N, M, i, Sigma, b, db_dL, xi)
    g2 = K(N, M, i, Sigma, b, db_dL, xi + dx/2.0d0)
    g3 = K(N, M, i, Sigma, b, db_dL, xi + dx)
    g_x(i+1) = g_x(i) + tmp * (g1 + 4.0d0 * g2 * g3) * dx / 6.0d0
  end do

end subroutine g_int

double precision function K(N, M, j, Sigma, b, db_dL, x)
  implicit none
  integer N, M, j
  double precision Sigma(0:2*N*M*240)
  double precision b(0:2*N*M*240)
  double precision db_dL(0:2*N*M*240)
  double precision x
  
  !K =  Sigma(j) * dsqrt( 1.0d0 - 1.0d0 / x ) / x**5.0d0

  b(j) = 1.0d0/x
  db_dL(j) = -1.0d0/x/x

  K =  - Sigma(j) * b(j) / x**2.0d0 * dsqrt(1.0d0 - b(j)) * db_dL(j)
  
  return
end function K



!begin calculation of Hill current

program plasma_velocity
  implicit none

  !parameter

  double precision, parameter :: pi = acos(-1.0d0)

  double precision x              !x
  double precision dx             !integration grid size : 1 / M
  double precision L              !distance normalized at RJ
  double precision, parameter :: L0 = 23.0d0 !corotation lag become significant at L0
  double precision f_u            !numerator of f(x)
  double precision f(20000)       !f = dw/wp : w/wp = 1 - f : corotation lag
  double precision f1             !g(x)
  double precision f2             !g(x+dx/2)
  double precision f3             !g(x+dx)
  double precision theta
  double precision Mdot, tmp, LT
  double precision, parameter :: Rp = 7.1492d+7

  integer, parameter :: N = 1   !grid number of L
  integer, parameter :: M = 100000  !grid number of x in 1 L
  integer, parameter :: tg = 10   !grid number of theta

  double precision g_x(0:2*N*M*240)
  double precision Sigma(0:2*N*M*240)
  double precision b(0:2*N*M*240)
  double precision db_dL(0:2*N*M*240)
  double precision L_theta(0:2*N*M*240)
  double precision Lb(0:2*N*M*240)

  double precision theta_to_L(181)


  integer i   !a grid number of L
  integer j   !a grid number of x
  integer ix, iy, iB, iL
  integer, parameter :: nx = 361
  integer, parameter :: ny = 181

  double precision hi_sigma_P(nx,ny), hi_sigma_H(nx,ny), hi_sigma_0(nx,ny)

  character(len=256) fname, num, char, dir_name

  !real(8) gamain, uicg_p1, uicg, alpha
  !real(8) Fe, be
  !integer ( kind = 4 ) ifault



  ! Parameter settings ----------------------------------------
  dir_name = "./no_metal_Hill"
  Mdot = 500.0d0
  LT   = 12.0d0



  ! Read mapping function ------------------------------
  print *, 'Reading b function and mapping data...'
  write(char,*) M
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
      read(11, *) tmp, L_theta(iL)
    end do
  close(11)

  print *, 'Finished reading b function and mapping data.'

  !stop



  ! read conductivity data ------------------------------------

  iB = 5
  write(char,*) iB

  fname = trim(ADJUSTL(dir_name))//'/output/conductivity/B_Jupiter_i='//trim(ADJUSTL(char))//'/hi_sigma_HP0.dat'
  open(11, file = fname, status = 'unknown' )
    do iy = 1, ny
    do ix = 1, nx
      read(11, *) tmp, tmp, hi_sigma_H(ix,iy), hi_sigma_P(ix,iy), hi_sigma_0(ix,iy)
    end do
    end do
  close(11)

  ! Pedersen conductivity as a function of L
  ix = nint(dble(nx-1)/24.0d0*LT+1.0d0)

  do j = 0, 2 * M * 240 * N - 1 ! L/dL = M
    L = 1.0d0 + ( dble(j)/dble(2*N*M) )
    theta = L_theta(j)
    iy = int(theta)
    Sigma(j) = hi_sigma_P(ix,iy+91)*(iy+1-theta) + hi_sigma_P(ix,iy+1+91)*(theta-iy)
    !print *, j, L, iy, theta, iy+1, Sigma(j)
  end do

  print *, 'Finished calculating Sigma as a function of L.'
  !stop




  ! Integration ------------------------------------------------
  print *, 'Integration started...'

  CALL g_int(N, M, Mdot, g_x, Sigma, b, db_dL)

  do i = 1, 100 * N
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
    print *, 'L = ',L, 'f = ', f(i)
  enddo
  print *, 'Integration finished.'

  ! output ----------------------------------------------------
  write(char,*) nint(Mdot)
  write(num,*) nint(LT)
  open(11, file = trim(ADJUSTL(dir_name))//'/output/plasmav/wi_Mdot'//trim(ADJUSTL(char))//'_LT'&
  & //trim(ADJUSTL(num))//'.dat', status = 'unknown' )
    do i = 1, 100 * N
      L = 1.0d0 + ( dble(i) ) / dble(N)
      write(11, *) L, 1.0d0 - f(i)
    enddo
  close(11)

  open(11, file = trim(ADJUSTL(dir_name))//'/output/plasmav/vi_Mdot'//trim(ADJUSTL(char))//'_LT'&
  & //trim(ADJUSTL(num))//'.dat', status = 'unknown' )
    do i = 1, 100 * N
      L =  1.0d0 + ( dble(i) ) / dble(N)
      write(11, *) L, (1.0d0 - f(i))*2.0d0*pi/35729.0d0 * L*Rp/1.0d3, 2.0d0*pi/35729.0d0 * L*Rp/1.0d3
    enddo
  close(11)



end program
