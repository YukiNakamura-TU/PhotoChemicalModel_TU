!integration in g(x) ; g(x) = exp( h(x) )　これサブルーチンにして、積分計算して呼び出すだけにしたいよね
subroutine g_int( N, M, Mdot, g_x, Sigma)
  implicit none
  integer N, M, i
  double precision g_x(0:2*N*M*240)
  double precision Sigma(0:2*N*M*240)
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
    g1 = K(N, M, i, Sigma, xi)
    g2 = K(N, M, i, Sigma, xi + dx/2.0d0)
    g3 = K(N, M, i, Sigma, xi + dx)
    g_x(i+1) = g_x(i) + tmp * (g1 + 4.0d0 * g2 * g3) * dx / 6.0d0
  end do

end subroutine g_int

double precision function K(N, M, j, Sigma, x)
  implicit none
  integer N, M, j
  double precision Sigma(0:2*N*M*240)
  double precision x
  
  K =  Sigma(j) * dsqrt( 1.0d0 - 1.0d0 / x )**0.5d0 / x**5.0d0
  
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
  double precision g1             !g(x)
  double precision g2             !g(x+dx/2)
  double precision g3             !g(x+dx)
  double precision theta
  double precision Mdot, tmp, LT
  double precision, parameter :: Rp = 7.1492d+7

  integer, parameter :: N = 1   !grid number of L
  integer, parameter :: M = 100000  !grid number of x in 1
  integer, parameter :: D = M * N !total number of grid x
  integer, parameter :: tg = 10   !grid number of theta

  double precision g_x(0:2*N*M*240)
  double precision Sigma(0:2*N*M*240)

  integer i   !a grid number of L
  integer j   !a grid number of x
  integer ix, iy, iB
  integer, parameter :: nx = 361
  integer, parameter :: ny = 181

  double precision hi_sigma_P(nx,ny), hi_sigma_H(nx,ny), hi_sigma_0(nx,ny)

  character(len=256) fname, num, char, dir_name



  ! Parameter settings ----------------------------------------
  dir_name = "./no_metal_Hill"
  Mdot = 500.0d0
  LT   = 12.0d0


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

  do j = 1, M * 240 * N ! L/dL = M
    L = 1.0d0 + ( dble(j)/dble(M) ) / dble(N)
    theta = dacos(dsqrt(1.0d0/L))*180.0d0/pi
    iy = int(theta)
    Sigma(j) = hi_sigma_P(ix,iy)*(iy+1-theta) + hi_sigma_P(ix,iy+1)*(theta-iy)
  end do

  print *, 'Finished calculating Sigma as a function of L.'
  !stop




  ! Integration ------------------------------------------------
  print *, 'Integration started...'

  CALL g_int(N,M,Mdot,g_x,Sigma)

  do i = 1, 100 * N
    L = 1.0d0 + ( dble(i) ) / dble(N)

    !integration
    f_u = 0.0d0
    do j = 1, M * i
      x = 1.0d0 + ( dble(j) - 1.0d0 ) / dble(D)
      dx = 1.0d0 / dble(D)
      g1 = 2.0d0 * x              * dexp( g_x(2*j)   - g_x(2*i*M) )
      g2 = 2.0d0 * (x + dx/2.0d0) * dexp( g_x(2*j+1) - g_x(2*i*M) )
      g3 = 2.0d0 * (x + dx)       * dexp( g_x(2*j+2) - g_x(2*i*M) )
      f_u = f_u + ( g1 + 4.0d0 * g2 + g3 ) * dx / 6.0d0
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
      L = 1.0d0 + ( dble(i) ) / dble(N)
      write(11, *) L, (1.0d0 - f(i))*2.0d0*pi/35729.0d0 * L*Rp/1.0d3, 2.0d0*pi/35729.0d0 * L*Rp/1.0d3
    enddo
  close(11)



end program
