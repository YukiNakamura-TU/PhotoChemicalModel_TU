!integration in g(x) ; g(x) = exp( h(x) )　これサブルーチンにして、積分計算して呼び出すだけにしたいよね
subroutine g_int( N, M, g_x, g_x2 )
  implicit none
  integer N, M, i
  double precision g_x(0:N*M*240)
  double precision g_x2(0:2*N*M*240)
  double precision K
  double precision dx, xi
  double precision Mdot, Rp, Bp, tmp
  double precision g1, g2, g3
  double precision, parameter :: pi = acos(-1.0d0)

  Rp = 7.1492d+7
  Bp = 4.2d-4
  Mdot = 500.0d0

  tmp = 4.0d0 * pi * Rp**2.0d0 * Bp**2.0d0 / Mdot

  g_x2 = 0.0d0
  do i = 0, 2 * N * M * 240 - 1
    xi = 1.0d0 + dble(i) / dble(2*N*M)
    dx = 1.0d0 / dble(2*N*M)
    g1 = K(xi)
    g2 = K(xi + dx/2.0d0)
    g3 = K(xi + dx)
    g_x2(i+1) = g_x2(i) + tmp * (g1 + 4.0d0 * g2 * g3) * dx / 6.0d0
  end do

end subroutine g_int

double precision function Sigma(x)
  implicit none
  double precision x
  ! Metal
  !Sigma = 0.35d0/x + 0.2d0 + 0.05d0 * dexp(-(x-19.0d0)*2.0d0/40.0d0)
  ! NO Metal
  Sigma = 0.2d0/x**2.0d0 + 0.03d0  

  return
end function Sigma

double precision function K(x)
  implicit none
  double precision x, Sigma
  
  K =  Sigma(x) * dsqrt( 1.0d0 - 1.0d0 / x )**0.5d0 / x**5.0d0
  
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
double precision L_f(10000)              !distance normalized at RJ
double precision, parameter :: L0 = 23.0d0 !corotation lag become significant at L0
double precision f_u            !numerator of f(x)
double precision f(20000)       !f = dw/wp : w/wp = 1 - f : corotation lag
double precision g1             !g(x)
double precision g2             !g(x+dx/2)
double precision g3             !g(x+dx)
double precision h, K       !integration in g(x) ; g(x) = exp( h(x) )
double precision P(10000)
double precision JII(10000)
double precision JII_sp(0:10000)
double precision theta
double precision Sum_j(0:10000)
double precision Num_j(0:10000)
double precision Mdot, Sigma, Sigma_i(10000), Rp, Bp, tmp

integer, parameter :: N = 1   !grid number of L
integer, parameter :: M = 100000  !grid number of x in 1
integer, parameter :: D = M * N !total number of grid x
integer, parameter :: tg = 10   !grid number of theta

double precision g_x(0:N*M*240)
double precision g_x2(0:2*N*M*240)

integer i   !a grid number of L
integer j   !a grid number of x
integer, parameter :: tn_j = 181

Rp = 7.1492d+7

CALL g_int(N,M,g_x,g_x2)

!calculation
do i = 1, 240 * N
  L = 1.0d0 + ( dble(i) ) / dble(N)

  !integration
  f_u = 0.0d0
  do j = 1, M * i
    x = 1.0d0 + ( dble(j) - 1.0d0 ) / dble(D)
    dx = 1.0d0 / dble(D)
    g1 = 2.0d0 * x              * dexp( g_x2(2*j)   - g_x2(2*i*M) )
    g2 = 2.0d0 * (x + dx/2.0d0) * dexp( g_x2(2*j+1) - g_x2(2*i*M) )
    g3 = 2.0d0 * (x + dx)       * dexp( g_x2(2*j+2) - g_x2(2*i*M) )
    f_u = f_u + ( g1 + 4.0d0 * g2 + g3 ) * dx / 6.0d0
  enddo
  f(i) = f_u / L / L
enddo

open(11, file = './wi_M500_NO_Metal.dat', status = 'unknown' )
  do i = 1, 240 * N
    L = 1.0d0 + ( dble(i) ) / dble(N)
    write(11, *) L, 1.0d0 - f(i)
  enddo
close(11)

open(11, file = './vi_M500_NO_Metal.dat', status = 'unknown' )
  do i = 1, 240 * N
    L = 1.0d0 + ( dble(i) ) / dble(N)
    write(11, *) L, (1.0d0 - f(i))*2.0d0*pi/35729.0d0 * L*Rp/1.0d3, 2.0d0*pi/35729.0d0 * L*Rp/1.0d3
  enddo
close(11)



end program
