reset
#set terminal aqua enhanced
set terminal postscript enhanced

dir1 = "./metal_Hill"
dir2 = "./no_metal_Hill"

Mdot = "500"
LT   = "12"

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,20"
set title font"Arial,20"
set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -3,0
#unset key
set key left

set xr [0:0.4]
set yr [0:10]

set xlabel "J_{||i} [{/Symbol m}A/m^2]"
set ylabel "{/Symbol S}_P [S]"


set out 'figures/j_sigmaP.eps'
plot 	dir1.'/output/conductivity/j_sigma_metal_x1.0.dat'   u 1:3   w l lw 5 title "metal x1.0", \
		dir1.'/output/conductivity/j_sigma_metal_x0.1.dat'   u 1:3   w l lw 5 title "metal x0.1", \
		dir2.'/output/conductivity/j_sigma_no_metal.dat'     u 1:3   w l lw 5 title "w/o metal"


#f(x) = a*x**3 + b*x**2 + c*x + d*log(1+x) + e

set yr [0.001:10]
set log y

f1(x) = a1 * log(b1+c1*x*x)

a1               = 2.23295       
b1               = 1.09797       
c1               = 345.821       

fit f1(x) dir1.'/output/conductivity/j_sigma_metal_x1.0.dat' u 1:3 via a1,b1,c1

set out 'figures/j_sigma_metal_x1.0_fit.eps'
plot dir1.'/output/conductivity/j_sigma_metal_x1.0.dat' u 1:3 w linespoints, f1(x)



f2(x) = a2 * log(b2+c2*x*x)

a2               = 1.44487   
b2               = 1.05008   
c2               = 242.836   

fit f2(x) dir1.'/output/conductivity/j_sigma_metal_x0.1.dat' u 1:3 via a2,b2,c2

set out 'figures/j_sigma_metal_x0.1_fit.eps'
plot dir1.'/output/conductivity/j_sigma_metal_x0.1.dat' u 1:3 w linespoints, f2(x)



f3(x) = a3 * log(b3+c3*x*x)

a3              = 1.28267        
b3              = 1.00612        
c3              = 173.301        

#fit f3(x) dir2.'/output/conductivity/j_sigma_no_metal.dat' u 1:3 via a3,b3,c3

set out 'figures/j_sigma_no_metal_fit.eps'
plot dir2.'/output/conductivity/j_sigma_no_metal.dat' u 1:3 w linespoints, f3(x)

unset log y
set yr [0:10]

set out 'figures/j_sigmaP_fit.eps'
plot 	f1(x)  lw 5 title "metal x1.0", \
		f2(x)  lw 5 title "metal x0.1", \
		f3(x)  lw 5 title "w/o metal"

