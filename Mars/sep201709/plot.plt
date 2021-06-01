reset

color_list="grey10 gray20 grey20 gray30 grey30 gray40 grey40 gray50 \
grey50 gray60 grey60 gray70 grey70 gray80 grey80 gray90 grey90 gray100 grey100 \
gray grey light-gray light-grey dark-gray dark-grey red light-red dark-red \
yellow light-yellow dark-yellow green light-green dark-green spring-green \
forest-green sea-green blue light-blue dark-blue midnight-blue navy \
medium-blue royalblue skyblue cyan light-cyan dark-cyan magenta light-magenta \
dark-magenta turquoise light-turquoise dark-turquoise pink light-pink \
dark-pink coral light-coral orange-red salmon light-salmon dark-salmon \
aquamarine khaki dark-khaki goldenrod light-goldenrod dark-goldenrod gold \
beige brown orange dark-orange violet dark-violet plum purple"

color_list_16="#1a1a1a #333333 #333333 #4d4d4d #4d4d4d #666666 #666666 #7f7f7f \
#7f7f7f #999999 #999999 #b3b3b3 #b3b3b3 #cccccc #cccccc #e5e5e5 #e5e5e5 \
#ffffff #ffffff #bebebe #bebebe #d3d3d3 #d3d3d3 #a9a9a9 #a9a9a9 #ff0000 \
#f03232 #8b0000 #ffff00 #ffffe0 #c8c800 #00ff00 #90ee90 #006400 #00ff7f \
#228b22 #2e8b57 #0000ff #add8e6 #00008b #191970 #000080 #0000cd #4169e1 \
#87ceeb #00ffff #e0ffff #008b8b #ff00ff #f055f0 #8b008b #40e0d0 #afeeee \
#00ced1 #ffc0cb #ffb6c1 #ff1493 #ff7f50 #f08080 #ff4500 #fa8072 #ffa07a \
#e9967a #7fffd4 #f0e68c #bdb76b #daa520 #eedd82 #b8860b #ffd700 #f5f5dc \
#a52a2a #ffa500 #ff8c00 #ee82ee #9400d3 #dda0dd #a020f0"

################################ output density #############################

set term postscript enhanced color
#set term aqua
set log x
set xlabel 'density [m^-^3]'
set ylabel 'altitude [km]'

#set xr [0.01:5]

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,15"
set title font"Arial,20"

set key spacing 1.5

set format x "10^{%L}"
#unset key

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -2,0

set out "output_density_neutral.eps"
set xr [1E5:1E24]
set yr [0:200]
plot "./output/density/num/CO2.dat"   u 2:1 w linespoints lw 3 title "CO2", \
     "./output/density/num/CO.dat"    u 2:1 w linespoints lw 3 title "CO", \
     "./output/density/num/O2.dat"    u 2:1 w linespoints lw 3 title "O2", \
     "./output/density/num/O.dat"     u 2:1 w linespoints lw 3 title "O", \
     "./output/density/num/O(1D).dat" u 2:1 w linespoints lw 3 title "O(^1D)", \
     "./output/density/num/H2O.dat"   u 2:1 w linespoints lw 3 title "H2O", \
     "./output/density/num/OH.dat"    u 2:1 w linespoints lw 3 title "OH", \
     "./output/density/num/HO2.dat"   u 2:1 w linespoints lw 3 title "HO2", \
     "./output/density/num/H2O2.dat"  u 2:1 w linespoints lw 3 title "H2O2", \
     "./output/density/num/O3.dat"    u 2:1 w linespoints lw 3 title "O3", \
     "./output/density/num/H2.dat"    u 2:1 w linespoints lw 3 title "H2", \
     "./output/density/num/H.dat"     u 2:1 w linespoints lw 3 title "H", \
     "./output/density/num/N2.dat"    u 2:1 w linespoints lw 3 title "N2"

set out "output_density_Chaffin.eps"
set xr [1E5:1E24]
set yr [0:200]
plot "./Chaffin_CO2.dat"   u ($2*1000000):($1-1) w linespoints lw 3 title "CO2", \
     "./Chaffin_CO.dat"    u ($2*1000000):($1-1) w linespoints lw 3 title "CO", \
     "./Chaffin_O2.dat"    u ($2*1000000):($1-1) w linespoints lw 3 title "O2", \
     "./Chaffin_O.dat"     u ($2*1000000):($1-1) w linespoints lw 3 title "O", \
     "./Chaffin_O1D.dat"   u ($2*1000000):($1-1) w linespoints lw 3 title "O(^1D)", \
     "./Chaffin_H2O.dat"   u ($2*1000000):($1-1) w linespoints lw 3 title "H2O", \
     "./Chaffin_OH.dat"    u ($2*1000000):($1-1) w linespoints lw 3 title "OH", \
     "./Chaffin_HO2.dat"   u ($2*1000000):($1-1) w linespoints lw 3 title "HO2", \
     "./Chaffin_H2O2.dat"  u ($2*1000000):($1-1) w linespoints lw 3 title "H2O2", \
     "./Chaffin_O3.dat"    u ($2*1000000):($1-1) w linespoints lw 3 title "O3", \
     "./Chaffin_H2.dat"    u ($2*1000000):($1-1) w linespoints lw 3 title "H2", \
     "./Chaffin_H.dat"     u ($2*1000000):($1-1) w linespoints lw 3 title "H"


set xr [1E-4:1E1]
set yr [0:200]
unset key
set out "CO_CO2.eps"
set xlabel 'CO/CO2 ratio'
plot '< paste output/density/num/CO2.dat output/density/num/CO.dat' u ($4/$2):1

set xr [1E-2:1E0]
set out "N2_CO2.eps"
set xlabel 'N2/CO2 ratio'
plot '< paste output/density/num/CO2.dat output/density/num/N2.dat' u ($4/$2):1
set key

#set out "output_density_ion.eps"
#set xr [1E5:1E13]
#set yr [0:200]
#plot "./output/density/num/e-.dat"    u 2:1 w l lw 3 title "e^-", \
#     "./output/density/num/CO2+.dat"  u 2:1 w l lw 3 title "CO_2^+", \
#     "./output/density/num/CO+.dat"   u 2:1 w l lw 3 title "CO^+", \
#     "./output/density/num/O2+.dat"   u 2:1 w l lw 3 title "O_2^+", \
#     "./output/density/num/C+.dat"    u 2:1 w l lw 3 title "C^+", \
#     "./output/density/num/OH+.dat"   u 2:1 w l lw 3 title "OH^+", \
#     "./output/density/num/NO+.dat"   u 2:1 w l lw 3 title "NO^+", \
#     "./output/density/num/H2+.dat"   u 2:1 w l lw 3 title "H_2^+", \
#     "./output/density/num/H+.dat"    u 2:1 w l lw 3 title "H^+", \
#     "./output/density/num/N2+.dat"   u 2:1 w l lw 3 title "N_2^+", \
#     "./output/density/num/N+.dat"    u 2:1 w l lw 3 title "N^+", \
#     "./output/density/num/O+(2D).dat" u 2:1 w l lw 3 title "O^+(^2D)", \
#     "./output/density/num/O+(2P).dat" u 2:1 w l lw 3 title "O^+(^2P)", \
#     "./output/density/num/O+(4S).dat" u 2:1 w l lw 3 title "O^+(^4S)"
#     #"./output/density/num/O(1D).dat" u 2:1 w linespoints lw 3 title "O(^1D)", \
#     #"./output/density/num/OH.dat"    u 2:1 w linespoints lw 3 title "OH", \
#     #"./output/density/num/NO.dat"    u 2:1 w linespoints lw 3 title "NO", \
#     #"./output/density/num/O3.dat"    u 2:1 w linespoints lw 3 title "O3", \
#     #"./output/density/num/CO2.dat"   u 2:1 w linespoints lw 3 title "CO2"


#set yr [0:200]
#set xr [1E-15:1]
#set xlabel 'volume mixing ratio [mol/mol]'
#set out "output_mixing_ratio.eps"
#
#plot "./output/density/vmr/vmr_CO2.dat"   u 2:1 w linespoints lw 3 title "CO2" , \
#     "./output/density/vmr/vmr_CO.dat"    u 2:1 w linespoints lw 3 title "CO" , \
#     "./output/density/vmr/vmr_O2.dat"    u 2:1 w linespoints lw 3 title "O2" , \
#     "./output/density/vmr/vmr_O.dat"     u 2:1 w linespoints lw 3 title "O" , \
#     "./output/density/vmr/vmr_OH.dat"    u 2:1 w linespoints lw 3 title "OH", \
#     "./output/density/vmr/vmr_HO2.dat"   u 2:1 w linespoints lw 3 title "HO2", \
#     "./output/density/vmr/vmr_NO.dat"    u 2:1 w linespoints lw 3 title "NO", \
#     "./output/density/vmr/vmr_O3.dat"    u 2:1 w linespoints lw 3 title "O3", \
#     "./output/density/vmr/vmr_e-.dat"    u 2:1 w l lw 3 title "e^-", \
#     "./output/density/vmr/vmr_CO2+.dat"  u 2:1 w l lw 3 title "CO_2^+", \
#     "./output/density/vmr/vmr_CO+.dat"   u 2:1 w l lw 3 title "CO^+", \
#     "./output/density/vmr/vmr_O2+.dat"   u 2:1 w l lw 3 title "O_2^+", \
#     "./output/density/vmr/vmr_C+.dat"    u 2:1 w l lw 3 title "C^+", \
#     "./output/density/vmr/vmr_O+.dat"    u 2:1 w l lw 3 title "O^+", \
#     "./output/density/vmr/vmr_OH+.dat"   u 2:1 w l lw 3 title "OH^+", \
#     "./output/density/vmr/vmr_NO+.dat"   u 2:1 w l lw 3 title "NO^+", \
#     "./output/density/vmr/vmr_H2+.dat"   u 2:1 w l lw 3 title "H_2^+", \
#     "./output/density/vmr/vmr_H+.dat"    u 2:1 w l lw 3 title "H^+", \
#     "./output/density/vmr/vmr_N2+.dat"   u 2:1 w l lw 3 title "N_2^+", \
#     "./output/density/vmr/vmr_N+.dat"    u 2:1 w l lw 3 title "N^+", \
#     "./output/density/vmr/vmr_O(1D).dat" u 2:1 w l lw 3 title "O(^1D)", \
#     "./output/density/vmr/vmr_O+(2D).dat" u 2:1 w l lw 3 title "O^+(^2D)", \
#     "./output/density/vmr/vmr_O+(2P).dat" u 2:1 w l lw 3 title "O^+(^2P)", \

################################ input density #############################
#set xr [1E5:1E25]
#set yr [0:300]
#set out "input_density.eps"
#
#plot "./input/density/CO2.dat"    u 2:1 w l lw 3 title "CO2", \
#     "./input/density/CO.dat"     u 2:1 w l lw 3 title "CO", \
#     "./input/density/O2.dat"     u 2:1 w l lw 3 title "O2", \
#     "./input/density/O.dat"      u 2:1 w l lw 3 title "O", \
#     "./input/density/N2.dat"     u 2:1 w l lw 3 title "N2", \
#     "./input/density/H2.dat"     u 2:1 w l lw 3 title "H2"
#

################################# Temperature ##############################
reset
set term postscript enhanced color
#set term aqua
set xlabel 'Temperature [K]'
set ylabel 'altitude [km]'
set xr [0:300]
set yr [0:500]
#set xr [0.01:5]

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,15"
set title font"Arial,20"

set key spacing 1.5

#unset key

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -2,0

set out "input_temperature.eps"

#plot "./input/Temperature/T_n_Mars.dat"    u 2:1 w l lw 3 title "neutral"


################################ output flux #############################

reset
set term postscript enhanced color
#set term aqua
set log x
set xlabel 'flux [cm^{-2}s^{-1}]'
set ylabel 'altitude [km]'
set xr [1E0:1E15]
set yr [0:200]
#set xr [0.01:5]

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,15"
set title font"Arial,20"

set key spacing 1.5

set format x "10^{%L}"
#unset key

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -2,0

set out "output_flux.eps"

plot "./output/flux/CO.dat"    u 2:1 w l lw 3           lc "red"    title "CO up  ", \
     "./output/flux/CO.dat"    u 3:1 w linespoints lw 3 lc "red"    title "CO down", \
     "./output/flux/O2.dat"    u 2:1 w l lw 3           lc "blue"   title "O2 up  ", \
     "./output/flux/O2.dat"    u 3:1 w linespoints lw 3 lc "blue"   title "O2 down", \
     "./output/flux/O.dat"     u 2:1 w l lw 3           lc "green"  title "O up  ", \
     "./output/flux/O.dat"     u 3:1 w linespoints lw 3 lc "green"  title "O down", \
     "./output/flux/H2.dat"    u 2:1 w l lw 3           lc "dark-green"  title "H2 up  ", \
     "./output/flux/H2.dat"    u 3:1 w linespoints lw 3 lc "dark-green"  title "H2 down", \
     "./output/flux/H.dat"     u 2:1 w l lw 3           lc "orange" title "H up  ", \
     "./output/flux/H.dat"     u 3:1 w linespoints lw 3 lc "orange" title "H down", \
     "./output/flux/OH.dat"    u 2:1 w l lw 3           lc "black"  title "OH up  ", \
     "./output/flux/OH.dat"    u 3:1 w linespoints lw 3 lc "black"  title "OH down", \
     "./output/flux/O3.dat"    u 2:1 w l lw 3           lc "purple" title "O3 up  ", \
     "./output/flux/O3.dat"    u 3:1 w linespoints lw 3 lc "purple" title "O3 down"

set xr [1E-10:1E5]
set xlabel 'velocity [m/s]'
set out "output_flux_velocity.eps"

plot "./output/flux/v_CO.dat"    u 2:1 w l lw 3           lc "red"    title "CO up  ", \
     "./output/flux/v_CO.dat"    u 3:1 w linespoints lw 3 lc "red"    title "CO down", \
     "./output/flux/v_O2.dat"    u 2:1 w l lw 3           lc "blue"   title "O2 up  ", \
     "./output/flux/v_O2.dat"    u 3:1 w linespoints lw 3 lc "blue"   title "O2 down", \
     "./output/flux/v_O.dat"     u 2:1 w l lw 3           lc "green"  title "O up  ", \
     "./output/flux/v_O.dat"     u 3:1 w linespoints lw 3 lc "green"  title "O down", \
     "./output/flux/v_H2.dat"    u 2:1 w l lw 3           lc "dark-green"  title "H2 up  ", \
     "./output/flux/v_H2.dat"    u 3:1 w linespoints lw 3 lc "dark-green"  title "H2 down", \
     "./output/flux/v_H.dat"     u 2:1 w l lw 3           lc "orange" title "H up  ", \
     "./output/flux/v_H.dat"     u 3:1 w linespoints lw 3 lc "orange" title "H down", \
     "./output/flux/v_OH.dat"    u 2:1 w l lw 3           lc "black"  title "OH up  ", \
     "./output/flux/v_OH.dat"    u 3:1 w linespoints lw 3 lc "black"  title "OH down", \
     "./output/flux/v_O3.dat"    u 2:1 w l lw 3           lc "purple" title "O3 up  ", \
     "./output/flux/v_O3.dat"    u 3:1 w linespoints lw 3 lc "purple" title "O3 down"


################################ output D, K #############################

reset
set term postscript enhanced color
#set term aqua
set log x
set xlabel 'diffusion coefficient [cm^2/s]'
set ylabel 'altitude [km]'
set xr [1E0:3E9]
set yr [0:200]
#set xr [0.01:5]

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,15"
set title font"Arial,20"

set key spacing 1.5

set format x "10^{%L}"
#unset key

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -2,0

set out "output_diffusion_coefficient.eps"

plot "./output/flux/K_eddy.dat"  u 2:1 w l lw 3 title "K", \
     "./output/flux/D_H.dat"     u 2:1 w l lw 3 title "D_H", \
     "./output/flux/D_H2.dat"    u 2:1 w l lw 3 title "D_{H_2}", \
     "./output/flux/D_O.dat"     u 2:1 w l lw 3 title "D_O", \


#reset
#set term aqua
#set pm3d 
#set pm3d map
#set ticslevel 0
#set log cb
#set cbrange[1E6:1E12]
#set xr [1:361]
#set yr [00:500]
#set zr [1E6:1E12]
#set format cb "10^{%L}"
##set palette defined (1E-10 "black", 1E-8 "blue" , 1E-6 "green", 1E-4 "red")
#splot "./output/density/rot_e-.dat" u 1:3:4 with pm3d 


##################################
reset
#set terminal aqua
#set view 135,270,1.2, 1.8
#set view 80,0,1.2, 1.5
#set mapping spherical
#set angle degrees
#set hidden3d front
#set nozeroaxis
#set noborder
#set nogrid
#set tics font "Arial,25"
#unset key
#unset xtics
#unset ytics
#unset ztics
#set log cb
#set format cb "10^{%L}"
#set cbrange[3E10:1E11]
##set cbtics 0.2
##set palette defined (8 "blue", 10 "green", 12 "red")
##set palette defined (-25 "blue", 0 "white", 10 "red")
#set palette rgbformulae 22, 13, -31
#set xyplane at -2
#set xrange [-1.2:1.2]
#set yrange [-1.2:1.2]
#set zrange [-1.2:1.2]
#set colorbox user vertical origin .8,.1 size .04,.8
#set parametric
#set isosamples 25
#set urange [0:360]
#set vrange [-90:90]
## Parametric functions for the sphere
#r = 1.01
#fx(v,u) = r*cos(v)*cos(u)
#fy(v,u) = r*cos(v)*sin(u)
#fz(v)   = r*sin(v)
#
##set term png size 480,360
##set out "global_e-_100km.png"

#splot "./output/density/global/global_65_O2+.dat" using 2:1:(1):3 w pm3d ,r*cos(v)*cos(u),r*cos(v)*sin(u),r*sin(v) w l ls 2 lc rgb "black"
