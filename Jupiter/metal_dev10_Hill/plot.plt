################################ output density #############################
reset 

set term postscript enhanced color
#set term aqua
set log x
set xlabel 'density [m^-^3]'
set ylabel 'altitude [km]'
set xr [1E6:3E12]
set yr [200:1000]
#set xr [0.01:5]

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,20"
set title font"Arial,20"

set key spacing 1.5

set format x "10^{%L}"
#unset key

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -2,0

set out "output_density_1000km.eps"

plot "./output/density/num/e-.dat"    u 2:1 w l lw 3 title "e^-", \
     "./output/density/num/H3+.dat"   u 2:1 w l lw 3 title "H_3^+", \
     "./output/density/num/H+.dat"    u 2:1 w l lw 3 title "H^+", \
     "./output/density/num/H2+.dat"   u 2:1 w l lw 3 title "H_2^+", \
     "./output/density/num/CH3+.dat"  u 2:1 w l lw 3 title "CH_3^+ ", \
     "./output/density/num/CH4+.dat"  u 2:1 w l lw 3 title "CH_4^+ ", \
     "./output/density/num/CH5+.dat"  u 2:1 w l lw 3 title "CH_5^+ ", \
     "./output/density/num/C2H2+.dat" u 2:1 w l lw 3 title "C_2H_2^+", \
     "./output/density/num/C2H3+.dat" u 2:1 w l lw 3 title "C_2H_3^+", \
     "./output/density/num/C2H5+.dat" u 2:1 w l lw 3 title "C_2H_5^+", \
     "./output/density/num/C3H5+.dat" u 2:1 w l lw 3 title "C_3H_5^+", \
     "./output/density/num/C4H7+.dat" u 2:1 w l lw 3 title "C_4H_7^+", \
     "./output/density/num/Fe+.dat"   u 2:1 w l lw 3 title "Fe^+", \
     "./output/density/num/Mg+.dat"   u 2:1 w l lw 3 title "Mg^+", \
     "./output/density/num/Na+.dat"   u 2:1 w l lw 3 title "Na^+", \
     "./output/density/num/Si+.dat"   u 2:1 w l lw 3 title "Si^+", \
     "./output/density/num/Fe.dat"    u 2:1 w l lw 3 title "Fe", \
     "./output/density/num/Mg.dat"    u 2:1 w l lw 3 title "Mg", \
     "./output/density/num/Na.dat"    u 2:1 w l lw 3 title "Na", \
     "./output/density/num/Si.dat"    u 2:1 w l lw 3 title "Si"

set yr [0:3000]
set out "output_density_3000km.eps"

plot "./output/density/num/e-.dat"    u 2:1 w l lw 3 title "e^-", \
     "./output/density/num/H3+.dat"   u 2:1 w l lw 3 title "H_3^+", \
     "./output/density/num/H+.dat"    u 2:1 w l lw 3 title "H^+", \
     "./output/density/num/H2+.dat"   u 2:1 w l lw 3 title "H_2^+", \
     "./output/density/num/CH3+.dat"  u 2:1 w l lw 3 title "CH_3^+ ", \
     "./output/density/num/CH4+.dat"  u 2:1 w l lw 3 title "CH_4^+ ", \
     "./output/density/num/CH5+.dat"  u 2:1 w l lw 3 title "CH_5^+ ", \
     "./output/density/num/C2H2+.dat" u 2:1 w l lw 3 title "C_2H_2^+", \
     "./output/density/num/C2H3+.dat" u 2:1 w l lw 3 title "C_2H_3^+", \
     "./output/density/num/C2H5+.dat" u 2:1 w l lw 3 title "C_2H_5^+", \
     "./output/density/num/C3H5+.dat" u 2:1 w l lw 3 title "C_3H_5^+", \
     "./output/density/num/C4H7+.dat" u 2:1 w l lw 3 title "C_4H_7^+", \
     "./output/density/num/Fe+.dat"   u 2:1 w l lw 3 title "Fe^+", \
     "./output/density/num/Mg+.dat"   u 2:1 w l lw 3 title "Mg^+", \
     "./output/density/num/Na+.dat"   u 2:1 w l lw 3 title "Na^+", \
     "./output/density/num/Si+.dat"   u 2:1 w l lw 3 title "Si^+", \
     "./output/density/num/Fe.dat"    u 2:1 w l lw 3 title "Fe", \
     "./output/density/num/Mg.dat"    u 2:1 w l lw 3 title "Mg", \
     "./output/density/num/Na.dat"    u 2:1 w l lw 3 title "Na", \
     "./output/density/num/Si.dat"    u 2:1 w l lw 3 title "Si"

################################ input density #############################
set xr [1E6:3E24]
set out "input_density.eps"

plot "./input/density/H2.dat"     u 2:1 w l lw 3 title "H_2", \
     "./input/density/He.dat"     u 2:1 w l lw 3 title "He", \
     "./input/density/CH4.dat"    u 2:1 w l lw 3 title "CH_4", \
     "./input/density/C2H2.dat"   u 2:1 w l lw 3 title "C_2H_2", \
     "./input/density/C2H4.dat"   u 2:1 w l lw 3 title "C_2H_4", \
     "./input/density/C2H6.dat"   u 2:1 w l lw 3 title "C_2H_6", \
     "./output/density/num/H.dat" u 2:1 w l lw 3 title "H"


################################# Temperature ##############################
reset 

set term postscript enhanced color
#set term aqua
set xlabel 'Temperature [K]'
set ylabel 'altitude [km]'
set xr [0:1200]
set yr [0:3000]
#set xr [0.01:5]

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,20"
set title font"Arial,20"

set key spacing 1.5

#unset key

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -2,0

#set out "input_temperature.eps"
#
#plot "./input/Temperature/T_n_Jupiter.dat"    u 2:1 w l lw 3 title "neutral"


################################ output flux #############################
reset 

set term postscript enhanced color
#set term aqua
set log x
set xlabel 'flux [/m^2/s]'
set ylabel 'altitude [km]'
set xr [1E-26:3E25]
set yr [0:3000]
#set xr [0.01:5]

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,20"
set title font"Arial,20"

set key spacing 1.5

set format x "10^{%L}"
#unset key

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -2,0

set out "output_flux.eps"

plot "./output/flux/H+.dat"    u 2:1 w l lw 3           title "H^+ up  ", \
     "./output/flux/H+.dat"    u 3:1 w linespoints lw 3 title "H^+ down", \
     "./output/flux/H2+.dat"   u 2:1 w l lw 3           title "H_2^+ up  ", \
     "./output/flux/H2+.dat"   u 3:1 w linespoints lw 3 title "H_2^+ down", \
     "./output/flux/H3+.dat"   u 2:1 w l lw 3           title "H_3^+ up  ", \
     "./output/flux/H3+.dat"   u 3:1 w linespoints lw 3 title "H_3^+ down", \
     "./output/flux/CH4+.dat"  u 2:1 w l lw 3           title "CH_4^+ up  ", \
     "./output/flux/CH4+.dat"  u 3:1 w linespoints lw 3 title "CH_4^+ down", \
     "./output/flux/H.dat"     u 2:1 w l lw 3           title "H up  ", \
     "./output/flux/H.dat"     u 3:1 w linespoints lw 3 title "H down"

set xlabel 'velocity[m/s]'
set xr [1E-10:3E5]
set out "output_fluxvel.eps"

plot "./output/flux/v_H+.dat"    u 2:1 w l lw 3           title "H^+ up  ", \
     "./output/flux/v_H+.dat"    u 3:1 w linespoints lw 3 title "H^+ down", \
     "./output/flux/v_H2+.dat"   u 2:1 w l lw 3           title "H_2^+ up  ", \
     "./output/flux/v_H2+.dat"   u 3:1 w linespoints lw 3 title "H_2^+ down", \
     "./output/flux/v_H3+.dat"   u 2:1 w l lw 3           title "H_3^+ up  ", \
     "./output/flux/v_H3+.dat"   u 3:1 w linespoints lw 3 title "H_3^+ down", \
     "./output/flux/v_CH4+.dat"  u 2:1 w l lw 3           title "CH_4^+ up  ", \
     "./output/flux/v_CH4+.dat"  u 3:1 w linespoints lw 3 title "CH_4^+ down", \
     "./output/flux/v_H.dat"     u 2:1 w l lw 3           title "H up  ", \
     "./output/flux/v_H.dat"     u 3:1 w linespoints lw 3 title "H down"

################################ output D, K #############################
reset 

set term postscript enhanced color
#set term aqua
set log x
set xlabel 'diffusion coefficient [m^2/s]'
set ylabel 'altitude [km]'
set xr [1E0:3E15]
set yr [0:3000]
#set xr [0.01:5]

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,20"
set title font"Arial,20"

set key spacing 1.5

set format x "10^{%L}"
#unset key

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -2,0

set out "output_diffusion_coefficient.eps"

plot "./output/flux/K_eddy.dat"  u 2:1 w l lw 3 title "K", \
     "./output/flux/D_H2+.dat"  u 2:1 w l lw 3 title "D_{H_2^+}"
    