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
set key font "Arial,25"
set title font"Arial,20"

set key spacing 1.5

set format x "10^{%L}"
#unset key

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -2,0

set out "ozone.eps"
set xr [1E10:1E16]
set yr [0:100]

plot "< paste ./stable201709/output/density/num/OH.dat ./stable201709/output/density/num/HO2.dat" using ($2+$4):1 w l lw 5 dt (10,5) lc "blue" title "HO_x (stable)", \
     "./stable201709/output/density/num/O3.dat"    u 2:1 w l lw 5 dt (10,5) lc "forest-green" title "O_3 (stable)", \
     "< paste ./sep201709/output/density/num/OH.dat ./sep201709/output/density/num/HO2.dat" using ($2+$4):1 w l lw 5 lc "blue" title "HO_x (SEP)", \
     "./sep201709/output/density/num/O3.dat"    u 2:1 w l lw 5 lc "forest-green" title "O_3 (SEP)"


