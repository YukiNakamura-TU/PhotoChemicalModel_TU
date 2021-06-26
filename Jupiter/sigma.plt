reset
set term postscript enhanced color
set pm3d map
set pm3d corners2color c1
#set pm3d interpolate 5,5


set log y

set log cb
#set format cb "10^{%L}"

set size 0.8, 0.8

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set cblabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,20"
set title font"Arial,20"

set xlabel 'Local time'
set ylabel 'B factor'
set cblabel 'Ratio'

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -3,0
set cblabel offset 5,0

unset key
set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')

set cbr [1:100]
set yr [0.1:10]


set xr [0:24]

set out 'figures/P_ratio.eps'
splot 'sigma_ratio.dat' u 1:2:4

set out 'figures/H_ratio.eps'
splot 'sigma_ratio.dat' u 1:2:3

