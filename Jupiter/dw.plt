reset
#set terminal aqua enhanced
set terminal postscript enhanced

dir1 = "./metal_Hill"
dir2 = "./no_metal_Hill"

Mdot = "300"
LT   = "12"

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,20"
set title font"Arial,20"
set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -5,0
#unset key

set xr [1:100]
set yr [0:1]

set xlabel "L"
set ylabel "{/Symbol w}_i / {/Symbol W}_J"
set out 'figures/wi_Mdot'.Mdot.'_LT'.LT.'.eps'
plot 	dir1.'/output/plasmav/wi_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:2  w l lw 5 title "Metal", \
		dir2.'/output/plasmav/wi_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:2  w l lw 5 title "NO Metal"

set yr [0:400]
set ylabel "v_{/Symbol f} [km/s]"
set out 'figures/vi_Mdot'.Mdot.'_LT'.LT.'.eps'
plot 	dir1.'/output/plasmav/vi_Mdot'.Mdot.'_LT'.LT.'.dat'    u 1:2 w l lw 5 title "Rigid corotation", \
		dir1.'/output/plasmav/vi_Mdot'.Mdot.'_LT'.LT.'.dat'    u 1:3 w l lw 5 title "Metal", \
		dir2.'/output/plasmav/vi_Mdot'.Mdot.'_LT'.LT.'.dat'    u 1:2 w l lw 5 title "NO Metal"

set yr [-3:3]
set ylabel "j_{||i} [{/Symbol m}A/m^2]"
set out 'figures/j_para_i_Mdot'.Mdot.'_LT'.LT.'.eps'
plot 	dir1.'/output/plasmav/j_para_i_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:($2*1e6) w l lw 5 title "Metal", \
		dir2.'/output/plasmav/j_para_i_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:($2*1e6) w l lw 5 title "NO Metal"


#unset key
#set xr [1:100]
#set yr [0:1]
#set out 'figures/wi_Mdot'.Mdot.'_LT'.LT.'.eps'
#plot 	dir2.'/output/plasmav/w_old.dat' w l lw 5 title "1/L", \
#		dir2.'/output/plasmav/w_new.dat' w l lw 5 title "b"



set terminal postscript enhanced color

set xr [1:100]
set yr [0:1]

set ylabel "{/Symbol w}_i / {/Symbol W}_J"
set out 'figures/wi_iter_Mdot'.Mdot.'_LT'.LT.'.eps'
plot 	dir1.'/output/plasmav/wi_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:2  w l lw 2 title "iter = 1", \
		dir1.'/output/plasmav/wi_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:3  w l lw 2 title "iter = 2", \
		dir1.'/output/plasmav/wi_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:4  w l lw 2 title "iter = 3", \
		dir1.'/output/plasmav/wi_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:5  w l lw 2 title "iter = 4"

set yr [-0.3:0.3]
set ylabel "j_{||i} [{/Symbol m}A/m^2]"
set out 'figures/j_iter_Mdot'.Mdot.'_LT'.LT.'.eps'
plot 	dir1.'/output/plasmav/j_para_i_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:($2*1e6) w l lw 2 title "iter = 1", \
		dir1.'/output/plasmav/j_para_i_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:($3*1e6) w l lw 2 title "iter = 2", \
		dir1.'/output/plasmav/j_para_i_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:($4*1e6) w l lw 2 title "iter = 3", \
		dir1.'/output/plasmav/j_para_i_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:($5*1e6) w l lw 2 title "iter = 4", \



set yr [0:5]
set ylabel "{/Symbol S}_P [S]"
set out 'figures/Sigma_P_Mdot'.Mdot.'_LT'.LT.'.eps'
plot 	dir1.'/output/plasmav/Sigma_P_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:2 w l lw 2 title "iter = 1", \
		dir1.'/output/plasmav/Sigma_P_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:3 w l lw 2 title "iter = 2", \
		dir1.'/output/plasmav/Sigma_P_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:4 w l lw 2 title "iter = 3", \
		dir1.'/output/plasmav/Sigma_P_Mdot'.Mdot.'_LT'.LT.'.dat'  u 1:5 w l lw 2 title "iter = 4", \