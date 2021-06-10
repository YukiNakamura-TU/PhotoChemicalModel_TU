
reset

dir = 'no_metal_Hill'


set terminal postscript enhanced color
set view 45,0,1.3, 0.8
set mapping spherical
set angle degrees
set hidden3d front
set nozeroaxis
set noborder
set nogrid


set tics font "Arial,25"
set cblabel font "Arial,25"
set cblabel offset 6, 0


unset key
unset xtics
unset ytics
unset ztics


set cbrange[-2:0.3]
set cbtics ("0.01" -2, "0.02" -1.69897, "0.03" -1.522878, "0.05" -1.30103, "0.07" -1.1549, "0.1" -1.0, "0.2" -0.69897, "0.3" -0.522878, "0.5" -0.30103, "0.7" -0.1549, "1.0" 0, "2.0" 0.3010, "3" 0.4771 )

set cblabel "{/Arial Pedersen conductance }{/Symbol S}{/Arial _P [S]}"
#set palette defined ( -1 '#000030', 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
set palette defined(0"#000030",0.8"#00008b",1.8"#2ca9e1",3"#008b00",4.2"#ffff00",5"#eb6101",5.5"#8b0000")
#set palette defined (8 "blue", 10 "green", 12 "red")
#set palette defined (-25 "blue", 0 "white", 10 "red")
#set palette rgbformulae 22, 13, -31


set xyplane at -2
set xrange [-1.2:1.2]
set yrange [-1.2:1.2]
set zrange [0:1.2]
set colorbox user vertical origin .83,.15 size .04,.6
set parametric
set isosamples 360
set urange [0:360]
set vrange [-90:90]


# Parametric functions for the sphere
r = 1.01
r = 1.003
fx(v,u) = r*cos(v)*cos(u)
fy(v,u) = r*cos(v)*sin(u)
fz(v)   = r*sin(v)
#set size 1200,900


set out './'.dir.'/output/conductivity/sigma_H.eps'
splot './'.dir.'/output/conductivity/hi_sigma_HP0_ij.dat' using 2:1:(1):3 w pm3d , \
       r*cos(0)   *cos(u), r*cos(0)   *sin(u), r*sin(0)    w l ls 1 lc rgb "black", \
       r*cos(7.5) *cos(u), r*cos(7.5) *sin(u), r*sin(7.5)  w l ls 1 lc rgb "black", \
       r*cos(15)  *cos(u), r*cos(15)  *sin(u), r*sin(15)   w l ls 1 lc rgb "black", \
       r*cos(22.5)*cos(u), r*cos(22.5)*sin(u), r*sin(22.5) w l ls 1 lc rgb "black", \
       r*cos(30)  *cos(u), r*cos(30)  *sin(u), r*sin(30)   w l ls 1 lc rgb "black", \
       r*cos(37.5)*cos(u), r*cos(37.5)*sin(u), r*sin(37.5) w l ls 1 lc rgb "black", \
       r*cos(45)  *cos(u), r*cos(45)  *sin(u), r*sin(45)   w l ls 1 lc rgb "black", \
       r*cos(52.5)*cos(u), r*cos(52.5)*sin(u), r*sin(52.5) w l ls 1 lc rgb "black", \
       r*cos(60)  *cos(u), r*cos(60)  *sin(u), r*sin(60)   w l ls 1 lc rgb "black", \
       r*cos(67.5)*cos(u), r*cos(67.5)*sin(u), r*sin(67.5) w l ls 1 lc rgb "black", \
       r*cos(75)  *cos(u), r*cos(75)  *sin(u), r*sin(75)   w l ls 1 lc rgb "black", \
       r*cos(82.5)*cos(u), r*cos(82.5)*sin(u), r*sin(82.5) w l ls 1 lc rgb "black", \
       r*cos(0)   *cos(u), r*sin(0)   *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(15)  *cos(u), r*sin(15)  *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(30)  *cos(u), r*sin(30)  *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(45)  *cos(u), r*sin(45)  *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(60)  *cos(u), r*sin(60)  *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(75)  *cos(u), r*sin(75)  *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(90)  *cos(u), r*sin(90)  *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(105) *cos(u), r*sin(105) *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(120) *cos(u), r*sin(120) *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(135) *cos(u), r*sin(135) *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(150) *cos(u), r*sin(150) *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(165) *cos(u), r*sin(165) *cos(u), r*sin(u)    w l ls 1 lc rgb "black"


set out './'.dir.'/output/conductivity/sigma_P.eps'
splot './'.dir.'/output/conductivity/hi_sigma_HP0_ij.dat' using 2:1:(1):4 w pm3d , \
       r*cos(0)   *cos(u), r*cos(0)   *sin(u), r*sin(0)    w l ls 1 lc rgb "black", \
       r*cos(7.5) *cos(u), r*cos(7.5) *sin(u), r*sin(7.5)  w l ls 1 lc rgb "black", \
       r*cos(15)  *cos(u), r*cos(15)  *sin(u), r*sin(15)   w l ls 1 lc rgb "black", \
       r*cos(22.5)*cos(u), r*cos(22.5)*sin(u), r*sin(22.5) w l ls 1 lc rgb "black", \
       r*cos(30)  *cos(u), r*cos(30)  *sin(u), r*sin(30)   w l ls 1 lc rgb "black", \
       r*cos(37.5)*cos(u), r*cos(37.5)*sin(u), r*sin(37.5) w l ls 1 lc rgb "black", \
       r*cos(45)  *cos(u), r*cos(45)  *sin(u), r*sin(45)   w l ls 1 lc rgb "black", \
       r*cos(52.5)*cos(u), r*cos(52.5)*sin(u), r*sin(52.5) w l ls 1 lc rgb "black", \
       r*cos(60)  *cos(u), r*cos(60)  *sin(u), r*sin(60)   w l ls 1 lc rgb "black", \
       r*cos(67.5)*cos(u), r*cos(67.5)*sin(u), r*sin(67.5) w l ls 1 lc rgb "black", \
       r*cos(75)  *cos(u), r*cos(75)  *sin(u), r*sin(75)   w l ls 1 lc rgb "black", \
       r*cos(82.5)*cos(u), r*cos(82.5)*sin(u), r*sin(82.5) w l ls 1 lc rgb "black", \
       r*cos(0)   *cos(u), r*sin(0)   *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(15)  *cos(u), r*sin(15)  *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(30)  *cos(u), r*sin(30)  *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(45)  *cos(u), r*sin(45)  *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(60)  *cos(u), r*sin(60)  *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(75)  *cos(u), r*sin(75)  *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(90)  *cos(u), r*sin(90)  *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(105) *cos(u), r*sin(105) *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(120) *cos(u), r*sin(120) *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(135) *cos(u), r*sin(135) *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(150) *cos(u), r*sin(150) *cos(u), r*sin(u)    w l ls 1 lc rgb "black", \
       r*cos(165) *cos(u), r*sin(165) *cos(u), r*sin(u)    w l ls 1 lc rgb "black"
