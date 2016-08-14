
set terminal postscript landscape enhanced color "Times-Roman" 15
#set term postscript eps enhanced color dashed lw 2 "Times-Roman" 22
#set term x11 enhanced font "times,25,italic"
set output "test.eps"

reset

unset key
set view 60, 40, 0.8, 1.2

set xrange [-20:20]
set yrange [-20:20]
set zrange [-20:20]

set datafile separator ","

set xlabel "Distance [ pc ]"
set ylabel "Distance [ pc ]"
set zlabel "Distance [ pc ]" rotate 

#14 mag
#17 color index
#set pm3d
#set palette @softpalette
#set palette @MATLAB
set palette @MATLAB1

#set colorbox
set colorbox vertical
set format cb ""
set cbrange [0.25:1.46]
set cbtics add ("M0" 1.41, "K0" 0.82, "G0" 0.59, "F0" 0.31, "A0" 0.0) scale 0.7

# x:y:z:size:color
splot "tgplot.csv" u (0):(0):(0):(cos($9)*cos($8*15)*$10):(cos($9)*sin($8*15)*$10):(sin($9)*$10) with vector nohead lt -1 lw 0.01  lc rgb "#DCDCDC", \
"tgplot.csv" u (cos($9)*cos($8*15)*$10):(cos($9)*sin($8*15)*$10):(sin($9)*$10):(10**(-0.06*$14)):17 w  points pt 7  ps var lt palette



#set format cb ""
#set cblabel "Fraction of Envelop Lost"
#set cbtics (""0.4,""0.6,""0.8,""1.0,""1.2,""1.4,""1.6,""1.8,""2.0,""2.2,""2.4)
#
#splot 'ref_red1e9.dat' u 51:5::($7<1000?1-$54:1/0) with points palette pt 7 lt 3 ps 0.3






#set format y "%.1t{/Symbol \264}10^%L"
#set format cb "10^{%L}"
#set xtics 0, 2e6  (mtic at 0, 2e6, 4e6,….)
#open circles: pt 65 filled cycles: pt 31 ps 0.5 for 3D maps

