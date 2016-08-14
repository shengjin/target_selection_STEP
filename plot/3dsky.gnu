reset

set terminal postscript landscape enhanced color "Times-Roman" 15
set output "test.eps"


set multiplot

set origin 0,0

set xrange [-15:15]
set yrange [-15:15]
set zrange [-15:15]

set dummy u,v
set angles degrees
unset key
set parametric
set view 60, 40, 0.8, 1.2
set samples 64, 64
#set samples 32, 32
#set isosamples 16, 16
set isosamples 12, 12
set mapping spherical
set yzeroaxis linetype 0 linewidth 1.000
set ticslevel 0
#set title "3D version using spherical coordinate system" 
set urange [ -90.0000 : 90.0000 ] noreverse nowriteback
set vrange [ 0.00000 : 360.000 ] noreverse nowriteback


splot cos(u)*cos(v)*15,cos(u)*sin(v)*15,sin(u)*15 with lines lt 0 lw 0.1 lc rgb "black" 


set datafile separator ","

unset mapping 
unset parametric
unset dummy
unset angles
unset samples
unset isosamples

set xlabel "Distance [ pc ]" offset 2,-0.6
set ylabel "Distance [ pc ]" offset -2,-0.6
set zlabel "Distance [ pc ]" rotate 

set xrange [-15:15]
set yrange [-15:15]
set zrange [-15:15]

set palette @MATLAB1
#set colorbox vertical
set colorbox vertical user origin 0.88,0.282 size 0.02,0.45
set format cb ""
set cbtics 0,10,100
set cbrange [0.25:1.46]
set cbtics add ("M0" 1.41, "K0" 0.82, "G0" 0.59, "F0" 0.31, "A0" 0.0) scale 0.7
splot "tgplot.csv" u (cos($9)*cos($8*15)*$10):(cos($9)*sin($8*15)*$10):(sin($9)*$10):(10**(-0.06*$14)):17 w  points pt 7  ps var lt palette



unset multiplot




