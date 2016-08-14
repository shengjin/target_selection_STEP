
reset
set terminal postscript landscape enhanced color "Times-Roman" 15
#set term postscript eps enhanced color dashed lw 2 "Times-Roman" 22
#set term x11 enhanced font "times,25,italic"
set output "test.eps"



unset key
set view 60, 40, 0.8, 1.2

set xrange [0:360]
set yrange [-90:90]

set ytics -90,30,90
set mytics 3
set xtics 0,30,360
set mxtics 3

set datafile separator ","

set xlabel "Right Ascension    ( Ep = J2000 )"
set ylabel "Declination    ( Ep = J2000 )"

#14 mag
#17 color index
set pm3d map
#set palette @softpalette
#set palette @MATLAB
set palette @MATLAB1
set grid ytics lt 0 lw 1 lc rgb "#999999"
set grid xtics lt 0 lw 1 lc rgb "#999999"
#set colorbox
#set colorbox vertical at 1.0,0.1 size 0.5,8
set colorbox vertical user origin 0.88,0.192 size 0.02,0.65
set cblabel "Distance  ( pc )"
#set format cb ""
#set cbrange [0.25:1.46]
#set cbtics add ("M0" 1.41, "K0" 0.82, "G0" 0.59, "F0" 0.31, "A0" 0.0) scale 0.7

# x:y:z:size:color
plot "tgplot.csv" u ($8*15):9 w  points pt 1 
#splot "tgplot.csv" u ($8*15):9::(10**(-0.06*$14)):10 w  points pt 7  ps var lt palette




