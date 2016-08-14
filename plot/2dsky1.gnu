
set term postscript eps enhanced color dashed lw 2 "Times-Roman" 22

set out "test.eps"

set xrange [25:28]
set yrange [63.7:64.2]

set size ratio 1

unset key

set datafile separator ","

set xlabel "Right Ascension    ( Ep = J2000 )"
set ylabel "Declination    ( Ep = J2000 )"

set pm3d map

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
splot "tgplot.csv" u ($8*15):9::(10**(-0.06*$14)):10 w  points pt 7  ps var lt palette




