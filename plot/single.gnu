# add 8 ref, and label

reset

set terminal postscript landscape enhanced color "Times-Roman" 15
set output "test.eps"


#set xrange [-20:20]
#set yrange [-20:20]

set size ratio 1.0

set xlabel "Right Ascension    ( Ep = J2000 )"
set ylabel "Declination    ( Ep = J2000 )"

set xtics format "%.1f"
set ytics format "%.1f"


set pm3d map
set palette @softpalette

set grid ytics lt 0 lw 1 lc rgb "#999999"
set grid xtics lt 0 lw 1 lc rgb "#999999"

set colorbox vertical user origin 0.78,0.192 size 0.02,0.62
set cblabel "Apparent Magnitude" offset 1,0
#set cbrange [0.25:1.46]
#set cbtics add ("M0" 1.41, "K0" 0.82, "G0" 0.59, "F0" 0.31, "A0" 0.0) scale 0.7

unset key

#set title "HD 69830"
set label "HD 69830" at graph 0.53,0.49 front

# x:y:z:size:color
#splot "hip2id_40693_ref" u 1:2::(5*10**(-0.06*$4)):(abs($9)+abs($10)) w  points pt 7  ps var lt palette
splot "hip2id_40693_ref" u 1:2::(5*10**(-0.06*$4)):4 w  points pt 7  ps var lt palette, \
"hip2id_40693" u ($1*15):2::(5*10**(-0.06*$3)):3 w  points pt 6  ps var lt palette, \
"hip2id_40693" u ($1*15):2::(10*10**(-0.06*$3)):3 w  points pt 3  ps var lt palette




