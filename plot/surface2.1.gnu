# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 500, 350 
# set output 'surface2.1.png'
set dummy u,v
set key bmargin center horizontal Right noreverse enhanced autotitles nobox
set parametric
set view 45, 50, 1, 1
set isosamples 50, 10
set hidden3d back offset 1 trianglepattern 3 undefined 1 altdiagonal bentover
#set ztics -1.00000,0.25,1.00000 norangelimit
set title "Parametric Sphere" 
set urange [ -1.57080 : 1.57080 ] noreverse nowriteback
set vrange [ 0.00000 : 6.28319 ] noreverse nowriteback
set datafile separator ","
#splot cos(u)*cos(v)*20,cos(u)*sin(v)*20,sin(u)*20, \
splot "tgplot.csv" u cos($9)*cos($8*15)*$10:cos($9)sin($8*15)*$10:sin($9)*$10 
