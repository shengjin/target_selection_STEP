reset

set terminal postscript portrait enhanced color "Times-Roman" 12
set output "test.eps"
#set term postscript eps enhanced color dashed lw 2 "Times-Roman" 22
#set term x11 enhanced font "times,25,italic"

set bmargin 0
set rmargin 0
set lmargin 0
set tmargin 0

nbin=100		#number of bins
max=17.0 		#max value
min=0.6 	        #min value
binwidth=(max-min)/nbin	#bin width

set log x
set log x2

set xrange [min:max]
set x2range [min:max]



#function used to map all the values into the middle of the bins
rint(x)=(x-int(x)>0.9999)?int(x)+1:int(x)
binmap(x,width)=width*rint(x/width)+width/2.0


set boxwidth binwidth

#plot 'ref_red1e8.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps , \
#'ref_red1e9.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps , \
#'ref_red1e10.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps

unset key

# linetype color for label
# heavy color
set style line 1  linetype 1 linecolor rgb "#1d4599"  linewidth 1.000 pointtype  3 pointsize default #blue
set style line 2  linetype 1 linecolor rgb "#11ad34"  linewidth 1.000 pointtype  3 pointsize default #green
set style line 3  linetype 1 linecolor rgb "#e62b17"  linewidth 1.000 pointtype  3 pointsize default #red
set style line 4  linetype 1 linecolor rgb "#000000"  linewidth 2.000 pointtype  3 pointsize default #red
set style line 5  linetype 1 linecolor rgb "#000000"  linewidth 2.000 pointtype  3 pointsize default #red
# soft color

xtotalsize = 1.1
ytotalsize = 0.35
Xxtotalsize = 1.2

set size Xxtotalsize,ytotalsize

# X Y ratio of each figure
# set size ratio 0.90
ytotratio=0.90

by = 1 
nby = 1	# number in the y direction
bx = 1	# number of blocks in the x-direction
nbx = 3	# plots per block
yintrablock = 0.0 #ysep = -0.32 # y seperation
xintrablock = 0.10	      # x seperation

set multiplot 

# begin x,y, point
ox = 0.11
oy1 = 0.04
# space left at the end of the multiple x,y plots
oyend = 0.05
oxend = -0.30

ysize = (ytotalsize-oy1-oyend-nby*yintrablock)/nby
xsize = (xtotalsize-ox-oxend-nbx*xintrablock)/nbx

set size xsize,ysize
dox = xsize*(1+xintrablock)
doy = ysize*(1+yintrablock)


# set ration of the two sub-figues in y direction
# ylower/(ytotal)
ysubratio=0.8
ysubratio1=0.8*1.15
# manually change
#ymanually=-0.03
#ymanually=-0.031
#ymanually=-0.0315
#ymanually=-0.0305 limit
#ymanually=-0.0329 # limit
ymanually=-0.090 # limit
ymanually1=-0.031
ymanually2=-0.031
#ymanually4=-0.0317225 # 31722-317226
ymanually4=ymanually

#set mxtics 5
set mytics 2
set mxtics 10
#set mx2tics 10

set ytics scale 0.6

#set style histogram rowstacked

# third row ############################################

unset xlabel
set format y #"10^{%L}"
#set xtics 0,2,13.5
set format x
set xtics  scale 0.65 nomirror
set  xtics add ("0.6" 0.6, "2" 2, "5" 5, "15" 15) scale 0.65
unset x2tics
set origin ox,oy1
unset label
unset label
set ylabel "N [ Planets ]" offset 1,1
set xlabel "Radius  [ R_{Earth} ]" offset -0
set size ratio ytotratio*ysubratio1
set yrange [0:300]
set ytics 0,50,275
set border 11 lt -1
plot 'AW_CD754_evap_xray_tau0p01/ref_red1e8.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 1, \
'AW_CD754_evap_xray_tau0p01/ref_red1e8_1.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 3, \
'AW_CD754_evap_xray_tau0p01/ref_red1e8_0p1.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 2
unset ylabel
set size ratio ytotratio*(1.0-ysubratio)
set origin ox,oy1+ysize*ysubratio+ymanually
set bmargin 0
set yrange [300:1000]
set ytics 300,300,900
set border 14 lt -1
unset xtics
set x2tics
#set mx2tics 10
#set  xtics add ("" 0.6, "" 2, "" 5, "" 15) scale 0.65
set x2tics format ""  scale 0.65 nomirror
set  x2tics add ("" 0.6, "" 2, "" 5, "" 15) scale 0.65
#set xtics format ""
#set x2tics format "" scale 0.65 nomirror
unset xlabel
plot 'AW_CD754_evap_xray_tau0p01/ref_red1e8.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 1, \
'AW_CD754_evap_xray_tau0p01/ref_red1e8_1.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 3, \
'AW_CD754_evap_xray_tau0p01/ref_red1e8_0p1.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 2
set size ratio ytotratio*(1.0-ysubratio)
set origin ox,oy1+ysize*ysubratio+ymanually4
set bmargin 0
set yrange [300:1000]
set ytics 300,300,900
set border 1 lt 2 lc -1
unset xtics
unset x2tics
unset ylabel
set ytics format ""
set label 1 "0.1 Gyr" at graph 0.92, graph 0.53 right font "Helvetica-Italic,12" textcolor rgb "#000000"
#set label 1 "0.1 Gyr" at graph 0.92, graph 0.53 right font "Helvetica-Italic,12" textcolor rgb "#000000"
#set x2tics 0,2,12
#set x2tics format ""
#set x2tics format "" scale 0.65 nomirror
plot  -1





set origin ox+dox,oy1
#set xtics 0,2,13.5
set format x
set xtics  scale 0.65 nomirror
set  xtics add ("0.6" 0.6, "2" 2, "5" 5, "15" 15) scale 0.65
unset x2tics
unset x2tics
set format y ""
unset label
set ylabel ""
unset label
set xlabel "Radius  [ R_{Earth} ]" offset -0
set size ratio ytotratio*ysubratio1
set yrange [0:300]
set ytics 0,50,275
set border 11 lt -1
plot 'AW_CD754_evap_xray_tau0p01/ref_red1e9.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 1, \
'AW_CD754_evap_xray_tau0p01/ref_red1e9_1.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 3, \
'AW_CD754_evap_xray_tau0p01/ref_red1e9_0p1.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 2
set size ratio ytotratio*(1.0-ysubratio)
set origin ox+dox,oy1+ysize*ysubratio+ymanually
set bmargin 0
set yrange [300:1000]
set ytics 300,300,900
set border 14 lt -1
unset xtics
set x2tics
set x2tics format "" scale 0.65 nomirror
set  x2tics add ("" 0.6, "" 2, "" 5, "" 15) scale 0.65
unset xlabel
plot 'AW_CD754_evap_xray_tau0p01/ref_red1e9.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 1, \
'AW_CD754_evap_xray_tau0p01/ref_red1e9_1.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 3, \
'AW_CD754_evap_xray_tau0p01/ref_red1e9_0p1.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 2
set size ratio ytotratio*(1.0-ysubratio)
set origin ox+dox,oy1+ysize*ysubratio+ymanually4
set bmargin 0
set yrange [300:1000]
set ytics 300,300,900
set border 1 lt 2 lc -1
unset xtics
unset x2tics
unset ylabel
set ytics format ""
set label 1 "1 Gyr" at graph 0.90, graph 0.53 right font "Helvetica-Italic,12" textcolor rgb "#000000"
#set x2tics 0,2,12
#set x2tics format ""
#set x2tics format "" scale 0.65 nomirror
plot  -1




set origin ox+2*dox,oy1
#set xtics 0,2,12
set format x
set xtics  scale 0.65 nomirror
set  xtics add ("0.6" 0.6, "2" 2, "5" 5, "15" 15) scale 0.65
unset x2tics
set format y ""
unset label
set ylabel ""
set xlabel "Radius  [ R_{Earth} ]" offset 0
set size ratio ytotratio*ysubratio1
set yrange [0:300]
set ytics 0,50,275
set border 11 lt -1
plot 'AW_CD754_evap_xray_tau0p01/ref_red5e9.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 1, \
'AW_CD754_evap_xray_tau0p01/ref_red5e9_1.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 3, \
'./kepler.dat' using (binmap($8,binwidth)):(1.0) smooth freq with histeps ls 4, \
'./kepler.txt' using (binmap($8,binwidth)):(1.0) smooth freq with histeps ls 5, \
'AW_CD754_evap_xray_tau0p01/ref_red5e9_0p1.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 2
set size ratio ytotratio*(1.0-ysubratio)
set origin ox+2*dox,oy1+ysize*ysubratio+ymanually
set bmargin 0
set yrange [300:1000]
set ytics 300,300,900
set border 14 lt -1
unset xtics
set x2tics format "" scale 0.65 nomirror
set  x2tics add ("" 0.6, "" 2, "" 5, "" 15) scale 0.65
unset xlabel
plot 'AW_CD754_evap_xray_tau0p01/ref_red5e9.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 1, \
'AW_CD754_evap_xray_tau0p01/ref_red5e9_1.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 3, \
'AW_CD754_evap_xray_tau0p01/ref_red5e9_0p1.dat' using (binmap($39,binwidth)):(1.0) smooth freq with histeps ls 2
set size ratio ytotratio*(1.0-ysubratio)
set origin ox+2*dox,oy1+ysize*ysubratio+ymanually4
set bmargin 0
set yrange [300:1000]
set ytics 300,300,900
set border 1 lt 2 lc -1
unset xtics
unset x2tics
unset ylabel
set ytics format ""
set label 1 "5 Gyr" at graph 0.90, graph 0.53 right font "Helvetica-Italic,12" textcolor rgb "#000000"
#set x2tics 0,2,12
#set x2tics format ""
#set x2tics format "" scale 0.65 nomirror
plot  -1




unset multiplot

