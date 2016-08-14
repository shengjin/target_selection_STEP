#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.6 patchlevel 5    last modified February 2014
#    	Build System: Linux x86_64
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2014
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
set terminal postscript landscape enhanced defaultplex \
   leveldefault color colortext \
   dashed dashlength 1.0 linewidth 1.0 butt noclip \
   nobackground \
   palfuncparam 2000,0.003 \
   "Times-Roman" 15  fontscale 1.0 
 set output 'test.eps'
unset clip points
set clip one
unset clip two
set bar 1.000000 front
set border 31 front linetype -1 linewidth 1.000
set timefmt z "%d/%m/%y,%H:%M"
set zdata 
set timefmt y "%d/%m/%y,%H:%M"
set ydata 
set timefmt x "%d/%m/%y,%H:%M"
set xdata 
set timefmt cb "%d/%m/%y,%H:%M"
set timefmt y2 "%d/%m/%y,%H:%M"
set y2data 
set timefmt x2 "%d/%m/%y,%H:%M"
set x2data 
set boxwidth 2.8 absolute
#set style fill  empty border
set style fill solid 0.50 noborder
set style rectangle back fc  lt -3 fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02, first 0, 0 
set style ellipse size graph 0.05, 0.03, first 0 angle 0 units xy
set dummy x,y
set format x ""
set format y "% g"
set format x2 "% g"
set format y2 "% g"
set format z "% g"
set format cb "% g"
set format r "% g"
set angles radians
unset grid
set raxis
set key title ""
set key inside right top vertical Right noreverse enhanced autotitles nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset key
unset label
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title  offset character 0, 0, 0
unset logscale
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
set view 60, 30, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
unset contour
set clabel '%8.3g'
set macros
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 1 1,1
set origin 0,0
set style data histograms
set style function lines
set xzeroaxis linetype -2 linewidth 1.000
set yzeroaxis linetype -2 linewidth 1.000
set zzeroaxis linetype -2 linewidth 1.000
set x2zeroaxis linetype -2 linewidth 1.000
set y2zeroaxis linetype -2 linewidth 1.000
set ticslevel 0.5
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set xtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set xtics -5.00000,10,100.000 norangelimit
set xtics add  ("F" 0.200000, "G" 1.20000, "K" 2.20000)
set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set ytics autofreq  norangelimit
set ztics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set ztics autofreq  norangelimit
set nox2tics
set noy2tics
set cbtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set cbtics autofreq  norangelimit
set rtics axis in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set rtics autofreq  norangelimit
set title "" 
set title  offset character 0, 0, 0 font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  offset character 0, 0, 0 font "" norotate
set rrange [ * : * ] noreverse nowriteback
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "Stellar Type" 
set xlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set x2label "" 
set x2label  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set xrange [ -0.500000 : 3.00000 ] noreverse nowriteback
set x2range [ * : * ] noreverse nowriteback
set ylabel "Number" 
set ylabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set y2label "" 
set y2label  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set yrange [ 0.00000 : 200.000 ] noreverse nowriteback
set y2range [ * : * ] noreverse nowriteback
set zlabel "" 
set zlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse nowriteback
set cblabel "" 
set cblabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set cbrange [ * : * ] noreverse nowriteback
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "en_US.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit noerrorvariables
GNUTERM = "qt"
MATLAB = "defined (0  0.0 0.0 0.5,                    1  0.0 0.0 1.0,                    2  0.0 0.5 1.0,                    3  0.0 1.0 1.0,                    4  0.5 1.0 0.5,                    5  1.0 1.0 0.0,                    6  1.0 0.5 0.0,                    7  1.0 0.0 0.0,                    8  0.5 0.0 0.0 )"
MATLAB1 = "defined (  0  0.0 0.0 1.0,                    1  0.0 0.5 1.0,                    2  0.0 1.0 1.0,                    3  0.5 1.0 0.5,                    4  1.0 1.0 0.0,                    5  1.0 0.5 0.0,                    6  1.0 0.0 0.0 )"
MATLAB2 = "defined (0  0.0 0.0 0.5,                    1  0.0 0.0 1.0,                    2  0.0 0.5 1.0,                    3  0.0 1.0 1.0,                    4  0.5 1.0 0.5,                    5  1.0 1.0 0.0,                    6  1.0 0.5 0.0,                    7  1.0 0.0 0.0)"
MATLABrv = "defined ( 0  0.5 0.0 0.0,                    1  1.0 0.0 0.0,                    2  1.0 0.5 0.0,                    3  1.0 1.0 0.0,                    4  0.5 1.0 0.5,                    5  0.0 1.0 1.0,                    6  0.0 0.5 1.0,                    7  0.0 0.0 1.0, \t           8  0.0 0.0 0.5 )"
MATLABrv1 = "defined ( 0  1.0 0.0 0.0,                    1  1.0 0.5 0.0,                    2  1.0 1.0 0.0,                    3  0.5 1.0 0.5,                    4  0.0 1.0 1.0,                    5  0.0 0.5 1.0, \t           6  0.0 0.0 0.5 )"
softpalette = " defined(0       0.2314  0.2980  0.7529,0.03125 0.2667  0.3529  0.8000,0.0625  0.3020  0.4078  0.8431,0.09375 0.3412  0.4588  0.8824,0.125   0.3843  0.5098  0.9176,0.15625 0.4235  0.5569  0.9451,0.1875  0.4667  0.6039  0.9686,0.21875 0.5098  0.6471  0.9843,0.25    0.5529  0.6902  0.9961,0.28125 0.5961  0.7255  1.0000,0.3125  0.6392  0.7608  1.0000,0.34375 0.6824  0.7882  0.9922,0.375   0.7216  0.8157  0.9765,0.40625 0.7608  0.8353  0.9569,0.4375  0.8000  0.8510  0.9333,0.46875 0.8353  0.8588  0.9020,0.5     0.8667  0.8667  0.8667,0.53125 0.8980  0.8471  0.8196,0.5625  0.9255  0.8275  0.7725,0.59375 0.9451  0.8000  0.7255,0.625   0.9608  0.7686  0.6784,0.65625 0.9686  0.7333  0.6275,0.6875  0.9686  0.6941  0.5804,0.71875 0.9686  0.6510  0.5294,0.75    0.9569  0.6039  0.4824,0.78125 0.9451  0.5529  0.4353,0.8125  0.9255  0.4980  0.3882,0.84375 0.8980  0.4392  0.3451,0.875   0.8706  0.3765  0.3020,0.90625 0.8353  0.3137  0.2588,0.9375  0.7961  0.2431  0.2196,0.96875 0.7529  0.1569  0.1843,1       0.7059  0.0157  0.1490)"
p "color.dat" u 3
#    EOF
