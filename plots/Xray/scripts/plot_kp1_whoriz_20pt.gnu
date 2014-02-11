set term postscript eps enhanced
#************************************************************
#change this line to whichever path you need to the output plot
#************************************************************
set output "HERA_II_compare_kp1_whoriz_20pt.eps"

#************************************************************
#redshift ticks and labels
#************************************************************
set xtics nomirror ( "5" 5,"6" 6,"7" 7,"8" 8,"9" 9,"10" 10,"11" 11,"12" 12,"13" 13,"14" 14,"15" 15,"16" 16,"17" 17,"18" 18,"19" 19,"20" 20,"" 21,"" 22,"" 23,"" 24,"25" 25,"" 26,"" 27,"" 28,"" 29,"30" 30,"" 31,"" 32,"" 33,"" 34,"35" 35) font ",20"
set xrange [7:19] noreverse
#************************************************************
#frequency ticks and labels
#************************************************************
set x2tics ("170" 7.35,"140" 9.143,"120" 10.83,"100" 13.2,"80" 16.75,"60" 22.67,"40" 35.51) font ",20"
set x2range [7:19] noreverse
set key bottom righ font ",20" spacing 3
#************************************************************
# axis labels
#************************************************************
set xlabel "Redshift" font ", 20"
set x2label "\nFrequency [MHz]\n" font ",20"
set ylabel "\n\n {/Symbol D}_{21}^2 [mK^2]" font ",20"
set ytics font ",20"

#************************************************************
#create parsons 2013 upper limit
#************************************************************
#set label 4 at 8.0,5e3
#set label 4 "Parsons 2013" font ",20"
#set arrow from 7.7,5e3 to 7.7,1e3

#************************************************************
#array labels
#************************************************************
set label 5 "HERA-331"
set label 5 at 7.5,0.05 rotate by 25 left
#set label 10 "HERA-568"
#set label 10 at 9.6,0.025 rotate by 25 left
#set label 10 font ",30"
set label 6 at 11,6000 right
set label 6 "LOFAR-Hi" font ",20"
set label 7 at 16.5,6000 right
set label 7 "MWA" font ",20"
set label 8 "LOFAR-Lo" font ",20"
set label 8 at 16.9,6000 left
set label 5 font ",30"
set label 9 at 13,6000 right
set label 9 "PAPER" font",20"

#************************************************************
#define line properties for each theory/instrument plot
#************************************************************

set linestyle 1 lt 3 lw 8 lc rgb "red"
set linestyle 2 lt 2 lw 8 lc rgb "#009900"
set linestyle 3 lt 3 lw 8 lc rgb "#0000FF"
set linestyle 4 lt 8 lw 8 lc rgb "#FF33CC"

set linestyle 30 lt 1 lw 4 lc rgb "#000000"
set linestyle 31 lt 1 lw 5 lc rgb "#000000"
set linestyle 32 lt 1 lw 7 lc rgb "#000000" 
set linestyle 33 lt 1 lw 9 lc rgb "#000000"
set linestyle 34 lt 1 lw 11 lc rgb "#000000"
set linestyle 35 lt 1 lw 6 lc rgb "#000000"
set style fill transparent solid 0.3

#************************************************************
#scale yaxis
#************************************************************
set yr [.01:10000]
set logscale y

#************************************************************
#plot noise and theory.
#************************************************************

plot "noise_kp1_mwa_whoriz_o03.dat" using 1:(exp($2*log(10))) notitle w lines ls 30,\
"noise_kp1_whoriz_o03.dat" using 1:(exp($3*log(10))) notitle w lines ls 31,\
"noise_kp1_whoriz_o03.dat" using 1:(exp($4*log(10))) notitle w lines ls 32,\
"noise_kp1_whoriz_o03.dat" using 1:(exp($5*log(10))) notitle w lines ls 33,\
"CDM.fiducial.no.0.1.ps.txt" using 1:(exp($4*log(10))) title 'Fiducial Model' ls 1 w lines, \
"CDM.extreme.no.0.1.ps.txt" using 1:(exp($4*log(10))) title 'Extreme X-ray Heating' ls 2 w lines,\
"CDM.fiducial.heating.0.1.ps.txt" using 1:(exp($4*log(10))) title 'Cold Dark Matter Annihilation' ls 3 w lines, \
"WDM.fiducial.no.0.1.ps.txt" using 1:(exp($4*log(10))) title 'Warm Dark Matter' ls 4 w lines 
#"parsons_line.dat" using 1:2 notitle w lines ls 30