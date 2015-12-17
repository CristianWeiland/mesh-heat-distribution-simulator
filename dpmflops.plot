set termoption dash
set grid
set logscale x 2
set logscale y 2
set xlabel 'Points'
set ylabel 'DP MFLOP/s'
set title 'Points x DP MFLOP/s'
set style line 1 lt 1 lc rgb '#9ecae1' lw 2
set style line 2 lt 1 lc rgb "#bdbdbd" lw 2
set style line 3 lt 2 lc rgb '#3182bd' lw 2
set style line 4 lt 2 lc rgb "#636363" lw 2
set key center bottom box

plot 'dpmflops.dat' u 1:2 w l t 'Residue - Version 1' ls 3\
   , 'dpmflops.dat' u 1:4 w l t 'SOR - Version 1' ls 1\
   , 'dpmflops.dat' u 1:3 w l t 'Residue - Version 2' ls 4\
   , 'dpmflops.dat' u 1:5 w l t 'SOR - Version 2' ls 2

pause -1