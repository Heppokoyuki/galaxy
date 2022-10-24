set xlabel "t"
set ylabel "E"
penergy="penergy.dat"
kenergy="kenergy.dat"
plot penergy using 1:2 ti "potential energy" w lp
replot kenergy using 1:2 ti "kinetic energy" w lp
replot "< paste penergy.dat kenergy.dat" using 1:($2+$4) ti "sum" w lp
