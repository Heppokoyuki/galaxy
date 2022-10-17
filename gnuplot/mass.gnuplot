set xlabel "r"
set ylabel "M"
set logscale x
plot "mass0.dat" using 1:2 ti "t=0" w lp
replot "mass1.dat" using 1:2 ti "t=0.003" w lp
replot "mass2.dat" using 1:2 ti "t=0.005" w lp
replot "mass3.dat" using 1:2 ti "t=0.008" w lp
replot "mass4.dat" using 1:2 ti "t=0.010" w lp
replot "mass5.dat" using 1:2 ti "t=0.013" w lp
replot "mass6.dat" using 1:2 ti "t=0.016" w lp
replot "mass7.dat" using 1:2 ti "t=0.019" w lp
replot "mass8.dat" using 1:2 ti "t=0.022" w lp
replot "mass9.dat" using 1:2 ti "t=0.024" w lp
