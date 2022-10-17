set xlabel "t"
set ylabel "E"
set logscale y
input="energy.dat"
plot input using 1:2 ti "4701 steps" w lp
