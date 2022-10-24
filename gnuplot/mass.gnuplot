set xlabel "r"
set ylabel "M"
set logscale x
set key outside
plot for [i=0:9] sprintf("mass%d.dat", i) every ::1 using 1:2 ti sprintf("t=%s", system("head -1 " . sprintf("mass%d.dat", i))) w lp
