set term png size 600,600
set logscale x
set format x "10^{%L}"
unset logscale y
unset format y
set xlabel "time"
set out "plots/N_0.png"
plot "./data/0/N_exact.dat" u 1:2 w p pt 7 ps 1 lc rgb "black" title "N_{0 exact}",\
    "./data/0/N_exact.dat" u 1:3 w p pt 7 ps 1 lc rgb "red" title "N_{1 exact}" ,\
    "./data/0/N_exact.dat" u 1:4 w p pt 7 ps 1 lc rgb "blue" title "N_{2 exact}" ,\
    "./data/0/N.dat" u 1:2 w l lc rgb "black" title "N_{0 num}",\
    "./data/0/N.dat" u 1:3 w l lc rgb "red" title "N_{1 num}",\
    "./data/0/N.dat" u 1:4 w l lc rgb "blue" title "N_{2 num}",\

set out

set term png size 600,600
set logscale x
set logscale y
set format y "10^{%L}"
set out "plots/ts0.png"
plot "./data/0/timestep.dat" u 1:2 w l title "TOL=0.0001" lc rgb "black"

set out





set term png size 600,600
set logscale x
unset logscale y
unset format y
set out "plots/N_1.png"
plot "./data/1/N_exact.dat" u 1:2 w p pt 7 ps 1 lc rgb "black" title "N_{0 exact}",\
    "./data/1/N_exact.dat" u 1:3 w p pt 7 ps 1 lc rgb "red" title "N_{1 exact}" ,\
    "./data/1/N_exact.dat" u 1:4 w p pt 7 ps 1 lc rgb "blue" title "N_{2 exact}" ,\
    "./data/1/N.dat" u 1:2 w l lc rgb "black" title "N_{0 num}",\
    "./data/1/N.dat" u 1:3 w l lc rgb "red" title "N_{1 num}",\
    "./data/1/N.dat" u 1:4 w l lc rgb "blue" title "N_{2 num}",\

set out

set term png size 600,600
set logscale x
set logscale y
set format y "10^{%L}"
set out "plots/ts1.png"
plot "./data/1/timestep.dat" u 1:2 w l lc rgb "black" title "TOL=0.000001"

set out





set term png size 600,600
set logscale x
unset logscale y
unset format y
set out "plots/N_2.png"
plot "./data/2/N_exact.dat" u 1:2 w p pt 7 ps 1 lc rgb "black" title "N_{0 exact}",\
    "./data/2/N_exact.dat" u 1:3 w p pt 7 ps 1 lc rgb "red" title "N_{1 exact}",\
    "./data/2/N_exact.dat" u 1:4 w p pt 7 ps 1 lc rgb "blue" title "N_{2 exact}",\
    "./data/2/N.dat" u 1:2 w l lc rgb "black" title "N_{0 num}",\
    "./data/2/N.dat" u 1:3 w l lc rgb "red" title "N_{1 num}",\
    "./data/2/N.dat" u 1:4 w l lc rgb "blue" title "N_{2 num}",\

set out

set term png size 600,600
set logscale x
set logscale y
set format y "10^{%L}"
set out "plots/ts2.png"
plot "./data/2/timestep.dat" u 1:2 w l lc rgb "black" title "TOL=0.001"

set out