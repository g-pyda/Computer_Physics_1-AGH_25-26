set term png size 1200,600
set yrange [-0.1:0.7]
set out "./plots/testing_correctness/energy.png"
plot "./data/testing_correctness/energy.dat" u 1:2 w l title "E_kin",\
    "./data/testing_correctness/energy.dat" u 1:3 w l title "E_pot",\
    "./data/testing_correctness/energy.dat" u 1:4 w l title "E_tot",\

set out

set term png size 600,600
unset yrange
set out "./plots/testing_correctness/v(t).png"
plot "./data/testing_correctness/v(t).dat" u 1:2 w l title "v(t)"

set out

set term png size 600,600
set out "./plots/testing_correctness/x(t).png"
plot "./data/testing_correctness/x(t).dat" u 1:2 w l title "x(t)"

set out

set term png size 600,600
set out "./plots/testing_correctness/x(v).png"
plot "./data/testing_correctness/x(v).dat" u 1:2 w l title "x(v)"

set out




set term png size 1200,600
set yrange [-0.1:0.7]
set out "./plots/testing_correctness-alpha0.1/energy.png"
plot "./data/testing_correctness-alpha0.1/energy.dat" u 1:2 w l title "E_kin",\
    "./data/testing_correctness-alpha0.1/energy.dat" u 1:3 w l title "E_pot",\
    "./data/testing_correctness-alpha0.1/energy.dat" u 1:4 w l title "E_tot",\

set out

set term png size 600,600
unset yrange
set out "./plots/testing_correctness-alpha0.1/x(t).png"
plot "./data/testing_correctness-alpha0.1/x(t).dat" u 1:2 w l title "num",\
    "./data/testing_correctness-alpha0.1/exact.dat" u 1:2 w p pointtype 6 title "exact"

set out

set term png size 600,600
set out "./plots/testing_correctness-alpha0.1/x(v).png"
plot "./data/testing_correctness-alpha0.1/x(v).dat" u 1:2 w l title "x(v)"

set out






set term png size 600,600
set out "./plots/energy_dissipation/x(t).png"
plot "./data/energy_dissipation-alpha/0.0001/x(t).dat" u 1:2 w l title "\alpha = 0.0001", \
    "./data/energy_dissipation-alpha/0.1/x(t).dat" u 1:2 w l title "\alpha = 0.1", \
    "./data/energy_dissipation-alpha/0.5/x(t).dat" u 1:2 w l title "\alpha = 0.5", \
    "./data/energy_dissipation-alpha/1.95/x(t).dat" u 1:2 w l title "\alpha = 1.95", \

set out





set term png size 6000,600
set out "./plots/external_force/x(t).png"
plot "./data/external_force/1.0/x(t).dat" u 1:2 w l title "x(t)"

set out





set term png size 600,600
set logscale y
set out "./plots/external_force-alpha/ampl.png"
plot "./data/external_force/alpha0.01.dat" u 1:2 w l title "alpha=0.01",\
    "./data/external_force/alpha0.1.dat" u 1:2 w l title "alpha=0.1",\
    "./data/external_force/alpha0.5.dat" u 1:2 w l title "alpha=0.5",\
    "./data/external_force/alpha1.0.dat" u 1:2 w l title "alpha=1.0",\

set out