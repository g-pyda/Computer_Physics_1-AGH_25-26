set term png size 600,600
set out "TR_no_drag.png"
plot "./data/no_drag_force/frajectory_n10.000000.dat" u 2:3 w l title "n=10",\
    "./data/no_drag_force/frajectory_n20.000000.dat" u 2:3 w l title "n=20",\
    "./data/no_drag_force/frajectory_n50.000000.dat" u 2:3 w l title "n=50",\
    "./data/no_drag_force/frajectory_n100.000000.dat" u 2:3 w l title "n=100",\
    "./data/no_drag_force/frajectory_n200.000000.dat" u 2:3 w l title "n=200",\
    "./data/no_drag_force/frajectory_n500.000000.dat" u 2:3 w l title "n=500"

set out

set term png size 600,600
set logscale x 10
set logscale y 10
set out "ERR_no_drag.png"
plot "./data/no_drag_force/errors.dat" u 2:1 with linespoints pointtype 7 pointsize 1.5 title "error" 

set out



set term png size 600,600
unset logscale x
unset logscale y
set out "TR_drag_no_altit.png"
plot "./data/drag_force_no_altitude/frajectory_D0.000000.dat" u 2:3 w l title "D=0.0",\
    "./data/drag_force_no_altitude/frajectory_D0.000100.dat" u 2:3 w l title "D=0.0001",\
    "./data/drag_force_no_altitude/frajectory_D0.000200.dat" u 2:3 w l title "D=0.0002",\
    "./data/drag_force_no_altitude/frajectory_D0.000500.dat" u 2:3 w l title "D=0.0005",\
    "./data/drag_force_no_altitude/frajectory_D0.001000.dat" u 2:3 w l title "D=0.0010"

set out



set term png size 600,600
set out "X_MAX_drag+angle_D0.0.png"
set xrange [15:65]
plot "./data/drag_force+angle/xmax_D0.000000.dat" u 1:2 w l title "D=0.0"

set out

set term png size 600,600
set out "X_MAX_drag+angle_D0.001.png"
set xrange [15:65]
plot "./data/drag_force+angle/xmax_D0.001000.dat" u 1:2 w l title "D=0.001"

set out

set term png size 600,600
set out "X_MAX_drag+angle_D0.002.png"
set xrange [15:65]
plot "./data/drag_force+angle/xmax_D0.002000.dat" u 1:2 w l title "D=0.002"

set out



set term png size 600,600
set out "TR_altitude_corr.png"
unset xrange
plot "./data/altitude_correction/trajectory_a0.000000theta35.dat" u 2:3 w l title "a=0.0 theta=35",\
    "./data/altitude_correction/trajectory_a0.000000theta45.dat" u 2:3 w l title "a=0.0 theta=45",\
    "./data/altitude_correction/trajectory_a0.006500theta35.dat" u 2:3 w l title "a=0.0065 theta=35",\
    "./data/altitude_correction/trajectory_a0.006500theta45.dat" u 2:3 w l title "a=0.0065 theta=45"

set out