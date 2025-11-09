set term png size 600,600
set out "./plots/tr_0,95T.png"
set xrange [-0.45:0.55]
set yrange [-0.45:0.55]
plot "./data/0,95T/trajectory.dat" u 1:2 w l notitle, \
'-' w p pt 7 ps 1.5
0 0
e

set out

set term png size 600,600
set out "./plots/tr_100T.png"
plot "./data/100T/trajectory.dat" u 1:2 w l notitle, \
'-' w p pt 7 ps 1.5
0 0
e

set out

set term png size 600,600
set output "./plots/tr+peri+ap.png"

plot \
    "./data/alpha0,01/trajectory.dat" u 1:2 w l title "trajectory", \
    '-' w p pt 7 ps 1.5 title "Sun", \
    '-' w l lt 2 lw 2 lc rgb "blue" title "aphelia", \
    '-' w l lt 2 lw 2 lc rgb "blue" notitle, \
    '-' w l lt 2 lw 2 lc rgb "blue" notitle, \
    '-' w l lt 2 lw 2 lc rgb "blue" notitle, \
    '-' w l lt 1 lw 2 lc rgb "red" title "perihelia", \
    '-' w l lt 1 lw 2 lc rgb "red" notitle, \
    '-' w l lt 1 lw 2 lc rgb "red" notitle

# --- Sun (point) ---
0 0
e

# --- Aphelia lines ---
-0.267574 -0.0721055
0 0
e
-0.194117 -0.197772
0 0
e
-0.0671027 -0.268872
0 0
e
0.0784252 -0.26579
0 0
e

# --- Perihelia lines ---
0.412831 0.242708
0 0
e
0.232876 0.418455
0 0
e
-0.0113265 0.478757
0 0
e

set output


a = 176.432
f(x) = a*x
set term png size 600,600
set format x "%1.1te%+02T"
set xtics rotate by -45
unset xrange
unset yrange
set out "./plots/prec_ang_vel.png"
plot "./data/prec_ang_vel/precession_angular_velocity.dat" u 1:2 w p pt 7 ps 1 title "numerical", \
    f(x) title "f(x) = a * x" with lines

set out