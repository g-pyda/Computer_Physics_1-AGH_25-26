# ------------------------------------------------------------1

set term pngcairo size 600,500
set output './plots/1/T10.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/1/tempDist_10.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 600,500
set output './plots/1/T100.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/1/tempDist_100.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 600,500
set output './plots/1/T1000.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/1/tempDist_1000.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 600,500
set output './plots/1/T10000.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/1/tempDist_10000.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term png size 600,500
set out "plots/1/T_evol.png"

unset xrange
unset yrange

plot "./data/1/tempDist.dat" u 1:2 w l title "T(t)"

set out





# ------------------------------------------------------------2

set term pngcairo size 600,500
set output './plots/2/T10.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/2/tempDist_10.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 600,500
set output './plots/2/T100.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/2/tempDist_100.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 600,500
set output './plots/2/T1000.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/2/tempDist_1000.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 600,500
set output './plots/2/T10000.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/2/tempDist_10000.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term png size 600,500
set out "plots/2/T_evol.png"

unset xrange
unset yrange

plot "./data/2/tempDist.dat" u 1:2 w l title "T(t)"

set out






# ------------------------------------------------------------3

set term pngcairo size 600,500
set output './plots/3/T10.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/3/tempDist_10.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 600,500
set output './plots/3/T100.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/3/tempDist_100.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 600,500
set output './plots/3/T1000.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/3/tempDist_1000.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 600,500
set output './plots/3/T10000.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/3/tempDist_10000.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term png size 600,500
set out "plots/3/T_evol.png"

unset xrange
unset yrange

plot "./data/3/tempDist.dat" u 1:2 w l title "T(t)"

set out







# ------------------------------------------------------------4

set term pngcairo size 600,500
set output './plots/4/T2500.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/4/tempDist_2500.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 600,500
set output './plots/4/T3000.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/4/tempDist_3000.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 600,500
set output './plots/4/T3500.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/4/tempDist_3500.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 600,500
set output './plots/4/T4000.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"

set xrange [0:10]
set yrange [0:10]
    
set palette defined (0 "light-blue", 0.25 "green", 0.5 "yellow", 0.75 "orange", 1 "red")

plot "./data/4/tempDist_4000.dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term png size 600,500
set out "plots/4/T_evol.png"

unset xrange
unset yrange

plot "./data/4/tempDist.dat" u 1:2 w l title "T(t)"

set out


# -------------------------------------------------------------------5

set terminal pngcairo size 1200,800 enhanced
set output "./plots/5/energy.png"

set xlabel "time"
set ylabel "Energy"
set yrange [0:50]

set grid
set key top left

plot "./data/5/tempDist.dat" using 1:3 with lines lw 2 title "E_{supplied}", \
     "./data/5/tempDist.dat" using 1:4 with lines lw 2 title "E_{lost}"

set out

# ---------------------------------------------------------------------6

set terminal pngcairo size 1200,800 enhanced
set output "./plots/6/energy.png"

set xlabel "time"
set ylabel "Energy"
set yrange [0:50]

set grid
set key top left

plot "./data/6/tempDist.dat" using 1:3 with lines lw 2 title "E_{supplied}", \
     "./data/6/tempDist.dat" using 1:4 with lines lw 2 title "E_{lost}"

set out