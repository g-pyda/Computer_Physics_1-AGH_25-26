# ------------------------------------------------------------initial_run

set term pngcairo size 1000,800
set output './plots/initial_run/potential.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"
set cblabel "V(x,y)"
set xrange [-4:4]
set yrange [-4:4]

set cbrange [-1:1]      
set palette defined (0 "blue", 0.5 "white", 1 "red")

plot "./data/initial_run/V(x,y).dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 1000,800
set output './plots/initial_run/error.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"
set cblabel "err(x,y)"

set cbrange [-0.025:0.025]
set palette defined (0 "blue", 0.5 "white", 1 "red")

plot "./data/initial_run/err(x,y).dat" u 1:2:3 with image

set out



# ------------------------------------------------------------initial_run_with_boundaries

set term pngcairo size 1000,800
set output './plots/initial_run_with_boundaries/potential.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"
set cblabel "V(x,y)"
set xrange [-4:4]
set yrange [-4:4]

set cbrange [-1:1]      
set palette defined (0 "blue", 0.5 "white", 1 "red")

plot "./data/initial_run_with_boundaries/V(x,y).dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 1000,800
set output './plots/initial_run_with_boundaries/error.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"
set cblabel "err(x,y)"

set cbrange [-0.025:0.025]
set palette defined (0 "blue", 0.5 "white", 1 "red")

plot "./data/initial_run_with_boundaries/err(x,y).dat" u 1:2:3 with image

set out



# ------------------------------------------------------------homogeneous_bc

set term pngcairo size 1000,800
set output './plots/homogeneous_bc/potential.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"
set cblabel "V(x,y)"
set xrange [-4:4]
set yrange [-4:4]

set cbrange [-1:1]      
set palette defined (0 "blue", 0.5 "white", 1 "red")

plot "./data/homogeneous_bc/V(x,y).dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 1000,800
set output './plots/homogeneous_bc/error.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"
set cblabel "err(x,y)"

set cbrange [-0.025:0.025]
set palette defined (0 "blue", 0.5 "white", 1 "red")

plot "./data/homogeneous_bc/err(x,y).dat" u 1:2:3 with image

set out



# ------------------------------------------------------------nohomogeneous_bc

set term pngcairo size 1000,800
set output './plots/nohomogeneous_bc/potential.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"
set cblabel "V(x,y)"
set xrange [-4:4]
set yrange [-4:4]

set cbrange [-1:1]      
set palette defined (0 "blue", 0.5 "white", 1 "red")

plot "./data/nohomogeneous_bc/V(x,y).dat" u 1:2:3 with image

set out

# -----------------------------------------------

set term pngcairo size 1000,800
set output './plots/nohomogeneous_bc/error.png'

set view map
set pm3d map
set pm3d interpolate 2,2     

set xlabel "x"
set ylabel "y"
set cblabel "err(x,y)"

set cbrange [-0.025:0.025]
set palette defined (0 "blue", 0.5 "white", 1 "red")

plot "./data/nohomogeneous_bc/err(x,y).dat" u 1:2:3 with image

set out