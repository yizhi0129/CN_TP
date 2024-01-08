set terminal png
set output './image/cvg_J.png'

set title "Convergence"
set xlabel "Iteration number"
set ylabel "Relative residual"

plot 'RESVEC.dat' with line title "Jacobi"