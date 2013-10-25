set terminal postscript eps enhanced monochrome 'Times-Roman' 24
set encoding iso_8859_1
set xlabel 'Eigenvalue number'
set ylabel '%'
set boxwidth 0.5
set xrange [0.2:LAST.5]
set output 'EIGENVALUES.eps'
plot 'eigenvalues.txt' using 1:3 title 'eigenvalues' with boxes
