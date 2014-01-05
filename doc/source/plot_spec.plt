set terminal postscript enhanced "Helvetica" 16 color eps
set output 'spec.eps'
# convert -density 100 spec.eps spec.png

set size 2,1

set multiplot layout 1,2

set size 1,1
set origin 0,0
set xrange [0.25:0.6]
set xtics 0.05
set yrange [0:1]
set ytics 0.2
set xlabel 'Normalized frequency (c/a)'
set ylabel 'Transmission coefficient'

set object 1 rect from 0.504,0 to 0.507,1 fs empty border rgb "red"

plot 'spec.dat' u 1:2 notitle w l ls 1 lw 3 lc rgb '#000000'


# make subplot
set origin 1,0
set size 1,1

unset arrow
unset object 1

set xrange [0.504:0.507]
set xtics 0.001
set yrange [0:1]
set ytics 0.2
unset ylabel

plot 'spec.dat' u 1:2 notitle w l ls 1 lw 3 lc rgb '#000000'
unset multiplot
set size 2,1