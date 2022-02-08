set term pngcairo

set output "sum.png"
set xlabel "it"
set ylabel "S"
set logscale x
set xrange[1:1000]
set title "S(it)"
plot 'sum.dat' i 0 u 1:2 w l ls 1 lw 3 t "k = 16", \
  'sum.dat' i 1 u 1:2 w l ls 2 lw 3 t "k = 8", \
  'sum.dat' i 2 u 1:2 w l ls 3 lw 3 t "k = 4", \
  'sum.dat' i 3 u 1:2 w lines ls 4 lw 3 t "k = 2", \
  'sum.dat' i 4 u 1:2 w l ls 5 lw 8 lc rgb "black" t "k = 1"

set nologscale x
set xlabel "x"
set ylabel "y"
set xrange[0:25.6]
set yrange[0:25.6]
set tics out
set tics nomirror

set output "k16.png"
set title "k = 16"
plot 'k16.dat' u 1:2:3 with image

set output "k8.png"
set title "k = 8"
plot 'k8.dat' u 1:2:3 with image

set output "k4.png"
set title "k = 4"
plot 'k4.dat' u 1:2:3 with image

set output "k2.png"
set title "k = 2"
plot 'k2.dat' u 1:2:3 with image

set output "k1.png"
set title "k = 1"
plot 'k1.dat' u 1:2:3 with image
