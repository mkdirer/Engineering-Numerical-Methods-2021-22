set term png
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1
set xlabel "x"
set ylabel "y"
set pm3d map

iterations = "100 200 500 1000 2000"

do for [it in iterations] {
  image_name = "images/d2T_it=".it.".png"
  set output image_name
  set title "liczba it=".it
  file = "data/out_d2T_".it.".dat"
  splot file i 0 u 1:2:3
} 
do for [it in iterations] {
  image_name = "images/T_it=".it.".png"
  set output image_name
  set title "T(x,y) liczba it=".it
  file = "data/out_T_".it.".dat"
  splot file i 0 u 1:2:3
} 
