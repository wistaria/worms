set xlabel "h" 
set yrange [0:]
set key left
pl 'diag-mag.dat' w l, 'diag-mag2.dat' w l, 'diag-smag2.dat' w l, 'worms-mag.dat' w e, 'worms-mag2.dat' w e, 'worms-smag2.dat' w e
