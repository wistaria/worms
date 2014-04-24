set xlabel "h" 
set key left bottom
pl 'diag-ene.dat' w l, 'diag-mag.dat' w l, 'diag-mag2.dat' w l, 'diag-smag2.dat' w l, 'worms-ene.dat' w e, 'worms-mag.dat' w e, 'worms-mag2.dat' w e, 'worms-smag2.dat' w e
