reset
set terminal png size 800,600
set xlabel 'time'
set output 'energyOMP.png'
plot 'energyOMP.dat' u 1:2 w line title 'kinetic', 'energyOMP.dat' u 1:3 w line title 'potential'

reset
set terminal png size 800,600
set xlabel 'time'
set ylabel 'Magnetization'
set output 'magnetOMP.png'
plot 'magnetOMP.dat' u 1:3 w line notitle

reset
set terminal png size 800,600
set xlabel 'position'
set ylabel 'momentum'
set output 'initialPhaseOMP.png'
plot 'initialPhaseOMP.dat' u 1:3 w dot notitle

reset
set terminal png size 800,600
set xlabel 'position'
set ylabel 'momentum'
set output 'finalPhaseOMP.png'
plot 'finalPhaseOMP.dat' u 1:3 w dot notitle
