set terminal png size 1024,800
set logscale x
set yrange [0:3.5]
g(x) = 3
f(x) = 1
set title 'HMF Identical Particles'
set xlabel 'time'
set output 'virialMoment.png'
plot 'energy.dat' u 1:(2*($2)/($3)) w line title '2<k>/<u>', f(x) notitle, 'moments.dat' title '<p^4>' w line, g(x) notitle

