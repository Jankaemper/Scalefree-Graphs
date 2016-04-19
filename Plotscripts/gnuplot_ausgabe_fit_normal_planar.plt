set term png
set output "../Results/wErrorbars_Normal_Planar.png"
set datafile missing '0.0000000000'
set multiplot layout 2,2
set key right reverse Right
set logscale y
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
set title "n=50"
set xrange [0:8]
fit f(x) "../OutputPlanar_Normed/histogram_N50.dat" using 1:3 via std, nu
plot f(x) title "", "../OutputPlanar_Normed/histogram_N50.dat" using 1:3:5 with yerrorbars title ""
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
set title "n=100"
set xrange [0:11]
fit f(x) "../OutputPlanar_Normed/histogram_N100.dat" using 1:3 via std, nu
plot f(x) title "", "../OutputPlanar_Normed/histogram_N100.dat" using 1:3:5 with yerrorbars title ""
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
set title "n=200"
set xrange [0:15]
fit f(x) "../OutputPlanar_Normed/histogram_N200.dat" using 1:3 via std, nu
plot f(x) title "", "../OutputPlanar_Normed/histogram_N200.dat" using 1:3:5 with yerrorbars title ""
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
set title "n=400"
set xrange [0:21]
fit f(x) "../OutputPlanar_Normed/histogram_N400.dat" using 1:3 via std, nu
plot f(x) title "", "../OutputPlanar_Normed/histogram_N400.dat" using 1:3:5 with yerrorbars title ""
