set term png
set output '../Results/wErrorbars_Normal_M1.png'
set datafile missing '0.0000000000'
set multiplot layout 2,2
set key right reverse Right
set logscale y
set xrange [0:15]
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
set title "n=50"
fit f(x) "../OutputScaleFree_Normed/histogram_N50_M1.dat" using 1:3 via std, nu
plot f(x) title "", "../OutputScaleFree_Normed/histogram_N50_M1.dat" using 1:3:5 with yerrorbars title ""
set xrange [0:17]
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
set title "n=100"
fit f(x) "../OutputScaleFree_Normed/histogram_N100_M1.dat" using 1:3 via std, nu
plot f(x) title "", "../OutputScaleFree_Normed/histogram_N100_M1.dat" using 1:3:5 with yerrorbars  title ""
set xrange [0:21]
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
set title "n=200"
fit f(x) "../OutputScaleFree_Normed/histogram_N200_M1.dat" using 1:3 via std, nu
plot f(x) title "", "../OutputScaleFree_Normed/histogram_N200_M1.dat" using 1:3:5 with yerrorbars  title ""
set xrange [0:21]
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
set title "n=400"
fit f(x) "../OutputScaleFree_Normed/histogram_N400_M1.dat" using 1:3 via std, nu
plot f(x) title "", "../OutputScaleFree_Normed/histogram_N400_M1.dat" using 1:3:5 with yerrorbars  title ""

