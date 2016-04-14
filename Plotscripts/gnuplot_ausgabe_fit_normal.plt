set datafile missing '0.0000000000'
set multiplot layout 3,1
set key right reverse Right
set logscale y
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
set title "n=50"
fit f(x) "Output_Normed/histogram_N50_M1.dat" using 1:3 via std, nu
plot f(x) title "", "Output_Normed/histogram_N50_M1.dat" using 1:3:5 with yerrorbars title ""
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
set title "n=100"
fit f(x) "Output_Normed/histogram_N100_M1.dat" using 1:3 via std, nu
plot f(x) title "", "Output_Normed/histogram_N100_M1.dat" using 1:3:5 with yerrorbars title ""
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
set title "n=200"
fit f(x) "Output_Normed/histogram_N200_M1.dat" using 1:3 via std, nu
plot f(x) title "", "Output_Normed/histogram_N200_M1.dat" using 1:3:5 with yerrorbars title ""
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
set title "n=400"
fit f(x) "Output_Normed/histogram_N400_M1.dat" using 1:3 via std, nu
plot f(x) title "", "Output_Normed/histogram_N400_M1.dat" using 1:3:5 with yerrorbars title ""
