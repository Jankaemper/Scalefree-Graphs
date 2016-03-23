set datafile missing '0.0000000000'
set multiplot layout 2,2
set key right reverse Right
set title "Normalverteilung"
set logscale y
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
fit f(x) "Output_Normed/histogram_N50_M1.dat" via std, nu
plot f(x), "Output_Normed/histogram_N50_M1.dat"
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
fit f(x) "Output_Normed/histogram_N100_M1.dat" via std, nu
plot f(x), "Output_Normed/histogram_N100_M1.dat"
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
fit f(x) "Output_Normed/histogram_N200_M1.dat" via std, nu
plot f(x), "Output_Normed/histogram_N200_M1.dat"
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
fit f(x) "Output_Normed/histogram_N400_M1.dat" via std, nu
plot f(x), "Output_Normed/histogram_N400_M1.dat"

