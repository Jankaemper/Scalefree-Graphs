set datafile missing '0.0000000000'
set multiplot layout 2,2
set key right reverse Right
set logscale y
set title "Normalverteilung m=2"
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
fit f(x) "../OutputScaleFree_Normed/histogram_N50_M2.dat" via std, nu
plot f(x), "../OutputScaleFree_Normed/histogram_N50_M2.dat" title 'n=50'
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 55
nu =1
fit f(x) "../OutputScaleFree_Normed/histogram_N100_M2.dat" via std, nu
plot f(x), "../OutputScaleFree_Normed/histogram_N100_M2.dat" title 'n=100'
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
fit f(x) "../OutputScaleFree_Normed/histogram_N200_M2.dat" via std, nu
plot f(x), "../OutputScaleFree_Normed/histogram_N200_M2.dat" title 'n=200'
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
fit f(x) "../OutputScaleFree_Normed/histogram_N400_M2.dat" via std, nu
plot f(x), "../OutputScaleFree_Normed/histogram_N400_M2.dat" title 'n=400'



