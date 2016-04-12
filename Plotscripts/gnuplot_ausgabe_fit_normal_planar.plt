set datafile missing '0.0000000000'
set multiplot layout 2,2
set key right reverse Right
set title "Normalverteilung"
set logscale y
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
fit f(x) "../OutputPlanar_Normed/histogram_N50.dat" via std, nu
plot f(x), "../OutputPlanar_Normed/histogram_N50.dat" title 'n=50'
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
fit f(x) "../OutputPlanar_Normed/histogram_N100.dat" via std, nu
plot f(x), "../OutputPlanar_Normed/histogram_N100.dat" title 'n=100'
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
fit f(x) "../OutputPlanar_Normed/histogram_N200.dat" via std, nu
plot f(x), "../OutputPlanar_Normed/histogram_N200.dat" title 'n=200'
f(x) = 1/((2*pi*std**2)**0.5) * exp(-(x-nu)**2/(2*std**2))
std = 50
nu = 6
fit f(x) "../OutputPlanar_Normed/histogram_N400.dat" via std, nu
plot f(x), "../OutputPlanar_Normed/histogram_N400.dat" title 'n=400'

