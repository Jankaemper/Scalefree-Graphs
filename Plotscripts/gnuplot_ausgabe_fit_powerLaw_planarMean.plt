set datafile missing '0.0000000000'
set key right reverse Right
set title "Mean Values"
set logscale y
f(x) = x**(B)
B = 0.9
fit f(x) "../OutputPlanar_Normed/meanValues.dat" using 1:2 via B
plot f(x) title "", "../OutputPlanar_Normed/meanValues.dat" using 1:2:3 with yerrorbars title ""
