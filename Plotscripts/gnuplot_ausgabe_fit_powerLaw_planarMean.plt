set datafile missing '0.0000000000'
set key right reverse Right
set title "Mean Values"
f(x) = x**(B)
B = 0.9
fit f(x) "../OutputPlanar_Normed/meanValues2.dat" via B
plot f(x), "../OutputPlanar_Normed/meanValues2.dat" title "Powerlaw"
