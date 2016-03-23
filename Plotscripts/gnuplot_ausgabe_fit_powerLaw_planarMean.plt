set datafile missing '0.0000000000'
set multiplot layout 1,2
set key right reverse Right
set title "Mean Values"
f(x) = x**(B)
B = 0.9
fit f(x) "OutputPlanar_Normed/meanValues.dat" via B
plot f(x), "OutputPlanar_Normed/meanValues.dat"
g(x) = log(x)
plot g(x), "OutputPlanar_Normed/meanValues.dat"
