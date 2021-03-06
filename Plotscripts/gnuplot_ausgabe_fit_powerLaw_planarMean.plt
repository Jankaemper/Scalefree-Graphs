set term png
set output '../Results/wErrorbars_Power_PlanarMean.png'
set datafile missing '0.0000000000'
set key right reverse Right
set logscale y
set title "Mean Values"
f(x) = x**(B)
B = 0.9
fit f(x) "../OutputPlanar_Normed/meanValues.dat" using 1:2 via B
plot f(x) title "", "../OutputPlanar_Normed/meanValues.dat" using 1:2:3 with yerrorbars title ""
