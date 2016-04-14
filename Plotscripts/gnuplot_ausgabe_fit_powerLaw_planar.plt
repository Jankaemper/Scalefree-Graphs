set datafile missing '0.000000'
set multiplot layout 2,2
set key right reverse Right
set logscale y
f(x) = k * (C*C/x)**(l*l*x)
k=1.2
l=0.0001
C=0.99
set title "n=50"
fit f(x) "../OutputPlanar_Normed/histogram_N50.dat" using 1:3 via k, l, C
plot f(x) title "", "../OutputPlanar_Normed/histogram_N50.dat" using 1:3:5 with yerrorbars title ""
g(x) = k_2 * (C_2/x)**(l_2*x)
k_2=0.1
l_2=1.1
C_2=10
set title "n=100"
fit g(x) "../OutputPlanar_Normed/histogram_N100.dat" using 1:3 via k_2, l_2, C_2
plot g(x) title "", "../OutputPlanar_Normed/histogram_N100.dat" using 1:3:5 with yerrorbars title ""
h(x) = k_3 * (C_3/x)**(l_3*x)
k_3=0.005
l_3=0.01
C_3=10
set title "n=200"
fit h(x) "../OutputPlanar_Normed/histogram_N200.dat" using 1:3 via k_3, l_3, C_3
plot h(x) title "", "../OutputPlanar_Normed/histogram_N200.dat" using 1:3:5 with yerrorbars title ""
i(x) = k_4 * (C_4/x)**(l_4*x)
k_4=8.6
l_4=1.2
C_4=13.1
set title "n=400"
fit i(x) "../OutputPlanar_Normed/histogram_N400.dat" using 1:3 via k_4, l_4, C_4
plot i(x) title "", "../OutputPlanar_Normed/histogram_N400.dat" using 1:3:5 with yerrorbars title ""
