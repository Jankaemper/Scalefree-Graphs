set datafile missing '0.000000'
set multiplot layout 2,2
set key right reverse Right
set title "nu=0.5"
set logscale y
f(x) = k * (C/x)**(l*x)
k=1.2
l=2.5
C=1.34
fit f(x) "OutputPlanar_Normed/histogram_N50.dat" via k, l, C
plot f(x), "OutputPlanar_Normed/histogram_N50.dat"
g(x) = k_2 * (C_2/x)**(l_2*x)
k_2=0.1
l_2=1.1
C_2=10
fit g(x) "OutputPlanar_Normed/histogram_N100.dat" via k_2, l_2, C_2
plot g(x), "OutputPlanar_Normed/histogram_N100.dat"
h(x) = k_3 * (C_3/x)**(l_3*x)
k_3=4.5
l_3=1.2
C_3=11
fit h(x) "OutputPlanar_Normed/histogram_N200.dat" via k_3, l_3, C_3
plot h(x), "OutputPlanar_Normed/histogram_N200.dat"
i(x) = k_4 * (C_4/x)**(l_4*x)
k_4=8.6
l_4=1.2
C_4=13.1
fit i(x) "OutputPlanar_Normed/histogram_N400.dat" via k_4, l_4, C_4
plot i(x), "OutputPlanar_Normed/histogram_N400.dat"
