set multiplot layout 2,2
set key right reverse Right
set title "nu=0.5"
set logscale y
set yrange [1:]
set xrange [0:40]
f(x) = k * (C/x)**(l*x)
k=1.2
l=0.01
C=1.34
fit f(x) "Output/histogram_N50_M1.dat" via k, l, C
plot f(x), "Output/histogram_N50_M1.dat"
g(x) = k_2 * (C_2/x)**(l_2*x)
k_2=1.2
l_2=0.01
C_2=1.34
fit g(x) "Output/histogram_N100_M1.dat" via k_2, l_2, C_2
plot g(x), "Output/histogram_N100_M1.dat"
h(x) = k_3 * (C_3/x)**(l_3*x)
k_3=1.2
l_3=0.01
C_3=1.34
fit h(x) "Output/histogram_N200_M1.dat" via k_3, l_3, C_3
plot h(x), "Output/histogram_N200_M1.dat"
i(x) = k_4 * (C_4/x)**(l_4*x)
k_4=1.2
l_4=0.0001
C_4=1.34
fit i(x) "Output/histogram_N400_M1.dat" via k_4, l_4, C_4
plot i(x), "Output/histogram_N400_M1.dat"
