set datafile missing '0.0000000000'
set multiplot layout 2,2
set key right reverse Right
set logscale y
f(x) = k * (C/x)**(l*x)
k=0.01
l=0.5
C=1
set title "n=50"
fit f(x) "../OutputScaleFree_Normed/histogram_N50_M1.dat" using 1:3 via k, l, C
plot f(x) title "", "../OutputScaleFree_Normed/histogram_N50_M1.dat" using 1:3:5 with yerrorbars title ""
g(x) = k_2 * (C_2/x)**(l_2*x)
k_2=0.01
l_2=0.5
C_2=1
set title "n=100"
fit g(x) "../OutputScaleFree_Normed/histogram_N100_M1.dat" using 1:3 via k_2, l_2, C_2
plot g(x) title "", "../OutputScaleFree_Normed/histogram_N100_M1.dat" using 1:3:5 with yerrorbars title ""
h(x) = k_3 * (C_3/x)**(l_3*x)
k_3=0.01
l_3=0.5
C_3=1
set title "n=200"
fit h(x) "../OutputScaleFree_Normed/histogram_N200_M1.dat" using 1:3 via k_3, l_3, C_3
plot h(x) title "", "../OutputScaleFree_Normed/histogram_N200_M1.dat" using 1:3:5 with yerrorbars title ""
i(x) = k_4 * (C_4/x)**(l_4*x)
k_4=0.01
l_4=0.5
C_4=1
set title "n=400"
fit i(x) "../OutputScaleFree_Normed/histogram_N400_M1.dat" using 1:3 via k_4, l_4, C_4
plot i(x) title "", "../OutputScaleFree_Normed/histogram_N400_M1.dat" using 1:3:5 with yerrorbars title ""
