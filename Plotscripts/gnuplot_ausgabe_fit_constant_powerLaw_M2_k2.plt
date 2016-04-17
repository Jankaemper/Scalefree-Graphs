set datafile missing '0.0000000000'
set term png
set output 'k2.png'
set multiplot layout 2,2
set logscale y
set key right reverse Right
f(x) = k * (C/x)**(l*x)
k=0.01
l=0.5
C=1
set title "n=50"
fit f(x) "Output_constant/histogram_N50_M2_k2.000000.dat" using 1:3 via k, l, C
plot f(x) title "", "Output_constant/histogram_N50_M2_k2.000000.dat" using 1:3:5 with yerrorbars title ""
g(x) = k_2 * (C_2/x)**(l_2*x)
k_2=0.1
l_2=1.1
C_2=1
set title "n=100"
fit g(x) "Output_constant/histogram_N100_M2_k2.000000.dat" using 1:3 via k_2, l_2, C_2
plot g(x) title "", "Output_constant/histogram_N100_M2_k2.000000.dat" using 1:3:5 with yerrorbars title ""
h(x) = k_3 * (C_3/x)**(l_3*x)
k_3=0.01
l_3=0.5
C_3=1
set title "n=200"
fit h(x) "Output_constant/histogram_N200_M2_k2.000000.dat" using 1:3 via k_3, l_3, C_3
plot h(x) title "", "Output_constant/histogram_N200_M2_k2.000000.dat" using 1:3:5 with yerrorbars title ""
i(x) = k_4 * (C_4/x)**(l_4*x)
k_4=0.01
l_4=0.5
C_4=1
set title "n=400"
fit i(x) "Output_constant/histogram_N400_M2_k2.000000.dat" using 1:3 via k_4, l_4, C_4
plot i(x) title "", "Output_constant/histogram_N400_M2_k2.000000.dat" using 1:3:5 with yerrorbars title ""
