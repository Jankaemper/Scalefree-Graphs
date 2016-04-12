set datafile missing '0.0000000000'
set multiplot layout 2,2
set key right reverse Right
set title "Powerlaw m=2"
f(x) = k * (C/x)**(l*x)
k=1.2
l=2.5
C=1.34
fit f(x) "../Output_constant/histogram_N50_M2_k0.500000.dat" via k, l, C
plot f(x), "../Output/histogram_N50_M2.dat" title 'n=50'
g(x) = k_2 * (C_2/x)**(l_2*x)
k_2=0.1
l_2=1.1
C_2=10
fit g(x) "../Output_constant/histogram_N100_M2_k0.500000.dat" via k_2, l_2, C_2
plot g(x), "../Output/histogram_N100_M2.dat" title 'n=100'
h(x) = k_3 * (C_3/x)**(l_3*x)
k_3=4.5
l_3=1.2
C_3=11
fit h(x) "../Output_constant/histogram_N200_M2_k0.500000.dat" via k_3, l_3, C_3
plot h(x), "../Output/histogram_N200_M2.dat" title 'n=200'
i(x) = k_4 * (C_4/x)**(l_4*x)
k_4=8.6
l_4=1.2
C_4=13.1
fit i(x) "../Output_constant/histogram_N400_M2_k0.500000.dat" via k_4, l_4, C_4
plot i(x), "../Output/histogram_N400_M2.dat" title 'n=400'
