set key right reverse Right
set title "nu=0.5"
set logscale y
f(x) = k * (C/x)**(l*x)
k=1.2
l=0.01
C=1.34
fit f(x) "Output/histogram_N50_M1.dat" via k, l, C
plot f(x), "Output/histogram_N50_M1.dat"
