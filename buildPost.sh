rm b.exe
rm *.mod
rm *.o
mpicxx -c postProcess.cpp -std=c++11 -mkl -O3 -mcmodel=large
mpicxx -o ./b.exe postProcess.cpp -std=c++11 -mkl -O3 -mcmodel=large
