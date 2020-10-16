rm a.exe
rm *.mod
rm *.o
mpicxx -c main.cpp -std=c++11 -mkl -O3 -mcmodel=medium
mpicxx -o ./a.exe main.cpp -std=c++11 -mkl -O3 -mcmodel=medium
