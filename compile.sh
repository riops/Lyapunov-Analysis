#!/bin/bash
module purge
module load comp/gcc/12.3.0
module load lib/openmpi/4.1.6
module load lib/gsl/2.7.0

# mpicxx -O3 -march=native -fopenmp -std=c++17 -I./include src/*.cpp -o builds/main
mpicxx -O3 -march=native -fopenmp -std=c++17 -I./include src/*.cpp -lgsl -lgslcblas -o builds/main
