#ifndef SYSTEMS_H
#define SYSTEMS_H

#include <vector>

// Functions to compute indices
int indexX(double a, double i, double l, double m, int N);

int indexXX(double a1, double i1, double l1, double m1, double a2, double i2,
            double l2, double m2, int N);

int indexXInt(int a, int i, int l, int m, int N);

int indexXXInt(int a1, int i1, int l1, int m1, int a2, int i2, int l2, int m2,
               int N);

std::vector<long double>
averagedEquationsPolarizationBasisSymmetryReducedParallel(
    const std::vector<long double> &allVectors);

std::vector<long double>
averagedEquationsPolarizationBasisSymmetryReducedParallelTrial(
    const std::vector<long double> &allVectors,
    const std::vector<std::vector<double>> &H);

std::vector<long double>
averagedEquationsPolarizationBasisSymmetryReducedParallelMpi(
    const std::vector<long double> &allVectors, int world_rank, int world_size);

std::vector<long double> GenerateSigmaAnsatzI(long double E, long double lambda,
                                              long double R, long double mu,
                                              int N);

std::vector<long double>
GenerateSigmaAnsatzII(long double E, long double lambda, long double R, int N);

std::vector<long double> GenerateInitialConditionAnsatzI(long double E, int N,
                                                         long double R,
                                                         long double mu,
                                                         long double lambda);

std::vector<long double> GenerateInitialConditionAnsatzII(long double E, int N,
                                                          long double R,
                                                          long double lambda);

std::vector<long double> GenerateInitialConditionReduced(int ansatz, int N,
                                                         long double E,
                                                         long double lambda,
                                                         long double R,
                                                         long double mu);

#endif // SYSTEMS_H
