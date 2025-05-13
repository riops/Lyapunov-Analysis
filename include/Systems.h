#ifndef SYSTEMS_H
#define SYSTEMS_H

#include <vector>

// Functions to compute indices
int indexX(double i, double l, double m, int N);

int indexXX(double i1, double l1, double m1, double i2, double l2, double m2,
            int N);

std::vector<long double>
averagedEquationsPolarizationBasisSymmetryReducedParallel(
    const std::vector<long double> &allVectors,
    const std::vector<std::vector<double>> &H, long double mu,
    long double lambda, long double R, int N);

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
