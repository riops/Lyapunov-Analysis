//
// Created by berk on 3/20/2023.
//
#include "../include/Systems.h"
#include "../include/ClebschGordanCoefficients.h"
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <random>
#include <vector>

// Functions to compute indices
int indexX(double i, double l, double m, int N) {
  int int_i = static_cast<int>(i);
  int int_l = static_cast<int>(l);
  int int_m = static_cast<int>(m);
  return ((i - 1) * (N * N) + (int_l * int_l + int_l + int_m));
}

int indexXX(double i1, double l1, double m1, double i2, double l2, double m2,
            int N) {
  int index1 = indexX(i1, l1, m1, N);
  int index2 = indexX(i2, l2, m2, N);
  int maxIndex = 2 * (N * N);        // Adjust as per your actual maximum index
  return index1 * maxIndex + index2; // Adjust based on how you store XX
}

std::vector<long double>
averagedEquationsPolarizationBasisSymmetryReducedParallel(
    const std::vector<long double> &allVectors,
    const std::vector<std::vector<double>> &H, long double mu,
    long double lambda, long double R, int N) {
  int Imax = 2 * (N * N);   // Number of averaged variables
  int CovMax = Imax * Imax; // Number of covariance terms

  // Extract covariance vectors (centered)
  std::vector<long double> XX(allVectors.begin(), allVectors.begin() + CovMax);
  std::vector<long double> PP(allVectors.begin() + CovMax,
                              allVectors.begin() + 2 * CovMax);
  std::vector<long double> XP(allVectors.begin() + 2 * CovMax,
                              allVectors.begin() + 3 * CovMax);

  std::vector<long double> result(3 * CovMax);

#pragma omp parallel for collapse(4)
  for (int lm1_int = 0; lm1_int < N * N; lm1_int++) {
    for (int i_int = 1; i_int <= 2; i_int++) {
      for (int lm2_int = 0; lm2_int < N * N; lm2_int++) {
        int l1_int = static_cast<int>(sqrt(lm1_int));
        int m1_int = lm1_int - (l1_int * l1_int + l1_int);
        double l1 = static_cast<double>(l1_int);
        double m1 = static_cast<double>(m1_int);
        double i = static_cast<double>(i_int);

        int l2_int = static_cast<int>(sqrt(lm2_int));
        int m2_int = lm2_int - (l2_int * l2_int + l2_int);
        double l2 = static_cast<double>(l2_int);
        double m2 = static_cast<double>(m2_int);
        double j = static_cast<double>(i_int);

        int idx_ej_l1m1_l2m2 = indexXX(i, l1, m1, j, l2, m2, N);

        result[idx_ej_l1m1_l2m2] = XP[idx_ej_l1m1_l2m2];

        result[CovMax + idx_ej_l1m1_l2m2] =
            -((std::pow(mu, 2) + l1 * (l1 + 1) / std::pow(R, 2)) *
              XP[idx_ej_l1m1_l2m2]);

        result[2 * CovMax + idx_ej_l1m1_l2m2] =
            -((std::pow(mu, 2) + l1 * (l1 + 1) / std::pow(R, 2)) *
              XX[idx_ej_l1m1_l2m2]) +
            PP[idx_ej_l1m1_l2m2];

        for (int h_index = 0; h_index < H.size(); h_index++) {
          for (double k = 1; k <= 2; k++) {

            double l3 = H[h_index][0];
            double m3 = H[h_index][1];
            double l4 = H[h_index][2];
            double m4 = H[h_index][3];
            double l5 = H[h_index][4];
            double m5 = H[h_index][5];
            double l1Orl2 = H[h_index][6];
            double m1Orm2 = H[h_index][7];

            int idx_di_l3m3_l4m4 = indexXX(k, l3, m3, i, l4, m4, N);
            int idx_di_l5m5_l4m4 = indexXX(k, l5, m5, i, l4, m4, N);
            int idx_bk_l3m3_l5m5 = indexXX(k, l3, m3, k, l5, m5, N);

            int idx_ej_l5m5_l2m2 = indexXX(k, l5, m5, j, l2, m2, N);
            int idx_ej_l3m3_l2m2 = indexXX(k, l3, m3, j, l2, m2, N);
            int idx_ej_l4m4_l2m2 = indexXX(i, l4, m4, j, l2, m2, N);
            // {i, a, l1, m1} <-> {j, e, l2, m2}
            int idx_ai_l5m5_l1m1 = indexXX(k, l5, m5, i, l1, m1, N);
            int idx_dj_l3m3_l4m4 = indexXX(k, l3, m3, j, l4, m4, N);
            int idx_ai_l3m3_l1m1 = indexXX(k, l3, m3, i, l1, m1, N);
            int idx_dj_l5m5_l4m4 = indexXX(k, l5, m5, j, l4, m4, N);
            int idx_ai_l4m4_l1m1 = indexXX(j, l4, m4, i, l1, m1, N);

            long double ppDotSubSum1 = 0;
            long double ppDotSubSum2 = 0;
            long double xpDotSubSum = 0;

            if (l1Orl2 == l1 && m1Orm2 == m1) {

              ppDotSubSum1 =
                  H[h_index][8] * (XP[idx_ej_l5m5_l2m2] * XX[idx_di_l3m3_l4m4] +
                                   XP[idx_ej_l3m3_l2m2] * XX[idx_di_l5m5_l4m4] +
                                   XP[idx_ej_l4m4_l2m2] * XX[idx_bk_l3m3_l5m5]);

              xpDotSubSum =
                  H[h_index][8] * (XX[idx_ej_l5m5_l2m2] * XX[idx_di_l3m3_l4m4] +
                                   XX[idx_ej_l3m3_l2m2] * XX[idx_di_l5m5_l4m4] +
                                   XX[idx_ej_l4m4_l2m2] * XX[idx_bk_l3m3_l5m5]);

            } else if (l1Orl2 == l2 && m1Orm2 == m2) {
              // HFunction(l3, l4, l5, l2, m3, m4, m5, m2, N)
              ppDotSubSum2 =
                  H[h_index][8] * (XP[idx_ai_l5m5_l1m1] * XX[idx_dj_l3m3_l4m4] +
                                   XP[idx_ai_l3m3_l1m1] * XX[idx_dj_l5m5_l4m4] +
                                   XP[idx_ai_l4m4_l1m1] * XX[idx_bk_l3m3_l5m5]);
            }

            result[CovMax + idx_ej_l1m1_l2m2] -=
                (lambda * (ppDotSubSum1 + ppDotSubSum2) / N);

            result[2 * CovMax + idx_ej_l1m1_l2m2] -= (lambda * xpDotSubSum / N);
          }
        }
      }
    }
  }
  return result;
}

std::vector<long double> GenerateSigmaAnsatzI(long double E, long double lambda,
                                              long double R, long double mu,
                                              int N) {
  long double sigmaXX =
      -((N * N - 1 + 2 * mu * mu * R * R) +
        std::sqrt(std::pow((N * N - 1 + 2 * mu * mu * R * R), 2) +
                  12 * lambda * E * R * R * R * R)) /
      (6 * lambda * N * N * R * R * R);
  long double sigmaPP = E * R / (N * N) -
                        lambda * N * N * R * R * R * sigmaXX * sigmaXX -
                        (N * N - 1 + 2 * mu * mu * R * R) * sigmaXX / 2;
  std::vector<long double> result = {sigmaXX, sigmaPP};
  return result;
}

std::vector<long double>
GenerateSigmaAnsatzII(long double E, long double lambda, long double R, int N) {
  long double sigmaXX =
      -((4 * N * N - 1) + std::sqrt(std::pow((4 * N * N - 1), 2) +
                                    864 * lambda * E * R * R * R * R)) /
      (144 * lambda * N * R * R * R);
  long double sigmaPP = E * R / (2 * N) - (4 * N * N - 1) * sigmaXX / 12 -
                        4 * N * R * R * R * lambda * sigmaXX * sigmaXX;
  std::vector<long double> result = {sigmaXX, sigmaPP};
  return result;
}

std::vector<long double> GenerateInitialConditionAnsatzI(long double E, int N,
                                                         long double R,
                                                         long double mu,
                                                         long double lambda) {
  int Imax = 2 * (N * N);   // Number of averaged variables
  int CovMax = Imax * Imax; // Number of covariance terms
  std::vector<long double> result(3 * CovMax);
  std::vector<long double> sigmas = GenerateSigmaAnsatzI(E, lambda, R, mu, N);
  long double sigmaXX = sigmas[0];
  long double sigmaPP = sigmas[1];
  for (double l1 = 0; l1 < N; l1++) {
    for (double m1 = -l1; m1 <= l1; m1++) {
      for (double i = 1; i <= 2; i++) {
        for (double l2 = 0; l2 < N; l2++) {
          for (double m2 = -l2; m2 <= l2; m2++) {
            for (double j = 1; j <= 2; j++) {
              int idx_ij_l1m1_l2m2 = indexXX(i, l1, m1, j, l2, m2, N);
              result[idx_ij_l1m1_l2m2] = R * sigmaXX * Delta(l1, l2) *
                                         Delta(i, j) * Delta(m1, -m2) *
                                         std::pow(-1, m2);
              result[CovMax + idx_ij_l1m1_l2m2] = sigmaPP * Delta(l1, l2) *
                                                  Delta(i, j) * Delta(m1, -m2) *
                                                  std::pow(-1, m2) / R;
              result[2 * CovMax + idx_ij_l1m1_l2m2] = 0;
            }
          }
        }
      }
    }
  }
  return result;
}

std::vector<long double> GenerateInitialConditionAnsatzII(long double E, int N,
                                                          long double R,
                                                          long double lambda) {
  int Imax = 2 * (N * N);   // Number of averaged variables
  int CovMax = Imax * Imax; // Number of covariance terms
  std::vector<long double> result(3 * CovMax);
  std::vector<long double> sigmas = GenerateSigmaAnsatzII(E, lambda, R, N);
  long double sigmaXX = sigmas[0];
  long double sigmaPP = sigmas[1];
  for (double l1 = 0; l1 < N; l1++) {
    for (double m1 = -l1; m1 <= l1; m1++) {
      for (double i = 1; i <= 2; i++) {
        for (double l2 = 0; l2 < N; l2++) {
          for (double m2 = -l2; m2 <= l2; m2++) {
            for (double j = 1; j <= 2; j++) {
              long double sigmaXXFactor = R * sigmaXX / (l1 + 0.5);
              long double sigmaPPFactor = sigmaPP * (l1 + 0.5) / R;
              int idx_ij_l1m1_l2m2 = indexXX(i, l1, m1, j, l2, m2, N);
              result[idx_ij_l1m1_l2m2] = sigmaXXFactor * Delta(l1, l2) *
                                         Delta(i, j) * Delta(m1, -m2) *
                                         std::pow(-1, m2);
              result[CovMax + idx_ij_l1m1_l2m2] =
                  sigmaPPFactor * Delta(l1, l2) * Delta(i, j) * Delta(m1, -m2) *
                  std::pow(-1, m2);
              result[2 * CovMax + idx_ij_l1m1_l2m2] = 0;
            }
          }
        }
      }
    }
  }
  return result;
}

std::vector<long double> GenerateInitialConditionReduced(int ansatz, int N,
                                                         long double E,
                                                         long double lambda,
                                                         long double R,
                                                         long double mu) {
  std::vector<long double> result;
  if (ansatz == 1) {
    result = GenerateInitialConditionAnsatzI(E, N, R, mu, lambda);
  } else if (ansatz == 2) {
    result = GenerateInitialConditionAnsatzII(E, N, R, lambda);
  } else {
    result = {0};
  }
  return result;
}

long double Sign(double m) {
  if (m < 0) {
    return -1.0;
  } else if (m == 0) {
    return 0.0;
  } else {
    return 1.0;
  }
}

long double R(double l, double m) {
  static constexpr long double tolerance = 1e-4L;

  // 1) Turn (l,m) into a deterministic 64-bit seed
  auto h1 = std::hash<double>{}(l);
  auto h2 = std::hash<double>{}(m);
  uint64_t seed = h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  // (the magic constant + bit-mix helps spread the bits)

  // 2) Make a generator seeded by that value
  std::mt19937_64 gen(seed);

  // 3) Draw uniformly in [0, tolerance)
  std::uniform_real_distribution<long double> dist(0.0L, tolerance);
  return dist(gen);
}

std::vector<long double>
ModifyInitialCondition(const std::vector<long double> &initialCondition,
                       int N) {
  int Imax = 2 * (N * N);
  int CovMax = Imax * Imax;
  std::vector<long double> result(3 * CovMax);
  for (double l1 = 0; l1 < N; l1++) {
    for (double m1 = -l1; m1 <= l1; m1++) {
      for (double i; i <= 2; i++) {
        int idx_ii_l1m1_l1m1 = indexXX(i, l1, m1, i, l1, m1, N);
        result[idx_ii_l1m1_l1m1] =
            initialCondition[idx_ii_l1m1_l1m1] + Sign(m1) * R(l1, std::abs(m1));
        result[CovMax + idx_ii_l1m1_l1m1] =
            initialCondition[CovMax + idx_ii_l1m1_l1m1] +
            Sign(m1) * R(l1, std::abs(m1));
        result[2 * CovMax + idx_ii_l1m1_l1m1] =
            initialCondition[2 * CovMax + idx_ii_l1m1_l1m1] +
            Sign(m1) * R(l1, std::abs(m1));
      }
    }
  }
  return result;
}
