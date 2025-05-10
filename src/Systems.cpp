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
int indexX(double a, double i, double l, double m, int N) {
  int int_a = static_cast<int>(a);
  int int_i = static_cast<int>(i);
  int int_l = static_cast<int>(l);
  int int_m = static_cast<int>(m);
  return ((2 * int_a + int_i - 3) * (N * N) + (int_l * int_l + int_l + int_m));
}

int indexXX(double a1, double i1, double l1, double m1, double a2, double i2,
            double l2, double m2, int N) {
  int index1 = indexX(a1, i1, l1, m1, N);
  int index2 = indexX(a2, i2, l2, m2, N);
  int maxIndex = 4 * (N * N);        // Adjust as per your actual maximum index
  return index1 * maxIndex + index2; // Adjust based on how you store XX
}

int indexXInt(int a, int i, int l, int m, int N) {
  return ((2 * a + i - 3) * (N * N) + (l * l + l + m));
}

int indexXXInt(int a1, int i1, int l1, int m1, int a2, int i2, int l2, int m2,
               int N) {
  int index1 = indexX(a1, i1, l1, m1, N);
  int index2 = indexX(a2, i2, l2, m2, N);
  int maxIndex = 4 * (N * N);        // Adjust as per your actual maximum index
  return index1 * maxIndex + index2; // Adjust based on how you store XX
}

std::vector<long double>
averagedEquationsPolarizationBasisSymmetryReducedParallel(
    const std::vector<long double> &allVectors) {
  long double mu = 0.25;    // Previous value was 1.0
  long double lambda = 1.0; // Previous value was 10.0
  long double R = 2.0;

  int N = 2;
  int Imax = 4 * (N * N);   // Number of averaged variables
  int CovMax = Imax * Imax; // Number of covariance terms
  // int currentIndex = 0;

  // Extract covariance vectors (centered)
  std::vector<long double> XX(allVectors.begin(), allVectors.begin() + CovMax);
  std::vector<long double> PP(allVectors.begin() + CovMax,
                              allVectors.begin() + 2 * CovMax);
  std::vector<long double> XP(allVectors.begin() + 2 * CovMax,
                              allVectors.begin() + 3 * CovMax);

  std::vector<long double> result(3 * CovMax);

#pragma omp parallel for collapse(9)
  for (int lm1_int = 0; lm1_int < N * N; lm1_int++) {
    for (int i_int = 1; i_int <= 2; i_int++) {
      for (int a_int = 1; a_int <= 2; a_int++) {
        for (int lm2_int = 0; lm2_int < N * N; lm2_int++) {
          for (int j_int = 1; j_int <= 2; j_int++) {
            for (int e_int = 1; e_int <= 2; e_int++) {
              for (int lm3_int = 0; lm3_int < N * N; lm3_int++) {
                for (int lm4_int = 0; lm4_int < N * N; lm4_int++) {
                  for (int lm5_int = 0; lm5_int < N * N; lm5_int++) {
                    for (double k = 1; k <= 2; k++) {
                      for (double b = 1; b <= 2; b++) {
                        for (double c = 1; c <= 2; c++) {
                          for (double d = 1; d <= 2; d++) {
                            int l1_int = static_cast<int>(sqrt(lm1_int));
                            int m1_int = lm1_int - (l1_int * l1_int + l1_int);
                            double l1 = static_cast<double>(l1_int);
                            double m1 = static_cast<double>(m1_int);
                            double i = static_cast<double>(i_int);
                            double a = static_cast<double>(a_int);

                            int l2_int = static_cast<int>(sqrt(lm2_int));
                            int m2_int = lm2_int - (l2_int * l2_int + l2_int);
                            double l2 = static_cast<double>(l2_int);
                            double m2 = static_cast<double>(m2_int);
                            double j = static_cast<double>(j_int);
                            double e = static_cast<double>(e_int);

                            int l3_int = static_cast<int>(sqrt(lm3_int));
                            int m3_int = lm3_int - (l3_int * l3_int + l3_int);
                            double l3 = static_cast<double>(l3_int);
                            double m3 = static_cast<double>(m3_int);

                            int l4_int = static_cast<int>(sqrt(lm4_int));
                            int m4_int = lm4_int - (l4_int * l4_int + l4_int);
                            double l4 = static_cast<double>(l4_int);
                            double m4 = static_cast<double>(m4_int);

                            int l5_int = static_cast<int>(sqrt(lm5_int));
                            int m5_int = lm5_int - (l5_int * l5_int + l5_int);
                            double l5 = static_cast<double>(l5_int);
                            double m5 = static_cast<double>(m5_int);

                            int idx_ai_l1m1 = indexX(a, i, l1, m1, N);
                            int idx_bk_l5m5 = indexX(b, k, l5, m5, N);
                            int idx_ck_di_l3m3_l4m4 =
                                indexXX(c, k, l3, m3, d, i, l4, m4, N);
                            int idx_ck_l3m3 = indexX(c, k, l3, m3, N);
                            int idx_bk_di_l5m5_l4m4 =
                                indexXX(b, k, l5, m5, d, i, l4, m4, N);
                            int idx_di_l4m4 = indexX(d, i, l4, m4, N);
                            int idx_ck_bk_l3m3_l5m5 =
                                indexXX(c, k, l3, m3, b, k, l5, m5, N);
                            // Compute the dot product of the XX, PP and
                            // XP terms.
                            int idx_ai_ej_l1m1_l2m2 =
                                indexXX(a, i, l1, m1, e, j, l2, m2, N);
                            // First term of the XXdot term
                            result[idx_ai_ej_l1m1_l2m2] =
                                XP[idx_ai_ej_l1m1_l2m2];
                            int idx_bk_ej_l5m5_l2m2 =
                                indexXX(b, k, l5, m5, e, j, l2, m2, N);
                            int idx_ck_ej_l3m3_l2m2 =
                                indexXX(c, k, l3, m3, e, j, l2, m2, N);
                            int idx_di_ej_l4m4_l2m2 =
                                indexXX(d, i, l4, m4, e, j, l2, m2, N);
                            // {i, a, l1, m1} <-> {j, e, l2, m2}
                            int idx_ej_ai_l2m2_l1m1 =
                                indexXX(e, j, l2, m2, a, i, l1, m1, N);
                            int idx_bk_ai_l5m5_l1m1 =
                                indexXX(b, k, l5, m5, a, i, l1, m1, N);
                            int idx_ck_dj_l3m3_l4m4 =
                                indexXX(c, k, l3, m3, d, j, l4, m4, N);
                            int idx_ck_ai_l3m3_l1m1 =
                                indexXX(c, k, l3, m3, a, i, l1, m1, N);
                            int idx_bk_dj_l5m5_l4m4 =
                                indexXX(b, k, l5, m5, d, j, l4, m4, N);
                            int idx_dj_ai_l4m4_l1m1 =
                                indexXX(d, j, l4, m4, a, i, l1, m1, N);
                            int idx_dj_l4m4 = indexX(d, j, l4, m4, N);
                            // Compute the dot product of the PP term.
                            // std::cout << "current index: " <<
                            // currentIndex << " / " << "262,144" <<
                            // "\n"; currentIndex += 1;
                            long double ppDotSubSum1 =
                                HFunction(l3, l4, l5, l1, m3, m4, m5, m1, N) *
                                G(a, b, c, d) *
                                (XP[idx_bk_ej_l5m5_l2m2] *
                                     XX[idx_ck_di_l3m3_l4m4] +
                                 XP[idx_ck_ej_l3m3_l2m2] *
                                     XX[idx_bk_di_l5m5_l4m4] +
                                 XP[idx_di_ej_l4m4_l2m2] *
                                     XX[idx_ck_bk_l3m3_l5m5]);
                            long double ppDotSubSum2 =
                                HFunction(l3, l4, l5, l2, m3, m4, m5, m2, N) *
                                G(e, b, c, d) *
                                (XP[idx_bk_ai_l5m5_l1m1] *
                                     XX[idx_ck_dj_l3m3_l4m4] +
                                 XP[idx_ck_ai_l3m3_l1m1] *
                                     XX[idx_bk_dj_l5m5_l4m4] +
                                 XP[idx_dj_ai_l4m4_l1m1] *
                                     XX[idx_ck_bk_l3m3_l5m5]);

                            result[CovMax + idx_ai_ej_l1m1_l2m2] =
                                -((mu + l1 * (l1 + 1) / std::pow(R, 2)) *
                                      XP[idx_ai_ej_l1m1_l2m2] -
                                  lambda * (ppDotSubSum1 + ppDotSubSum2) / N);

                            // Compute the dot product of the XP term.
                            long double xpDotSubSum =
                                HFunction(l3, l4, l5, l1, m3, m4, m5, m1, N) *
                                G(a, b, c, d) *
                                (XX[idx_bk_ej_l5m5_l2m2] *
                                     XX[idx_ck_di_l3m3_l4m4] +
                                 XX[idx_ck_ej_l3m3_l2m2] *
                                     XX[idx_bk_di_l5m5_l4m4] +
                                 XX[idx_di_ej_l4m4_l2m2] *
                                     XX[idx_ck_bk_l3m3_l5m5]);

                            result[2 * CovMax + idx_ai_ej_l1m1_l2m2] =
                                -((mu + l1 * (l1 + 1) / std::pow(R, 2)) *
                                      XX[idx_ai_ej_l1m1_l2m2] -
                                  lambda * xpDotSubSum / N) +
                                PP[idx_ai_ej_l1m1_l2m2];
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return result;
}

std::vector<long double>
averagedEquationsPolarizationBasisSymmetryReducedParallelTrial(
    const std::vector<long double> &allVectors,
    const std::vector<std::vector<double>> &H) {
  long double mu = 0.25;    // Previous value was 1.0
  long double lambda = 1.0; // Previous value was 10.0
  long double R = 2.0;

  int N = 2;
  int Imax = 4 * (N * N);   // Number of averaged variables
  int CovMax = Imax * Imax; // Number of covariance terms

  // Extract covariance vectors (centered)
  std::vector<long double> XX(allVectors.begin(), allVectors.begin() + CovMax);
  std::vector<long double> PP(allVectors.begin() + CovMax,
                              allVectors.begin() + 2 * CovMax);
  std::vector<long double> XP(allVectors.begin() + 2 * CovMax,
                              allVectors.begin() + 3 * CovMax);

  std::vector<long double> result(3 * CovMax);

  // #pragma omp parallel for collapse(4)
  for (int lm1_int = 0; lm1_int < N * N; lm1_int++) {
    for (int i_int = 1; i_int <= 2; i_int++) {
      for (int a_int = 1; a_int <= 2; a_int++) {
        for (int lm2_int = 0; lm2_int < N * N; lm2_int++) {
          int l1_int = static_cast<int>(sqrt(lm1_int));
          int m1_int = lm1_int - (l1_int * l1_int + l1_int);
          double l1 = static_cast<double>(l1_int);
          double m1 = static_cast<double>(m1_int);
          double i = static_cast<double>(i_int);
          double a = static_cast<double>(a_int);

          int l2_int = static_cast<int>(sqrt(lm2_int));
          int m2_int = lm2_int - (l2_int * l2_int + l2_int);
          double l2 = static_cast<double>(l2_int);
          double m2 = static_cast<double>(m2_int);
          double j = static_cast<double>(i_int);
          double e = static_cast<double>(a_int);

          int idx_ai_ej_l1m1_l2m2 = indexXX(a, i, l1, m1, e, j, l2, m2, N);

          // First term of the XXdot term
          result[idx_ai_ej_l1m1_l2m2] = XP[idx_ai_ej_l1m1_l2m2];

          result[CovMax + idx_ai_ej_l1m1_l2m2] =
              -((std::pow(mu, 2) + l1 * (l1 + 1) / std::pow(R, 2)) *
                XP[idx_ai_ej_l1m1_l2m2]);

          result[2 * CovMax + idx_ai_ej_l1m1_l2m2] =
              -((std::pow(mu, 2) + l1 * (l1 + 1) / std::pow(R, 2)) *
                XX[idx_ai_ej_l1m1_l2m2]) +
              PP[idx_ai_ej_l1m1_l2m2];

          for (int h_index = 0; h_index < H.size(); h_index++) {
            for (double k = 1; k <= 2; k++) {
              for (double b = 1; b <= 2; b++) {
                for (double c = 1; c <= 2; c++) {
                  for (double d = 1; d <= 2; d++) {

                    double l3 = H[h_index][0];
                    double m3 = H[h_index][1];
                    double l4 = H[h_index][2];
                    double m4 = H[h_index][3];
                    double l5 = H[h_index][4];
                    double m5 = H[h_index][5];
                    double l1Orl2 = H[h_index][6];
                    double m1Orm2 = H[h_index][7];

                    int idx_ai_l1m1 = indexX(a, i, l1, m1, N);
                    int idx_bk_l5m5 = indexX(b, k, l5, m5, N);
                    int idx_ck_di_l3m3_l4m4 =
                        indexXX(c, k, l3, m3, d, i, l4, m4, N);
                    int idx_ck_l3m3 = indexX(c, k, l3, m3, N);
                    int idx_bk_di_l5m5_l4m4 =
                        indexXX(b, k, l5, m5, d, i, l4, m4, N);
                    int idx_di_l4m4 = indexX(d, i, l4, m4, N);
                    int idx_ck_bk_l3m3_l5m5 =
                        indexXX(c, k, l3, m3, b, k, l5, m5, N);
                    // Compute the dot product of the XX, PP and
                    // XP terms.

                    // First term of the XXdot term
                    int idx_bk_ej_l5m5_l2m2 =
                        indexXX(b, k, l5, m5, e, j, l2, m2, N);
                    int idx_ck_ej_l3m3_l2m2 =
                        indexXX(c, k, l3, m3, e, j, l2, m2, N);
                    int idx_di_ej_l4m4_l2m2 =
                        indexXX(d, i, l4, m4, e, j, l2, m2, N);
                    // {i, a, l1, m1} <-> {j, e, l2, m2}
                    int idx_ej_ai_l2m2_l1m1 =
                        indexXX(e, j, l2, m2, a, i, l1, m1, N);
                    int idx_bk_ai_l5m5_l1m1 =
                        indexXX(b, k, l5, m5, a, i, l1, m1, N);
                    int idx_ck_dj_l3m3_l4m4 =
                        indexXX(c, k, l3, m3, d, j, l4, m4, N);
                    int idx_ck_ai_l3m3_l1m1 =
                        indexXX(c, k, l3, m3, a, i, l1, m1, N);
                    int idx_bk_dj_l5m5_l4m4 =
                        indexXX(b, k, l5, m5, d, j, l4, m4, N);
                    int idx_dj_ai_l4m4_l1m1 =
                        indexXX(d, j, l4, m4, a, i, l1, m1, N);
                    int idx_dj_l4m4 = indexX(d, j, l4, m4, N);

                    long double ppDotSubSum1 = 0;
                    long double ppDotSubSum2 = 0;
                    long double xpDotSubSum = 0;

                    if (l1Orl2 == l1 && m1Orm2 == m1) {

                      // HFunction(l3, l4, l5, l1, m3, m4, m5, m1, N) *
                      ppDotSubSum1 =
                          H[h_index][8] * G(a, b, c, d) *
                          (XP[idx_bk_ej_l5m5_l2m2] * XX[idx_ck_di_l3m3_l4m4] +
                           XP[idx_ck_ej_l3m3_l2m2] * XX[idx_bk_di_l5m5_l4m4] +
                           XP[idx_di_ej_l4m4_l2m2] * XX[idx_ck_bk_l3m3_l5m5]);

                      // Compute the dot product of the XP term.
                      // HFunction(l3, l4, l5, l1, m3, m4, m5, m1, N)
                      xpDotSubSum =
                          H[h_index][8] * G(a, b, c, d) *
                          (XX[idx_bk_ej_l5m5_l2m2] * XX[idx_ck_di_l3m3_l4m4] +
                           XX[idx_ck_ej_l3m3_l2m2] * XX[idx_bk_di_l5m5_l4m4] +
                           XX[idx_di_ej_l4m4_l2m2] * XX[idx_ck_bk_l3m3_l5m5]);

                      // if (ppDotSubSum1 != 0 || xpDotSubSum != 0) {
                      //   std::cout << "ppDotSubSum1 = " << ppDotSubSum1
                      //             << ", xpDotSubSum = " << xpDotSubSum
                      //             << ", Hsymbol(" << l3 << ", " << l4 << ", "
                      //             << l5 << ", " << l1 << ", " << m3 << ", "
                      //             << m4 << ", " << m5 << ", " << m1
                      //             << ") = " << H[h_index][8] << ", "
                      //             << G(a, b, c, d) << "\n";
                      // }

                    } else if (l1Orl2 == l2 && m1Orm2 == m2) {
                      // HFunction(l3, l4, l5, l2, m3, m4, m5, m2, N)
                      ppDotSubSum2 =
                          H[h_index][8] * G(e, b, c, d) *
                          (XP[idx_bk_ai_l5m5_l1m1] * XX[idx_ck_dj_l3m3_l4m4] +
                           XP[idx_ck_ai_l3m3_l1m1] * XX[idx_bk_dj_l5m5_l4m4] +
                           XP[idx_dj_ai_l4m4_l1m1] * XX[idx_ck_bk_l3m3_l5m5]);
                      // if (ppDotSubSum2 != 0) {
                      //   std::cout << "ppDotSubSum2 = " << ppDotSubSum2
                      //             << ", Hsymbol(" << l3 << ", " << l4 << ", "
                      //             << l5 << ", " << l1 << ", " << m3 << ", "
                      //             << m4 << ", " << m5 << ", " << m1
                      //             << ") = " << H[h_index][8] << ", "
                      //             << G(a, b, c, d) << "\n";
                      // }
                    }

                    result[CovMax + idx_ai_ej_l1m1_l2m2] -=
                        (lambda * (ppDotSubSum1 + ppDotSubSum2) / N);

                    result[2 * CovMax + idx_ai_ej_l1m1_l2m2] -=
                        (lambda * xpDotSubSum / N);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return result;
}

std::vector<long double>
averagedEquationsPolarizationBasisSymmetryReducedParallelMpi(
    const std::vector<long double> &allVectors, int world_rank,
    int world_size) {

  long double mu = 1.0;
  long double lambda = 5.0;
  long double R = 2.0;
  int N = 2;
  int Imax = 4 * (N * N);
  int CovMax = Imax * Imax;

  std::vector<long double> XX(allVectors.begin(), allVectors.begin() + CovMax);
  std::vector<long double> PP(allVectors.begin() + CovMax,
                              allVectors.begin() + 2 * CovMax);
  std::vector<long double> XP(allVectors.begin() + 2 * CovMax,
                              allVectors.begin() + 3 * CovMax);

  std::vector<long double> local_result(3 * CovMax, 0.0);

  int total = 4 * N * N;
  int chunk = (total + world_size - 1) / world_size;
  int start = world_rank * chunk;
  int end = std::min(start + chunk, total);

#pragma omp parallel for collapse(7)
  for (int lm1ia_int = start; lm1ia_int < end; ++lm1ia_int) {
    for (int lm2_int = 0; lm2_int < N * N; lm2_int++) {
      for (int j_int = 1; j_int <= 2; j_int++) {
        for (int e_int = 1; e_int <= 2; e_int++) {
          for (int lm3_int = 0; lm3_int < N * N; lm3_int++) {
            for (int lm4_int = 0; lm4_int < N * N; lm4_int++) {
              for (int lm5_int = 0; lm5_int < N * N; lm5_int++) {
                for (double k = 1; k <= 2; k++) {
                  for (double b = 1; b <= 2; b++) {
                    for (double c = 1; c <= 2; c++) {
                      for (double d = 1; d <= 2; d++) {
                        // Unflatten the first loop
                        int a_int = lm1ia_int % 2 + 1;
                        int i_int = (lm1ia_int / 2) % 2 + 1;
                        int lm1_int = lm1ia_int / 4;

                        int l1_int = static_cast<int>(sqrt(lm1_int));
                        int m1_int = lm1_int - (l1_int * l1_int + l1_int);
                        double l1 = static_cast<double>(l1_int);
                        double m1 = static_cast<double>(m1_int);
                        double i = i_int, a = a_int;

                        int l2_int = static_cast<int>(sqrt(lm2_int));
                        int m2_int = lm2_int - (l2_int * l2_int + l2_int);
                        double l2 = static_cast<double>(l2_int);
                        double m2 = static_cast<double>(m2_int);
                        double j = j_int, e = e_int;

                        int l3_int = static_cast<int>(sqrt(lm3_int));
                        int m3_int = lm3_int - (l3_int * l3_int + l3_int);
                        double l3 = static_cast<double>(l3_int);
                        double m3 = static_cast<double>(m3_int);

                        int l4_int = static_cast<int>(sqrt(lm4_int));
                        int m4_int = lm4_int - (l4_int * l4_int + l4_int);
                        double l4 = static_cast<double>(l4_int);
                        double m4 = static_cast<double>(m4_int);

                        int l5_int = static_cast<int>(sqrt(lm5_int));
                        int m5_int = lm5_int - (l5_int * l5_int + l5_int);
                        double l5 = static_cast<double>(l5_int);
                        double m5 = static_cast<double>(m5_int);

                        int idx_ai_l1m1 = indexX(a, i, l1, m1, N);
                        int idx_ai_ej_l1m1_l2m2 =
                            indexXX(a, i, l1, m1, e, j, l2, m2, N);
                        int idx_bk_l5m5 = indexX(b, k, l5, m5, N);
                        int idx_ck_di_l3m3_l4m4 =
                            indexXX(c, k, l3, m3, d, i, l4, m4, N);
                        int idx_bk_di_l5m5_l4m4 =
                            indexXX(b, k, l5, m5, d, i, l4, m4, N);
                        int idx_ck_bk_l3m3_l5m5 =
                            indexXX(c, k, l3, m3, b, k, l5, m5, N);
                        int idx_bk_ej_l5m5_l2m2 =
                            indexXX(b, k, l5, m5, e, j, l2, m2, N);
                        int idx_ck_ej_l3m3_l2m2 =
                            indexXX(c, k, l3, m3, e, j, l2, m2, N);
                        int idx_di_ej_l4m4_l2m2 =
                            indexXX(d, i, l4, m4, e, j, l2, m2, N);
                        int idx_bk_ai_l5m5_l1m1 =
                            indexXX(b, k, l5, m5, a, i, l1, m1, N);
                        int idx_ck_dj_l3m3_l4m4 =
                            indexXX(c, k, l3, m3, d, j, l4, m4, N);
                        int idx_ck_ai_l3m3_l1m1 =
                            indexXX(c, k, l3, m3, a, i, l1, m1, N);
                        int idx_bk_dj_l5m5_l4m4 =
                            indexXX(b, k, l5, m5, d, j, l4, m4, N);
                        int idx_dj_ai_l4m4_l1m1 =
                            indexXX(d, j, l4, m4, a, i, l1, m1, N);
                        int idx_dj_l4m4 = indexX(d, j, l4, m4, N);

                        local_result[idx_ai_ej_l1m1_l2m2] =
                            XP[idx_ai_ej_l1m1_l2m2];

                        long double ppDotSubSum1 =
                            HFunction(l3, l4, l5, l1, m3, m4, m5, m1, N) *
                            G(a, b, c, d) *
                            (XP[idx_bk_ej_l5m5_l2m2] * XX[idx_ck_di_l3m3_l4m4] +
                             XP[idx_ck_ej_l3m3_l2m2] * XX[idx_bk_di_l5m5_l4m4] +
                             XP[idx_di_ej_l4m4_l2m2] * XX[idx_ck_bk_l3m3_l5m5]);

                        long double ppDotSubSum2 =
                            HFunction(l3, l4, l5, l2, m3, m4, m5, m2, N) *
                            G(e, b, c, d) *
                            (XP[idx_bk_ai_l5m5_l1m1] * XX[idx_ck_dj_l3m3_l4m4] +
                             XP[idx_ck_ai_l3m3_l1m1] * XX[idx_bk_dj_l5m5_l4m4] +
                             XP[idx_dj_ai_l4m4_l1m1] * XX[idx_ck_bk_l3m3_l5m5]);

                        local_result[CovMax + idx_ai_ej_l1m1_l2m2] =
                            -((mu + l1 * (l1 + 1) / (R * R)) *
                                  XP[idx_ai_ej_l1m1_l2m2] -
                              lambda * (ppDotSubSum1 + ppDotSubSum2) / N);

                        long double xpDotSubSum =
                            HFunction(l3, l4, l5, l1, m3, m4, m5, m1, N) *
                            G(a, b, c, d) *
                            (XX[idx_bk_ej_l5m5_l2m2] * XX[idx_ck_di_l3m3_l4m4] +
                             XX[idx_ck_ej_l3m3_l2m2] * XX[idx_bk_di_l5m5_l4m4] +
                             XX[idx_di_ej_l4m4_l2m2] * XX[idx_ck_bk_l3m3_l5m5]);

                        local_result[2 * CovMax + idx_ai_ej_l1m1_l2m2] =
                            -((mu + l1 * (l1 + 1) / (R * R)) *
                                  XX[idx_ai_ej_l1m1_l2m2] -
                              lambda * xpDotSubSum / N) +
                            PP[idx_ai_ej_l1m1_l2m2];
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  std::vector<long double> global_result(3 * CovMax, 0.0);
  MPI_Reduce(local_result.data(), global_result.data(), 3 * CovMax,
             MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return global_result;
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
  std::vector<long double> result(48 * N * N * N * N);
  int Imax = 4 * (N * N);   // Number of averaged variables
  int CovMax = Imax * Imax; // Number of covariance terms
  std::vector<long double> sigmas = GenerateSigmaAnsatzI(E, lambda, R, mu, N);
  long double sigmaXX = sigmas[0];
  long double sigmaPP = sigmas[1];
  for (double l1 = 0; l1 < N; l1++) {
    for (double m1 = -l1; m1 <= l1; m1++) {
      for (double a = 1; a <= 2; a++) {
        for (double i = 1; i <= 2; i++) {
          for (double l2 = 0; l2 < N; l2++) {
            for (double m2 = -l2; m2 <= l2; m2++) {
              for (double e = 1; e <= 2; e++) {
                for (double j = 1; j <= 2; j++) {
                  int idx_ai_l1m1 = indexX(a, i, l1, m1, N);
                  int idx_ai_ej_l1m1_l2m2 =
                      indexXX(a, i, l1, m1, e, j, l2, m2, N);
                  result[idx_ai_ej_l1m1_l2m2] =
                      R * sigmaXX * Delta(l1, l2) * Delta(i, j) * Delta(a, e) *
                      (Delta(m1, m2) +
                       std::pow(-1, m2 + a + 1) * Delta(m1, -m2));
                  result[CovMax + idx_ai_ej_l1m1_l2m2] =
                      sigmaPP * Delta(l1, l2) * Delta(i, j) * Delta(a, e) *
                      (Delta(m1, m2) +
                       std::pow(-1, m2 + a + 1) * Delta(m1, -m2)) /
                      R;
                  result[2 * CovMax + idx_ai_ej_l1m1_l2m2] = 0;
                }
              }
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
  std::vector<long double> result(48 * N * N * N * N);
  int Imax = 4 * (N * N);   // Number of averaged variables
  int CovMax = Imax * Imax; // Number of covariance terms
  std::vector<long double> sigmas = GenerateSigmaAnsatzII(E, lambda, R, N);
  long double sigmaXX = sigmas[0];
  long double sigmaPP = sigmas[1];
  for (double l1 = 0; l1 < N; l1++) {
    for (double m1 = -l1; m1 <= l1; m1++) {
      for (double a = 1; a <= 2; a++) {
        for (double i = 1; i <= 2; i++) {
          for (double l2 = 0; l2 < N; l2++) {
            for (double m2 = -l2; m2 <= l2; m2++) {
              for (double e = 1; e <= 2; e++) {
                for (double j = 1; j <= 2; j++) {
                  long double sigmaXXFactor = R * sigmaXX / (l1 + 0.5);
                  long double sigmaPPFactor = sigmaPP * (l1 + 0.5) / R;
                  int idx_ai_l1m1 = indexX(a, i, l1, m1, N);
                  int idx_ai_ej_l1m1_l2m2 =
                      indexXX(a, i, l1, m1, e, j, l2, m2, N);
                  result[idx_ai_ej_l1m1_l2m2] =
                      sigmaXXFactor * Delta(l1, l2) * Delta(i, j) *
                      Delta(a, e) *
                      (Delta(m1, m2) +
                       std::pow(-1, m2 + a + 1) * Delta(m1, -m2));
                  result[CovMax + idx_ai_ej_l1m1_l2m2] =
                      sigmaPPFactor * Delta(l1, l2) * Delta(i, j) *
                      Delta(a, e) *
                      (Delta(m1, m2) +
                       std::pow(-1, m2 + a + 1) * Delta(m1, -m2));
                  result[2 * CovMax + idx_ai_ej_l1m1_l2m2] = 0;
                }
              }
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
  int Imax = 4 * (N * N);
  int CovMax = Imax * Imax;
  std::vector<long double> result(3 * CovMax);
  for (double l1 = 0; l1 < N; l1++) {
    for (double m1 = -l1; m1 <= l1; m1++) {
      for (double i; i <= 2; i++) {
        for (double a = 1; a <= 2; a++) {
          int idx_ai_ai_l1m1_l1m1 = indexXX(a, i, l1, m1, a, i, l1, m1, N);
          result[idx_ai_ai_l1m1_l1m1] = initialCondition[idx_ai_ai_l1m1_l1m1] +
                                        Sign(m1) * R(l1, std::abs(m1));
          result[CovMax + idx_ai_ai_l1m1_l1m1] =
              initialCondition[CovMax + idx_ai_ai_l1m1_l1m1] +
              Sign(m1) * R(l1, std::abs(m1));
          result[2 * CovMax + idx_ai_ai_l1m1_l1m1] =
              initialCondition[2 * CovMax + idx_ai_ai_l1m1_l1m1] +
              Sign(m1) * R(l1, std::abs(m1));
        }
      }
    }
  }
  return result;
}
