//
// Created by berko on 3/20/2023.
//
#include "../include/ClebschGordanCoefficients.h"
#include "../include/MatrixOperations.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

long double Delta(double i, double j) {
  long double result;
  if (i == j) {
    result = 1.0;
  } else {
    result = 0.0;
  }
  return result;
}

long double Epsilon(double i, double j) {
  long double result;
  if (i == 1 && j == 2) {
    result = 1.0;
  } else if (i == 2 && j == 1) {
    result = -1.0;
  } else {
    result = 0.0;
  }
  return result;
}

long double SigmaX(double i, double j) {
  long double result;
  if (i == 1 && j == 2) {
    result = 1.0;
  } else if (i == 2 && j == 1) {
    result = 1.0;
  } else {
    result = 0.0;
  }
  return result;
}

long double SigmaZ(double i, double j) {
  long double result;
  if (i == 1 && j == 1) {
    result = 1.0;
  } else if (i == 2 && j == 2) {
    result = -1.0;
  } else {
    result = 0.0;
  }
  return result;
}

bool InVector(const std::vector<double> &v, double element) {
  bool isIn = false;
  for (size_t i = 0; i < v.size(); i++) {
    if (v[i] == element) {
      isIn = true;
    }
  }
  return isIn;
}

bool IsIntegerOrHalfInteger(double x) {
  return std::floor(x) == x || std::floor(x) + 0.5 == x;
}

bool IsInteger(double x) { return std::floor(x) == x; }

bool IsHalfInteger(double x) { return std::floor(x) + 0.5 == x; }

bool IsSumZero(double x) {
  if (x == 0.0) {
    return true;
  } else {
    return false;
  }
}

bool IsEvenInteger(double number) {
  // Check if the number is an integer
  if (std::floor(number) != number) {
    return false;
  }

  // Check if the integer is even
  int intPart = static_cast<int>(number);
  return intPart % 2 == 0;
}

long double FactorialDelta(double j1, double j2, double j3) {
  long double result = 0.0;
  result =
      std::sqrt((std::tgamma(1 + j1 + j2 - j3) * std::tgamma(1 + j1 - j2 + j3) *
                 std::tgamma(1 + -j1 + j2 + j3)) /
                (std::tgamma(1 + j1 + j2 + j3 + 1)));
  return result;
}

bool Wigner3jSelectionRule(double j1, double j2, double j3, double m1,
                           double m2, double m3) {
  std::vector<double> m1Values;
  for (double i = -j1; i <= j1; i++) {
    m1Values.push_back(i);
  }
  std::vector<double> m2Values;
  for (double i = -j2; i <= j2; i++) {
    m2Values.push_back(i);
  }
  std::vector<double> m3Values;
  for (double i = -j3; i <= j3; i++) {
    m3Values.push_back(i);
  }
  bool isSumZero = IsSumZero(m1 + m2 + m3);

  if (!(InVector(m1Values, m1) && InVector(m2Values, m2) &&
        InVector(m3Values, m3))) {
    return false;
  } else if (m1 + m2 + m3 != 0) {
    return false;
  } else if ((std::abs(j1 - j2) > j3) || (j3 > j1 + j2)) {
    return false;
  } else if (!(IsEvenInteger(j1 + j2 + j3))) {
    return false;
  } else if (isSumZero) {
    if (!(IsEvenInteger(j1 + j2 + j3))) {
      return false;
    }
  }
  return true;
}

long double Wigner3j(double j1, double j2, double j3, double m1, double m2,
                     double m3) {
  long double result;
  double maxSum = std::min(std::min(j1 + j2 - j3, j1 - m1), j2 + m2);
  double minSum = std::max(std::max(0.0, j2 - j3 - m1), j1 - j3 + m2);
  double sumResult = 0;
  if (Wigner3jSelectionRule(j1, j2, j3, m1, m2, m3)) {
    for (double k = minSum; k <= maxSum; k++) {
      sumResult +=
          std::pow(-1, k) /
          (std::tgamma(k + 1) * std::tgamma(j1 + j2 - j3 - k + 1) *
           std::tgamma(j1 - m1 - k + 1) * std::tgamma(j2 + m2 - k + 1) *
           std::tgamma(j3 - j2 + m1 + k + 1) *
           std::tgamma(j3 - j1 - m2 + k + 1));
    }
    result = Delta(m1 + m2 + m3, 0) * std::pow(-1, j1 - j2 - m3);
    result *= std::sqrt(
        (std::tgamma(j1 + j2 - j3 + 1) * std::tgamma(j1 - j2 + j3 + 1) *
         std::tgamma(-j1 + j2 + j3 + 1)) /
        std::tgamma(j1 + j2 + j3 + 2));
    result *= std::sqrt(std::tgamma(j1 - m1 + 1) * std::tgamma(j1 + m1 + 1) *
                        std::tgamma(j2 - m2 + 1) * std::tgamma(j2 + m2 + 1) *
                        std::tgamma(j3 - m3 + 1) * std::tgamma(j3 + m3 + 1));
    result *= sumResult;
  } else {
    result = 0.0;
  }
  return result;
}

bool TriangularInequality(double j1, double j2, double j3) {
  return (j3 <= j1 + j2 && j3 >= std::abs(j1 - j2) && j2 <= j1 + j3 &&
          j2 >= std::abs(j1 - j3) && j1 <= j3 + j2 && j1 >= std::abs(j3 - j2));
}

bool Wigner6jSelectionRule(double j1, double j2, double j3, double j4,
                           double j5, double j6) {
  return (
      TriangularInequality(j1, j2, j3) && TriangularInequality(j3, j4, j5) &&
      TriangularInequality(j1, j5, j6) && TriangularInequality(j2, j4, j6) &&
      IsInteger(j1 + j2 + j3) && IsInteger(j3 + j4 + j5) &&
      IsInteger(j1 + j5 + j6) && IsInteger(j2 + j4 + j6));
}

long double Wigner6j(double j1, double j2, double j3, double j4, double j5,
                     double j6) {
  long double result;
  if (Wigner6jSelectionRule(j1, j2, j3, j4, j5, j6)) {
    long double term1 = FactorialDelta(j1, j2, j3) *
                        FactorialDelta(j3, j4, j5) *
                        FactorialDelta(j1, j5, j6) * FactorialDelta(j2, j4, j6);
    long double term2 = 0;
    double startValue =
        std::max(j1 + j2 + j3,
                 std::max(j3 + j4 + j5, std::max(j1 + j5 + j6, j2 + j4 + j6)));
    double endValue = std::min(j1 + j2 + j4 + j5,
                               std::min(j1 + j3 + j4 + j6, j2 + j3 + j5 + j6));
    for (double n = startValue; n <= endValue; n++) {
      long double subTerm = (std::pow(-1, n) * std::tgamma(2 + n)) /
                            (std::tgamma(1 + j1 + j3 + j4 + j6 - n) *
                             std::tgamma(1 + j2 + j3 + j5 + j6 - n) *
                             std::tgamma(1 + n - j1 - j2 - j3) *
                             std::tgamma(1 + n - j3 - j4 - j5) *
                             std::tgamma(1 + n - j1 - j5 - j6) *
                             std::tgamma(1 + n - j2 - j4 - j6) *
                             std::tgamma(1 + j1 + j2 + j4 + j5 - n));
      term2 += subTerm;
    }
    result = term1 * term2;
  } else {
    result = 0;
  }
  return result;
}

long double ClesbchGordan(double j1, double j2, double j3, double m1, double m2,
                          double m3) {
  long double result;
  long double result1;
  long double result2;
  long double result3;
  result1 = std::pow(-1.0, -j1 + j2 - m3);
  result2 = std::sqrt(2.0 * j3 + 1.0);
  result3 = Wigner3j(j1, j2, j3, m1, m2, -m3);
  result = result1 * result2 * result3;
  return result;
}

// long double G(double l1, double l2, double l3, double l, double m1, double
// m2, double m3, double m, double N){
//   long double result = std::pow(-1, l);
//   long double summationTerm = 0;
//   for (double l4 = 0; l4 <= N - 1; l4++){
//     for (double m4 = -l4; m4 <= l4; m4++){
//       summationTerm += std::pow(-1, l4) * std::sqrt((2 * l1 + 1) * (2 * l2 +
//       1) * (2 * l3 + 1) * (2 * l4 + 1)) * Wigner6j(l1, l2, l4, (N - 1)/2, (N
//       - 1)/2, (N - 1)/2) * Wigner6j(l4, l3, l, (N - 1)/2, (N - 1)/2, (N -
//       1)/2) * ClesbchGordan(l4, l1, l2, m4, m1, m2) * ClesbchGordan(l, l4,
//       l3, m, m4, m3);
//     }
//   }
//   return result;
// }

// long double K(double l1, double l2, double l3, double l4, double m1, double
// m2, double m3, double m4, double N){
//   long double result;
//   if (m1 + m2 + m3 + m4 != 0){
//     result = 0;
//   }
//   else {
//     result = std::sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) * (2 * l4 +
//     1)); long double summationTerm = 0; for (double l = 0; l <= N - 1; l++){
//       for (double m = -l; m <= l; m++){
//         summationTerm += ClesbchGordan((N - 1)/2, l1, (N - 1)/2, m + m1, m1,
//         m) * ClesbchGordan((N - 1)/2, l2, (N - 1)/2, m + m1 + m2, m2, m + m1)
//         * ClesbchGordan((N - 1)/2, l3, (N - 1)/2, m + m1 + m2 + m3, m3, m +
//         m1 + m2) * ClesbchGordan((N - 1)/2, l4, (N - 1)/2, m, m4, m + m1 + m2
//         + m3);
//       }
//     }
//   }
//   return result;
// }

long double HFunction(double l1, double l2, double l3, double l, double m1,
                      double m2, double m3, double m, int N) {
  long double result = 0;
  for (double l4 = 0; l4 <= N; l4++) {
    for (double m4 = -l4; m4 <= l4; m4++) {
      result +=
          std::pow(-1, l1 + l2 + l4) * std::pow(-1, l3 + l4 + l) *
          std::sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) * (2 * l4 + 1)) *
          Wigner6j(l1, l2, l4, (N - 1) / 2.0, (N - 1) / 2.0, (N - 1) / 2.0) *
          Wigner6j(l3, l4, l, (N - 1) / 2.0, (N - 1) / 2.0, (N - 1) / 2.0) *
          ClesbchGordan(l4, l1, l2, m4, m1, m2) *
          ClesbchGordan(l, l3, l4, m, m3, m4);
    }
  }
  return result;
}

long double G(double a, double b, double c, double d) {
  long double result =
      -Epsilon(a, b) * SigmaX(c, d) + Delta(a, b) * SigmaZ(c, d);
  return result;
}

long double HFunctionInt(int l1_int, int l2_int, int l3_int, int l_int,
                         int m1_int, int m2_int, int m3_int, int m_int, int N) {
  double l1 = static_cast<double>(l1_int);
  double l2 = static_cast<double>(l2_int);
  double l3 = static_cast<double>(l3_int);
  double l = static_cast<double>(l_int);
  double m1 = static_cast<double>(m1_int);
  double m2 = static_cast<double>(m2_int);
  double m3 = static_cast<double>(m3_int);
  double m = static_cast<double>(m_int);

  long double result = 0;
  for (double l4 = 0; l4 <= N; l4++) {
    for (double m4 = -l4; m4 <= l4; m4++) {
      result +=
          std::pow(-1, l1 + l2 + l4) * std::pow(-1, l3 + l4 + l) *
          std::sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) * (2 * l4 + 1)) *
          Wigner6j(l1, l2, l4, (N - 1) / 2.0, (N - 1) / 2.0, (N - 1) / 2.0) *
          Wigner6j(l3, l4, l, (N - 1) / 2.0, (N - 1) / 2.0, (N - 1) / 2.0) *
          ClesbchGordan(l4, l1, l2, m4, m1, m2) *
          ClesbchGordan(l, l3, l4, m, m3, m4);
    }
  }
  return result;
}

long double GInt(int a_int, int b_int, int c_int, int d_int) {
  double a = static_cast<double>(a_int);
  double b = static_cast<double>(b_int);
  double c = static_cast<double>(c_int);
  double d = static_cast<double>(d_int);
  long double result =
      -Epsilon(a, b) * SigmaX(c, d) + Delta(a, b) * SigmaZ(c, d);
  return result;
}
