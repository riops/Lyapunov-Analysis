#ifndef CLEBSCH_GORDAN_COEFFICIENTS_H
#define CLEBSCH_GORDAN_COEFFICIENTS_H

#include <vector>

long double Delta(double i, double j);
long double Epsilon(double i, double j);
long double SigmaX(double i, double j);
long double SigmaZ(double i, double j);
bool InVector(const std::vector<double>& v, double element);
bool IsIntegerOrHalfInteger(double x);
bool IsInteger(double x);
bool IsHalfInteger(double x);
bool IsSumZero(double x);
bool IsEvenInteger(double number);
long double FactorialDelta(double j1, double j2, double j3);
bool Wigner3jSelectionRule(double j1, double j2, double j3, double m1, double m2, double m3);
long double Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3);
bool TriangularInequality(double j1, double j2, double j3);
bool Wigner6jSelectionRule(double j1, double j2, double j3, double j4, double j5, double j6);
long double Wigner6j(double j1, double j2, double j3, double j4, double j5, double j6);
long double ClesbchGordan(double j1, double j2, double j3, double m1, double m2, double m3);
long double HFunction(double l1, double l2, double l3, double l, double m1, double m2, double m3, double m, int N);
long double G(double a, double b, double c, double d);
long double HFunctionInt(int l1, int l2, int l3, int l, int m1, int m2, int m3, int m, int N);
long double GInt(int a, int b, int c, int d);

#endif // CLEBSCH_GORDAN_COEFFICIENTS_H
