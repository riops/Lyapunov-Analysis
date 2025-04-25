#ifndef LYAPUNOV_ANALYSIS_H
#define LYAPUNOV_ANALYSIS_H

#include <vector>
#include <functional>

std::vector<long double> JacobianSystem(std::function<std::vector<long double>(const std::vector<long double>&)> f, const std::vector<long double>& point);

std::vector<long double> RungeKuttaJacobian4(std::function<std::vector<long double>(const std::vector<long double>&)> f, const std::vector<long double>& point, long double dt);

std::vector<long double> RungeKuttaJacobian5(std::function<std::vector<long double>(const std::vector<long double>&)> f, const std::vector<long double>& point, long double dt);

std::vector<long double> IntegrateJacobianSystem4(std::function<std::vector<long double>(const std::vector<long double>&)> f, const std::vector<long double>& point, long double dt, int numSteps);

std::vector<long double> IntegrateJacobianSystem5(std::function<std::vector<long double>(const std::vector<long double>&)> f, const std::vector<long double>& point, long double dt, int numSteps);

std::vector<long double> IntegrateJacobianSystem45(std::function<std::vector<long double>(const std::vector<long double>&)> f, const std::vector<long double>& point, long double dt, int numSteps);

std::vector<long double> GetPhi(const std::vector<long double>& point, int dimension);

std::vector<long double> LogNorms(const std::vector<long double>& vector1);

std::vector<long double> Normalize(const std::vector<long double>& vector1);

std::vector<long double> ReplacePhi(const std::vector<long double>& point, const std::vector<long double>& newPhi);

std::vector<std::vector<long double>> SumLogNorms(const std::vector<std::vector<long double>>& logNorms, long double dt, int numSteps);

std::vector<std::vector<long double>> LyapunovSpectrum4(std::function<std::vector<long double>(const std::vector<long double>&)> f, const std::vector<long double>& point, long double dt, int numSteps, int numIterations);

std::vector<std::vector<long double>> LyapunovSpectrum5(std::function<std::vector<long double>(const std::vector<long double>&)> f, const std::vector<long double>& point, long double dt, int numSteps, int numIterations);

std::vector<std::vector<long double>> LyapunovSpectrum45(std::function<std::vector<long double>(const std::vector<long double>&)> f, const std::vector<long double>& point, long double dt, int numSteps, int numIterations);

long double CalculateLargestLyapunov(std::function<std::vector<long double>(const std::vector<long double>&)> f, std::vector<long double> point, long double dt, int numSteps, int numIterations);

#endif // LYAPUNOV_ANALYSIS_H
