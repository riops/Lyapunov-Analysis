#ifndef FUNCTIONAL_ANALYSIS_H
#define FUNCTIONAL_ANALYSIS_H

#include <functional>
#include <vector>

std::vector<long double> Jacobian(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point);

std::vector<long double> RungeKutta4(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt);

std::vector<long double> RungeKutta5(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt);

long double AdjustStepSize(long double dt, long double error,
                           long double tolerance);

std::vector<std::vector<long double>> RungeKutta45(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt);
std::vector<std::vector<long double>>
RungeKutta45Mpi(const std::function<std::vector<long double>(
                    const std::vector<long double> &)> &f,
                const std::vector<long double> &point_full, long double dt);

std::vector<std::vector<long double>> IntegrateSystem45(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps,
    bool allValues);

std::vector<std::vector<long double>>
IntegrateSystem45Mpi(const std::function<std::vector<long double>(
                         const std::vector<long double> &)> &f,
                     const std::vector<long double> &init, long double dt,
                     int numSteps, bool allValues);

std::vector<std::vector<long double>> IntegrateSystem4(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps,
    bool allValues);

std::vector<std::vector<long double>> IntegrateSystem5(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps,
    bool allValues);

std::vector<std::vector<long double>> IntegrateSystem45_MpiRoot(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &initialPoint, long double dt, int numSteps,
    bool allValues, int world_rank, int world_size);

#endif // FUNCTIONAL_ANALYSIS_H
