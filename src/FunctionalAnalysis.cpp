//
// Created by berko on 3/20/2023.
//
#include "../include/FunctionalAnalysis.h"
#include "../include/MatrixOperations.h"
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <mpi.h>
#include <vector>

// This function takes in a vector valued function and a point in the phase
// space, and returns the jacobian matrix of the function at that point as a
// vector.
std::vector<long double> Jacobian(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point) {
  std::vector<long double> result;
  long double epsilon = 0.001;
  int dimension = f(point).size();
  result.resize(dimension * dimension);
  std::vector<long double> point1 = point;
  std::vector<long double> point2 = point;

  for (int i = 0; i < dimension; i++) {
    point1[i] += epsilon;
    point2[i] -= epsilon;
    std::vector<long double> f1 = f(point1);
    std::vector<long double> f2 = f(point2);
    for (int j = 0; j < dimension; j++) {
      result[j * dimension + i] = (f1[j] - f2[j]) / (2 * epsilon);
    }
    point1[i] -= epsilon;
    point2[i] += epsilon;
  }
  return result;
}

// This function takes in a vector valued function and a point in the phase
// space, then integrates it using the Runge-Kutta method of order 4. The
// function returns the value of the function at the next time step.
std::vector<long double> RungeKutta4(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt) {
  std::vector<long double> k1 = f(point);
  std::vector<long double> k2 =
      f(VectorAddition(point, ScalarVectorMultiplication(dt / 2, k1)));
  std::vector<long double> k3 =
      f(VectorAddition(point, ScalarVectorMultiplication(dt / 2, k2)));
  std::vector<long double> k4 =
      f(VectorAddition(point, ScalarVectorMultiplication(dt, k3)));
  std::vector<long double> result = VectorAddition(
      point,
      ScalarVectorMultiplication(
          dt / 6,
          VectorAddition(
              k1, VectorAddition(
                      ScalarVectorMultiplication(2, k2),
                      VectorAddition(ScalarVectorMultiplication(2, k3), k4)))));
  return result;
}

// This function takes in a vector valued function and a point in the phase
// space, then integrates it using the Runge-Kutta method of order 4. The
// function returns the value of the function at the next time step.
std::vector<long double> RungeKutta5(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt) {
  std::vector<long double> k1 = f(point);
  std::vector<long double> k2 =
      f(VectorAddition(point, ScalarVectorMultiplication(dt / 4, k1)));
  std::vector<long double> k3 = f(VectorAddition(
      point, VectorAddition(ScalarVectorMultiplication(dt / 8, k1),
                            ScalarVectorMultiplication(dt / 8, k2))));
  std::vector<long double> k4 = f(VectorAddition(
      point, VectorAddition(ScalarVectorMultiplication(-dt / 2, k2),
                            ScalarVectorMultiplication(dt, k3))));
  std::vector<long double> k5 = f(VectorAddition(
      point, VectorAddition(ScalarVectorMultiplication(3 * dt / 16, k1),
                            ScalarVectorMultiplication(9 * dt / 16, k4))));
  std::vector<long double> k6 = f(VectorAddition(
      point, VectorAddition(
                 ScalarVectorMultiplication(-3 * dt / 7, k1),
                 VectorAddition(
                     ScalarVectorMultiplication(2 * dt / 7, k2),
                     VectorAddition(
                         ScalarVectorMultiplication(12 * dt / 7, k3),
                         VectorAddition(
                             ScalarVectorMultiplication(-12 * dt / 7, k4),
                             ScalarVectorMultiplication(8 * dt / 7, k5)))))));

  // Adjusting the addition to work with a function that only adds two vectors
  // at a time
  std::vector<long double> result = VectorAddition(
      point, VectorAddition(
                 ScalarVectorMultiplication(7 * dt / 90, k1),
                 VectorAddition(
                     ScalarVectorMultiplication(32 * dt / 90, k3),
                     VectorAddition(
                         ScalarVectorMultiplication(12 * dt / 90, k4),
                         VectorAddition(
                             ScalarVectorMultiplication(32 * dt / 90, k5),
                             ScalarVectorMultiplication(7 * dt / 90, k6))))));
  return result;
}

long double AdjustStepSize(long double dt, long double error, long double tol) {
  const long double safety = 0.9L;
  const long double order = 5.0L; // for RK45 use exponent 1/(order)
  const long double p = 1.0L / order;

  if (error <= 0.0L) {
    // no error ⇒ crank it up
    return dt * 2.0L;
  }

  long double scale = safety * std::pow(tol / error, p);
  // clamp to [0.1, 5.0]
  scale = std::max(0.1L, std::min(5.0L, scale));
  return dt * scale;
}

// Returns {y4, y5}: 4th‐order and 5th‐order estimates at point+dt
std::vector<std::vector<long double>> RungeKutta45(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt) {
  // Fehlberg coefficients
  const long double a21 = 1.0L / 4;
  const long double a31 = 3.0L / 32, a32 = 9.0L / 32;
  const long double a41 = 1932.0L / 2197, a42 = -7200.0L / 2197,
                    a43 = 7296.0L / 2197;
  const long double a51 = 439.0L / 216, a52 = -8.0L, a53 = 3680.0L / 513,
                    a54 = -845.0L / 4104;
  const long double a61 = -8.0L / 27, a62 = 2.0L, a63 = -3544.0L / 2565,
                    a64 = 1859.0L / 4104, a65 = -11.0L / 40;

  const long double b4_1 = 25.0L / 216, b4_3 = 1408.0L / 2565,
                    b4_4 = 2197.0L / 4104, b4_5 = -1.0L / 5;
  const long double b5_1 = 16.0L / 135, b5_3 = 6656.0L / 12825,
                    b5_4 = 28561.0L / 56430, b5_5 = -9.0L / 50,
                    b5_6 = 2.0L / 55;

  const long double tolerance = 1e-5L;

  // Stage 1
  auto k1 = f(point);

  // Stage 2
  auto y2 = VectorAddition(point, ScalarVectorMultiplication(dt * a21, k1));
  auto k2 = f(y2);

  // Stage 3
  auto tmp3 = VectorAddition(ScalarVectorMultiplication(dt * a31, k1),
                             ScalarVectorMultiplication(dt * a32, k2));
  auto y3 = VectorAddition(point, tmp3);
  auto k3 = f(y3);

  // Stage 4
  auto tmp4 =
      VectorAddition(ScalarVectorMultiplication(dt * a41, k1),
                     VectorAddition(ScalarVectorMultiplication(dt * a42, k2),
                                    ScalarVectorMultiplication(dt * a43, k3)));
  auto y4 = VectorAddition(point, tmp4);
  auto k4 = f(y4);

  // Stage 5
  auto tmp5 = VectorAddition(
      ScalarVectorMultiplication(dt * a51, k1),
      VectorAddition(ScalarVectorMultiplication(dt * a52, k2),
                     VectorAddition(ScalarVectorMultiplication(dt * a53, k3),
                                    ScalarVectorMultiplication(dt * a54, k4))));
  auto y5 = VectorAddition(point, tmp5);
  auto k5 = f(y5);

  // Stage 6
  auto tmp6 = VectorAddition(
      ScalarVectorMultiplication(dt * a61, k1),
      VectorAddition(
          ScalarVectorMultiplication(dt * a62, k2),
          VectorAddition(
              ScalarVectorMultiplication(dt * a63, k3),
              VectorAddition(ScalarVectorMultiplication(dt * a64, k4),
                             ScalarVectorMultiplication(dt * a65, k5)))));
  auto y6 = VectorAddition(point, tmp6);
  auto k6 = f(y6);

  // Combine for 4th- and 5th-order results
  size_t M = point.size();
  std::vector<long double> res4(M), res5(M);
  for (size_t i = 0; i < M; ++i) {
    res4[i] = point[i] +
              dt * (b4_1 * k1[i] + b4_3 * k3[i] + b4_4 * k4[i] + b4_5 * k5[i]);
    res5[i] = point[i] + dt * (b5_1 * k1[i] + b5_3 * k3[i] + b5_4 * k4[i] +
                               b5_5 * k5[i] + b5_6 * k6[i]);
  }

  long double error2 = 0.0L;
  for (size_t i = 0; i < point.size(); ++i) {
    long double d = res5[i] - res4[i];
    error2 += d * d;
  }
  long double error = std::sqrt(error2);
  long double dt_new = AdjustStepSize(dt, error, tolerance);
  std::vector<long double> dtReturn = {dt_new};
  return {res4, dtReturn};
}

// assume VectorAddition, ScalarVectorMultiplication and AdjustStepSize are
// available

std::vector<std::vector<long double>> RungeKutta45Mpi(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const size_t M = point.size();
  const long double tolerance = 1e-5L;

  // Fehlberg coefficients (as in your original)
  const long double a21 = 1.0L / 4;
  const long double a31 = 3.0L / 32, a32 = 9.0L / 32;
  const long double a41 = 1932.0L / 2197, a42 = -7200.0L / 2197,
                    a43 = 7296.0L / 2197;
  const long double a51 = 439.0L / 216, a52 = -8.0L, a53 = 3680.0L / 513,
                    a54 = -845.0L / 4104;
  const long double a61 = -8.0L / 27, a62 = 2.0L, a63 = -3544.0L / 2565,
                    a64 = 1859.0L / 4104, a65 = -11.0L / 40;

  const long double b4_1 = 25.0L / 216, b4_3 = 1408.0L / 2565,
                    b4_4 = 2197.0L / 4104, b4_5 = -1.0L / 5;
  const long double b5_1 = 16.0L / 135, b5_3 = 6656.0L / 12825,
                    b5_4 = 28561.0L / 56430, b5_5 = -9.0L / 50,
                    b5_6 = 2.0L / 55;

  // storage for the six k‑vectors
  std::vector<std::vector<long double>> k(6, std::vector<long double>(M));

  // helper to compute and broadcast stage i
  auto do_stage = [&](int stage) {
    int root = (stage - 1) % size; // round‑robin assignment
    std::vector<long double> y(M);

    // build the stage input y_i = point + dt * sum_j a_{ij} k_j
    switch (stage) {
    case 1:
      y = point;
      break;
    case 2:
      y = VectorAddition(point, ScalarVectorMultiplication(dt * a21, k[0]));
      break;
    case 3: {
      auto t = VectorAddition(ScalarVectorMultiplication(dt * a31, k[0]),
                              ScalarVectorMultiplication(dt * a32, k[1]));
      y = VectorAddition(point, t);
    } break;
    case 4: {
      auto t = VectorAddition(
          ScalarVectorMultiplication(dt * a41, k[0]),
          VectorAddition(ScalarVectorMultiplication(dt * a42, k[1]),
                         ScalarVectorMultiplication(dt * a43, k[2])));
      y = VectorAddition(point, t);
    } break;
    case 5: {
      auto t = VectorAddition(
          ScalarVectorMultiplication(dt * a51, k[0]),
          VectorAddition(
              ScalarVectorMultiplication(dt * a52, k[1]),
              VectorAddition(ScalarVectorMultiplication(dt * a53, k[2]),
                             ScalarVectorMultiplication(dt * a54, k[3]))));
      y = VectorAddition(point, t);
    } break;
    case 6: {
      auto t = VectorAddition(
          ScalarVectorMultiplication(dt * a61, k[0]),
          VectorAddition(
              ScalarVectorMultiplication(dt * a62, k[1]),
              VectorAddition(
                  ScalarVectorMultiplication(dt * a63, k[2]),
                  VectorAddition(ScalarVectorMultiplication(dt * a64, k[3]),
                                 ScalarVectorMultiplication(dt * a65, k[4])))));
      y = VectorAddition(point, t);
    } break;
    }

    // only the assigned root computes f(y)
    if (rank == root) {
      k[stage - 1] = f(y);
    }
    // then broadcast k[stage-1] from that root to all
    MPI_Bcast(k[stage - 1].data(), M, MPI_LONG_DOUBLE, root, MPI_COMM_WORLD);
  };

  // do all six stages
  for (int i = 1; i <= 6; ++i) {
    do_stage(i);
  }

  // now rank 0 assembles the 4th‐order result and new dt
  std::vector<long double> res4(M), dt_ret(1);
  if (rank == 0) {
    // compute 4th‐ and 5th‐order estimates
    std::vector<long double> res5(M);
    for (size_t i = 0; i < M; ++i) {
      res4[i] = point[i] + dt * (b4_1 * k[0][i] + b4_3 * k[2][i] +
                                 b4_4 * k[3][i] + b4_5 * k[4][i]);
      res5[i] =
          point[i] + dt * (b5_1 * k[0][i] + b5_3 * k[2][i] + b5_4 * k[3][i] +
                           b5_5 * k[4][i] + b5_6 * k[5][i]);
    }
    // error estimate and adjust dt
    long double err2 = 0;
    for (size_t i = 0; i < M; ++i) {
      long double d = res5[i] - res4[i];
      err2 += d * d;
    }
    long double err = std::sqrt(err2);
    dt_ret[0] = AdjustStepSize(dt, err, tolerance);
  }

  // broadcast the final res4 and new dt back to everyone
  MPI_Bcast(res4.data(), M, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(dt_ret.data(), 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);

  return {res4, dt_ret};
}

std::vector<std::vector<long double>> IntegrateSystem45(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps,
    bool allValues) {
  int dimension = f(point).size();

  std::vector<std::vector<long double>> result;
  std::vector<long double> currentPoint = point;
  std::vector<long double> dtVector = {dt};

  for (int i = 0; i < numSteps; i++) {
    using Clock = std::chrono::high_resolution_clock;
    auto start = Clock::now();

    std::vector<std::vector<long double>> rk45Result =
        RungeKutta45(f, currentPoint, dt);
    currentPoint = rk45Result[0];
    dt = rk45Result[1][0];
    dtVector.push_back(dt + dtVector[i]);

    auto end = Clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Elapsed time: " << duration.count()
              << ". Estimated Time: " << duration.count() * (numSteps - i)
              << "ms\n";

    if (allValues) {
      result.push_back(currentPoint);
    }
  }
  if (!allValues) {
    result.push_back(currentPoint);
  }

  result.push_back(dtVector);
  return result;
}

std::vector<std::vector<long double>> IntegrateSystem45Mpi(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps,
    bool allValues) {
  int dimension = f(point).size();

  std::vector<std::vector<long double>> result;
  std::vector<long double> currentPoint = point;
  std::vector<long double> dtVector = {dt};

  for (int i = 0; i < numSteps; i++) {
    using Clock = std::chrono::high_resolution_clock;
    auto start = Clock::now();

    std::vector<std::vector<long double>> rk45Result =
        RungeKutta45Mpi(f, currentPoint, dt);
    currentPoint = rk45Result[0];
    dt = rk45Result[1][0];
    dtVector.push_back(dt + dtVector[i]);

    auto end = Clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Elapsed time: " << duration.count()
              << ". Estimated Time: " << duration.count() * (numSteps - i)
              << "ms\n";

    if (allValues) {
      result.push_back(currentPoint);
    }
  }
  if (!allValues) {
    result.push_back(currentPoint);
  }

  result.push_back(dtVector);
  return result;
}

std::vector<std::vector<long double>> IntegrateSystem4(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps,
    bool allValues) {
  int dimension = f(point).size();

  std::vector<std::vector<long double>> result;
  std::vector<long double> currentPoint = point;
  std::vector<long double> dtVector = {dt};

  for (int i = 0; i < numSteps; i++) {
    currentPoint = RungeKutta4(f, currentPoint, dt);
    if (allValues) {
      result.push_back(currentPoint);
    }
  }
  if (!allValues) {
    result.push_back(currentPoint);
  }

  result.push_back(dtVector);
  return result;
}

std::vector<std::vector<long double>> IntegrateSystem5(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps,
    bool allValues) {
  int dimension = f(point).size();

  std::vector<std::vector<long double>> result;
  std::vector<long double> currentPoint = point;
  std::vector<long double> dtVector = {dt};

  for (int i = 0; i < numSteps; i++) {
    currentPoint = RungeKutta5(f, currentPoint, dt);
    if (allValues) {
      result.push_back(currentPoint);
    }
  }
  if (!allValues) {
    result.push_back(currentPoint);
  }

  result.push_back(dtVector);
  return result;
}

std::vector<std::vector<long double>> IntegrateSystem45_MpiRoot(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &initialPoint, long double dt, int numSteps,
    bool allValues, int world_rank, int world_size) {
  // dimension of your phase‑space
  int dimension = static_cast<int>(initialPoint.size());

  // the “shared” variables
  std::vector<long double> currentPoint = initialPoint;
  std::vector<std::vector<long double>> result;
  std::vector<long double> dtHistory;
  dtHistory.reserve(numSteps + 1);
  dtHistory.push_back(dt);

  // Record the initial data.
  if (allValues && world_rank == 0) {
    result.push_back(currentPoint);
  }

  for (int step = 0; step < numSteps; ++step) {
    // ——— 1) collective RK45: every rank calls this,
    //                      each sub‑call to f(...) does an MPI_Reduce
    //                      internally
    using Clock = std::chrono::high_resolution_clock;
    auto start = Clock::now();

    auto rk45res = RungeKutta45(f, currentPoint, dt);
    auto &nextPoint = rk45res[0];
    long double newDt = rk45res[1][0];

    auto end = Clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // ——— 2) only rank 0 updates & records
    if (world_rank == 0) {
      std::cout << "Elapsed time: " << duration.count()
                << ". Estimated Time: " << duration.count() * (numSteps - step)
                << "ms\n";
      currentPoint = nextPoint;
      dt = newDt;

      if (allValues) {
        result.push_back(currentPoint);
      }
      dtHistory.push_back(dtHistory.back() + dt);
    }

    // ——— 3) broadcast the new state/step‐size to all ranks
    MPI_Bcast(currentPoint.data(), dimension, MPI_LONG_DOUBLE, 0,
              MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
  }

  // at the end, only rank 0 returns the full trajectory + dtHistory
  if (world_rank == 0) {
    result.push_back(dtHistory);
    return result;
  } else {
    return {}; // dummy on non‑roots
  }
}
