//
// Created by berko on 3/20/2023.
//
#include "../include/LyapunovAnalysis.h"
#include "../include/FunctionalAnalysis.h"
#include "../include/MatrixOperations.h"
#include <cmath>
#include <iostream>
#include <vector>

// This function takes a system as an input (system should be a vector valued
// function) and a point in the phase space and returns the combined system with
// the Jacobian matrix of the system at that point.
std::vector<long double> JacobianSystem(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point) {
  std::vector<long double> result;

  std::vector<long double> f1 = f(point);
  int dimension = f1.size();

  std::vector<long double> initialPhi;
  initialPhi.resize(dimension * dimension);
  for (int i = 0; i < dimension * dimension; i++) {
    initialPhi[i] = point[i + dimension];
  }

  std::vector<long double> jacobian = Jacobian(f, point);
  std::vector<long double> phi = MatrixProduct(initialPhi, Transpose(jacobian));
  result = ExtendVector(f1, phi);
  return result;
}

std::vector<long double> RungeKuttaJacobian4(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt) {
  std::vector<long double> k1 = JacobianSystem(f, point);
  std::vector<long double> k2 = JacobianSystem(
      f, VectorAddition(point, ScalarVectorMultiplication(dt / 2, k1)));
  std::vector<long double> k3 = JacobianSystem(
      f, VectorAddition(point, ScalarVectorMultiplication(dt / 2, k2)));
  std::vector<long double> k4 = JacobianSystem(
      f, VectorAddition(point, ScalarVectorMultiplication(dt, k3)));
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

std::vector<long double> RungeKuttaJacobian5(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt) {
  std::vector<long double> k1 = JacobianSystem(f, point);
  std::vector<long double> k2 = JacobianSystem(
      f, VectorAddition(point, ScalarVectorMultiplication(dt / 4.0, k1)));
  std::vector<long double> k3 = JacobianSystem(
      f, VectorAddition(
             point, VectorAddition(ScalarVectorMultiplication(dt / 8.0, k1),
                                   ScalarVectorMultiplication(dt / 8.0, k2))));
  std::vector<long double> k4 = JacobianSystem(
      f, VectorAddition(
             point, VectorAddition(ScalarVectorMultiplication(-dt / 2.0, k2),
                                   ScalarVectorMultiplication(dt, k3))));
  std::vector<long double> k5 = JacobianSystem(
      f, VectorAddition(
             point,
             VectorAddition(ScalarVectorMultiplication(3.0 * dt / 16.0, k1),
                            ScalarVectorMultiplication(9.0 * dt / 16, k4))));
  std::vector<long double> k6 = JacobianSystem(
      f,
      VectorAddition(
          point,
          VectorAddition(
              ScalarVectorMultiplication(-3.0 * dt / 7.0, k1),
              VectorAddition(
                  ScalarVectorMultiplication(2.0 * dt / 7.0, k2),
                  VectorAddition(
                      ScalarVectorMultiplication(12.0 * dt / 7.0, k3),
                      VectorAddition(
                          ScalarVectorMultiplication(-12.0 * dt / 7.0, k4),
                          ScalarVectorMultiplication(8.0 * dt / 7.0, k5)))))));

  // Adjusting the addition to work with a function that only adds two vectors
  // at a time
  std::vector<long double> result = VectorAddition(
      point,
      VectorAddition(
          ScalarVectorMultiplication(7.0 * dt / 90.0, k1),
          VectorAddition(
              ScalarVectorMultiplication(32.0 * dt / 90.0, k3),
              VectorAddition(
                  ScalarVectorMultiplication(12.0 * dt / 90.0, k4),
                  VectorAddition(
                      ScalarVectorMultiplication(32.0 * dt / 90.0, k5),
                      ScalarVectorMultiplication(7.0 * dt / 90.0, k6))))));
  return result;
}

// This function integrates the Jacobian system for a given number of steps and
// returns the final point.
std::vector<long double> IntegrateJacobianSystem4(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps) {
  std::vector<long double> currentPoint = point;
  for (int i = 0; i < numSteps; i++) {
    currentPoint = RungeKuttaJacobian4(f, currentPoint, dt);
  }
  return currentPoint;
}

std::vector<long double> IntegrateJacobianSystem5(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps) {
  std::vector<long double> currentPoint = point;
  for (int i = 0; i < numSteps; i++) {
    currentPoint = RungeKuttaJacobian5(f, currentPoint, dt);
  }
  return currentPoint;
}

std::vector<long double> IntegrateJacobianSystem45(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps) {
  std::vector<long double> currentPoint = point;
  for (int i = 0; i < numSteps; i++) {
    currentPoint = RungeKuttaJacobian5(f, currentPoint, dt);
  }
  return currentPoint;
}

// This function takes a point in the phase space and returns the phi matrix.
std::vector<long double> GetPhi(const std::vector<long double> &point,
                                int dimension) {
  std::vector<long double> phi;
  phi.resize(dimension * dimension);
  for (int i = 0; i < dimension * dimension; i++) {
    phi[i] = point[i + dimension];
  }
  return phi;
}

// This function takes in a vector and returns the natural logarithm of the
// norms of the vectors.
std::vector<long double> LogNorms(const std::vector<long double> &vector1) {
  std::vector<long double> result;
  int dimension = sqrt(vector1.size());
  for (int i = 0; i < dimension; i++) {
    std::vector<long double> currentVector;
    currentVector.resize(dimension);
    for (int j = 0; j < dimension; j++) {
      currentVector[j] = vector1[i * dimension + j];
    }
    result.push_back(log(Norm(currentVector)));
  }
  return result;
}

// This function takes a vector and normalizes it.
std::vector<long double> Normalize(const std::vector<long double> &vector1) {
  std::vector<long double> result;
  int dimension = sqrt(vector1.size());
  for (int i = 0; i < dimension; i++) {
    std::vector<long double> currentVector;
    currentVector.resize(dimension);
    for (int j = 0; j < dimension; j++) {
      currentVector[j] = vector1[i * dimension + j];
    }
    currentVector =
        ScalarVectorMultiplication(1 / Norm(currentVector), currentVector);
    for (int j = 0; j < dimension; j++) {
      result.push_back(currentVector[j]);
    }
  }
  return result;
}

// This function takes a point in the phase space and a vector and replaces the
// phi matrix that is given as vector with the given vector.
std::vector<long double> ReplacePhi(const std::vector<long double> &point,
                                    const std::vector<long double> &newPhi) {
  std::vector<long double> result;
  int dimension = sqrt(newPhi.size());
  for (int i = 0; i < dimension; i++) {
    result.push_back(point[i]);
  }
  for (int i = 0; i < dimension * dimension; i++) {
    result.push_back(newPhi[i]);
  }
  return result;
}

// This function takes a vector of vectors and sums the vectors.
std::vector<std::vector<long double>>
SumLogNorms(const std::vector<std::vector<long double>> &logNorms,
            long double dt, int numSteps) {
  std::vector<std::vector<long double>> result;
  int dimension = logNorms[0].size();
  result.resize(logNorms.size());
  for (int i = 0; i < logNorms.size(); i++) {
    result[i].resize(dimension);
  }

  for (int i = 0; i < dimension; i++) {
    long double lyapunov = 0.0;
    for (int j = 0; j < logNorms.size(); j++) {
      lyapunov += logNorms[j][i];
      result[j][i] = lyapunov / ((j + 1) * numSteps * dt);
    }
  }
  return result;
}

std::vector<std::vector<long double>> LyapunovSpectrum4(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps,
    int numIterations) {
  std::vector<std::vector<long double>> result;
  std::vector<std::vector<long double>> logNorms;
  int dimension = f(point).size();
  std::cout << "Starting Lyapunov spectrum calculation...\n";

  std::vector<long double> currentPoint = point;
  for (int i = 0; i < numIterations; i++) {
    // cout << "\033[A\33[2K\rIteration number: " << i+1 << " / " <<
    // numIterations << "\n";
    currentPoint = IntegrateJacobianSystem4(f, currentPoint, dt, numSteps);
    std::vector<long double> phi = GetPhi(currentPoint, dimension);
    phi = Orthogonalize(phi);

    logNorms.push_back(LogNorms(phi));

    phi = Normalize(phi);

    currentPoint = ReplacePhi(currentPoint, phi);
  }

  result = SumLogNorms(logNorms, dt, numSteps);
  return result;
}

std::vector<std::vector<long double>> LyapunovSpectrum5(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps,
    int numIterations) {
  std::vector<std::vector<long double>> result;
  std::vector<std::vector<long double>> logNorms;
  int dimension = f(point).size();
  std::cout << "Starting Lyapunov spectrum calculation...\n";

  std::vector<long double> currentPoint = point;
  for (int i = 0; i < numIterations; i++) {
    // cout << "\033[A\33[2K\rIteration number: " << i+1 << " / " <<
    // numIterations << "\n";
    currentPoint = IntegrateJacobianSystem5(f, currentPoint, dt, numSteps);
    std::vector<long double> phi = GetPhi(currentPoint, dimension);
    phi = Orthogonalize(phi);

    logNorms.push_back(LogNorms(phi));

    phi = Normalize(phi);

    currentPoint = ReplacePhi(currentPoint, phi);
  }

  result = SumLogNorms(logNorms, dt, numSteps);
  return result;
}

std::vector<std::vector<long double>> LyapunovSpectrum45(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    const std::vector<long double> &point, long double dt, int numSteps,
    int numIterations) {
  std::vector<std::vector<long double>> result;
  std::vector<std::vector<long double>> logNorms;
  int dimension = f(point).size();
  std::cout << "Starting Lyapunov spectrum calculation...\n";

  std::vector<long double> currentPoint = point;
  for (int i = 0; i < numIterations; i++) {
    // cout << "\033[A\33[2K\rIteration number: " << i+1 << " / " <<
    // numIterations << "\n";
    currentPoint = IntegrateJacobianSystem45(f, currentPoint, dt, numSteps);
    std::vector<long double> phi = GetPhi(currentPoint, dimension);
    phi = Orthogonalize(phi);

    logNorms.push_back(LogNorms(phi));

    phi = Normalize(phi);

    currentPoint = ReplacePhi(currentPoint, phi);
  }

  result = SumLogNorms(logNorms, dt, numSteps);
  return result;
}

long double CalculateLargestLyapunov(
    std::function<std::vector<long double>(const std::vector<long double> &)> f,
    std::vector<long double> point, long double dt, int numSteps,
    int numIterations) {
  int dimension = point.size();

  std::vector<long double> epsilon;
  epsilon.resize(dimension);

  long double seperationSize = 1e-8;

  for (int i = 0; i < dimension; i++) {
    int random = rand();
    long double randomLongDouble = random * seperationSize;
    epsilon[i] = randomLongDouble;
  }

  long double initialNorm = Norm(epsilon);
  epsilon = ScalarVectorMultiplication(seperationSize / initialNorm, epsilon);

  std::vector<long double> adjustedPoint = VectorAddition(point, epsilon);

  std::vector<long double> logNorms;

  long double totalTime = 0.0;

  for (int i = 0; i < numIterations; i++) {
    std::vector<std::vector<long double>> finalResult =
        LyapunovSpectrum45(f, point, dt, numSteps, true);

    std::vector<long double> finalPoint = finalResult[numSteps - 1];
    std::vector<long double> finalTimeList = finalResult[numSteps];
    PrintVector(finalTimeList);

    std::vector<long double> adjustedFinalPoint =
        IntegrateSystem4(f, adjustedPoint, dt, numSteps, true)[numSteps - 1];

    totalTime += SumValues(finalTimeList);

    long double finalPointNorm = Norm(finalPoint);
    long double adjustedFinalPointNorm = Norm(adjustedFinalPoint);

    std::vector<long double> delta =
        VectorSubtraction(adjustedFinalPoint, finalPoint);
    long double normDelta = Norm(delta);

    long double logNorm = log(normDelta / seperationSize);

    logNorms.push_back(logNorm);
    // cout << "Logarithm of the iteration no " << i << " is " << logNorm << "
    // norm delta is " << normDelta << endl;

    std::vector<long double> normalizedDelta =
        ScalarVectorMultiplication(seperationSize / normDelta, delta);

    point = finalPoint;

    adjustedPoint = VectorAddition(finalPoint, normalizedDelta);
  }

  long double largestLyapunov = 0.0;

  for (int i = 0; i < numIterations; i++) {
    largestLyapunov += logNorms[i];
  }

  return largestLyapunov / totalTime;
}
