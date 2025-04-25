#include "../include/Analysis.h"
#include "../include/FunctionalAnalysis.h"
#include "../include/Systems.h"
#include <chrono>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <sstream>
#include <vector>

struct HybridFunction {
  int rank, size;
  std::vector<long double>
  operator()(const std::vector<long double> &input) const {
    return averagedEquationsPolarizationBasisSymmetryReducedParallelMpi(
        input, rank, size);
  }
};

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int world_size, world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  omp_set_num_threads(omp_get_max_threads());

  int matrixDimension = 2;
  int dimension = 48 * matrixDimension * matrixDimension * matrixDimension *
                  matrixDimension;

  int numIterations = 10;
  long double dt = 1e-2L;

  // build your initial condition
  std::vector<long double> sigmaXX = {1.22577L, 0.603123L};
  std::vector<long double> sigmaPP = {19.8276L, 10.0575L};
  auto initialCondition =
      GenerateInitialConditionReduced(sigmaXX, sigmaPP, matrixDimension);

  if (world_rank == 0) {
    std::cout << "Initial condition generated. Starting integration…\n";
  }

  // — every rank MUST call this —
  auto spectrum = IntegrateSystem45Mpi(
      averagedEquationsPolarizationBasisSymmetryReducedParallel,
      initialCondition, dt, numIterations, true);

  // — only rank 0 writes the CSV and traced file —
  if (world_rank == 0) {
    std::cout << "\nDone!\n\nWriting to file...\n\n";

    // timestamp for filename
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    std::string datetime = oss.str();

    std::string base = "./data/csv/EOMPolarizationBasisMpi_" + datetime;
    std::ofstream fs(base + ".csv");
    // write the trajectory (first numIterations rows × dimension cols)
    for (int j = 0; j < numIterations; ++j) {
      for (int k = 0; k < dimension; ++k) {
        fs << spectrum[j][k] << (k + 1 < dimension ? "," : "");
      }
      fs << "\n";
    }
    // write the dtHistory (last row)
    auto &dtRow = spectrum.back();
    for (size_t i = 0; i < dtRow.size(); ++i) {
      fs << dtRow[i] << (i + 1 < dtRow.size() ? "," : "");
    }
    fs << "\n";
    fs.close();

    // compute traced values
    CalculateTracedValues(base + ".csv", base + "_Traced.csv", matrixDimension);
  }

  MPI_Finalize();
  return 0;
}
