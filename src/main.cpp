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

  // Let each rank spawn 110 threads in the RHS
  omp_set_num_threads(110);

  // --- 1) Build your list of initial conditions (length ≥ world_size) ---
  int matrixDimension = 2;
  std::vector<std::vector<long double>> allICs;
  // Example: here we just make 6 slight variations; replace with your real data
  for (int i = 0; i < world_size; ++i) {
    long double scale = 1.0L + 0.1L * i;
    std::vector<long double> sigmaXX = {1.22577L * scale, 0.603123L * scale};
    std::vector<long double> sigmaPP = {19.8276L * scale, 10.0575L * scale};
    allICs.push_back(
        GenerateInitialConditionReduced(sigmaXX, sigmaPP, matrixDimension));
  }

  // Make sure we have at least one IC per rank:
  if ((int)allICs.size() < world_size) {
    if (world_rank == 0)
      std::cerr << "Need at least " << world_size << " initial conditions!\n";
    MPI_Finalize();
    return 1;
  }

  // --- 2) Pick this rank’s initial condition ---
  auto myIC = allICs[world_rank];

  // Parameters
  long double dt = 1e-2L;
  int numSteps = 10;

  if (world_rank == 0) {
    std::cout << "Rank 0: starting " << world_size
              << " trajectories in parallel...\n";
  }

  // --- 3) Do the serial RK45 integrator on your one trajectory ---
  auto traj = IntegrateSystem45(
      averagedEquationsPolarizationBasisSymmetryReducedParallel, myIC, dt,
      numSteps, /*allValues=*/true);

  // --- 4) Write out your own CSV and traced file ---
  {
    std::ostringstream fn;
    fn << "./data/csv/traject_rank" << world_rank << ".csv";
    std::ofstream fs(fn.str());
    for (int step = 0; step < numSteps; ++step) {
      auto &row = traj[step];
      for (size_t j = 0; j < row.size(); ++j) {
        fs << row[j] << (j + 1 < row.size() ? "," : "");
      }
      fs << "\n";
    }
    fs.close();
  }
  {
    std::ostringstream fn;
    fn << "./data/csv/traject_rank" << world_rank << "_traced.csv";
    CalculateTracedValues(std::string("traject_rank") +
                              std::to_string(world_rank) + ".csv",
                          fn.str(), matrixDimension);
  }

  if (world_rank == 0) {
    std::cout << "All ranks done. Check data/csv/traject_rank*.csv\n";
  }

  MPI_Finalize();
  return 0;
}
