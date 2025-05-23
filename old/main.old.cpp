#include "../include/Analysis.h"
#include "../include/FunctionalAnalysis.h" // for IntegrateSystem45
#include "../include/HTableLoader.h"
#include "../include/Systems.h" // for averagedEquations…Parallel
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <sstream>
#include <vector>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // 1) Make one ISO‐style timestamp on rank 0
  char datetime_cstr[64];
  if (world_rank == 0) {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&t);
    // e.g. "2025-04-25_19-34-02"
    std::strftime(datetime_cstr, sizeof(datetime_cstr), "%Y-%m-%d_%H-%M-%S",
                  &tm);
  }
  // 2) Broadcast that timestamp to everyone
  MPI_Bcast(datetime_cstr, sizeof(datetime_cstr), MPI_CHAR, 0, MPI_COMM_WORLD);
  std::string datetime{datetime_cstr};

  omp_set_num_threads(110);

  // 3) Build initial conditions (one per rank)
  std::vector<long double> EnergyValues = {5, 10, 15, 20, 25, 30};
  int matrixDimension = 3;
  long double mu = 0.025;
  long double lambda = 10.0;
  long double R = 0.2;
  auto HTable = LoadHTable(matrixDimension);
  std::vector<std::vector<long double>> allICs(world_size);
  for (int r = 0; r < world_size; ++r) {
    allICs[r] = GenerateInitialConditionReduced(1, matrixDimension,
                                                EnergyValues[r], lambda, R, mu);
  }
  auto myIC = allICs[world_rank];

  long double dt = 1e-3L;
  int numSteps = 16000;

  if (world_rank == 0) {
    std::cout << "Running " << world_size << " trajectories in parallel...\n";
  }

  auto AveragedFunction = [&HTable](const std::vector<long double> &allVectors)
      -> std::vector<long double> {
    return averagedEquationsPolarizationBasisSymmetryReducedParallelTrial(
        allVectors, HTable);
  };

  // 4) Integrate serial-RK45 with OpenMP RHS
  MPI_Barrier(MPI_COMM_WORLD);
  double t0 = MPI_Wtime();
  auto traj = IntegrateSystem4(AveragedFunction, myIC, dt, numSteps,
                               /*allValues=*/true);
  MPI_Barrier(MPI_COMM_WORLD);
  double t1 = MPI_Wtime();
  if (world_rank == 0) {
    std::cout << "Pure compute time: " << (t1 - t0) << " s\n";
  }

  // 5) Write out trajectory CSV:
  //    EOMPolarizationBasis_<datetime>_<rank>.csv
  {
    std::ostringstream fn;
    fn << "./data/csv/EOMPolarizationBasis_" << datetime << "_"
       << EnergyValues[world_rank] << ".csv";
    std::ofstream fs(fn.str());
    for (size_t step = 0; step < traj.size(); ++step) {
      auto &row = traj[step];
      for (size_t j = 0; j < row.size(); ++j)
        fs << row[j] << (j + 1 < row.size() ? "," : "");
      fs << "\n";
    }
  }

  {
    long double E = EnergyValues[world_rank];
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(0) << E;
    std::string energy_str = oss.str();
    CalculateTracedValues(std::string("EOMPolarizationBasis_") + datetime +
                              "_" + energy_str + ".csv",
                          std::string("EOMPolarizationBasis_") + datetime +
                              "_" + energy_str + "_Traced.csv",
                          matrixDimension);
  }

  if (world_rank == 0) {
    std::cout << "All done; files are EOMPolarizationBasis_" << datetime
              << "_<rank>.csv\n";
  }

  MPI_Finalize();
  return 0;
}
