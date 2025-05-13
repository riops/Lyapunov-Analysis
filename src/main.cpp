// main.cpp
#include "../include/Analysis.h"
#include "../include/HTableLoader.h"
#include "../include/Systems.h"
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <sstream>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

// extend GSLParams to carry the extra parameters:
struct GSLParams {
  decltype(LoadHTable(0)) *HTablePtr;
  size_t dim;
  long double mu;
  long double lambda;
  long double R;
  int N;
};

// Updated gsl_rhs that unpacks and passes mu, lambda, R, N:
int gsl_rhs(double /*t*/, const double y[], double dydt[], void *p) {
  auto *params = static_cast<GSLParams *>(p);

  auto &HTable = *params->HTablePtr;
  size_t dim = params->dim;
  long double mu = params->mu;
  long double lambda = params->lambda;
  long double R = params->R;
  int N = params->N;

  // copy into vector<long double>
  std::vector<long double> x(dim);
  for (size_t i = 0; i < dim; ++i) {
    x[i] = static_cast<long double>(y[i]);
  }

  // call the new signature
  auto dx = averagedEquationsPolarizationBasisSymmetryReducedParallel(
      x, HTable, mu, lambda, R, N);

  // copy back to GSL array
  for (size_t i = 0; i < dim; ++i) {
    dydt[i] = static_cast<double>(dx[i]);
  }
  return GSL_SUCCESS;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // 1) ISO‐style timestamp on rank 0
  char datetime_cstr[64];
  if (world_rank == 0) {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&t);
    std::strftime(datetime_cstr, sizeof(datetime_cstr), "%Y-%m-%d_%H-%M-%S",
                  &tm);
  }
  // 2) Broadcast timestamp
  MPI_Bcast(datetime_cstr, sizeof(datetime_cstr), MPI_CHAR, 0, MPI_COMM_WORLD);
  std::string datetime{datetime_cstr};

  omp_set_num_threads(110);

  // 3) Build initial conditions
  std::vector<long double> EnergyValues = {8, 10, 15, 20, 25, 30};
  int matrixDimension = 3;
  long double mu = 10.0;
  long double lambda = 0.5;
  long double R = 1.0;
  auto HTable = LoadHTable(matrixDimension);

  std::vector<std::vector<long double>> allICs(world_size);
  for (int r = 0; r < world_size; ++r) {
    allICs[r] = GenerateInitialConditionReduced(1, matrixDimension,
                                                EnergyValues[r], lambda, R, mu);
  }
  auto myIC = allICs[world_rank];

  long double dt = 1e-3L;
  int numSteps = 50000;

  if (world_rank == 0) {
    std::cout << "Running " << world_size << " trajectories in parallel...\n";
  }

  // 4) Set up GSL integrator
  size_t dim = myIC.size();
  GSLParams params{&HTable, dim, mu, lambda, R, matrixDimension};

  gsl_odeiv2_system sys;
  sys.function = gsl_rhs;
  sys.jacobian = nullptr;
  sys.dimension = static_cast<int>(dim);
  sys.params = &params;

  gsl_odeiv2_driver *driver =
      gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
                                    static_cast<double>(dt), // initial step
                                    1e-3,                    // abs tol
                                    1e-3                     // rel tol
      );

  // storage for states and times
  std::vector<std::vector<long double>> traj;
  std::vector<double> times;
  traj.reserve(numSteps + 1);
  times.reserve(numSteps + 1);

  // initial
  traj.push_back(myIC);
  times.push_back(0.0);

  std::vector<double> y(dim);
  for (size_t i = 0; i < dim; ++i)
    y[i] = static_cast<double>(myIC[i]);

  MPI_Barrier(MPI_COMM_WORLD);
  double t0 = MPI_Wtime();

  double t = 0.0;
  for (int step = 1; step <= numSteps; ++step) {
    double ti = step * static_cast<double>(dt);
    int status = gsl_odeiv2_driver_apply(driver, &t, ti, y.data());
    if (status != GSL_SUCCESS) {
      std::cerr << "GSL error at step " << step << ": " << gsl_strerror(status)
                << "\n";
      break;
    }
    // record state
    std::vector<long double> row(dim);
    for (size_t i = 0; i < dim; ++i)
      row[i] = static_cast<long double>(y[i]);
    traj.push_back(std::move(row));
    // record time
    times.push_back(t);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double t1 = MPI_Wtime();
  if (world_rank == 0) {
    std::cout << "Pure compute time: " << (t1 - t0) << " s\n";
  }

  gsl_odeiv2_driver_free(driver);

  // 5) Write out CSV
  {
    std::ostringstream fn;
    fn << "./data/csv/EOMPolarizationBasis_" << datetime << "_"
       << EnergyValues[world_rank] << ".csv";
    std::ofstream fs(fn.str());

    // states
    for (auto &row : traj) {
      for (size_t j = 0; j < row.size(); ++j)
        fs << row[j] << (j + 1 < row.size() ? "," : "");
      fs << "\n";
    }
    // final line: times
    for (size_t i = 0; i < times.size(); ++i) {
      fs << times[i] << (i + 1 < times.size() ? "," : "");
    }
    fs << "\n";
  }

  // 6) Traced values
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
