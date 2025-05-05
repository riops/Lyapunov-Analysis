//
// Created by berko on 3/20/2023.
//
#include "../include/MatrixOperations.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

// This function prints a matrix that is given as a vector.
void PrintMatrix(const std::vector<long double> &matrix) {
  int dimension = sqrt(matrix.size());
  for (int i = 0; i < dimension; i++) {
    std::cout << "| ";
    for (int j = 0; j < dimension; j++) {
      std::cout << matrix[i * dimension + j] << " ";
    }
    std::cout << "|";
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

// This function writes matrix to the given ofstream.
void WriteMatrix(const std::vector<long double> &matrix) {
  std::ofstream file("IntegrationResults.txt", std::ios_base::app);
  int dimension = sqrt(matrix.size());
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      file << matrix[i * dimension + j] << " ";
    }
    file << std::endl;
  }
  file << std::endl << std::endl;
}

// This function prints a vector.
void PrintVector(const std::vector<long double> &vector) {
  int dimension = vector.size();
  for (int i = 0; i < dimension; i++) {
    std::cout << vector[i] << " ";
  }
  std::cout << std::endl;
}

std::vector<long double> ZeroMatrix(int dimension) {
  std::vector<long double> result;
  result.resize(dimension * dimension);
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      result[i * dimension + j] = 0;
    }
  }
  return result;
}

// This function returns the identity matrix of a given dimension as a vector.
std::vector<long double> IdentityMatrix(int dimension) {
  std::vector<long double> result;
  result.resize(dimension * dimension);
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      if (i == j) {
        result[i * dimension + j] = 1;
      } else {
        result[i * dimension + j] = 0;
      }
    }
  }
  return result;
}

// This function extends a vector with another vector.
std::vector<long double> ExtendVector(const std::vector<long double> &vector1,
                                      const std::vector<long double> &vector2) {
  std::vector<long double> result;
  int dimension1 = vector1.size();
  int dimension2 = vector2.size();
  result.resize(dimension1 + dimension2);
  for (int i = 0; i < dimension1; i++) {
    result[i] = vector1[i];
  }
  for (int i = 0; i < dimension2; i++) {
    result[dimension1 + i] = vector2[i];
  }
  return result;
}

// This function adds two vectors/ or matrices that are given as vectors.
std::vector<long double>
VectorAddition(const std::vector<long double> &vector1,
               const std::vector<long double> &vector2) {
  std::vector<long double> result;
  int dimension = vector1.size();
  result.resize(dimension);
  for (int i = 0; i < dimension; i++) {
    result[i] = vector1[i] + vector2[i];
  }
  return result;
}

// This function subtracts two vectors/ or matrices that are given as vectors.
std::vector<long double>
VectorSubtraction(const std::vector<long double> &vector1,
                  const std::vector<long double> &vector2) {
  std::vector<long double> result;
  int dimension = vector1.size();
  result.resize(dimension);
  for (int i = 0; i < dimension; i++) {
    result[i] = vector1[i] - vector2[i];
  }
  return result;
}

// This function calculates the transpose of a matrix that is given as a vector.
std::vector<long double> Transpose(const std::vector<long double> &matrix) {
  std::vector<long double> result;
  int dimension = sqrt(matrix.size());
  result.resize(dimension * dimension);
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      result[i * dimension + j] = matrix[j * dimension + i];
    }
  }
  return result;
}

// This function calculates the product of two matrices that are given as
// vectors.
std::vector<long double>
MatrixProduct(const std::vector<long double> &matrix1,
              const std::vector<long double> &matrix2) {
  std::vector<long double> result;
  int dimension = sqrt(matrix1.size());
  result.resize(dimension * dimension);
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      for (int k = 0; k < dimension; k++) {
        result[i * dimension + j] +=
            matrix1[i * dimension + k] * matrix2[k * dimension + j];
      }
    }
  }
  return result;
}

// This function calculates the product of a number and a vector.
std::vector<long double>
ScalarVectorMultiplication(long double scalar,
                           const std::vector<long double> &vector1) {
  std::vector<long double> result;
  int dimension = vector1.size();
  result.resize(dimension);

  for (int i = 0; i < dimension; i++) {
    result[i] = vector1[i] * scalar;
  }
  return result;
}

// This function calculates the inner product of two vectors.
long double InnerProduct(const std::vector<long double> &vector1,
                         const std::vector<long double> &vector2) {
  long double result = 0.0;
  int dimension = vector1.size();
  for (int i = 0; i < dimension; i++) {
    result += vector1[i] * vector2[i];
  }
  return result;
}

// This function calculates the norm of a given vector.
long double Norm(const std::vector<long double> &vector1) {
  long double result;
  result = InnerProduct(vector1, vector1);
  return sqrt(result);
}

// This function calculates the projection of a vector onto another vector.
std::vector<long double> Projection(const std::vector<long double> &vector1,
                                    const std::vector<long double> &vector2) {
  std::vector<long double> result;
  result.resize(vector1.size());
  if (Norm(vector2) == 0) {
    return result;
  } else {
    long double scalar =
        InnerProduct(vector1, vector2) / InnerProduct(vector2, vector2);
    result = ScalarVectorMultiplication(scalar, vector2);
    return result;
  }
}

// This function orthonormalizes given vectors. Vectors are given as a vector,
// if the given vector is n dimensional then there are sqrt(n) vectors that are
// sqrt(n) dimensional. It uses the Gram-Schmidt process.

std::vector<long double>
Orthogonalize(const std::vector<long double> &vectors) {
  std::vector<long double> result;
  int dimension = sqrt(vectors.size());
  result.resize(dimension * dimension);
  for (int i = 0; i < dimension; i++) {
    std::vector<long double> currentVector;
    currentVector.resize(dimension);
    for (int j = 0; j < dimension; j++) {
      currentVector[j] = vectors[i * dimension + j];
    }
    // PrintVector(currentVector);
    for (int j = 0; j < i; j++) {
      std::vector<long double> previousVector;
      previousVector.resize(dimension);
      for (int k = 0; k < dimension; k++) {
        previousVector[k] = result[j * dimension + k];
      }
      // PrintVector(previousVector);
      currentVector = VectorSubtraction(
          currentVector, Projection(currentVector, previousVector));
    }
    for (int j = 0; j < dimension; j++) {
      result[i * dimension + j] = currentVector[j];
    }
  }
  return result;
}

// This function will be used to rewrite a block of a matrix that is given as a
// vector with another matrix that is also given as a vector.
std::vector<long double> WriteBlock(std::vector<long double> matrix,
                                    const std::vector<long double> block, int i,
                                    int j, int numVariables) {
  int dimension = sqrt(matrix.size());
  for (int k = 0; k < dimension / numVariables; k++) {
    for (int l = 0; l < dimension / numVariables; l++) {
      matrix[(k + i * dimension / numVariables) * dimension +
             (l + j * dimension / numVariables)] =
          block[k * dimension / numVariables + l];
    }
  }
  return matrix;
}

long double SumValues(const std::vector<long double> &vector) {
  long double result = 0;
  for (int i = 0; i < vector.size(); i++) {
    result += vector[i];
  }
  return result;
}

std::vector<std::vector<long double>>
SliceVectorVector(const std::vector<std::vector<long double>> &inputVector,
                  int start, int end) {
  std::vector<std::vector<long double>> result;
  for (int i = start; i < end && i < inputVector.size(); i++) {
    result.push_back(inputVector[i]);
  }
  return result;
}
