#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <vector>

void PrintMatrix(const std::vector<long double>& matrix);

void WriteMatrix(const std::vector<long double>& matrix);

void PrintVector(const std::vector<long double>& vector);

std::vector<long double> ZeroMatrix(int dimension);

std::vector<long double> IdentityMatrix(int dimension);

std::vector<long double> ExtendVector(const std::vector<long double>& vector1, const std::vector<long double>& vector2);

std::vector<long double> VectorAddition(const std::vector<long double>& vector1, const std::vector<long double>& vector2);

std::vector<long double> VectorSubtraction(const std::vector<long double>& vector1, const std::vector<long double>& vector2);

std::vector<long double> Transpose(const std::vector<long double>& matrix);

std::vector<long double> MatrixProduct(const std::vector<long double>& matrix1, const std::vector<long double>& matrix2);

std::vector<long double> ScalarVectorMultiplication(long double scalar, const std::vector<long double>& vector1);

long double InnerProduct(const std::vector<long double>& vector1, const std::vector<long double>& vector2);

long double Norm(const std::vector<long double>& vector1);

std::vector<long double> Projection(const std::vector<long double>& vector1, const std::vector<long double>& vector2);

std::vector<long double> Orthogonalize(const std::vector<long double>& vectors);

std::vector<long double> WriteBlock(std::vector<long double> matrix, const std::vector<long double>& block, int i, int j, int numVariables);

long double SumValues(const std::vector<long double>& vector);

std::vector<std::vector<long double>> SliceVectorVector(const std::vector<std::vector<long double>>& inputVector, int start, int end);


#endif // MATRIX_OPERATIONS_H
