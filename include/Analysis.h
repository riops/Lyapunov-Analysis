#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <vector>
#include <string>

long double CalculateAverage(const std::vector<long double>& numbers);
long double CalculateStandardDeviationError(const std::vector<long double>& values);
std::vector<long double> Sort(const std::vector<long double>& numbers);
std::vector<std::vector<long double>> SeparateAsMinStd(const std::vector<long double>& numbers);
void CalculateTracedValues(const std::string& inputFileName, const std::string& outputFileName, int N);
#endif // ANALYSIS_H
