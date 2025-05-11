//
// Created by berko on 3/20/2023.
//
// #include <iostream>
#include "../include/Analysis.h"
#include "../include/Systems.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Function to calculate the average of a vector of numbers
long double CalculateAverage(const std::vector<long double> &numbers) {
  long double sum = 0.0;
  for (const auto &num : numbers) {
    sum += num;
  }
  return static_cast<double>(sum) / numbers.size();
}

// Function to calculate the standard deviation of a vector of numbers
long double
CalculateStandardDeviationError(const std::vector<long double> &values) {
  // Step 1: Calculate the mean
  long double sum = 0.0;
  for (const auto &value : values) {
    sum += value;
  }
  long double mean = sum / values.size();

  // Step 2: Calculate the sum of squared differences
  long double squaredDifferencesSum = 0.0;
  for (const auto &value : values) {
    squaredDifferencesSum += (value - mean) * (value - mean);
  }

  // Step 3: Divide the sum of squared differences by the number of elements
  long double variance = squaredDifferencesSum / values.size();

  // Step 4: Take the square root of the variance to obtain the standard
  // deviation
  long double standardDeviationError = std::sqrt(variance);

  return standardDeviationError;
}

std::vector<long double> Sort(const std::vector<long double> &numbers) {
  std::vector<long double> sortedNumbers;
  sortedNumbers.resize(numbers.size());
  for (int i = 0; i < numbers.size(); i++) {
    sortedNumbers[i] = numbers[i];
  }
  for (int i = 0; i < sortedNumbers.size(); i++) {
    for (int j = i + 1; j < sortedNumbers.size(); j++) {
      if (sortedNumbers[i] > sortedNumbers[j]) {
        long double temp = sortedNumbers[i];
        sortedNumbers[i] = sortedNumbers[j];
        sortedNumbers[j] = temp;
      }
    }
  }
  return sortedNumbers;
}

std::vector<std::vector<long double>>
SeparateAsMinStd(const std::vector<long double> &numbers) {
  std::vector<long double> sortedNumbers = Sort(numbers);
  int numberSize = sortedNumbers.size();

  std::vector<long double> stdErrors;
  stdErrors.resize(numberSize - 1);
  std::vector<long double> stdErrors1;
  stdErrors1.resize(numberSize - 1);
  std::vector<long double> stdErrors2;
  stdErrors2.resize(numberSize - 1);
  std::vector<long double> averages1;
  averages1.resize(numberSize - 1);
  std::vector<long double> averages2;
  averages2.resize(numberSize - 1);

  for (int separator = 1; separator < numberSize; separator++) {
    int firstGroupSize = separator;
    int secondGroupSize = numberSize - separator;

    std::vector<long double> group1;
    group1.resize(firstGroupSize);
    std::vector<long double> group2;
    group2.resize(secondGroupSize);

    for (int i = 0; i < firstGroupSize; i++) {
      group1[i] = sortedNumbers[i];
    }
    for (int i = 0; i < secondGroupSize; i++) {
      group2[i] = sortedNumbers[i + firstGroupSize];
    }
    long double stdError1 = CalculateStandardDeviationError(group1);
    stdErrors1[separator - 1] = stdError1;
    long double stdError2 = CalculateStandardDeviationError(group2);
    stdErrors2[separator - 1] = stdError2;
    long double average1 = CalculateAverage(group1);
    averages1[separator - 1] = average1;
    long double average2 = CalculateAverage(group2);
    averages2[separator - 1] = average2;
    stdErrors[separator - 1] = stdError1 + stdError2;
  }

  long double currentMin = stdErrors[0];
  int separationPoint = 1;
  for (int i = 0; i < numberSize - 1; i++) {
    if (stdErrors[i] < currentMin) {
      currentMin = stdErrors[i];
      separationPoint = i + 1;
    }
  }

  std::vector<long double> group1;
  group1.resize(separationPoint);
  std::vector<long double> group2;
  group2.resize(numberSize - separationPoint);

  for (int i = 0; i < separationPoint; i++) {
    group1[i] = sortedNumbers[i];
  }
  for (int i = 0; i < numberSize - separationPoint; i++) {
    group2[i] = sortedNumbers[i + separationPoint];
  }

  long double realAverage1 = CalculateAverage(group1);
  long double realAverage2 = CalculateAverage(group2);

  long double realStdError1 = CalculateStandardDeviationError(group1);
  long double realStdError2 = CalculateStandardDeviationError(group2);

  return {{realAverage1, realStdError1, realAverage2, realStdError2},
          group1,
          group2};
}

void CalculateTracedValues(const std::string &inputFileName,
                           const std::string &outputFileName, int N) {
  // Open the input CSV file
  std::cout << inputFileName << "\n";
  std::cout << outputFileName << "\n";
  std::ifstream inputFile("./data/csv/" + inputFileName);
  if (!inputFile.is_open()) {
    std::cerr << "Error opening input file" << std::endl;
    return;
  }

  // Open the output CSV file
  std::ofstream outputFile("./data/csv/" + outputFileName);
  if (!outputFile.is_open()) {
    std::cerr << "Error opening output file" << std::endl;
    return;
  }

  int Imax = 4 * N * N;
  int CovMax = Imax * Imax;

  std::string lineStr;
  std::vector<std::string> lines;
  // Read all lines into a vector
  while (std::getline(inputFile, lineStr)) {
    lines.push_back(lineStr);
  }

  // Process all lines except the last one
  // std::cout << lines.size() << '\n';
  for (size_t lineNum = 0; lineNum < lines.size(); ++lineNum) {
    lineStr = lines[lineNum];

    // Check if this is the last line
    if (lineNum == lines.size() - 1) {
      // Copy the last line as is
      outputFile << lineStr << "\n";
      continue;
    }

    // Parse the line into a vector of long doubles
    std::vector<long double> lineValues;
    std::stringstream ss(lineStr);
    std::string valueStr;
    while (std::getline(ss, valueStr, ',')) {
      try {
        lineValues.push_back(std::stold(valueStr));
      } catch (const std::invalid_argument &e) {
        std::cerr << "Invalid number in line " << lineNum + 1 << ": "
                  << valueStr << std::endl;
        lineValues.push_back(0.0); // Default value
      }
    }

    // Ensure that the line has enough elements
    if (lineValues.size() < static_cast<size_t>(3 * CovMax)) {
      std::cerr << "Line " << lineNum + 1 << " does not have enough elements."
                << std::endl;
      continue;
    }

    // Initialize trace variables
    long double XTrace = 0.0;
    // long double XTrace = 0.0;
    long double PTrace = 0.0;
    // long double PTrace = 0.0;

    // Perform the calculations
    for (int l1 = 0; l1 < N; l1++) {
      for (int m1 = -l1; m1 <= l1; m1++) {
        for (int a = 1; a <= 2; a++) {
          int idx_a1_a1_l1m1_l1m1 = indexXX(a, 1, l1, m1, a, 1, l1, m1, N);
          int idx_a2_a2_l1m1_l1m1 = indexXX(a, 2, l1, m1, a, 2, l1, m1, N);
          XTrace += lineValues[idx_a1_a1_l1m1_l1m1];
          XTrace += lineValues[idx_a2_a2_l1m1_l1m1];
          PTrace += lineValues[CovMax + idx_a1_a1_l1m1_l1m1];
          PTrace += lineValues[CovMax + idx_a2_a2_l1m1_l1m1];
        }
      }
    }

    // Write the results to the output file
    outputFile << XTrace / N << "," << PTrace / N << "\n";
  }

  // Close the files
  inputFile.close();
  outputFile.close();

  std::cout << "Processing complete. Results written to " << outputFileName
            << std::endl;
}
