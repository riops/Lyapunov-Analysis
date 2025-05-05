#include "json.hpp"
#include <fstream>
#include <iostream>

using json = nlohmann::json;

std::vector<std::vector<double>> LoadHTable(int N);
