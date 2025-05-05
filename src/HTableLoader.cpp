#include "../include/HTableLoader.h"
#include "../include/json.hpp"
#include <fstream>
#include <iostream>

using json = nlohmann::json;

std::vector<std::vector<double>> LoadHTable(int N) {
  const std::string filename = "data/json/HSymbol_nonzero.json";
  std::ifstream in{filename};
  if (!in)
    throw std::runtime_error("cannot open " + filename);

  json j;
  in >> j;

  auto key = std::to_string(N);
  if (!j.contains(key) || !j[key].is_array())
    throw std::runtime_error("no array for N=" + key);

  std::vector<std::vector<double>> table;
  for (auto &elem : j[key]) {
    table.emplace_back(std::vector<double>{
        elem.at("l1").get<double>(), elem.at("m1").get<double>(),
        elem.at("l2").get<double>(), elem.at("m2").get<double>(),
        elem.at("l3").get<double>(), elem.at("m3").get<double>(),
        elem.at("l4").get<double>(), elem.at("m4").get<double>(),
        elem.at("value").get<double>()});
  }
  return table;
}
