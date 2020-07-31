#include "zipf.h"
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <map>
#include <vector>


void generate_build(const std::string &name, int build_rows) {
  std::ofstream table;
  table.open(name + "_build.txt");

  for (int row = 1; row <= build_rows; ++row) {
    table << row << "|" << row << std::endl;
  }

  table.close();
}

template <class R>
void generate_probe(const std::string &name, int build_rows,
                    int probe_build_ratio, R rng) {

}

int main(int argc, char **argv) {

  if (argc != 4) {
    std::cout << "Wrong Args" << std::endl;
  }

  int data_size = std::stoi(argv[1]);
  int data_accesses = std::stoi(argv[2]);
  float skew = std::stof(argv[3]);

  std::map<int, int> m;


  std::string name = std::to_string(data_size) + "_" + std::to_string(data_accesses) + "_" + std::to_string(skew) + ".csv";
  std::ofstream table;
  table.open(name);

  if (skew == 0.0) {
    std::uniform_int_distribution<int> uniform(1, data_size);
    std::random_device rd;
    std::mt19937 mt(rd());

    for (int row = 1; row <= data_accesses; ++row) {
      if (row % 1'000'000 == 0) { std::cout << row << std::endl;}
      int c = uniform(mt);
      table << c << "\n";
      if (m.find(c) != m.end())  { m[c] += 1; } else {m[c] =1 ;}
    }
  } else {
    ZipfDistribution zipf(data_size, skew);
    std::random_device rd;
    std::mt19937 mt(rd());

    for (int row = 1; row <= data_accesses; ++row) {
      if (row % 1'000'000 == 0) { std::cout << row << std::endl;}
      int c = zipf(mt);
      table << c << "\n";
      if (m.find(c) != m.end()) { m[c] += 1; } else {m[c] =1 ;}
    }
  }

  table.close();

  std::vector<int> vec;

  std::map<int, int> :: iterator it;
  for (it=m.begin(); it!=m.end(); it++)
  {
    vec.push_back(it->second);
  }

  sort(vec.begin(), vec.end());

  std::ofstream f;
  std::string fname = std::to_string(data_size) + "_" + std::to_string(data_accesses) + "_" + std::to_string(skew) + "_sorted_freq.csv";
  f.open(fname);

  for (int x : vec) {
    f << x << "\n";
  }
  f.close();


  return 0;
}
