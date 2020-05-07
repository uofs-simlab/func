#include <iostream>
#include <random>
#include "RngInterface.hpp"
#pragma once

// TODO implement smart pointers
// #include <memory>
template <class DistType>
class StdRng : public RngInterface<> {
  std::unique_ptr<DistType> mp_distribution=nullptr;
  unsigned int m_seed;
  std::unique_ptr<std::mt19937> mp_generator=nullptr;

  public:
    StdRng(DistType *dist) : mp_distribution(std::move(dist)) {}

    void init(unsigned int seed)
    {
      m_seed = seed;
      mp_generator.reset(new std::mt19937(seed));
    }

    unsigned int seed(){return m_seed;}

    double getPt(){ return (*mp_distribution)(*mp_generator); }
    
    ~StdRng(){}
};

