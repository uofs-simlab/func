/* 
  An implementation of RngInterface for sampling from 
  the distributions defined in std::random.

  Usage example:
    auto rng = new StdRng<std::uniform_real_distribution>(new std::uniform_real_distribution<double>(0.0,1.0));
    rng.init(2020);
    rng.getPt();
 */
#pragma once
#include "RngInterface.hpp"
#include <iostream>
#include <random>

template <class DistType>
class StdRng : public RngInterface<> 
{
  std::unique_ptr<DistType> mp_distribution;
  unsigned int m_seed;
  std::unique_ptr<std::mt19937> mp_generator;

  public:
    // create a new StdRng based on the pointer to random distribution passed in
    StdRng(DistType *dist) : mp_distribution(std::move(dist)){}

    // set the seed and generator of the distribution
    void init(unsigned int seed)
    {
      m_seed = seed;
      mp_generator.reset(new std::mt19937(seed));
    }

    // return the seed
    unsigned int seed(){ return m_seed; }

    // get a random point from the given distribution
    double getPt(){ return (*mp_distribution)(*mp_generator); }

    ~StdRng(){}
};

