/* 
  An implementation of RngInterface, intended to be used with the 
  generators/ distributions defined in std::random, but it should
  be compatible with every generator constructed with a single 
  unsigned int seed (eg. linear_congruential_engine,
  subtract_with_carry_engine, mersenne_twister_engine) and every
  probability distribution compatible with those generators.

  Usage example:
    // uniform_real_distribution<double> generated from std::mt19937 
    // within the range 0.0 to 1.0
    StdRng<double> rng(0.0,1.0);
    rng->init(2020);
    rng->get_point();

    // or if you want to get fancy, a normal distribution with mean 0.0
    // and standard deviation 1.0, using points generated from minstd_rand0
    StdRng<float,std::normal_distribution<float>,minstd_rand0> rng2(0.0,1.0));
 */
#pragma once
#include "RngInterface.hpp"
#include <iostream>
#include <random>

template <typename POINT_TYPE,
         class DIST_TYPE = std::uniform_real_distribution<POINT_TYPE>,
         class RNG_TYPE  = std::mt19937>
class StdRng : public RngInterface<POINT_TYPE>
{
  std::unique_ptr<DIST_TYPE> mp_distribution;
  std::unique_ptr<RNG_TYPE> mp_generator;
  unsigned int m_seed;

  public:
    // create a new StdRng based on the pointer to probability distribution passed in
    StdRng(DIST_TYPE *dist) : mp_distribution(std::move(dist)){}

    // create a new StdRng based on the pointer to probability distribution passed in
    template <typename ... DIST_TYPE_ARGS>
    StdRng(DIST_TYPE_ARGS ... args) :
      mp_distribution(std::unique_ptr<DIST_TYPE>(new DIST_TYPE(args ...))){}

    // set the seed and generator of the distribution
    void init(unsigned int seed)
    {
      m_seed = seed;
      mp_generator.reset(new RNG_TYPE(seed));
    }

    // return the seed
    unsigned int seed(){ return m_seed; }

    // get a random point from the given distribution
    POINT_TYPE get_point(){ return (*mp_distribution)(*mp_generator); }

    ~StdRng(){}
};
