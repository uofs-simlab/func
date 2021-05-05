/* 
  An implementation of RngInterface, intended to be used with the 
  generators/ distributions defined in std::random.
  Notes:
  - Will take ownership if given a probability distribution
  - Will build its own probability distribution corresponding to
  DIST_TYPE if given the correct constructor args.
  - Builds its own number generator when init(seed) is called

  Usage example:
    // uniform_real_distribution<double> generated from std::mt19937 
    // within the range 0.0 to 1.0
    StdRng<double> rng(0.0,1.0); // build a std::uniform_real_distribution<double>
    rng.init(2020);             // build a std::mt19937
    rng.get_point();            // return a random number

    // or if you want to get fancy and use a faster generator, here's
    // a normal distribution with mean 0.0 and standard deviation 1.0,
    // using points generated from minstd_rand0
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
  std::unique_ptr<RNG_TYPE>  mp_generator;
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
