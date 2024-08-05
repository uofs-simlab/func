#pragma once
namespace func {

/**
  Abstract interface for interfacing with classes
  that can generate random numbers
*/
template <typename POINT_TYPE>
class RngInterface {
  public:
    // build a random generator for sampling from a distribution
    virtual void init(unsigned int seed)=0;

    // return the current seed
    virtual unsigned int seed()=0;

    // get a random point from the random distribution
    virtual POINT_TYPE get_point()=0;

    virtual ~RngInterface(){}
};
}
