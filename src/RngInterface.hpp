#pragma once

/* Abstract interface for interfacing with random distributions */
template <typename OutType=double>
class RngInterface{
  public:
    // build a random generator for sampling from a distribution
    virtual void init(unsigned int seed)=0;

    // return the current seed
    virtual unsigned int seed()=0;

    // get a random point from the random distribution
    virtual OutType getPt()=0;

    virtual ~RngInterface(){}
};

