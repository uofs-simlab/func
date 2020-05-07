#pragma once

/* Abstract interface for interfacing with random distributions */
template <typename OutType=double>
class RngInterface{
  public:
    virtual void init(unsigned int seed)=0;
    virtual unsigned int seed()=0;
    virtual OutType getPt()=0;
    virtual ~RngInterface(){}
};

