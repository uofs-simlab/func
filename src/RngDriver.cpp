#include <iostream>

#include "RngInterface.hpp"
#include "StdRng.hpp"

// a dummy version of ImplementationComparator. Just samples points based on the 
// type of RngInterface it's given
class Driver{
  RngInterface<double> *mp_sampler;
  public:
  Driver(RngInterface<double> *inRng=nullptr, unsigned int seed=2020) {
    mp_sampler=inRng;
    if(mp_sampler==nullptr)
      mp_sampler=std::move(new StdRng<std::uniform_real_distribution<double>>
        (new std::uniform_real_distribution<double>(0.0, 1.0)));
    mp_sampler->init(seed);
  }

  double samplePoints(){return mp_sampler->getPt();}

  ~Driver(){}
};


int main(){
  // sample from a normal distribution
  StdRng<std::normal_distribution<double>> standardRNG(new std::normal_distribution<double>(-5.0, 1.0));
  Driver driver1(&standardRNG);
  std::cout << driver1.samplePoints() << std::endl;

  // sample from the default uniform distribution
  Driver driver2;
  std::cout << driver2.samplePoints() << std::endl;
}

