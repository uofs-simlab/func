#include <iostream>
#include <string>
#include <cmath>
#include "func.hpp"

template <typename T>
T f(T x){
  return std::cos(x);
}

template <typename T>
T f1(T x){
  return -std::sin(x);
}

template <typename T>
T f2(T x){
  return -std::cos(x);
}

template <typename T>
T f3(T x){
  return std::sin(x);
}

#define MIN -1
#define MAX 1

#define TYPE double
#define STEP 0.05
#define N 1000.0
using namespace func;

int main(){
  auto LUT = func::UniformEqSpaceInterpTable<3,TYPE>({f<TYPE>},{MIN,MAX,STEP});

  //for(unsigned int i=0; i<=N; i++){
  //  auto x = MIN + (MAX-MIN)*i/N;
  //  std::cout << x << " " << f1(x) << " " << LUT.diff(1,x) << std::endl;
  //}

  //for(unsigned int i=0; i<=N; i++){
  //  auto x = MIN + (MAX-MIN)*i/N;
  //  std::cout << x << " " << f2(x) << " " << LUT.diff(2,x) << std::endl;
  //}

  for(unsigned int i=0; i<=N; i++){
    auto x = MIN + (MAX-MIN)*i/N;
    std::cout << x << " " << f3(x) << " " << LUT.diff(3,x) << std::endl;
  }
  
  //for(unsigned int i=0; i<=N; i++){
  //  auto x = MIN + (MAX-MIN)*i/N;
  //  std::cout << x << " " << f(x) << " " << LUT.diff(4,x) << std::endl;
  //}
}
