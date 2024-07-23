#include <iostream>
#include <string>
#include <cmath>
#include "func.hpp"

template <typename T>
T f(T x, T y){
  return x*x*x*x*y*y*y;
}

template <typename T>
T f10(T x, T y){
  return 4*x*x*x*y*y*y;
}

template <typename T>
T f11(T x, T y){
  return 12*x*x*x*y*y;
}



#define MIN -1
#define MAX 1

#define TYPE double
#define STEP 0.25
#define N 50.0
using namespace func;

template <class T1, class T2>
using MyLUT = func::UniformExactInterpTable<3,T1,T2>;

int main(){
  std::vector<func::LookupTableParameters<TYPE>> params = {{MIN,MAX,STEP}, {MIN,MAX,STEP}};
  auto LUT = func::ndimLUT<2,TYPE,TYPE,MyLUT>(f<TYPE>, params);

  // heatmap of a LUT of f
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<N; j++){
      auto x = MIN + (MAX - MIN)*i/N;
      auto y = MIN + (MAX - MIN)*j/N;
      //std::cout << std::abs((LUT.diff(1,x))(y) - f10(x,y)) << " ";
      //std::cout << std::abs(LUT.diff(1,x,1,y) - f11(x,y)) << " ";
      //std::cout << LUT.diff(1,x,1,y) << " ";
      //std::cout << LUT.diff(1,x,0,y) << " ";
      //std::cout << LUT(x).diff(1,y) << " ";
      //std::cout << LUT(x).diff(1,y) << " ";
      std::cout << LUT.diff(1,x,2,y) << " ";
      //std::cout << std::abs(LUT(x,y) - f(x,y)) << " ";
    }
    std::cout << "\n";
  }
}
