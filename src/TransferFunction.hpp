/* interface for a class that builds FunC transfer function pairs
   g and g^{-1}. Used by the NonUniformLookupTables to map a 
   uniform grid in [0,1] to a non-uniform grid in [0,1]. 
   g must satisfy the following conditions
   g(0) = 0
   g(1) = 1
   x <= y implies g(x) <= g(y)

  Notes:
    - g^{-1} must be quick to evaluate
    - g should equally distribute an interpolants error throughout [a,b] for any given f
*/
#pragma once
#include <functional> // std::function
#include <utility> // std::pair

class TransferFunction
{
  protected:
    /* --- Member variables --- */
    double m_minArg, m_maxArg;
    std::pair<std::function<double(double)>,std::function<double(double)>> mp_g_and_g_inv;
    std::function<double(double)> m_base_function; // might want autodiff functionality here?

  public:
    // build the function pair
    TransferFunction(std::function<double(double)> f, double minArg, double maxArg) :
      m_minArg(minArg), m_maxArg(maxArg), m_base_function(f) {}
    virtual ~TransferFunction(){}

    // public access to private vars
    std::pair<std::function<double(double)>,std::function<double(double)>> function_pair(){ return mp_g_and_g_inv; }
    std::pair<double,double> arg_bounds_of_interval(){ return std::make_pair(m_minArg, m_maxArg); }
    std::function<double(double)> function(){ return m_base_function; }
};
