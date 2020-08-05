/* 
   Interface for a class that builds and contains FunC transfer
   function pairs g and g^{-1}. Used by the NonUniformLookupTables
   to map a uniform grid in [0,1] to a non-uniform grid in [0,1].
   this new grid will ideally do a better job of distributing error
   when used for interpolation points.
   g must satisfy the following conditions
   g(0) = 0
   g(1) = 1
   x <= y implies g(x) <= g(y) (ie, g must be monotone increasing)

   Since more grid points will exist where g(x/(b-a)) changes the slowest,
   in order to better distribute error, we want
   g'(x/(b-a)) to be similar to 1/f' on [a,b]

  Notes:
    - g^{-1} must be quick to evaluate to see any speedup compared to uniform lookup tables
*/
#pragma once
#include <functional> // std::function
#include <utility> // std::pair

#define INHERIT_TRANSFER_FUNCTION(IN_TYPE) \
  using TransferFunction<IN_TYPE>::m_minArg; \
  using TransferFunction<IN_TYPE>::m_maxArg; \
  using TransferFunction<IN_TYPE>::g; \
  using TransferFunction<IN_TYPE>::g_inv

template <typename IN_TYPE>
class TransferFunction
{
  protected:
    IN_TYPE m_minArg, m_maxArg;

  public:
    // build the function pair
    TransferFunction(IN_TYPE minArg, IN_TYPE maxArg) : m_minArg(minArg), m_maxArg(maxArg){}
    virtual ~TransferFunction(){}

    virtual void print_details(std::ostream& out){};

    // the main functionality of the class. TODO profile with a class
    // more lightweight than std::function
    std::function<IN_TYPE(IN_TYPE)> g;
    std::function<IN_TYPE(IN_TYPE)> g_inv;

    // public access to private vars
    std::pair<IN_TYPE,IN_TYPE> arg_bounds_of_interval(){ return std::make_pair(m_minArg, m_maxArg); }
};
