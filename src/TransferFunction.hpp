/* 
   Interface for a class that builds and contains FunC transfer
   function pairs g and g^{-1}. These are the backbones of
   NonUniformLookupTables and are used to map a uniform grid
   in [a,b] to a non-uniform grid in [a,b].
   The new grid will ideally do a better job of distributing error
   when used for grid points.
   g must satisfy the following conditions:
   g(a) = a;
   g(b) = b;
   x <= y implies g(x) <= g(y) (ie, g must be monotone increasing).

  Notes:
    - g^{-1} must be quick to evaluate to see any speedup compared to uniform
    lookup tables
    - Since more grid points will exist where g changes the slowest we
    want g' to grow like 1/f' in order to distribute error.
*/
#pragma once
#include <utility> // std::pair
#include <memory> // std::unique_ptr

/* inheritance macro */
#define INHERIT_TRANSFER_FUNCTION(IN_TYPE) \
  using TransferFunction<IN_TYPE>::m_minArg; \
  using TransferFunction<IN_TYPE>::m_maxArg; \
  using TransferFunction<IN_TYPE>::mp_g; \
  using TransferFunction<IN_TYPE>::mp_g_inv

/* a lightweight functor used to speedup evaluations of g_inv 
   Unfortunate that we'll have to do so much pointer chasing
   here but this implementation is convenient and flexible */
template <typename IN_TYPE>
struct LightweightFunctor
{
  virtual IN_TYPE operator()(IN_TYPE x)=0;
  virtual ~LightweightFunctor(){};
};

template <typename IN_TYPE>
class TransferFunction
{
protected:
  /* This min must be the same as the corresponding table's min, and the max
   * must be equal to the tables max arg (not necessarily the actual max arg) */
  IN_TYPE m_minArg, m_maxArg;
  IN_TYPE m_stepSize;

  // the main functionality of the class
  std::unique_ptr<LightweightFunctor<IN_TYPE>> mp_g;
  std::unique_ptr<LightweightFunctor<IN_TYPE>> mp_g_inv;

public:
  // build the function pair. The FunctionContainer is unused but it's nice as an interface
  template <typename OUT_TYPE>
  TransferFunction(FunctionContainer<IN_TYPE,OUT_TYPE> *fc, IN_TYPE minArg, IN_TYPE maxArg, IN_TYPE stepSize) :
    m_minArg(minArg), m_maxArg(maxArg), m_stepSize(stepSize) {}
  virtual ~TransferFunction(){}

  virtual void print_details(std::ostream& out){};

  // public access to private vars
  std::pair<IN_TYPE,IN_TYPE> arg_bounds_of_interval(){ return std::make_pair(m_minArg, m_maxArg); }

  // getters used to avoid having to use pointer syntax
  IN_TYPE g(IN_TYPE x){ return (*mp_g)(x); }
  IN_TYPE g_inv(IN_TYPE x){ return (*mp_g_inv)(x); }
};
