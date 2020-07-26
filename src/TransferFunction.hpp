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

class TransferFunction
{
  protected:
    /* --- Member variables --- */
    double m_minArg, m_maxArg;
    std::pair<std::function<double(double)>,std::function<double(double)>> m_function_pair;
    std::function<double(double)> m_base_function;

    /* simplest example of a usable g, but since it doesn't distribute points
     * based on where f changes, it makes nonuniform lookup tables into 
     * expensive uniform lookup tables */
    std::function<double(double)> identity(){ return [](double x){ return x; }; }

  public:
    // build the function pair
    TransferFunction(std::function<double(double)> f, double minArg, double maxArg) :
      m_minArg(minArg), m_maxArg(maxArg), m_base_function(f) {}
    virtual ~TransferFunction(){}

    virtual void print_details(std::ostream& out){};

    // the main functionality of the class.
    std::function<double(double)> g;
    std::function<double(double)> g_inv;

    // public access to private vars
    std::pair<std::function<double(double)>,std::function<double(double)>> function_pair(){ return m_function_pair; }
    std::pair<double,double> arg_bounds_of_interval(){ return std::make_pair(m_minArg, m_maxArg); }
    std::function<double(double)> function(){ return m_base_function; }
};
