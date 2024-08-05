#pragma once
#include "MetaTable.hpp"

/* nobody wants 'static_cast<TIN>(...)' all over the code. That is a
 * completely illegible way to write numeric literals. Programming math
 * in C++ is embarrassing */
#define S(x) static_cast<TIN>(x)

namespace func {

/**
  \brief Interpolation over Chebyshev nodes of the second kind. The
  inverse Vandermonde matrix is hard-coded. Allows for full type
  generality, but numerical output is not quite as good as
  ChebyInterpTable because Armadillo can do iterative refinement.
  \ingroup MetaTable

  \code{.cpp}
  // ExactInterpTable works with an untemplated function
  double foo(double x){ return x; }
 
  int main(){
    double min = 0.0, max = 10.0, step = 0.0001;
    UniformExactInterpTable<double>    L({foo}, {min, max, step}); // uniform partition
    NonUniformExactInterpTable<double> L({foo}, {min, max, step}); // nonuniform partition
    auto val = L(0.87354);
  }
  \endcode


  Notes:
  - this LUT precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
*/
template <unsigned int N, typename TIN, typename TOUT, GridTypes GT=GridTypes::UNIFORM>
class ExactInterpTable final : public MetaTable<N+1,TIN,TOUT,GT>
{
  INHERIT_META(N+1,TIN,TOUT,GT);

  //template <std::size_t K>
  //inline TOUT dot(std::array<TIN,K> x, std::array<TOUT,K> y){
  //  TOUT sum = x[0]*y[0];
  //  for(std::size_t k = 1; k<K; k++) sum += x[k]*y[k];
  //  return sum;
  //}

public:
  ExactInterpTable() = default;
  ExactInterpTable(const MetaTable<N+1,TIN,TOUT,GT>& L): MetaTable<N+1,TIN,TOUT,GT>(L) {}

  // Either build the LUT from scratch or read data from json
  ExactInterpTable(const FunctionContainer<TIN,TOUT>& fun_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<N+1,TIN,TOUT,GT>(fun_container, par, jsonStats)
  {
    if(!jsonStats.empty())
      return; // we already read all the LUT data from json

    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "ExactInterpTable<" + std::to_string(N) + ">";
    m_order = N+1;
    m_numTableEntries = m_numIntervals+1; // this is what it should be
    m_dataSize = m_numTableEntries * sizeof(m_table[0]);

    auto fun = fun_container.standard_fun;
    if(fun == nullptr)
      throw std::invalid_argument("Error in func::ExactInterpTable: Given an invalid FunctionContainer");

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,N+1>[m_numTableEntries]);
    FUNC_BUILDPAR
    for(unsigned int ii=0; ii<m_numTableEntries-1; ++ii) {
      TIN x = m_minArg + ii*m_stepSize;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      FUNC_IF_CONSTEXPR(GT != GridTypes::UNIFORM){
        x = m_transferFunction(x);
        h = m_transferFunction(m_minArg + (ii+1)*m_stepSize) - x;
      }

      // TIN needs an overload for pi and sqrt
      
      // TODO make a more general pi function for TIN
      using boost::math::constants::pi;
      std::array<TOUT,N+1> y;
      if (N == 0)
        y[0] = fun(static_cast<TIN>(x + h/2.0));
      else
        for(unsigned int k=0; k<N+1; k++)
          y[k] = fun(static_cast<TIN>(x + h/2.0 - h*cos(pi<TIN>()*k/N)/2.0));

      /* Hardcoded solutions to the Vandermonde system V(x)*c=f(x).
       * Using the macro S makes the code more legible. This is the most generic way I could write these C++
       * literals and it is embarrassing */
      /* Another option is to use Armadillo's `inv` here, but that produces
       * objectively worse results! Armadillo (and every other high performance
       * linear algebra library Shawn has looked into) does not have a solver for general vector spaces.
       * That is, most libraries cannot solve Ax=b where the entries of A are
       * elements of some field (not necessarily float or double) and the
       * entries of b aren't the same type as the entries of A.*/

      using std::sqrt;
      switch(N){
        case 0:
m_table[ii].coefs[0]=y[0];
          break;
        case 1: {
m_table[ii].coefs[0]=y[0];
m_table[ii].coefs[1]=y[1]-y[0];
          break;
        } case 2: {
m_table[ii].coefs[0] = y[0];
m_table[ii].coefs[1] = S(-3)*y[0]+S( 4)*y[1]+S(-1)*y[2];
m_table[ii].coefs[2] = S( 2)*y[0]+S(-4)*y[1]+S( 2)*y[2];
          break;
        } case 3: {
m_table[ii].coefs[0] = y[0];
m_table[ii].coefs[1] = S(-19)*y[0]/S(3)+S(  8)*y[1]     +S(- 8)*y[2]/S(3)+S(  1)*y[3];
m_table[ii].coefs[2] = S( 32)*y[0]/S(3)+S(-56)*y[1]/S(3)+S( 40)*y[2]/S(3)+S(-16)*y[3]/S(3);
m_table[ii].coefs[3] = S(-16)*y[0]/S(3)+S( 32)*y[1]/S(3)+S(-32)*y[2]/S(3)+S( 16)*y[3]/S(3);
          break;
        } case 4: {
TIN sq = sqrt(S(2));
m_table[ii].coefs[0]=y[0];
m_table[ii].coefs[1]=S(-11)*y[0]+(S(  8)+S(  4)*sq)*y[1]+S(- 4)*y[2]+(S(  8)+S(-4)*sq)*y[3]+S( -1)*y[4];
m_table[ii].coefs[2]=S( 34)*y[0]+(S(-40)+S(-12)*sq)*y[1]+S( 36)*y[2]+(S(-40)+S(12)*sq)*y[3]+S( 10)*y[4];
m_table[ii].coefs[3]=S(-40)*y[0]+(S( 64)+S(  8)*sq)*y[1]+S(-64)*y[2]+(S( 64)+S(-8)*sq)*y[3]+S(-24)*y[4];
m_table[ii].coefs[4]=S( 16)*y[0]+            S(-32)*y[1]+S( 32)*y[2]+           S(-32)*y[3]+S( 16)*y[4];
          break;
        } case 5: {
TIN sq = sqrt(S(5));
m_table[ii].coefs[0]=y[0];
m_table[ii].coefs[1]=S( -17)*y[0]     +(S(   12)+S(  4)*sq)*y[1]     +(S(  -4)+S(  -4)/sq)*y[2]     +(S(   12)+S(  -4)*sq)*y[3]     +(S(  -4)+S(   4)/sq)*y[4]     +S(   1)*y[5];
m_table[ii].coefs[2]=S( 416)*y[0]/S(5)+(S(  -92)+S(-20)*sq)*y[1]     +(S( 292)+S(  52)*sq)*y[2]/S(5)+(S(  -92)+S(  20)*sq)*y[3]     +(S( 292)+S( -52)*sq)*y[4]/S(5)+S( -16)*y[5];
m_table[ii].coefs[3]=S(-848)*y[0]/S(5)+(S( 1232)+S(144)*sq)*y[1]/S(5)+(S(-976)+S(-112)*sq)*y[2]/S(5)+(S( 1232)+S(-144)*sq)*y[3]/S(5)+(S(-976)+S( 112)*sq)*y[4]/S(5)+S( 336)*y[5]/S(5);
m_table[ii].coefs[4]=S( 768)*y[0]/S(5)+(S(-1344)+S(-64)*sq)*y[1]/S(5)+(S(1216)+S(  64)*sq)*y[2]/S(5)+(S(-1344)+S(  64)*sq)*y[3]/S(5)+(S(1216)+S( -64)*sq)*y[4]/S(5)+S(-512)*y[5]/S(5);
m_table[ii].coefs[5]=S(-256)*y[0]/S(5)+              S(512)*y[1]/S(5)+             S(-512)*y[2]/S(5)+               S(512)*y[3]/S(5)+             S(-512)*y[4]/S(5)+S( 256)*y[5]/S(5);
          break;
        } case 6: {
TIN sq = sqrt(S(3));
m_table[ii].coefs[0]=y[0];
m_table[ii].coefs[1]=S(  -73)*y[0]/S(3)+(S(   16)+S(   8)*sq)*y[1]     +S(   -8)*y[2]     +S(    4)*y[3]     +S(   -8)*y[4]/S(3)+(S(   16)+S(  -8)*sq)*y[5]     +S(   -1)*y[6];
m_table[ii].coefs[2]=S(  518)*y[0]/S(3)+(S( -496)+S(-200)*sq)*y[1]/S(3)+S(  488)*y[2]/S(3)+S( -268)*y[3]/S(3)+S(  184)*y[4]/S(3)+(S( -496)+S( 200)*sq)*y[5]/S(3)+S(   70)*y[6]/S(3);
m_table[ii].coefs[3]=S(-1600)*y[0]/S(3)+(S(  640)+S( 560)/sq)*y[1]     +S(-2192)*y[2]/S(3)+S(  512)*y[3]     +S(-1136)*y[4]/S(3)+(S(  640)+S(-560)/sq)*y[5]     +S( -448)*y[6]/S(3);
m_table[ii].coefs[4]=S( 2432)*y[0]/S(3)+(S(-3520)+S(-640)*sq)*y[1]/S(3)+S( 1344)*y[2]     +S(-3328)*y[3]/S(3)+S( 2752)*y[4]/S(3)+(S(-3520)+S( 640)*sq)*y[5]/S(3)+S(  384)*y[6];
m_table[ii].coefs[5]=S(-1792)*y[0]/S(3)+(S( 1024)+S( 256)/sq)*y[1]     +S(-3328)*y[2]/S(3)+S( 1024)*y[3]     +S(-2816)*y[4]/S(3)+(S( 1024)+S(-256)/sq)*y[5]     +S(-1280)*y[6]/S(3);
m_table[ii].coefs[6]=S(  512)*y[0]/S(3)+             S(-1024)*y[1]/S(3)+S( 1024)*y[2]/S(3)+S(-1024)*y[3]/S(3)+S( 1024)*y[4]/S(3)+             S(-1024)*y[5]/S(3)+S(  512)*y[6]/S(3);
          break;
        } // TODO case 7: in Maple? or Mathematica? Must involve expressions like cos(pi/7) and cos(pi/14)
          // which _might_ simplify in the final formula. Writing the formula in terms of cos might work too...
        default: { throw std::invalid_argument(std::string("ExactInterpTables<N> only support N=0,1,2,3,4,5,6 but given N=") + std::to_string(N)); }
      }

      FUNC_IF_CONSTEXPR(GT == GridTypes::NONUNIFORM)
        m_table[ii] = taylor_shift(m_table[ii], static_cast<TIN>(0), static_cast<TIN>(1), x, x+h);
    }
    /* special case to make lut(tableMaxArg) work. Move the second last polynomial into the last interval (shifting args) */
    FUNC_IF_CONSTEXPR(GT == GridTypes::UNIFORM)
      m_table[m_numTableEntries-1] = taylor_shift(m_table[m_numTableEntries-2], static_cast<TIN>(1), static_cast<TIN>(2), static_cast<TIN>(0), static_cast<TIN>(1));
    else
      m_table[m_numTableEntries-1] = m_table[m_numTableEntries-2];
  }

  // operator() is in MetaTable
};

// define friendlier names
template <unsigned int N, typename TIN, typename TOUT=TIN>
using UniformExactInterpTable = ExactInterpTable<N,TIN,TOUT,GridTypes::UNIFORM>;
template <unsigned int N, typename TIN, typename TOUT=TIN>
using NonUniformExactInterpTable = ExactInterpTable<N,TIN,TOUT,GridTypes::NONUNIFORM>;

} // namespace func
