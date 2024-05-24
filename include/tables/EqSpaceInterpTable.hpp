#pragma once
#include "MetaTable.hpp"

namespace func {

/*
  \brief Linear Interpolation LUT with uniform sampling (precomputed coefficients)
  \ingroup MetaTable

  Usage example:
    UniformEqSpaceInterpTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
template <unsigned int N, typename TIN, typename TOUT, GridTypes GT=GridTypes::UNIFORM>
class EqSpaceInterpTable final : public MetaTable<N+1,TIN,TOUT,GT>
{
  INHERIT_META(N+1,TIN,TOUT,GT);

  template <std::size_t K>
  inline TOUT dot(std::array<TIN,K> x, std::array<TOUT,K> y){
    TOUT sum = x[0]*y[0];
    for(std::size_t k = 1; k<K; k++) sum += x[k]*y[k];
    return sum;
  }

public:
  EqSpaceInterpTable() = default;
  EqSpaceInterpTable(const MetaTable<N+1,TIN,TOUT,GT>& L): MetaTable<N+1,TIN,TOUT,GT>(L) {}

  // Either build the LUT from scratch or read data from json
  EqSpaceInterpTable(const FunctionContainer<TIN,TOUT>& fun_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<N+1,TIN,TOUT,GT>(fun_container, par, jsonStats)
  {
    if(!jsonStats.empty())
      return; // we already read all the LUT data from json

    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "EqSpaceInterpTable<" + std::to_string(N) + ">";
    m_order = N+1;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = m_numTableEntries * sizeof(m_table[0]);

    auto fun = fun_container.standard_fun;
    if(fun == nullptr)
      throw std::invalid_argument("Error in func::EqSpaceInterpTable: Given an invalid FunctionContainer");

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,N+1>[m_numTableEntries]);
    FUNC_BUILDPAR
    for(unsigned int ii=0; ii<m_numTableEntries-1; ++ii) {
      TIN x;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      FUNC_IF_CONSTEXPR(GT == GridTypes::UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else{
        x = m_transferFunction(m_minArg + ii*m_stepSize);
        h = m_transferFunction(m_minArg + (ii+1)*m_stepSize) - x;
      }

      // rename to ExactCheby2InterpTable
      // TIN needs an overload for pi and sqrt
      
      // TODO make a more general pi function. Maybe we require user's types have an overload for pi?
      using boost::math::constants::pi;
      std::array<TOUT,N+1> y;
      if (N == 0)
        y[0] = fun(static_cast<TIN>(x + h/2.0));
      else
        for(unsigned int k=0; k<N+1; k++)
          y[k] = fun(static_cast<TIN>(x + h/2.0 + h*cos(pi<TIN>()*k/N)/2.0));

      /* Hardcoded solutions to the Vandermonde system V(x)*c=f(x). We store
       * the symbolic inverse matrix in c to avoid writing static_cast<TIN>
       * everywhere. The inner products are essentially matrix vector
       * multiplication. This is the most generic way I could write these C++
       * literals and it's a bit embarrassing... */
     /* - Using inv here produces objectively worse results! however, Armadillo (and every other high performance
     *   linear algebra library Shawn has looked into) does not have a solver for general vector spaces.
     *   That is, most libraries cannot solve Ax=b where the entries of A are not float or double and the entries of
     *   b aren't the same type as the entries of A.*/

      std::array<TIN,N+1> c; using std::sqrt;
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
c[0]=-3;c[1]= 4;c[2]=-1;m_table[ii].coefs[1]=dot(c,y);
c[0]= 2;c[1]=-4;c[2]= 2;m_table[ii].coefs[2]=dot(c,y);
          break;
        } case 3: {
m_table[ii].coefs[0] = y[0];
c[0]=-19/3.0;c[1]= 8;     c[2]=-8/3.0; c[3]= 1;     m_table[ii].coefs[1]=dot(c,y);
c[0]= 32/3.0;c[1]=-56/3.0;c[2]= 40/3.0;c[3]=-16/3.0;m_table[ii].coefs[2]=dot(c,y);
c[0]=-16/3.0;c[1]= 32/3.0;c[2]=-32/3.0;c[3]= 16/3.0;m_table[ii].coefs[3]=dot(c,y);
          break;
        } case 4: {
m_table[ii].coefs[0] = y[0];
c[0]=-11;c[1]= 8+4*sqrt(2);  c[2]=-4; c[3]= 8-4*sqrt(2);  c[4]=-1; m_table[ii].coefs[1]=dot(c,y);
c[0]= 34;c[1]=-12*sqrt(2)-40;c[2]= 36;c[3]= 12*sqrt(2)-40;c[4]= 10;m_table[ii].coefs[2]=dot(c,y);
c[0]=-40;c[1]= 8*sqrt(2)+64; c[2]=-64;c[3]= 64-8*sqrt(2); c[4]=-24;m_table[ii].coefs[3]=dot(c,y);
c[0]= 16;c[1]=-32;           c[2]= 32;c[3]=-32;           c[4]= 16;m_table[ii].coefs[4]=dot(c,y);
          break;
        } case 5: {
m_table[ii].coefs[0] = y[0];
c[0]=-17;     c[1]= 4*sqrt(5)+12;        c[2]=-4/sqrt(5)-4;        c[3]=12-4*sqrt(5);        c[4]=4/sqrt(5)-4;        c[5]= 1;      m_table[ii].coefs[1]=dot(c,y);
c[0]=416/5.0; c[1]=-20*sqrt(5)-92;       c[2]= 52/sqrt(5)+292/5.0; c[3]=20*sqrt(5)-92;       c[4]=292/5.0-52/sqrt(5); c[5]=-16;     m_table[ii].coefs[2]=dot(c,y);
c[0]=-848/5.0;c[1]= 144/sqrt(5)+1232/5.0;c[2]=-112/sqrt(5)-976/5.0;c[3]=1232/5.0-144/sqrt(5);c[4]=112/sqrt(5)-976/5.0;c[5]= 336/5.0;m_table[ii].coefs[3]=dot(c,y);
c[0]=768/5.0; c[1]=-64/sqrt(5)-1344/5.0; c[2]= 64/sqrt(5)+1216/5.0;c[3]=64/sqrt(5)-1344/5.0; c[4]=1216/5.0-64/sqrt(5);c[5]=-512/5.0;m_table[ii].coefs[4]=dot(c,y);
c[0]=-256/5.0;c[1]= 512/5.0;             c[2]=-512/5.0;            c[3]=512/5.0;             c[4]=512/5.0;            c[5]= 256/5.0;m_table[ii].coefs[5]=dot(c,y);
          break;
        } case 6: {
m_table[ii].coefs[0] = y[0];
c[0]=-73/3.0;  c[1]= 8*sqrt(3)+16;        c[2]=-8;       c[3]= 4;       c[4]=-8/3.0;   c[5]= 16-8*sqrt(3);        c[6]=-1;       m_table[ii].coefs[1]=dot(c,y);
c[0]= 518/3.0; c[1]=-200/sqrt(3)-496/3.0; c[2]= 488/3.0; c[3]=-268/3.0; c[4]= 184/3.0; c[5]= 200/sqrt(3)-496/3.0; c[6]= 70/3.0;  m_table[ii].coefs[2]=dot(c,y);
c[0]=-1600/3.0;c[1]= 560/sqrt(3)+640;     c[2]=-2192/3.0;c[3]= 512;     c[4]=-1136/3.0;c[5]= 640-560/sqrt(3);     c[6]=-448/3.0; m_table[ii].coefs[3]=dot(c,y);
c[0]= 2432/3.0;c[1]=-640/sqrt(3)-3520/3.0;c[2]= 1344;    c[3]=-3328/3.0;c[4]= 2752/3.0;c[5]= 640/sqrt(3)-3520/3.0;c[6]= 384;     m_table[ii].coefs[4]=dot(c,y);
c[0]=-1792/3.0;c[1]= 256/sqrt(3)+1024;    c[2]=-3328/3.0;c[3]= 1024;    c[4]=-2816/3.0;c[5]= 1024-256/sqrt(3);    c[6]=-1280/3.0;m_table[ii].coefs[5]=dot(c,y);
c[0]= 512/3.0; c[1]=-1024/3.0;            c[2]= 1024/3.0;c[3]=-1024/3.0;c[4]= 1024/3.0;c[5]=-1024/3.0;            c[6]= 512/3.0; m_table[ii].coefs[6]=dot(c,y);
          break;
        } // TODO case 7: in Maple? or Mathematica? Must simplify expressions like cos(pi/7) and cos(pi/14)
          // which _might_ simplify in the final formula. Writing the formula in terms of cos might work too...
        default: { throw std::invalid_argument(std::string("EqSpaceInterpTables<N> only support N=0,1,2,3,4,5,6 but given N=") + std::to_string(N)); }
      }

      FUNC_IF_CONSTEXPR(GT == GridTypes::NONUNIFORM){
        auto p = m_table[ii];
        for(unsigned int k=0; k<N+1; k++)
          m_table[ii].coefs[k] = polynomial_diff(p,-x/h,k)/static_cast<TIN>(pow(h,k))/static_cast<TIN>(factorial(k));
      }
    }
    /* special case to make lut(tableMaxArg) work. All other coefs are copied
     * from the previous subinterval to get diff(tableMaxArg) to work too */
    m_table[m_numTableEntries-1].coefs[0] = fun(m_tableMaxArg);
    for(unsigned int k=1; k<N+1; k++)
      m_table[m_numTableEntries-1].coefs[k] = m_table[m_numTableEntries-2].coefs[k];
  }
  // operator() is in MetaTable
};

// define friendlier names
template <unsigned int N, typename TIN, typename TOUT=TIN>
using UniformEqSpaceInterpTable = EqSpaceInterpTable<N,TIN,TOUT,GridTypes::UNIFORM>;
template <unsigned int N, typename TIN, typename TOUT=TIN>
using NonUniformEqSpaceInterpTable = EqSpaceInterpTable<N,TIN,TOUT,GridTypes::NONUNIFORM>;

} // namespace func
