/*
   Attempts to factor out common differences between table types
   (eg hash type, grid type, setup/reading polynomial coefficients)
   which can be pieced together based on template parameters.

   N = number of coefficients used in underlying piecewise polynomials
   Provided Horner's method which is the most common table evaluation method in FunC

   UNIFORM: Distance between each subinterval is always the same so we
   can use a faster hash in exchange for (probably) higher error.
   NONUNIFORM: Use a transfer function to create a nonuniform grid
   with a quicker hash.
   NONUNIFORM_PSEUDO: same as NONUNIFORM but uses a faster, less accurate
   hash.

   This template style approach greatly reduces code redundancy
*/
#pragma once
#include "LookupTable.hpp"
#include <array>
#include <stdexcept>

#define INHERIT_META(TIN,TOUT,N,GT) \
  using MetaTable<TIN,TOUT,N,GT>::m_table; \
  using MetaTable<TIN,TOUT,N,GT>::m_transferFunction

/* Parallelization macro.
 * Play around with this to see which OpenMP option for parallelizing a for loop is best
 * Might be nice to have simd table generation so something like the LUT generator can use actual
 * parallelism. We also know the alignment of m_table so that might give some speedup. */
//_Pragma("omp simd aligned(m_table:sizeof(TOUT))")
//#pragma omp simd aligned(m_table:sizeof(TOUT)) // needs the constructor to be declared simd
// assuming each iteration will take about the same amount of time
//#pragma omp parallel for schedule(static)

//namespace func {
// TODO do we actually need/want to assign numbers to this enum? Use an enum class instead?
enum GridTypes {UNIFORM = 0, NONUNIFORM = 1, NONUNIFORM_PSEUDO = 2};

template <GridTypes GT>
std::string grid_type_to_string() {
  switch(GT){
    case UNIFORM:
      return "Uniform";
    case NONUNIFORM:
      return "NonUniform";
    case NONUNIFORM_PSEUDO:
      return "NonUniformPseudo";
    default: { throw std::invalid_argument("Broken switch case in func::MetaTable"); }
  } 
}

template <typename TIN, typename TOUT, unsigned int N, GridTypes GT=UNIFORM>
class MetaTable : public LookupTable<TIN,TOUT>
{
protected:
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);

  __attribute__((aligned)) std::unique_ptr<polynomial<TOUT,N>[]> m_table;
  TransferFunctionSinh<TIN> m_transferFunction; // used to make nonuniform grids

  TOUT get_table_entry(unsigned int i, unsigned int j) override { return m_table[i].coefs[j]; }
  unsigned int get_num_coefs() override { return N; }
  std::array<TIN,4> get_transfer_function_coefs() override { return m_transferFunction.get_coefs(); }

public:
  // std::unique_ptr m_array implicitly deletes the copy ctor so we have to explicitly
  // ask for the default copy ctor
  MetaTable(MetaTable&&) = default;

  MetaTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par) :
    LookupTable<TIN,TOUT>(func_container, par),
    m_transferFunction(TransferFunctionSinh<TIN>(m_minArg,m_tableMaxArg,m_stepSize))
  {
    // initialize the transfer function to something useful
    if(GT != UNIFORM)
      m_transferFunction = TransferFunctionSinh<TIN>(func_container,m_minArg,m_minArg + m_numIntervals*m_stepSize,m_stepSize);
    //m_transferFunction.print_details(std::cout);
  }

  /* build this table from a file. Everything other than m_table is built by LookupTable */
  MetaTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename, std::string tablename) :
    LookupTable<TIN,TOUT>(func_container, filename),
    m_transferFunction(TransferFunctionSinh<TIN>(m_minArg,m_tableMaxArg,m_stepSize))
  {
    nlohmann::json jsonStats;
    std::ifstream(filename) >> jsonStats;

    // check that the names match
    m_name = jsonStats["name"].get<std::string>();
    if(m_name != tablename)
      throw std::invalid_argument("Error while building " + tablename + " : " + filename +
          " contains data for building a " + m_name + " which is not compatible");

    m_table.reset(new polynomial<TOUT,N>[m_numTableEntries]);
    for(unsigned int i=0; i<m_numTableEntries; i++)
      for(unsigned int j=0; j<m_table[i].num_coefs; j++)
        m_table[i].coefs[j] = jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)].get<TOUT>();

    // rebuild the transfer function
    std::array<TIN,4> inv_coefs = jsonStats["transfer_function_coefs"].get<std::array<TIN,4>>();
    m_transferFunction = TransferFunctionSinh<TIN>(m_minArg,m_tableMaxArg,m_stepSize,inv_coefs);
  }

  /* Provide the most common hash. The compiler should simplify this method when templates are instantiated */
  TOUT operator()(TIN x) override
  {
    TOUT dx;
    unsigned int x0;
    switch(GT){
      case UNIFORM:
        {
        // nondimensionalized x position, scaled by step size
        dx  = (TOUT) m_stepSize_inv*(x-m_minArg);
        // index of previous table entry
        x0  = (unsigned) dx;
        // value of table entries around x position
        dx -= x0;
        break;
        }
      case NONUNIFORM:
        {
        // find the subinterval x lives in
        x0 = m_transferFunction.g_inv(x);
        // find where x is within that interval
        TIN h = m_grid[x0+1] - m_grid[x0];
        dx    = (x - m_grid[x0])/h;
        break;
        }
      case NONUNIFORM_PSEUDO:
        {
        // find the subinterval x lives in
        dx  = m_transferFunction.g_inv(x);
        // just take the fractional part of dx as x's location in this interval
        x0  = (unsigned) dx;
        dx -= x0;
        break;
        }
    }
    
    // general degree horners method, evaluated from the inside out.
    TOUT sum = 0;
    for(int k=N-1; k>0; k--)
      sum = dx*(m_table[x0].coefs[k] + sum);
    return m_table[x0].coefs[0]+sum;
  }
};

/* TODO
void to_json(json& j, const person& p) {
  j = json{{"name", p.name}, {"address", p.address}, {"age", p.age}};
}

void from_json(const json& j, person& p) {
  j.at("name").get_to(p.name);
  j.at("address").get_to(p.address);
  j.at("age").get_to(p.age);
}

} // namespace func */
