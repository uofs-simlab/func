/*
  LUT using [M/N] pade approximants with uniform sampling. Polynomial coefficients are calculated using
  Armadillo. 

  Usage example using [4/3] approximants:
    UniformPadeTable<4,3> look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - Available template values are all M,N such that 0 < N <= M and M+N<=7
  - Template values where M < N are not supported
*/
#pragma once
#include "UniformLookupTable.hpp"
#include "config.hpp"

#ifndef FUNC_USE_BOOST_AUTODIFF
#error "UniformPadeTable needs boost version >= 1.71"
#endif

#ifndef FUNC_USE_ARMADILLO
#error "UniformPadeTable needs Armadillo"
#endif

#include <armadillo>
#include <iostream>
#include <stdexcept>
#include <cmath> //isinfite

static double const fact[] = {1,1,2,6,24,120,720,5040};

template <typename IN_TYPE, typename OUT_TYPE, unsigned int M, unsigned int N> 
class UniformPadeTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  FUNC_REGISTER_LUT(UniformPadeTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE, M+N+1>[]> m_table;
  std::function<adVar<OUT_TYPE,M+N>(adVar<OUT_TYPE,M+N>)> mp_boost_func;
  OUT_TYPE get_table_entry(unsigned int i, unsigned int j) override { return m_table[i].coefs[j]; }
  unsigned int get_num_coefs() override { return m_table[0].num_coefs; }

public:
  UniformPadeTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    using boost::math::differentiation::make_fvar;

    /* Base class default variables */
    m_name = "UniformPadeTable<" + std::to_string(M) + "," + std::to_string(N) + ">";
    m_order = M+N+1;  
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

    __IS_NULL(func_container->template get_nth_func<M+N>()); // the least descriptive exception
    mp_boost_func = func_container->template get_nth_func<M+N>();

    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,M+N+1>[m_numTableEntries]);
    for (int ii=0;ii<m_numIntervals;++ii) {
      const IN_TYPE x = m_minArg + ii*m_stepSize;
      // grid points
      m_grid[ii] = x;
      
      // build the matrix of taylor coefficients
      arma::mat T = arma::zeros(M+N+1, N+1);
      const auto derivs = (mp_boost_func)(make_fvar<IN_TYPE,M+N>(x));
      for(unsigned int i=0; i<M+N+1; i++)
        T(i,0) = derivs.derivative(i)/(fact[i]);

      // copy a shifted column 1 of T into the rest of the matrix
      for(unsigned int i=0; i<N+1; i++)
        T(arma::span(i,N+M), i) = T(arma::span(0,N+M-i), 0);

      // find the coefficients of Q.
      arma::mat Q = arma::null(T.rows(M+1, M+N));
      if(Q.n_rows > 1)
        throw std::range_error(m_name + " is too poorly conditioned");

      // scale Q such that its first entry equals 1.
      Q=Q/Q[0];
      for(unsigned int i=1; i<M+N+1; i++)
        if(!std::isfinite(Q[i]))
          throw std::range_error(m_name + " is too poorly conditioned");

      // find the coefficients of P
      arma::vec P = T.rows(0,M)*Q;

      /* Check if the denominator Q has any roots
         within the subinterval [-m_stepSize/2,m_stepSize/2).
         If any roots exist, then lower the degree of the denominator.
        
         We'll check for the existence of a root by building a bracket,
         using Q(0)=1 as our positive endpoint. Thus, we just need
         to find a point where Q is negative. Hence this helper function: */
      auto Q_is_negative = [this, &Q, &ii](IN_TYPE x) -> bool {
        // Tell us if this point is within this subinterval's range
        if(((ii == 0 && x < 0.0) || (ii == m_numIntervals - 1 && x > 0.0)))
          return false;

        // compute Q(x) using horners, evaluating from the inside out
        OUT_TYPE sum = x*Q[N];
        for (int k=N-1; k>0; k--)
          sum = x*(Q[k] + sum);
        sum += 1;
        // Tell us if this point is negative
        return sum < 0.0;
      };

      for(unsigned int k=N; k>0; k--){
        // check Q at the subinterval endpoints
        bool Q_has_root = Q_is_negative(-m_stepSize/2.0) || Q_is_negative(m_stepSize/2.0);

        // Check Q for negativity at any of its vertexes
        double desc = 0.0;
        if(!Q_has_root)
          switch(k){
            case 1:
              break;
            case 2:
              Q_has_root = Q_is_negative(-Q[1]/(2.0*Q[2]));
              break;
            case 3:
              desc = Q[2]*Q[2]-3*Q[1]*Q[3];
              Q_has_root = desc > 0.0 && 
                (Q_is_negative(-Q[2]+sqrt(desc)/(3*Q[3])) || Q_is_negative(-Q[2]+sqrt(desc)/(3*Q[3])));
              break;
          }

        // switch to using the [M,k-1] pade approximant on this table interval
        if(Q_has_root){
          if(k == 1){
            // our familiar Taylor series
            Q = arma::zeros(N+1);
            Q[0] = 1.0;
            P = T(arma::span(0,M),0);
          }else{
            Q[k] = 0.0;
            Q.rows(0,k-1) = arma::null(T(arma::span(M+1, M+k-1), arma::span(0,k-1)));
            Q = Q/Q[0];
            P = T.rows(0,M)*Q;
          }
        }else // Q is free to go if it has no roots here
          break;
      }

      // move these coefs into m_table
      for (unsigned int k=0; k<M+1; k++)
        m_table[ii].coefs[k] = P[k];

      for (unsigned int k=0; k<N; k++)
        m_table[ii].coefs[M+1+k] = Q[k+1]; // ignore the first coef of Q b/c it's always 1.
    }
  }

  /* build this table from a file. Everything other than m_table is built by UniformLookupTable */
  UniformPadeTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, std::string filename) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, filename)
  {
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;

    // double check the names match
    std::string temp_name = jsonStats["name"].get<std::string>();
    if(temp_name != "UniformPadeTable<" + std::to_string(M) + "," + std::to_string(N) + ">")
      throw std::invalid_argument("Error while reading " + filename + ": "
          "Cannot build a " + temp_name + " from a UniformPadeTable<" +
          std::to_string(M) + "," + std::to_string(N) + ">");

    m_table.reset(new polynomial<OUT_TYPE,M+N+1>[m_numTableEntries]);
    for(unsigned int i=0; i<m_numTableEntries; i++)
      for(unsigned int j=0; j<m_table[i].num_coefs; j++)
        m_table[i].coefs[j] = jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)].get<OUT_TYPE>();
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // nondimensionalized x position
    OUT_TYPE dx  = (x-m_minArg);
    OUT_TYPE x1r = dx/m_stepSize+0.5;
    // index of previous table entry
    unsigned x1 = ((unsigned) x1r);
    dx -= x1*m_stepSize;
    
    // general degree horners method, evaluated from the inside out.
    OUT_TYPE P = dx*m_table[x1].coefs[M];
    for (int k=M-1; k>0; k--)
      P = dx*(m_table[x1].coefs[k] + P);
    P = P+m_table[x1].coefs[0];

    OUT_TYPE Q = dx*m_table[x1].coefs[M+N];
    for (int k=N-1; k>0; k--)
      Q = dx*(m_table[x1].coefs[M+k] + Q);
    Q = 1+Q;  // the constant term in Q will always be 1
    return P/Q;
  }

  std::function<adVar<OUT_TYPE,M+N>(adVar<OUT_TYPE,M+N>)> boost_function(){ return mp_boost_func; }
};
