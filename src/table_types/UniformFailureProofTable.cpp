#include "UniformFailureProofTable.hpp"

#ifdef NDEBUG
  #define RECORD_ARG(x)
  #define PRINT_ARGS(out)
#else
  #include <iostream>
  #include <vector>
  #include <mutex>

  // make sure we don't swallow the semicolon
  #define RECORD_ARG(x)                   \
    do{                                   \
      const std::lock_guard<std::mutex>   \
        lock(m_args_mutex);               \
      m_args.push_back((x));              \
    } while(0)
  #define PRINT_ARGS(out)                                 \
    do {                                                  \
      if(!m_args.empty()){                                \
        out << "args outside table range:" << std::endl;  \
        out << m_args.front();                            \
        m_args.erase(m_args.begin());                     \
        for(auto x : m_args)                              \
          out << ", " << x;                               \
        out << std::endl;                                 \
      }                                                   \
    } while(0)
#endif

// copy everything from the given LUT
UniformFailureProofTable::UniformFailureProofTable(std::unique_ptr<UniformLookupTable> LUT) :
  mp_LUT(std::move(LUT))
{
  mp_func    = mp_LUT->function();
  m_minArg   = mp_LUT->min_arg();
  m_maxArg   = mp_LUT->max_arg();
  m_order    = mp_LUT->order();
  m_name     = mp_LUT->name();
  m_dataSize = mp_LUT->size();
}

double UniformFailureProofTable::operator()(double x)
{
  // check if x is in the range of the table
  if(x<m_minArg || x>m_maxArg){
    RECORD_ARG(x);
    return (*mp_func)(x);
  }
  return (*mp_LUT)(x);
}

void UniformFailureProofTable::print_details(std::ostream &out)
{
  mp_LUT->print_details(out);
}

UniformFailureProofTable::~UniformFailureProofTable()
{
  // is a null op if iostream isn't included
  PRINT_ARGS(std::cerr);
}
