#include "CompositeLookupTable.hpp"
#include <stdexcept> // domain_error
#include <string> // to_string

// Currently incompatible with failure proof table b/c it doesn't have a function
CompositeLookupTable::CompositeLookupTable(std::initializer_list<std::shared_ptr<UniformLookupTable>> LUT_list) :
  mv_LUT(LUT_list), EvaluationImplementation(NULL, "")
{
  // sort the vector based on the table's min args with insertion sort
  for(int i=1; i<mv_LUT.size(); i++){
    if(mv_LUT[i]->min_arg()<mv_LUT[i-1]->min_arg()){
      int j = i;
      while(j > 0 && mv_LUT[j]->min_arg()<mv_LUT[j-1]->min_arg()){
        auto temp       = mv_LUT[j-1];
        mv_LUT[j-1] = mv_LUT[j];
        mv_LUT[j]   = temp;
        j -= 1;
      }
    }
  }
  
  // set the global min/max
  m_minArg = mv_LUT.front()->min_arg();
  m_maxArg = mv_LUT.back()->max_arg();
  for(int i=1; i<mv_LUT.size(); i++){
    // assert the table intervals do not overlap
    double space = mv_LUT[i]->min_arg() - mv_LUT[i-1]->max_arg();
    assert(space >= 0.0);

    // find any holes in the range
    if(space > std::numeric_limits<double>::epsilon()){
      mv_discontinuities.push_back(std::make_pair(mv_LUT[i-1]->max_arg(),mv_LUT[i]->min_arg()));
    }
  }

  // set other member vars cumulatively
  m_order = 0;
  m_dataSize = 0;
  m_name  = "CompoundLookupTable";
  for(auto lut : mv_LUT){
    // get every special points
    for(auto pt : lut->special_points())
      this->m_special_points.push_back(pt);
    m_order    += lut->order();
    m_dataSize += lut->size();
  }

  // init most recently used index to the midpoint
  mostRecentlyUsed_idx = (int) mv_LUT.size()/2;
}

double CompositeLookupTable::binarySearch(double x, int i, int min_idx, int max_idx)
{
  // Binary search for the correct interval starting with i (most
  // recently used table index). Best for seemly random table evaluations.
  // Note: Currently doesn't check to see if the arg is less/greater than
  // overall table min/max
  if(x < mv_LUT[i]->min_arg()){
    if(i == min_idx)
      throw std::domain_error(std::string("Composite table undefined for x=") +
          std::to_string(x));
    return binarySearch(x, (int)(i+min_idx)/2, min_idx, i);
  }
  else if(x > mv_LUT[i]->max_arg()){
    if(i == max_idx)
      throw std::domain_error(std::string("Composite table undefined for x=") +
          std::to_string(x));
    return binarySearch(x, (int)(i+max_idx)/2, i, max_idx);
  }
  return (*mv_LUT[i])(x);
}

// TODO template or enum specifying binary vs linear search
double CompositeLookupTable::operator()(double x)
{
  return binarySearch(x, mostRecentlyUsed_idx,0,mv_LUT.size()-1);
}

void CompositeLookupTable::print_details(std::ostream &out)
{
  out << m_name;
  for(auto lut : mv_LUT)
    lut->print_details(out);
}

CompositeLookupTable::~CompositeLookupTable(){}
