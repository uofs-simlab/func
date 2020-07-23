#include <utility> // std::pair

// provide more information about function's domain and behaviour
class SpecialPoint
{
  std::pair<double,double> m_point; // x,y coordinate

  // specify why this point is special
  enum DiscontType { None=-1, Discont=0, FirstDiscont=1, SecondDiscont=2, ThirdDiscont=3 };
  enum LimitType { Equals, Approaches, Inf };
  DiscontType m_discType;
  LimitType m_limType;

public:
  SpecialPoint(double x, double y, DiscontType dt, LimitType lt) : m_point(std::make_pair(x,y)), m_discType(dt), m_limType(lt) {}
  SpecialPoint(std::pair<double,double> pt, DiscontType dt, LimitType lt) : m_point(pt), m_discType(dt), m_limType(lt) {}
};
