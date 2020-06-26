#include <boost/math/differentiation/autodiff.hpp>

using boost::math::differentiation::autodiff_fvar;
typedef autodiff_fvar<double,1> fvar1;
typedef autodiff_fvar<double,2> fvar2;
typedef autodiff_fvar<double,3> fvar3;
typedef autodiff_fvar<double,4> fvar4;
typedef autodiff_fvar<double,5> fvar5;
typedef autodiff_fvar<double,6> fvar6;
typedef autodiff_fvar<double,7> fvar7;

struct FunctionHelper{
  EvaluationFunctor<double,double> *double_func;
  EvaluationFunctor<fvar1,fvar1> *fvar1_func;
  EvaluationFunctor<fvar2,fvar2> *fvar2_func;
  EvaluationFunctor<fvar3,fvar3> *fvar3_func;
  EvaluationFunctor<fvar4,fvar4> *fvar4_func;
  EvaluationFunctor<fvar5,fvar5> *fvar5_func;
  EvaluationFunctor<fvar6,fvar6> *fvar6_func;
  EvaluationFunctor<fvar7,fvar7> *fvar7_func;
}
