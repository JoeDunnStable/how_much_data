//
///  \file myFloat.h
/// Declares four different types of floats to be used:
/// 1. double, which is the default
/// 2. mpreal, when the preprocessor variable MPREAL is defined
/// 3. mprf_float, the boost multiprecision wrapper around mpfr when MPFR_FLOAT_50 is defined
/// 4. cpp_bin_float, the boost multiprecion cpp_bin_float when CPP_BIN_FLOAT is defined
///
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef myFloat_h
#define myFloat_h

#undef BOOST_MATH_OVERFLOW_ERROR_POLICY
#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <Eigen/Eigenvalues>

#ifdef CPP_BIN_FLOAT
#include <boost/multiprecision/cpp_bin_float.hpp>
using CppBinFloat = boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<28>, boost::multiprecision::et_off>;
using BigCppBinFloat = boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<38>, boost::multiprecision::et_off>;
#endif

#if defined(MPFR_FLOAT)
#include <boost/multiprecision/mpfr.hpp>

using MpfrFloat = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<28>, boost::multiprecision::et_off>;
using BigMpfrFloat = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<38>, boost::multiprecision::et_off>;
#endif

#ifdef MPREAL
#include <mpreal.h>
using mpfr::mpreal;
using mpfr::digamma;
using mpfr::const_pi;
using mpfr::const_euler;

#include <boost/math/tools/real_cast.hpp>
#include <boost/math/special_functions/trunc.hpp>
#endif //MPREAL

// These are for double

#include <cmath>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/sinc.hpp>
//using fmax = boost::multiprecision::max;
//using fmin = boost::multiprecision::min;

/*
using boost::math::isfinite;
using boost::math::isnan;
using boost::math::isinf;
using boost::math::expm1;
using boost::math::log1p;
using boost::math::erf;
using boost::math::erfc;
using boost::math::tgamma;
using boost::math::lgamma;
*/
using boost::math::digamma;
using boost::math::tgamma_ratio;
using boost::math::binomial_coefficient;
using boost::math::factorial;
using boost::math::zeta;
using boost::math::sinc_pi;

template<typename myFloat>
inline myFloat erf_inv(myFloat p) {
  return boost::math::erf_inv(p);
}

template<typename myFloat>
inline myFloat erfc_inv(myFloat p) {
  return boost::math::erfc_inv(p);
}

template<typename myFloat>
inline myFloat my_erfc(myFloat x) {
  return erfc(x);
}

#ifdef CPP_BIN_FLOAT

#include <boost/math/constants/constants.hpp>

template<> CppBinFloat my_erfc<CppBinFloat>(CppBinFloat x) {
  // https://dlmf.nist.gov/7.12.E1
  using namespace boost::math::constants;
  if (x<=100) return erfc(x);
  CppBinFloat x_sq{x*x};
  CppBinFloat fac = exp(-x_sq)/(x*root_pi<CppBinFloat>());
  CppBinFloat sum = 0;
  CppBinFloat term = 1;
  int n = 0;
  for (n=0; n<10; n++) {
    sum += term;
    term *= (-1)*(static_cast<CppBinFloat>(.5)+n)/x_sq;
    if (fabs(term) < std::numeric_limits<CppBinFloat>::epsilon()) break;
  }
  return fac * sum;
}

#endif

#ifdef MPREAL

#include <boost/multiprecision/mpfr.hpp>
using MpfrProxy = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<29>, boost::multiprecision::et_off>;


// the boost version of sinc_pi, tgamma_ratio, erf_inv, erfc_inv is broken
// for variable precision mpreal

mpreal sinc_pi(mpreal x) {
  MpfrProxy x_mpfr{x.mpfr_ptr()};
  mpreal ret = static_cast<mpreal>(boost::math::sinc_pi(x_mpfr).backend().data());
  ret.set_prec(mpreal::get_default_prec());
  return ret;
}

mpreal tgamma_ratio(mpreal num, mpreal denom) {
  MpfrProxy num_mpfr{num.mpfr_ptr()}, denom_mpfr{denom.mpfr_ptr()};
  mpreal ret = static_cast<mpreal>(tgamma_ratio(num_mpfr, denom_mpfr).backend().data());
  ret.set_prec(mpreal::get_default_prec());
  return ret;
}

mpreal erf_inv(mpreal x) {
  MpfrProxy x_mpfr{x.mpfr_ptr()};
  mpreal ret = static_cast<mpreal>(boost::math::erf_inv(x_mpfr).backend().data());
  ret.set_prec(mpreal::get_default_prec());
  return ret;
}

mpreal erfc_inv(mpreal x) {
  MpfrProxy x_mpfr{x.mpfr_ptr()};
  mpreal ret = static_cast<mpreal>(boost::math::erfc_inv(x_mpfr).backend().data());
  ret.set_prec(mpreal::get_default_prec());
  return ret;
}

namespace mpfr {
inline long long lltrunc(mpfr::mpreal const& x)
{
  MpfrProxy x_mpfr{x.mpfr_ptr()};
  return boost::math::lltrunc(x_mpfr);
}
}


#endif

/*
template<typename myFloat>
myFloat pow(myFloat x, myFloat y) {
  return pow(x, y);
}

#ifdef CPP_BIN_FLOAT
template<> CppBinFloat pow<CppBinFloat>(CppBinFloat x, CppBinFloat y) {
  if (x<0)
    return NAN;
  else if (x == 0 && y == 0)
    return NAN;
  else if (x == 0 && y < 0)
    return std::numeric_limits<CppBinFloat>::infinity();
  else if (x == 0 && y>0)
    return 0;
  else
    return pow(x, y);
}
#endif

#ifdef MPFR_FLOAT
template<> MpfrFloat pow<MpfrFloat>(MpfrFloat x, MpfrFloat y) {
  if (x<0)
    return NAN;
  else if (x == 0 && y == 0)
    return NAN;
  else if (x == 0 && y < 0)
    return std::numeric_limits<MpfrFloat>::infinity();
  else if (x == 0 && y>0)
    return 0;
  else
    return pow(x, y);
}
#endif
 
 */
#include <boost/math/constants/constants.hpp>
template<typename myFloat>
inline myFloat const_pi() {
  return boost::math::constants::pi<myFloat>();
}

#ifdef MPREAL
template<>
inline mpreal const_pi<mpreal>() {return const_pi();}
#endif

template<typename myFloat>
inline myFloat const_euler() {
  return boost::math::constants::euler<myFloat>();
}

#ifdef MPREAL
template<>
inline mpreal const_euler<mpreal>() {return const_euler();}
#endif

template<typename myFloat>
void reset_prec(myFloat& x) {}

#ifdef MPREAL

template<>
void reset_prec<mpreal>(mpreal& x) {
  x.set_prec(mpreal::get_default_prec());
}

#endif

template<typename T>
class Fmt {
public:
  int digits10;
  int width;
  Fmt() {
    digits10=static_cast<int>(-log(std::numeric_limits<T>::epsilon())/log(10));
    width=digits10+8;
  }
  friend
  std::ostream& operator<< (std::ostream& os, const Fmt<T>& fmt) {
    os << std::setw(fmt.width)
    << std::setprecision(fmt.digits10)
    << std::scientific;
    return os;
  }
};

template<typename myFloat>
myFloat exp_m_xi(myFloat alpha, myFloat xB) {
  myFloat xi = (alpha==1)
                ? exp(xB-1)
                : fabs(1-alpha) * pow(xB/alpha, alpha/(alpha-1));
  return exp(-xi);
}

#include <boost/multiprecision/cpp_bin_float.hpp>
template<>
double exp_m_xi<double> (double alpha, double xB) {
  using namespace boost::multiprecision;
  using BigFloat = number<backends::cpp_bin_float<20>, et_off>;
  BigFloat a = alpha;
  BigFloat x = xB;
  BigFloat xi = (a==1)
                 ? exp(x-1)
                 : fabs(1-a) * pow(x/a, a/(a-1));
  return static_cast<double>(exp(-xi));
}

#ifdef MPREAL
template<>
mpreal exp_m_xi<mpreal> (mpreal alpha, mpreal xB) {
  mpreal a = alpha; a.set_prec(128);
  mpreal x = xB; x.set_prec(128);
  mpreal xi = (a==1)
               ? exp(x-1)
               : fabs(1-a) * pow(x/a, a/(a-1));
  mpreal ret = exp(-xi); ret.set_prec(mpreal::get_default_prec());
  return ret;
}
#endif

#ifdef MPFR_FLOAT
template<>
MpfrFloat exp_m_xi<MpfrFloat> (MpfrFloat alpha, MpfrFloat xB) {
  BigMpfrFloat a = alpha;
  BigMpfrFloat x = xB;
  BigMpfrFloat xi = (a==1)
                 ? exp(x-1)
                 : fabs(1-a) * pow(x/a, a/(a-1));
  return static_cast<MpfrFloat>(exp(-xi));
}
#endif

#ifdef CPP_BIN_FLOAT


template<>
CppBinFloat exp_m_xi<CppBinFloat> (CppBinFloat alpha, CppBinFloat xB) {
  BigCppBinFloat a = alpha;
  BigCppBinFloat x = xB;
  BigCppBinFloat xi = (a==1)
  ? exp(x-1)
  : fabs(1-a) * pow(x/a, a/(a-1));
  return static_cast<CppBinFloat>(exp(-xi));
}
#endif

#include <string>
using std::string;

namespace Eigen {
  
#ifdef CPP_BIN_FLOAT
  template<> struct NumTraits<CppBinFloat> : GenericNumTraits<CppBinFloat>
  {
    static inline Real dummy_precision() { return epsilon()*static_cast<CppBinFloat>(1024); }
  };
#endif
    
#ifdef MPFR_FLOAT
    template<> struct NumTraits<MpfrFloat> : GenericNumTraits<MpfrFloat>
  {
    static inline Real dummy_precision() { return epsilon()*static_cast<MpfrFloat>(1024); }
  };
#endif
  
#ifdef MPREAL
  template<> struct NumTraits<mpreal> : GenericNumTraits<mpreal>
  {
    static inline Real dummy_precision() { return epsilon()*static_cast<mpreal>(1024); }
  };
#endif
  
} // namespace eigen

#endif /* myFloat_h */

