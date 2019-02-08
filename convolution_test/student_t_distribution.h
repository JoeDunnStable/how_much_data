//
/// \file student_t_distribution.h
/// \package how_much_data
//
/// \author  Created by Joseph Dunn on 1/3/19.
/// \copyright © 2019 Joseph Dunn. All rights reserved.
//

#ifndef student_t_distribution_h
#define student_t_distribution_h

#include <boost/random.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/beta.hpp>
using boost::math::beta;
#include <boost/math/special_functions/bessel.hpp>
using boost::math::cyl_bessel_k;

/** Instances of class student_t_distribution give random variates
   for a student t distribution with parameter n = alpha
 */
template<class RealType = double>
struct student_t_distribution : boost::random::student_t_distribution<RealType> {
  /// construct an instnace give alpha
  student_t_distribution(RealType alpha) : boost::random::student_t_distribution<RealType>(alpha),
  t_dist(alpha) {}
  
  /// return the cdf or the complement of the cdf
  RealType cdf(RealType x, bool lower_tail = true) const{
    return lower_tail ? boost::math::cdf(t_dist, x)
                      : boost::math::cdf(boost::math::complement(t_dist,x));
  }
  
  /// return the probability density function at x
  RealType pdf(RealType x) const{
    return boost::math::pdf(t_dist, x);
  }
  
  /// return the quantile for probability p
  RealType quantile(RealType p) const {
    return boost::math::quantile(t_dist, p);
  }
  
  /// return the alpha = n of the distribution
  RealType alpha() const {return boost::random::student_t_distribution<RealType>::n();}
  
  /** return the alpha of the asymptotic stable distribution */
  RealType alpha_stable() const {return min(RealType(2.),alpha());}
  
  /** Return mean of distribution */
  RealType mean() const {return RealType(0);}
  
  /** Return mean average deviation of the distribution */
  RealType mad() const {
    return 2 * sqrt(alpha())/((alpha()-1)*beta(alpha()/2,.5));}
  
  /** Return the mad of the square of the distribution */
  RealType mad2() const {
    double a = alpha();
    double num =(pow(2,3 - a) * sqrt(a) * tgamma(.5 + a/2) * tgamma(-.5 + a));
    double denom =((-1 + a) * pow(tgamma(a/2),3));
    double ret = num/denom;
    return ret;
  }
  
  /** Return the confidence interval of the distribution */
  RealType ci(RealType level = RealType(.05)) const
  {
    return boost::math::quantile(t_dist,1-level/2)-boost::math::quantile(t_dist,level/2);
  }
  
  /// return the characteristic function of the distribution at omega
  RealType characteristic_function(RealType omega) const{
    RealType a=alpha();
    if (omega == 0) return 1;
    else {
      RealType besselk = cyl_bessel_k(a/2.,sqrt(a)*abs(omega));
      return (pow(2,1 - a/2.)*pow(a,a/4.)*pow(abs(omega),a/2.)*besselk)/tgamma(a/2.);
    }
  }
  /**
   * Write distribution to std::ostream
   */
  template <class charT, class traits>
  friend std::basic_ostream<charT,traits>& operator<< ( std::basic_ostream<charT,traits>& os,
                                                       const student_t_distribution& dist ) {
    os << "Student t distribution with n = " << dist.n() << endl;
    os << "Mean   = " << dist.mean() << endl
    << "MAD    = " << dist.mad() << endl
    << "MAD2   = " << dist.mad2() << endl
    << "kappa2 = " << 2 - log(2)/(log(dist.mad2())-log(dist.mad())) << endl
    << "95% CI = " << dist.ci() << endl;
    return os;
  }
private:
  boost::math::students_t_distribution<RealType> t_dist;
  
};


#endif /* student_t_distribution_h */
