//
//  student_t_distribution.h
//  pareto_test
//
//  Created by Joseph Dunn on 1/3/19.
//  Copyright © 2019 Joseph Dunn. All rights reserved.
//

#ifndef student_t_distribution_h
#define student_t_distribution_h

#include <boost/random.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/beta.hpp>
using boost::math::beta;

template<class RealType = double>
struct student_t_distribution : boost::random::student_t_distribution<RealType> {
  student_t_distribution(RealType alpha) : boost::random::student_t_distribution<RealType>(alpha),
  t_dist(alpha) {}
  
  RealType cdf(RealType x) const{
    return boost::math::cdf(t_dist, x);
  }
  
  RealType quantile(RealType p) const {
    return boost::math::quantile(t_dist, p);
  }
  
  RealType alpha() const {return boost::random::student_t_distribution<RealType>::n();}
  
  /** Return mean of distribution */
  RealType mean() const {return RealType(0);}
  
  /** Return mean average deviation of the distribution */
  RealType mad() const {
    return 2 * sqrt(alpha())/((alpha()-1)*beta(alpha()/2,.5));}
  
  /** Return the confidence interval of the distribution */
  RealType ci(RealType level = RealType(.05)) const
  {
    return boost::math::quantile(t_dist,1-level/2)-boost::math::quantile(t_dist,level/2);
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
    << "95% CI = " << dist.ci() << endl;
    return os;
  }
private:
  boost::math::students_t_distribution<RealType> t_dist;
  
};


#endif /* student_t_distribution_h */
