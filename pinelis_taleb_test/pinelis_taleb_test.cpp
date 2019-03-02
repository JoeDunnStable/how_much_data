//
/// \file  pinelis_taleb_test.cpp
/// \package how_much_data
//
/// \author Created by Joseph Dunn on 12/31/18.
/// \copyright  © 2018, 2019 Joseph Dunn. All rights reserved.
//

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
#include <iomanip>
using std::setw;
using std::setprecision;
using std::right;
using std::fixed;
using std::scientific;
using std::defaultfloat;
#include <string>
using std::string;
using std::to_string;
#include <sstream>
using std::istringstream;
using std::ostringstream;
#include <fstream>
using std::ofstream;
#include <vector>
using std::vector;
#include <array>
using std::array;
#include <algorithm>
using std::max;
using std::min;
using std::max_element;
using std::sort;
#include <mutex>
using std::mutex;
using std::unique_lock;
#include <thread>
using std::thread;

#include <complex>
using dcomplex = std::complex<double>;
using std::real;

#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;
using boost::timer::cpu_timer;
#include <boost/filesystem.hpp>
using boost::filesystem::path;
using boost::filesystem::is_directory;
using boost::filesystem::create_directory;
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;
using boost::math::constants::one_div_root_two_pi;
using boost::math::constants::one_div_two_pi;

#include "adaptive_integration.h"
using namespace adaptive_integration;

Kronrod<double> k_big(10);
int noext = 0;
double epsabs_double = 0;
double epsrel_double = 64 * std::numeric_limits<double>::epsilon();
int limit = 2000;
int verbose_integration = 0;

/// Controller of Pinelis Taleb integral.
static IntegrationController<double> ctl_double(noext, k_big,
                                         epsabs_double, epsrel_double,
                                         limit, verbose_integration);

/// Controller of beta function integrals.
static IntegrationController<double> ctl_beta(noext, k_big,
                                                epsabs_double, epsrel_double,
                                                limit, verbose_integration);

/// Controller of mad2 integrals in pareto.
static IntegrationController<double> ctl_mad2(noext, k_big,
                                              epsabs_double, epsrel_double,
                                              limit, verbose_integration);

#include "pareto_distribution.h"
template<>
IntegrationController<double>& pareto_distribution<double>::_ctl_beta{ctl_beta};
template<>
IntegrationController<double>& pareto_distribution<double>::_ctl_mad2{ctl_mad2};

#include "student_t_distribution.h"
#include "exponential_distribution.h"
#include "lognormal_distribution.h"
#include "taleb_results.h"

/// sturcture passed by imbue to ostreams to use commas in numbers
struct Numpunct: public std::numpunct<char>{
protected:
  virtual char do_thousands_sep() const{return ',';}
  virtual std::string do_grouping() const{return "\03";}
};

/// return the relative differnece between two numbers
template<typename RealType>
RealType rel_err(RealType a, RealType b) {
  return fabs(a-b)/std::max(a,b);
}

/// structure holding the results of a run for a single parameter value
struct KappaResult {
  /// constructor
  KappaResult (vector<int> ns  ///< the durations of the output
               ) : ns(ns) {}
  double param;        ///< the parameter for the run
  double mad_rel_err;  ///< the relative error of mad vs theory
  vector<int> ns;      ///< the durations saved
  vector<double> mad;  ///< the mean absolute deviation by duration
  vector<double> kappa_mad; /// the kappa_mad by duration
  /// calculate the kappa from the mad and ci variables
  void calc_kappa() {
    for (size_t i = 1; i<ns.size(); ++i) {
      kappa_mad.push_back(2- (log(ns.at(i))-log(ns.at(0)))/
                          (log(mad.at(i))-log(mad.at(0))));
    }
  }
};

 /// output the result for one parameter value to an ostream
ostream& operator<< (ostream& os, const KappaResult& k) {
  os << fixed << setw(7) << setprecision(2) << k.param
  << setw(12) << setprecision(2) << k.mad_rel_err*100 << "%";
  for (size_t i=1; i<k.mad.size(); ++i) {
    os << setw(13) << setprecision(3) << k.kappa_mad.at(i-1);
  }
  os << endl;
  return os;
}

/// structure holding the results of all parameter values
struct KappaResults {
  /// constructor
  KappaResults(const vector<int>& ns,  ///< the durations
               size_t n_params,       ///< the # of params
               string param_label,     ///< the param_label
               size_t taleb_offset     ///< the column offset into Talebs table
               )
  : ns(ns), param_label(param_label), kr(n_params, ns),
    taleb_offset(taleb_offset) {}
  /// return reference to a particular result
  KappaResult& at(size_t i) {
    unique_lock<mutex> lock(kr_mutex);
    return this->kr.at(i);
  }
  
  void dump_results(ostream& os, string dist_name) {
    for (auto k : kr ) {
      for (size_t j=0; j<ns.size()-1 ; ++j) {
        os << setw(15) << dist_name << ","
           << setw(15) << param_label << ","
           << setw(20) << setprecision(8) << k.param << ","
           << setw(10) << ns.at(j+1) << ","
           << setw(20) << setprecision(8) << k.kappa_mad.at(j)
           << endl;
      }
    }
  }
  vector<int> ns;           ///< the durations saved
  string param_label;       ///< the name of the parameter
  vector<KappaResult> kr;   ///< the results of the runs by duration
  size_t taleb_offset;      ///< the column offset into Taleb's table
  mutex kr_mutex;           ///< a mutex for writing results
};

/// output the results for all parameter values to an ostream
ostream& operator<< (ostream& os, KappaResults& ks) {
  os << setw(7) << right << ks.param_label
     << setw(13) << right << "mad_rel_err";
  for (int i=1; i<ks.ns.size(); ++i)
    os << setw(13) << right << "kappa_mad";
  os << endl << setw(20) << " ";
  for (int i=1; i<ks.ns.size(); ++i)
    os << setw(13) << ks.ns.at(i);
  os << endl << endl;
  for (auto kr : ks.kr)
    os << kr;
  os << endl;
  if (ks.taleb_offset >0) {
    os << setw(54) << right << "Error vs Taleb's Results" << endl << endl;
    for (size_t i=0; i<ks.kr.size(); ++i) {
      os << setw(7) << setprecision(2) << ks.kr.at(i).param << setw(13) << " ";
      for (size_t j=0; j<3; ++j) {
        double kappa_mad = ks.kr.at(i).kappa_mad.at(j);
        double err = kappa_mad - taleb_results.at(i).at(j+ks.taleb_offset);
        os << setw(13) << setprecision(3) << err;
      }
      os << endl;
    }
    os << endl;
  }
  return os;
}

/// Functor used in the integratoin of the Pinelis Taleb result
template<typename Dist>
class PinelisTalebIntegrand {
  Dist& dist;
  int n;
public:
  
  /// construct the integrand
  PinelisTalebIntegrand(Dist& dist, ///< [in] the distribution
                        const int n        ///< [in] the duration
                       ) : dist(dist), n(n) {}
  
  /// return the value of the Pinelis Taleb integrand
  double operator() (double x)
  {
    double exponent = 3.;
    double omega = pow(x/(1-x),exponent);
    double domega = exponent * pow(x/(1-x),exponent-1)*(1/pow(1-x,2));
    if (!isfinite(omega)) return 0.;
    double mean = dist.mean();
    dcomplex fac = exp(dcomplex(0,-mean * omega));
    dcomplex adj_cf = fac * dist.characteristic_function(omega);
    dcomplex adj_cfp = fac* dist.characteristic_function_prime(omega)-dcomplex(0,mean)*adj_cf;
    dcomplex pow_adj_cf = (n==1) ? 1 : pow(adj_cf,n-1);
    double ret =  n * real(pow_adj_cf * adj_cfp) * domega/omega;
    if (isnan(ret) || !isfinite(ret)) {
      ostringstream oss;
      oss << " PinelisTalebIntegrand: ret is nan or infinite" << endl
          << scientific
          << "omega        = " << omega  << endl
          << "adj_cf       =  " << adj_cf.real() << " + i*" << adj_cf.imag()  << endl
          << "adj_cf^(n-1) = " << pow_adj_cf.real() << "+ i*" << pow_adj_cf.imag()  << endl
      << "adj_cfp      = " << adj_cfp.real() << " + i*" << adj_cfp.imag()  << endl;
      throw std::runtime_error(oss.str());
    } else
      return ret;
  }
  
  void operator() (std::vector<double>& xs)
  {
    for (double& x : xs)
      x = operator() (x);
  }
};

/// Functor used in the integration of the Pinelis Taleb result for pareto dist.
template<>
class PinelisTalebIntegrand<pareto_distribution<> > {
  pareto_distribution<>& dist;
  int n;
public:
  
  /// construct the integrand
  PinelisTalebIntegrand(pareto_distribution<>& dist, ///< [in] the distribution
                        const int n        ///< [in] the duration
  ) : dist(dist), n(n) {}
  
  /// return the value of the Pinelis Taleb integrand for pareto distribution
  /// Uses an adjusted contour of integration
  double operator() (double x) const
  {
    const dcomplex i{0.,1.};
    double exponent = 2.;
    dcomplex omega = (x<=.5) ? pow(2*x,exponent) : 1. -i * (x - .5)/(1.-x);
    dcomplex domega =(x<=.5) ? exponent * pow(2,exponent) * pow(x,exponent-1)
                             : - .5 * i / pow(1-x,2);
    if (!isfinite(omega.real()) || !isfinite(omega.imag()) ) return 0.;
    double mean = dist.mean();
    dcomplex fac = exp(-i * mean * omega);
    dcomplex adj_cf = fac * dist.characteristic_function(omega);
    dcomplex adj_cfp = fac* dist.characteristic_function_prime(omega)- i * mean * adj_cf;
    dcomplex pow_adj_cf = (n==1) ? 1 : pow(adj_cf,n-1);
    double ret =  n * real(pow_adj_cf * adj_cfp * domega/omega);
    if (isnan(ret) || !isfinite(ret)) {
      ostringstream oss;
      oss << " PinelisTalebIntegrand: ret is nan or infinite" << endl
          << "omega        = " << omega.real() << " + i*" << omega.imag()  << endl
          << "adj_cf       =  " << adj_cf.real() << " + i*" << adj_cf.imag()  << endl
          << "adj_cf^(n-1) = " << pow_adj_cf.real() << "+ i*" << pow_adj_cf.imag()  << endl
          << "adj_cfp      = " << adj_cfp.real() << " + i*" << adj_cfp.imag()  << endl;
      throw std::runtime_error(oss.str());
    } else
      return ret;
  }

  void operator() (std::vector<double>& xs) const
  {
    for (double& x : xs)
      x = operator() (x);
  }
};

/// Functor used in the integration of the Pinelis Taleb result for exp. dist.
template<>
class PinelisTalebIntegrand<exponential_distribution<> > {
  const exponential_distribution<>& dist;
  int n;
public:
  
  /// construct the integrand
  PinelisTalebIntegrand(const exponential_distribution<>& dist, ///< [in] the distribution
                        const int n        ///< [in] the duration
  ) : dist(dist), n(n) {}
  
  /// return the value of the Pinelis Taleb integrand for exp distribution
  /// Uses an adjusted contour of integration
  double operator() (double x) const
  {
    const dcomplex i{0.,1.};
    const double exponent = 2.;
    dcomplex omega = (x<=.5) ? pow(2*x,exponent) : 1. -i * (x - .5)/(1.-x);
    dcomplex domega =(x<=.5) ? exponent * pow(2,exponent) * pow(x,exponent-1)
                             : - .5 * i / pow(1-x,2);
    if (!isfinite(omega.real()) || !isfinite(omega.imag()) ) return 0.;
    double mean = dist.mean();
    dcomplex fac = exp(-i * mean * omega);
    dcomplex adj_cf = fac * dist.characteristic_function(omega);
    dcomplex adj_cfp = fac* dist.characteristic_function_prime(omega)- i * mean * adj_cf;
    dcomplex pow_adj_cf = (n==1) ? 1 : pow(adj_cf,n-1);
    double ret =  n * real(pow_adj_cf * adj_cfp * domega/omega);
    if (isnan(ret) || !isfinite(ret)) {
      ostringstream oss;
      oss << " PinelisTalebIntegrand: ret is nan or infinite" << endl
          << "omega        = " << omega.real() << " + i*" << omega.imag()  << endl
          << "adj_cf       =  " << adj_cf.real() << " + i*" << adj_cf.imag()  << endl
          << "adj_cf^(n-1) = " << pow_adj_cf.real() << "+ i*" << pow_adj_cf.imag()  << endl
          << "adj_cfp      = " << adj_cfp.real() << " + i*" << adj_cfp.imag()  << endl;
      throw std::runtime_error(oss.str());
    } else
      return ret;
  }

  void operator() (std::vector<double>& xs) const
  {
    for (double& x : xs)
      x = operator() (x);
  }
};
  
/// set up integrand and calculate the results for one parameter value
template<typename Dist>
void calculate_kappa(
                     vector<int> ns,  ///< [in] the durations to calculate
                     Dist dist,       ///< [in] the distribution
                     KappaResult& k,  ///< [out] the results
                     int verbose = 0 ///< [in] flag for trace
                     )
{
  if (verbose)
    cout << dist << endl;
  for (int n : ns) {
    PinelisTalebIntegrand<Dist> pti(dist, n);
    const double lower_limit = 0.;
    const double upper_limit = 1;
    vector<double> points;
    double result;
    double error;
    int neval;
    IntegrationController<double>::TerminationCode t_c;
    int last;
    double mad;
  
    do {
      if (points.size()==0)
        points = {lower_limit, upper_limit};
      else {
        // the last iteration failed to find a non-zero value
        // split first interval in two.
        // this happens when n gets very large
        points.insert(points.begin(),.5*(points.at(0)+points.at(1)));
      }
      ctl_double.integrate(pti, points, result, error, neval, t_c, last);
      if (verbose > 3) {
        Fmt<double> fmt;
        double rsum=0, esum=0;
        for (int i=0; i<last; i++) {
          rsum += ctl_double.subs.at(i).r;
          esum += ctl_double.subs.at(i).e;
        }
        
        if (t_c > 0)
          cout << ctl_double.tc_msg(t_c) << ":" << endl;
        cout << "Integral from " << fmt << points.front()
        << " to " << fmt << points.back()
        << " = " << fmt << result
        << ", with absolute error = " << fmt << error << endl
        << "Number of evaluations = " << neval
        << ", subintervals = " << last << endl
        << "rsum = " << fmt << rsum << ", esum = " << fmt << esum << endl;
        print_subs_summary(cout, ctl_double.subs, last, points);
      }
      if (verbose>=4){
        print_subs(cout, ctl_double.subs, last, points);

      }
      mad = (-2./pi<double>())*result;
    } while (false);
    if (verbose)
      cout << "n = " << n
           << ", mad = " << mad
           << ", error = " << error
           << ", neval = " << neval
           << ", t_c = " << ctl_double.tc_msg(t_c) << endl;
    k.mad.push_back(mad);
  }
  
  k.calc_kappa();
  k.mad_rel_err = rel_err(k.mad.at(0), dist.mad());
}

/// show proper usage after improper command line argument
void show_usage(path p     ///< [in] the path of the executable
                ) {
  cerr << "Usage: " << p.filename().string() << endl;
}


/// main program for pinelis-taleb run
int main(int argc, const char * argv[]) {
  path p(argv[0]);
  if ( argc !=1  ) {
    show_usage(p);
    return 1;
  }
/*
  {
    Kronrod<double> k_big(10);
    int noext = 0;
    double epsabs_double = std::numeric_limits<double>::epsilon()/128;
    double epsrel_double = boost::math::tools::root_epsilon<double>();
    int limit = 1000;
    int verbose_integration = 0;
    IntegrationController<double> cf_ctl(noext, k_big,
                                         epsabs_double, epsrel_double,
                                         limit, verbose_integration);
    cout << endl;
    double sigma = 1;

    int type = 2;    // Use special contour
    lognormal_distribution<> lnd(0,sigma, cf_ctl, type);
    vector<int> ns{1000000};
    KappaResult kr(ns);
    kr.param = sigma;
    calculate_kappa(ns, lnd, kr, 4);
    ofstream out("../output/pinelis_taleb_trace.out");
    ns ={1, 2, 30, 100, 1000, 1000000};
    out << setw(16) << "sigma,"
        << setw(11) << "n,"
        << setw(16) << "s,"
    << setw(20) << "integrand" << endl;
    for (auto n : ns) {
      PinelisTalebIntegrand<lognormal_distribution<> > pti(lnd, n);
      for (int i=0; i<1000; i++) {
        double s = (i+.5)/1000;
        out << setw(15) << sigma << ","
            << setw(10) << n << ","
            << setw(15) << s << ","
            << setw(20) << pti(s) << endl;
      }
    }


  }
*/
  


  // Add commas to the cout format
  cout.imbue({std::locale(), new Numpunct});
  
  // for pareto and student t distribution use the same alphas as Taleb
  vector<double> alphas;
  for (size_t i=0; i<taleb_results.size(); ++i)
    alphas.push_back(taleb_results.at(i).at(0));
  vector<double> lambdas{.5, 1, 2};
  
  // The durations for which MAD will be calculated
  vector<int> ns{1, 2, 30, 100, 1000, 1000000};
  
  string out_dir{"../output"};
  if (!is_directory(out_dir))
    create_directory(out_dir);
  ostringstream oss;
  oss << out_dir << "/pinelis_taleb_test.out";
  cout << "Writing to file " << oss.str() << endl << endl;
  ofstream out{oss.str()};
  out.imbue({std::locale(), new Numpunct});
  
  string dump_file_name = out_dir + "/pinelis_taleb_dump.out";
  cout << "Writing to file " << dump_file_name << endl << endl;
  cout.flush();
  ofstream dump{dump_file_name};
  dump << setw(15) << "distribution" << ","
  << setw(15) << "param_name" << ","
  << setw(20) << "param_value" << ","
  << setw(10) << "n" << ","
  << setw(20) << "kappa_mad"
  << endl;


  auto_cpu_timer t(out);
  {
    KappaResults ks_pareto(ns, alphas.size(), "alpha", 1);
    for (size_t i=0; i<alphas.size(); ++i) {
      double alpha=alphas.at(i);
      pareto_distribution<> pd(alpha);
      KappaResult kr(ns);
      kr.param = alpha;
      calculate_kappa(ns, pd, kr, true);
      ks_pareto.at(i) = kr;
    }
    out << "Pareto Distribution" << endl << endl;
    out << ks_pareto;
    ks_pareto.dump_results(dump, "pareto");
    out.flush();
  }

  {
    KappaResults ks_student(ns, alphas.size(), "alpha",4);
    for (size_t i=0; i<alphas.size(); ++i) {
      double alpha = alphas.at(i);
      student_t_distribution<> td(alpha);
      KappaResult kr(ns);
      kr.param = alpha;
      calculate_kappa(ns, td, kr, true);
      ks_student.at(i)=kr;
    }
    out << "Student t Distribution" << endl << endl;
    out << ks_student;
    ks_student.dump_results(dump, "student_t");
  }

  {
    KappaResults ks_exponential(ns, lambdas.size(), "lambda", 0);
    for (size_t i=0; i<lambdas.size(); ++i) {
      double lambda = lambdas.at(i);
      exponential_distribution<> ed(lambda);
      KappaResult kr(ns);
      kr.param = lambda;
      calculate_kappa(ns, ed, kr, true);
      ks_exponential.at(i) = kr;
    }
    out << "Exponential Distribution" << endl << endl;
    out << ks_exponential;
    ks_exponential.dump_results(dump, "exponential");
  }

  {
    Kronrod<double> k_big(10);
    int noext = 0;
    double epsabs_double = std::numeric_limits<double>::epsilon()/128;
    double epsrel_double = boost::math::tools::root_epsilon<double>();
    int limit = 1000;
    int verbose_integration = 0;
    IntegrationController<double> cf_ctl(noext, k_big,
                                         epsabs_double, epsrel_double,
                                         limit, verbose_integration);
    cout << endl;
    vector<double> sigmas = {.01, .1, .5, 1, 2, 5};
    KappaResults ks_lognormal(ns, sigmas.size(), "sigma",0);
    for (size_t i=0; i<sigmas.size(); ++i) {
      double sigma = sigmas.at(i);
      int type = 2;    // Use special contour
      lognormal_distribution<> lnd(0,sigma, cf_ctl, type);
      KappaResult kr(ns);
      kr.param = sigma;
      calculate_kappa(ns, lnd, kr, true);
      ks_lognormal.at(i)=kr;
    }
    out << "Lognormal Distribution" << endl << endl;
    out << ks_lognormal;
    ks_lognormal.dump_results(dump, "lognormal");
  }
  
  return 0;
}
