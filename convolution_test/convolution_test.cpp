//
/// \file  convolution_test.cpp
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
#include <utility>
using std::pair;
#include <mutex>
using std::mutex;
using std::unique_lock;
#include <thread>
using std::thread;

#include <complex>
using dcomplex = std::complex<double>;
#include <Eigen/Dense>
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::Map;

/// cause Eigen/FFT to use fftw.  delete to avoid GPL
#define EIGEN_FFTW_DEFAULT
#include <unsupported/Eigen/FFT>

Eigen::FFT<double> fft_eng;

#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;
using boost::timer::cpu_timer;
#include <boost/filesystem.hpp>
using boost::filesystem::path;
#include <boost/math/tools/roots.hpp>
using boost::math::tools::bracket_and_solve_root;
using boost::math::tools::eps_tolerance;
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;
using boost::math::constants::one_div_root_two_pi;
using boost::math::constants::one_div_two_pi;

#include "pareto_distribution.h"
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

/// structure holding the results of a run for a single alpha
struct KappaResult {
  /// constructor
  KappaResult (vector<int> ns  ///< the durations of the output
               ) : ns(ns) {}
  double param;        ///< the parameter for the run
  int nsize;           ///< the size of the vector pased to fft
  double mad_rel_err;  ///< the relative error of mad vs theory
  vector<int> ns;      ///< the durations saved
  vector<double> mad;  ///< the mean absolute deviation by duration
  vector<double> kappa_mad; /// the kappa_mad by duration
  double ci_rel_err;   ///< the relative error of conf. int. vs theory
  vector<double> ci;   ///< the conficence interval by duration
  vector<double> kappa_ci; ///< the kappa ci by duration
  /// calculate the kappa from the mad and ci variables
  void calc_kappa() {
    for (size_t i = 1; i<ns.size(); ++i) {
      kappa_mad.push_back(2- (log(ns.at(i))-log(ns.at(0)))/
                          (log(mad.at(i))-log(mad.at(0))));
      kappa_ci.push_back(2- (log(ns.at(i))-log(ns.at(0)))/
                         (log(ci.at(i))-log(ci.at(0))));
    }
  }
};

 /// output the results to an ostream
ostream& operator<< (ostream& os, const KappaResult& k) {
  os << fixed << setw(7) << setprecision(2) << k.param
  << setw(14) << k.nsize
  << setw(12) << setprecision(2) << k.mad_rel_err*100 << "%";
  for (size_t i=1; i<k.mad.size(); ++i) {
    os << setw(13) << setprecision(3) << k.kappa_mad.at(i-1);
  }
  os << setw(12) << setprecision(2) << k.ci_rel_err*100 << "%";
  for (size_t i=1; i<k.ci.size(); ++i) {
    os << setw(13) << setprecision(3) << k.kappa_ci.at(i-1);
  }
  os << endl;
  return os;
}

/// structure holding the results of all runs
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
  vector<int> ns;           ///< the durations saved
  string param_label;       ///< the name of the parameter
  vector<KappaResult> kr;   ///< the results of the runs by duration
  size_t taleb_offset;      ///< the column offset into Taleb's table
  mutex kr_mutex;           ///< a mutex for writing results
};

/// output the results of all runs to an ostream
ostream& operator<< (ostream& os, KappaResults& ks) {
  os << setw(7) << right << ks.param_label
     << setw(14) << right << "nsize"
     << setw(13) << right << "mad_rel_err";
  for (int i=1; i<ks.ns.size(); ++i)
    os << setw(13) << right << "kappa"+to_string(ks.ns.at(i))+"_mad";
  os << setw(13) << right << "ci_rel_err";
  for (int i=1; i<ks.ns.size(); ++i)
    os << setw(13) << right << "kappa"+to_string(ks.ns.at(i))+"_ci";
  os << endl << endl;
  for (auto kr : ks.kr)
    os << kr;
  os << endl;
  if (ks.taleb_offset >0) {
    os << setw(68) << right << "Error vs Taleb's Results" << endl << endl;
    for (size_t i=0; i<ks.kr.size(); ++i) {
      os << setw(7) << setprecision(2) << ks.kr.at(i).param << setw(27) << " ";
      for (size_t j=0; j<ks.ns.size()-1; ++j) {
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

/// modulo calculation return nonnegative number less than the modulus
int mod (int a, int b)
{
  int ret = a % b;
  if(ret < 0)
    ret+=b;
  return ret;
}

/// calcuate confidence interval given pdf at points in range
double confidence_interval(int nmin,    ///< [in] the minimum index
                           int nmax,    ///< [in] the maximum index
                           double delta,///< [in] the spacing of the x's
                           const Map<Matrix<double, Dynamic, 1> >& pdf, ///< [in] the calculated pdf's
                           double ci_level  ///< [in] the confidence level
                           ) {
  int nsize = static_cast<int>(pdf.size());
  double ci_lower_left{nmin * delta};
  double ci_lower_right=ci_lower_left;
  double cdf_lower_left=.5*pdf(mod(nmin,nsize))*delta;
  double cdf_lower_right = cdf_lower_left;
  for (int j=nmin; j<=nmax ; ++j) {
    if (cdf_lower_right > ci_level/2) break;
    ci_lower_left = ci_lower_right;
    cdf_lower_left = cdf_lower_right;
    ci_lower_right= (j+1)*delta;
    cdf_lower_right = cdf_lower_left+.5*(pdf(mod(j, nsize))+pdf(mod(j+1,nsize)))*delta;
  }
  double ci_lower;
  if (ci_lower_left == ci_lower_right)
    ci_lower = ci_lower_right;
  else {
    double t = (ci_level/2-cdf_lower_left)/(cdf_lower_right-cdf_lower_left);
    ci_lower = (1-t) * ci_lower_left + t * ci_lower_right;
  }
  
  double ci_upper_left{nmax * delta};
  double ci_upper_right=ci_upper_left;
  double cdf_upper_left{.5*pdf(mod(nmax,nsize))*delta};
  double cdf_upper_right=cdf_upper_left;
  for (int j=nmax; j>=nmin ; --j) {
    if (cdf_upper_left > ci_level/2) break;
    ci_upper_right = ci_upper_left;
    cdf_upper_right = cdf_upper_left;
    ci_upper_left= (j-1) * delta;
    cdf_upper_left = cdf_upper_right+.5*(pdf(mod(j,nsize))+pdf(mod(j-1, nsize)))*delta;
  }
  double ci_upper;
  if (ci_upper_right == ci_upper_left)
    ci_upper = ci_upper_left;
  else {
    double t = (ci_level/2-cdf_upper_right)/(cdf_upper_left-cdf_upper_right);
    ci_upper = (1-t) * ci_upper_right + t * ci_upper_left;
  }
  return ci_upper-ci_lower;
}

/// functor for determining upper limit for target mad
template<typename Dist>
class Upper {
public:
  /// constructor
  Upper (double delta,  ///< [in] target mad arising from upper tail
         Dist& dist     ///< [in] reference to distribution
         ) : delta(delta), dist(dist) {}
  /// return excess of estimated mad in tail over target
  double operator() (double x   ///< [in] the quantile
                     ) {
    double a = dist.alpha_stable()/(dist.alpha_stable()-1);
    double ret = max(1.,fabs(x-dist.mean()))* a * (dist.cdf(x, false))-delta;
    return ret;
  }
private:
  double delta;
  Dist& dist;
};

/// functor for determining lower limit for target mad
template<typename Dist>
class Lower {
public:
  /// constructor
  Lower (double delta,    ///< [in] the target mad arising from lower tail
         Dist& dist       ///< [in] the distribution
         ) : delta(delta), dist(dist) {}
  /// return excess of estimated mad in tail over target
  double operator() (double x    ///< [in] the trial position
                     ) {
    double a = dist.alpha_stable()/(dist.alpha_stable()-1);
    double ret = max(1.,fabs(x-dist.mean())) * a * (dist.cdf(x, true))-delta;
    return ret;
  }
private:
  double delta;
  Dist& dist;
};

/// Check whether the factorization of n contains only listed primes
bool factor_check(int n,              ///< [in] the number to check
                  vector<int> primes  ///< [in] the array of candidate primes
                  )
{
  int n_current = n;
  for (auto p : primes) {
    while (n_current % p == 0 ) {
      n_current = n_current / p;
    }
  }
  return n_current == 1;
}

/// instances] of Job used to parcel characteristic funciion calculation to threads
struct Job {
  int nmin;          ///< the overall minimum index
  int nmax;          ///< the overall maximum index
  int ncurrent;      ///< the next index to parcel out
  int nchunk;        ///< the size of the range to parcel out
  mutex job_mutex;   ///< to prevent multiple accesses
  /// consturctor
  Job(int nmin,      ///< [in] the overall minimum index
      int nmax,      ///< [in] the overall maximum index
      int nchunk     ///< [in] the size of the range to parcel out
      ) : nmin(nmin), nmax(nmax), ncurrent(nmin),
                                        nchunk(nchunk) {}
  /// get a range and return true if okay
  bool get_next(int& nstart,   ///< [out] start of the assigned range
                int& nend      ///< [out] next index after assigned range
                ) {
    unique_lock<mutex> lock;
    if (ncurrent > nmax) return false;
    nstart = ncurrent;
    ncurrent = min(nmax+1,ncurrent+nchunk);
    nend = ncurrent;
    return true;
  }
};

/// calculate the characteristic function for a range assigned to trhead
template<typename Dist>
void calc_characteristic_function(
                                const Dist& dist, ///< [in] the distribution
                                const double mean,  ///< [in] the mean to remove
                                const int n,        ///< [in] the duration
                                const double delta_omega, ///< [in] step for omega
                                Job *job,           ///< [in,out]ptr to job assigner
                                Matrix<dcomplex, Dynamic, 1> *adj_cf ///< [out] the characeristic
                                  ///< function of the normalized distrbibution
                                  ///> for non-negative frequencies
                                  )
{
  int n_start, n_end;
  while (job->get_next(n_start, n_end)) {
    for (int i=n_start; i<n_end; ++i) {
      double omega = delta_omega*i/pow(n,1/dist.alpha_stable());
      // remove the location parameter from the characteristic function
      dcomplex fac = exp(dcomplex(0,-mean*omega));
      dcomplex cf = pow(fac*dist.characteristic_function(omega),n);
      adj_cf->operator()(i) = cf;
    }
  } // while job
}

/// set up ranges and step sizes for one alpham pass calculation of cf to threads
/// and use fft to estimate distribution
template<typename Dist>
void calculate_kappa(double delta,    ///< [in] the step size in x / dist.mad
                     double delta2,   ///< [in] cap on % mad from the tail
                     int m,           ///< [in] cap on maximum index
                     vector<int> ns,  ///< [in] the durations to calculate
                     Dist dist,       ///< [in] the distribution
                     double ci_level, ///< [in] the confidence level for kappa_ci
                     KappaResult& k,  ///< [out] the results
                     bool verbose = false ///< [in] flag for trace
                     )
{
  double mean = dist.mean();
  delta *= dist.mad();
  delta2 *= dist.mad();
  if (verbose) {
    cout << scientific
         << "Rescaled delta  = " << setw(12) << setprecision(3) << delta << endl
         << "Rescaled delta2 = " << setw(12) << setprecision(3) << delta2 << endl
         << defaultfloat;
  }
  double alpha_stable = dist.alpha_stable();
  if (verbose)
    cout << dist << endl;
  
  // Calculate xmax0 and xmin0
  // so that the contriubution of the tails to MAD is small

  Upper<Dist> upper(delta2, dist);
  double guess = dist.quantile(1-delta2/2);
  double factor = 2;
  bool rising = false;
  eps_tolerance<double> tol;
  boost::uintmax_t max_iter{1000};
  pair<double, double> root = bracket_and_solve_root(upper, guess, factor,
                                                     rising, tol, max_iter);
  double xmax0 = root.second - mean;
  
  Lower<Dist> lower(delta2, dist);
  guess = dist.quantile(delta2/2);
  rising = true;
  max_iter = 1000;
  root = bracket_and_solve_root(lower, guess, factor, rising, tol, max_iter);
  double xmin0 = root.first - mean;
  
  // Since we're centering on the mean we'll set abs(xmax)=abs(xmin)
  double xmax = max(10.,max(fabs(xmax0),fabs(xmin0)));
  // Cap the result to the m input on the command line
  int nmax = static_cast<int>(min(double(m), xmax/delta));
  if (nmax < 1000000) {
    // Add extra points spliting between increasing the range and the density
    delta = delta/sqrt(1000000./nmax);
    nmax = 1000000;
    if (verbose)
      cout << "Revising delta to " << delta << endl;
  }
  int nsize = 2*nmax+1;
  
  // fftw works best for array sizes that are products of small primes
  vector<int> primes{3,5,7,11};
  for (; !factor_check(nsize, primes); nsize+=2) {}
  nmax = ((nsize-1)/2);
  int nmin = -nmax;

  xmax = nmax * delta;
  // The amounts of mad left in the tail after all of the adjustments
  double uppermad = delta2 + upper(xmax);
  double lowermad = delta2 + lower(-xmax);
  k.nsize = nsize;
  if (verbose) {
    cout << scientific
         << setw(12) << "xmin0 = " << setw(15) << setprecision(2) << xmin0 << endl
         << setw(12) << "xmax0 = " << setw(15) << setprecision(2) << xmax0 << endl
         << defaultfloat
         << setw(12) << "nmin = " << setw(15) << nmin << endl
         << setw(12) << "nmax = " << setw(15) << nmax << endl
         << setw(12) << "xmin = " << setw(15) << setprecision(2) << -xmax << endl
         << setw(12) << "xmax = " << setw(15) << setprecision(2) << xmax << endl
         << scientific
         << setw(12) << "mad = " << setw(15) << setprecision(3) << dist.mad() << endl
         << setw(12) << "uppermad = " << setw(15) << setprecision(3) << uppermad << endl
         << setw(12) << "lowermad = " << setw(15) << setprecision(3) << lowermad << endl
         << defaultfloat
         << setw(12) << "nsize = " << setw(15) << nsize << endl << endl;
  }
  // the maximum absolute value of he angular frequency
  double omega_max = pi<double>()/delta;
  double delta_omega = (2*omega_max/nsize);
  // the adjusted characteristic function for the non-negative half spectrum
  // use Eigen::Matrix to assure alignment to 16 butye boundary for fftw
  Matrix<dcomplex,Dynamic,1> adj_cf(nmax+1);
  
  if (verbose) {
    cout << setw(5) << " "
    << setw(13) << right << "x.at(nmin)"
    << setw(13) << right << "x.at(-1)"
    << setw(13) << right << "x.at(0)"
    << setw(13) << right << "x.at(1)"
    << setw(13) << right << "x.at(nmax)"
    << endl;
  
    cout << setw(5) << " "
    << setw(13) << fixed << setprecision(3) << nmin * delta
    << setw(13) << fixed << setprecision(3) << -1. * delta
    << setw(13) << fixed << setprecision(3) << 0. * delta
    << setw(13) << fixed << setprecision(3) << 1. * delta
    << setw(13) << fixed << setprecision(3) << nmax * delta
    << endl << endl;
    cout << setw(5) << right << "n"
    << setw(13) << right << "pdf.at(nmin)"
    << setw(13) << right << "pdf.at(-1)"
    << setw(13) << right << "pdf.at(0)"
    << setw(13) << right << "pdf.at(1)"
    << setw(13) << right << "pdf.at(nmax)"
    << setw(13) << right << "p_total"
    << setw(12) << right << "mad"
    << setw(12) << right << "madlower"
    << setw(12) << right << "madupper"
    << endl;
  }
  for (size_t j=0; j<ns.size(); ++j){
    int n = ns.at(j);
    array<thread,8> threads;
    int nchunk = 1000;
    // Just do the positive half spectrum
    Job job(0, nmax, nchunk);
    for (int i=0; i<8; ++i) {
      threads.at(i)=thread{calc_characteristic_function<Dist>, dist,
                           mean, n, delta_omega,
                           &job, &adj_cf};
    }
    for (int i=0; i<8; ++i)
      threads.at(i).join();
    // fft_eng does an in-place DFT if the pointers to the src & dst match
    // In our case the source is adj_cf and the dst is pdf
    Map<Matrix<double, Dynamic, 1> >pdf(reinterpret_cast<double*>(&adj_cf(0)), nsize);
    fft_eng.inv(&pdf(0), &adj_cf(0), nsize);
    
    for (int i=-nmax; i<=nmax; ++i)
      pdf(mod(i,nsize)) *= delta_omega * nsize * one_div_two_pi<double>();
    double mad = 0;
    for (int i=-nmax; i<=nmax; ++i) {
      mad += pdf(mod(i,nsize)) * fabs(i * delta)*delta;
    }
    // Assume that the relative weight in the tail is the same as the original distribution
    double mad_upper_tail = mad*uppermad/dist.mad();
    double mad_lower_tail = mad*lowermad/dist.mad();;
    mad += mad_upper_tail+mad_lower_tail;
    // Scal mad.  This reverse the scaling that was done in the adj_cf calculation
    mad *= pow(n,1/alpha_stable);
    k.mad.push_back(mad);
    double ci = confidence_interval(nmin, nmax, delta, pdf, ci_level);
    ci *= pow(n,1/alpha_stable);
    k.ci.push_back(ci);
    if (verbose) {
      double p_total = pdf.sum()*delta;
      cout << setw(5) << n
      << setw(13) << fixed << setprecision(8) << pdf(mod(nmin,nsize))
      << setw(13) << fixed << setprecision(8) << pdf(mod(-1,nsize))
      << setw(13) << fixed << setprecision(8) << pdf(mod(0,nsize))
      << setw(13) << fixed << setprecision(8) << pdf(mod(1,nsize))
      << setw(13) << fixed << setprecision(8) << pdf(mod(nmax,nsize))
      << setw(13) << fixed << setprecision(8) << p_total
      << setw(12) << scientific<< setprecision(3) << mad
      << setw(12) << scientific << setprecision(3) << mad_lower_tail*pow(n,1/alpha_stable)
      << setw(12) << scientific << setprecision(3) << mad_upper_tail*pow(n,1/alpha_stable)
      << endl;
    }
  } // j over ns
  k.calc_kappa();
  k.mad_rel_err = rel_err(k.mad.at(0), dist.mad());
  k.ci_rel_err = rel_err(k.ci.at(0), dist.ci(ci_level));
}

/// show proper usage after improper command line argument
void show_usage(path p     ///< [in] the path of the executable
                ) {
  cerr << "Usage: " << p.filename().string() << " delta [delta2=.001] [m=1e7]" << endl
       << "  where delta is the interval of integration/mad" << endl
       << "  and delta2 is the tail contribution to mad" << endl;
}


/// main program for convolution_testl
int main(int argc, const char * argv[]) {
  path p(argv[0]);
  if ( (argc < 2) || (argc > 4) ) {
    show_usage(p);
    return 1;
  }
  istringstream iss1{string(argv[1])};
  double delta;
  if ( !(iss1 >> delta) || ( delta <=0 ))
    show_usage(p);
  
  double delta2 = .001;
  if (argc >= 3) {
    istringstream iss2{string(argv[2])};
    if ( !(iss2 >> delta2) || (delta2 <=0 ) )
      show_usage(p);
  }
  int m = 1e7;
  if (argc >= 4) {
    istringstream iss3{string(argv[3])};
    if ( !(iss3 >> m) || (m <=0 ) )
      show_usage(p);
  }
  // Add commas to the cout format
  cout.imbue({std::locale(), new Numpunct});
  
  // for pareto and student t distribution use the same alphas as Taleb
  vector<double> alphas;
  for (size_t i=0; i<taleb_results.size(); ++i)
    alphas.push_back(taleb_results.at(i).at(0));
  
  // The durations for which MAD will be calculated
  vector<int> ns{1, 2, 30, 100};
  
  // In addtion to MAD to measure scale we'll use the length of a conficence
  // interval with ci_level confidence
  double ci_level = .05;
  
  string out_dir{"../output"};
  ostringstream oss;
  oss << out_dir << "/convolution_test_"
      << delta << "_"
      << delta2 << "_"
      << m
      << ".out";
  cout << "Writing to file " << oss.str() << endl << endl;
  cout.flush();
  
  ofstream out{oss.str()};
  out.imbue({std::locale(), new Numpunct});
  
  auto_cpu_timer t(out);
  
#ifdef EIGEN_FFTW_DEFAULT
  if (fftw_init_threads() == 0)
    throw std::runtime_error("fftw_init_threads failed");
  fftw_plan_with_nthreads(8);
#endif
  
  {
    KappaResults ks_pareto(ns, alphas.size(), "alpha", 1);
    for (size_t i=0; i<alphas.size(); ++i) {
      double alpha=alphas.at(i);
      pareto_distribution<> pd(alpha);
      KappaResult kr(ns);
      kr.param = alpha;
      calculate_kappa(delta, delta2, m, ns, pd, ci_level, kr, true);
      ks_pareto.at(i) = kr;
    }
    out << "Pareto Distribution" << endl << endl;
    out << ks_pareto;
    out.flush();
  }
 
  {
    KappaResults ks_student(ns, alphas.size(), "alpha",4);
    for (size_t i=0; i<alphas.size(); ++i) {
      double alpha = alphas.at(i);
      student_t_distribution<> td(alpha);
      KappaResult kr(ns);
      kr.param = alpha;
      calculate_kappa(delta, delta2, m, ns, td, ci_level, kr, true);
      ks_student.at(i)=kr;
    }
    out << "Student t Distribution" << endl << endl;
    out << ks_student;
  }

  
  {
    KappaResults ks_exponential(ns, 1, "lambda", 0);
    exponential_distribution<> ed(1);
    KappaResult kr(ns);
    kr.param = 1;
    calculate_kappa(delta, delta2, m, ns, ed, ci_level, kr, true);
    ks_exponential.at(0) = kr;
    out << ed << endl << endl;
    out << ks_exponential;
  }

  {
    cout << endl;
    vector<double> sigmas = {.01, .1, .5, 1, 2, 5};
    KappaResults ks_lognormal(ns, sigmas.size(), "sigma",0);
    for (size_t i=0; i<sigmas.size(); ++i) {
      double sigma = sigmas.at(i);
      int type = 3;    // Use the Lambert W approximation
      lognormal_distribution<> lnd(0,sigma, type);
      KappaResult kr(ns);
      kr.param = sigma;
      calculate_kappa(delta, delta2, m, ns, lnd, ci_level, kr, true);
      ks_lognormal.at(i)=kr;
    }
    out << "Lognormal Distribution" << endl << endl;
    out << ks_lognormal;
  }
  fftw_cleanup_threads();
  
  return 0;
}
