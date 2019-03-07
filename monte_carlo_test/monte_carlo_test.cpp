//
/// \file monte_carlo_test.cpp
/// \package how_much_data
//
/// \author Created by Joseph Dunn on 12/31/18.
/// \copyright © 2018, 2019 Joseph Dunn. All rights reserved.
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
#include <string>
using std::string;
using std::to_string;
#include <sstream>
using std::istringstream;
using std::ostringstream;
#include <fstream>
using std::ofstream;
#include <random>
using std::mt19937;
#include <vector>
using std::vector;
#include <array>
using std::array;
#include <algorithm>
using std::max;
using std::min;
using std::max_element;
using std::sort;
#include <numeric>
using std::accumulate;
#include <mutex>
using std::mutex;
using std::unique_lock;
#include <thread>
using std::thread;

#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;
using boost::timer::cpu_timer;
using boost::timer::cpu_times;
using boost::timer::nanosecond_type;
#include <boost/filesystem.hpp>
using boost::filesystem::path;
#include <boost/math/special_functions/beta.hpp>
using boost::math::beta;
#include <boost/math/distributions/students_t.hpp>
using boost::math::students_t;

#include "adaptive_integration.h"
using namespace adaptive_integration;
Kronrod<double> k_big(10);
int noext = 0;
double epsabs_double = 0;
double epsrel_double = 64 * std::numeric_limits<double>::epsilon();
int limit = 2000;
int verbose_integration = 0;

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
#include "normal_switch_mean.h"
#include "normal_switch_stddev.h"
#include "taleb_results.h"

/// return the relative error between two numbers
template<typename RealType>
RealType rel_err(RealType a,   ///< [in] the first number
                 RealType b    ///< [in] the second number
                 ) {
  return fabs(a-b)/std::max(a,b);
}

/// return quantiles of the ensemble of trials
template<typename RealType>
vector<RealType> quantile(const vector<RealType> &x,   ///< [in] the result by trial
                         const vector<RealType>& probs ///< [in] the desired probabilities
                         )
{
  RealType eps = 100 * std::numeric_limits<RealType>::epsilon();
  long np = probs.size();
  for (auto p : probs) {
    if (p < -eps || p > 1 + eps)
      throw std::range_error("quantile: 'probs' outside [0,1]");
  }
  vector<RealType> qs(np);
  long n = x.size();
  if (n > 0 && np > 0) {
    vector<RealType> x_sort = x;
    sort(x_sort.data(), x_sort.data()+n);
    for (int j=0; j<np; ++j) {
      RealType index = (n - 1) * probs.at(j);
      int lo = static_cast<int>(floor(index));
      int hi = static_cast<int>(ceil(index));
      RealType h = index - lo;
      qs.at(j) = (1-h) * x_sort.at(lo) + h * x_sort.at(hi);
    }
    return qs;
  } else {
    throw std::range_error("quantile: Both x and prob must be of length greater than 0");
  }
}

/// a mutex for writing to KappaResults
mutex kr_mutex;
/// a mutex for writing to cout
mutex cout_mutex;

/// the struct holding the resutls of a single alpha
struct KappaResult {
  double param;         ///< the parameter for the run
  double mad_rel_err;   ///< the relative error of the mad vs theory
  size_t m = 0;         ///< the number of trials run for mad
  size_t m_ci = 0;      ///< the number of trials for the conf. interval
  vector<int> ns;       ///< the durations calculated
  /// initialize the structure
  void initialize(double param_in,         ///< the alpha of the run
                  const vector<int>& ns_in ///< the durations of the run
                  ) {
    param = param_in;
    ns = ns_in;
    sum_dev.resize(ns_in.size(),0);
    sum_abs_dev.resize(ns_in.size(),0);
    conf_int.resize(ns_in.size());
    x_ci.resize(ns_in.size());
  }
  /// update the deviations w result from one thread
  void update_dev(size_t m,                 ///< the number of trials for mad
                  size_t m_ci,              ///< the number of trials for ci
                  const vector<double>& dev, ///< sum of deviation by duration
                  const vector<double>& abs_dev, ///<  sum of abs dev by duration
                  const vector<vector<double> >& x_in ///< the variate by trial
                  ) {
    unique_lock<mutex> lock(kr_mutex);
    this->m += m;
    this->m_ci += m_ci;
    for (size_t i=0; i<ns.size(); ++i){
      sum_dev.at(i) += dev.at(i);
      sum_abs_dev.at(i) += abs_dev.at(i);
      x_ci.at(i).insert(x_ci.at(i).end(), x_in.at(i).begin(), x_in.at(i).end());
    }
  }
  /// calculate the confidence interval for all trials
  void update_conf_int(double ci_level     ///< the confidence level to use
                       ) {
    unique_lock<mutex> lock(kr_mutex);
    vector<double> p{ci_level/2,1-ci_level/2};
    for (size_t j=0; j<ns.size(); ++j) {
      vector<double> tmp = quantile(x_ci.at(j),p);
      conf_int.at(j) = tmp.at(1)-tmp.at(0);
      x_ci.at(j).clear();
    }
  }
  vector<double> sum_dev;        ///< sum or raw deviation by duration
  vector<double> sum_abs_dev;    ///< sum of abs deviations by duration
  vector<vector<double> > x_ci;   ///< the results by trial and duration

  double ci_rel_err;              ///< the relative error of the ci vs theory
  vector<double> conf_int;        ///< the calculated ci by duration
  cpu_times elapsed_time;         ///< the elapsed time for the run
};

/// output the results from a single run
ostream& operator<< (ostream& os,          ///< [in,out] the output stream
                     const KappaResult& k  ///< [in] the sturct with results
                     ) {
  os << fixed << setw(7) << setprecision(2) << k.param
     << setw(13) << k.m
  << setw(12) << setprecision(2) << k.mad_rel_err*100 << "%";
  double mad0 = static_cast<double>(k.sum_abs_dev.at(0)/k.m);
  for (size_t i=1; i<k.ns.size(); ++i) {
    double madn = static_cast<double>(k.sum_abs_dev.at(i)/k.m);
    double kappa_mad = 2 - (log(k.ns.at(i))-log(k.ns.at(0)))
                           /(log(madn)-log(mad0));
    os << setw(13) << setprecision(3) << kappa_mad;
  }
  os << setw(10) << k.m_ci
     << setw(11) << setprecision(2) << k.ci_rel_err*100 << "%";
  double ci0 = k.conf_int.at(0);
  for (size_t i=1; i<k.ns.size(); ++i){
    double cin = k.conf_int.at(i);
    double kappa_ci = 2 -(log(k.ns.at(i))-log(k.ns.at(0)))/
                         (log(cin) - log(ci0));
    os << setw(13) << setprecision(3) << kappa_ci;
  }
  os << format(k.elapsed_time,1,"   %ws wall, %ts CPU (%p%)");
  os << endl;
  return os;
}

/// sturcture holding the results from all runs
struct KappaResults {
  /// constructor
  KappaResults(const vector<int>& ns,       ///< the durations calculated
               const vector<double> &params,///< the params for each run
               const string param_label,    ///< the label for the param
               size_t taleb_offset          ///< the column offset into the table
                                            ///< of taleg's results. 0 for none
               )
  : ns(ns), param_label(param_label), kr(params.size()),
    taleb_offset(taleb_offset){
    for (size_t i=0; i<params.size(); ++i) {
      kr.at(i).initialize(params.at(i),ns);
    }
  }
  /// the durations calculated
  vector<int> ns;
  /// a vector holding the results of each run
  /// the label for the parameter
  string param_label;
  vector<KappaResult> kr;
  /// the offset into Taleb's table of results.  =0 for none available.
  size_t taleb_offset;
};

/// output the results for all runs
ostream& operator<< (ostream& os, KappaResults& ks) {
  os << setw(7) << right << ks.param_label
  << setw(13) << right << "m_mad"
  << setw(13) << right << "mad_rel_err";
  for (size_t j=1; j<ks.ns.size(); ++j)
    os << setw(13) << right << "kappa"+to_string(ks.ns.at(j))+"_mad";
  os << setw(10) << right << "m_ci"
     << setw(12) << right << "ci_rel_err";
  for (size_t j=1; j<ks.ns.size(); ++j)
    os << setw(13) << right << "kappa"+to_string(ks.ns.at(j))+"_ci";
  os << endl << endl;
  for (auto& kr : ks.kr)
    os << kr;
  os << endl;
  if (ks.taleb_offset > 0) {
    os << setw(68) << right << "Error vs Taleb's Results" << endl << endl;
    for (size_t i=0; i<ks.kr.size(); ++i) {
      os << setw(7) << setprecision(2) << ks.kr.at(i).param << setw(26) << " ";
      double mad0 = static_cast<double>(ks.kr.at(i).sum_abs_dev.at(0)/ks.kr.at(i).m);
      for (size_t j=1; j<ks.ns.size(); ++j) {
        double madn = static_cast<double>(ks.kr.at(i).sum_abs_dev.at(j)/ks.kr.at(i).m);
        double kappa_mad = 2 - (log(ks.kr.at(i).ns.at(j))-log(ks.kr.at(i).ns.at(0)))
        /(log(madn)-log(mad0));
        double err = kappa_mad - taleb_results.at(i).at(j-1+ks.taleb_offset);
        os << setw(13) << setprecision(3) << err;
      }
      os << endl;
    }
    os << endl;
  }
  return os;
}

/// the per thread cacluaiton engine
template<typename Dist>
void calc_kappa(unsigned int thread_id,    ///< [in] the number of the thread
                                           ///< used as seed for urng
                size_t m,                  ///< [in] the maximum # of trials for mad
                size_t m_ci_limit,         ///< [in] the maximum # of trials for ci
                vector<int> ns,            ///< [in] the durations to save
                Dist dist,                 ///< [in] the distribution
                double ci_level,           ///< [in] the confidence level for kappa_ci
                KappaResult* kp,           ///< [out]the results
                bool verbose = false       ///< [in] flag for trace infomation
                ) {
  if (verbose) {
    unique_lock<mutex> lock(cout_mutex);
    cout << "Starting thread " << thread_id
         << ", with m = " << m
         << ", and m_ci_limit = " << m_ci_limit << endl;
  }
  mt19937 eng(thread_id);
  double mean = dist.mean();
  
  // Up to m_ci_limit we'll save the individual run results calculate
  // confidence intervals
  vector<vector<double> > x_ci(ns.size(), vector<double>(m_ci_limit));
  
  vector<double> dev(ns.size(),0.);
  vector<double> abs_dev(ns.size(),0.);

  for (size_t i=0; i<m; ++i) {
    int n_old = 0;
    double sum = 0;
    for (size_t j=0; j<ns.size(); j++) {
      int n = ns.at(j);
      for (int jj=n_old; jj<n; ++jj)
        sum += dist(eng);
      if ( i<m_ci_limit )
        x_ci.at(j).at(i) =sum;
      dev.at(j) += sum;
      abs_dev.at(j) += fabs(sum-n*mean);
      n_old = n;
    }
  } // i
  kp->update_dev(m, m_ci_limit, dev, abs_dev, x_ci);

  return;
}

/// calculate kappa for sums of iid variables at specified durations
template<typename Dist>
void calculate_kappa(size_t m,               ///< [in] the number of scenarios
                     vector<int> ns,         ///< [in] the durations to save
                     Dist dist,              ///< [in] the distribution
                     double ci_level,        ///< [in] the confidence level to use
                     KappaResult* kp,        ///< [out] a ptr to the results
                     bool verbose = false    ///< [in] a flag to generate trace
                     )
{
  cpu_timer timer;
  
  size_t m_ci_limit = min(static_cast<size_t>(2000000/ci_level),m);
                       
  // Subdivide m between eight threads
  array<thread, 8> threads;
  size_t m_unassigned = m;
  size_t m_ci_unassigned = m_ci_limit;
  for (int i = 0; i<static_cast<int>(threads.size()); ++i) {
    int threads_to_go = static_cast<int>(threads.size())-i;
    size_t m_assigned = m_unassigned/threads_to_go;
    size_t m_ci_assigned = m_ci_unassigned/threads_to_go;
    threads.at(i)=thread(calc_kappa<Dist>, i+1, m_assigned, m_ci_assigned,
                         ns, dist, ci_level, kp, verbose);
    m_unassigned -= m_assigned;
    m_ci_unassigned -= m_ci_assigned;
  }
  for (size_t i=0; i<threads.size(); ++i)
    threads.at(i).join();
  kp->update_conf_int(ci_level);

  kp->mad_rel_err = rel_err(dist.mad(), static_cast<double>(kp->sum_abs_dev.at(0)/kp->m));
  kp->ci_rel_err = rel_err(dist.ci(ci_level), kp->conf_int.at(0));
  kp->elapsed_time = timer.elapsed();
}

/// show the usage.  called when the wrong # of arguments is used
void show_usage(path p) {
  cerr << "Usage: " << p.filename().string() << " m" << endl
  << "  where m is number of trials" << endl;
}


/// main program for convolution test with one input parameter the # of trials
int main(int argc, const char * argv[]) {
  path p(argv[0]);
  if ( argc != 2) {
    show_usage(p);
    return 1;
  }
  size_t m;
  istringstream iss1{argv[1]};
  if (!(iss1 >> m) || m < 1) {
    show_usage(p);
    return 1;
  }

  // Needed for the lognormal characteristic function
  Kronrod<double> k_big(10);
  int noext = 0;
  double epsabs_double = std::numeric_limits<double>::epsilon()/128;
  double epsrel_double = boost::math::tools::root_epsilon<double>();
  int limit = 1000;
  int verbose_integration = 0;
  IntegrationController<double> cf_ctl(noext, k_big,
                                       epsabs_double, epsrel_double,
                                       limit, verbose_integration);

  vector<double> alphas;
  
  for (double alpha=1.25; alpha<=4; alpha+=.25)
    alphas.push_back(alpha);
  
  vector<int> ns{1, 2, 30, 100};
  
  double ci_level = .05;
  size_t m_ci_limit = min(static_cast<size_t>(2000000/ci_level), m);

  string out_dir{"../output"};
  ostringstream oss;
  oss << out_dir << "/monte_carlo_test_"  << m << ".out";
  cout << "Writing to file " << oss.str() << endl << endl;
  cout.flush();
  
  ofstream out{oss.str()};
  auto_cpu_timer t(out);
  out << "m =     " << m << endl
      << "ns =     ";
  for (auto n : ns)
    out << n << " ";
  out << endl<< endl;
  out.flush();
  bool verbose = true;
  {
    KappaResults ks_pareto(ns, alphas, "alpha", 1);
    for (size_t i=0; i<alphas.size(); ++i) {
      pareto_distribution<> pd(alphas.at(i));
      // m_alpha is the number of samples needed to reduce the
      // uncertainty in the average of m_alpha stable variates to
      // .001 of the uncertainty of a single sample
      double alpha_cap = log(2)/(log(pd.mad2())-log(pd.mad()));
      size_t m_alpha = static_cast<size_t>(pow(.001,alpha_cap/(1-alpha_cap)));
      m_alpha = min(max(m_ci_limit, m_alpha), m);
      if (verbose)
        cout << pd << endl;
      KappaResult kr;
      kr.initialize(alphas.at(i), ns);
      calculate_kappa(m_alpha, ns, pd, ci_level, &kr, verbose);
      ks_pareto.kr.at(i) = kr;
    }
    out << "Pareto Distribution" << endl << endl;
    out << ks_pareto;
  }
  
  {
    KappaResults ks_student(ns, alphas, "alpha", 4);
    for (size_t i=0; i<alphas.size(); ++i) {
      student_t_distribution<> td(alphas.at(i));
      // m_alpha is the number of samples needed to reduce the
      // uncertainty in the average of m_alpha stable variates to
      // .001 of the uncertainty of a single sample
      double alpha_cap = log(2)/(log(td.mad2())-log(td.mad()));
      size_t m_alpha = static_cast<size_t>(pow(.001,alpha_cap/(1-alpha_cap)));
      m_alpha = min(max(m_ci_limit, m_alpha), m);
      if (verbose)
        cout << td << endl;
      KappaResult kr;
      kr.initialize(alphas.at(i),ns);
      calculate_kappa(m_alpha, ns, td, ci_level, &kr, verbose);
      ks_student.kr.at(i) = kr;
    }
    out << "Student Distribution" << endl << endl;
    out << ks_student;
  }

  {
    vector<double> lambdas = {1.};
    KappaResults ks_exponential(ns, lambdas, "lambda" , 0);
    exponential_distribution<> ed(1);
    // m_alpha is the number of samples needed to reduce the
    // uncertainty in the average of m_alpha stable variates to
    // .001 of the uncertainty of a single sample
    double alpha_cap = log(2)/(log(ed.mad2())-log(ed.mad()));
    size_t m_alpha = static_cast<size_t>(pow(.001,alpha_cap/(1-alpha_cap)));
    m_alpha = min(max(m_ci_limit, m_alpha), m);
    if (verbose)
      cout << ed << endl;
    KappaResult kr;
    kr.initialize(1,ns);
    calculate_kappa(m_alpha, ns, ed, ci_level, &kr, verbose);
    ks_exponential.kr.at(0) = kr;
    out << ed << endl;
    out << ks_exponential << endl;
  }

  {
    vector<double> sigmas = {.1, .2, 1, 5, 10};
    KappaResults ks_lognormal(ns, sigmas, "sigma", 0);
    for (size_t i=0; i<sigmas.size(); ++i) {
      lognormal_distribution<> lnd(0,sigmas.at(i), cf_ctl);
      // m_alpha is the number of samples needed to reduce the
      // uncertainty in the average of m_alpha stable variates to
      // .001 of the uncertainty of a single sample
      double alpha_cap = log(2)/(log(lnd.mad2())-log(lnd.mad()));
      size_t m_alpha = static_cast<size_t>(pow(.001,alpha_cap/(1-alpha_cap)));
      m_alpha = min(max(m_ci_limit, m_alpha), m);
      if (verbose)
        cout << lnd << endl;
      KappaResult kr;
      kr.initialize(sigmas.at(i),ns);
      calculate_kappa(m_alpha, ns, lnd, ci_level, &kr, verbose);
      ks_lognormal.kr.at(i) = kr;
    }
    out << "Lognormal Distribution" << endl << endl;
    out << ks_lognormal;
  }
  
  {
    vector<double> ds = {0, 1, 2, 3, 4, 5};
    KappaResults ks_normal_switch_mean(ns, ds, "d" , 0);
    for(size_t j=0; j<ds.size(); ++j) {
      double d = ds.at(j);
      normal_switch_mean<> nsm(d);
      // m_alpha is the number of samples needed to reduce the
      // uncertainty in the average of m_alpha stable variates to
      // .001 of the uncertainty of a single sample
      double alpha_cap = log(2)/(log(nsm.mad2())-log(nsm.mad()));
      size_t m_alpha = static_cast<size_t>(pow(.001,alpha_cap/(1-alpha_cap)));
      m_alpha = min(max(m_ci_limit, m_alpha), m);
      if (verbose)
        cout << nsm << endl;
      KappaResult kr;
      kr.initialize(ds.at(j),ns);
      calculate_kappa(m_alpha, ns, nsm, ci_level, &kr, verbose);
      ks_normal_switch_mean.kr.at(j) = kr;
    }
    out << "Normal with Switching Mean But Common sigma" << endl << endl;
    out << ks_normal_switch_mean << endl;
  }
  
  {
    vector<double> as = {1, 2, 3, 4, 5};
    KappaResults ks_normal_switch_stddev(ns, as, "a" , 0);
    for(size_t j=0; j<as.size(); ++j) {
      double a = as.at(j);
      normal_switch_stddev<> nss(a, .1);
      // m_alpha is the number of samples needed to reduce the
      // uncertainty in the average of m_alpha stable variates to
      // .001 of the uncertainty of a single sample
      double alpha_cap = log(2)/(log(nss.mad2())-log(nss.mad()));
      size_t m_alpha = static_cast<size_t>(pow(.001,alpha_cap/(1-alpha_cap)));
      m_alpha = min(max(m_ci_limit, m_alpha), m);
      if (verbose)
        cout << nss << endl;
      KappaResult kr;
      kr.initialize(as.at(j),ns);
      calculate_kappa(m_alpha, ns, nss, ci_level, &kr, verbose);
      ks_normal_switch_stddev.kr.at(j) = kr;
    }
    out << "Normal with Switching Stddev w Mean = 0" << endl << endl;
    out << ks_normal_switch_stddev << endl;
  }
  
  return 0;
}
