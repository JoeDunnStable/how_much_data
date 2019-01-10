//
//  monte_carlo_test.cpp
//  pareto_test
//
//  Created by Joseph Dunn on 12/31/18.
//  Copyright © 2018 Joseph Dunn. All rights reserved.
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
/*
#include <boost/multiprecision/mpfr.hpp>
using myFloat = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<30> >;
 */
using myFloat = double;

#include "pareto_distribution.h"
#include "student_t_distribution.h"

template<typename RealType>
RealType rel_err(RealType a, RealType b) {
  return fabs(a-b)/std::max(a,b);
}

template<typename myFloat>
vector<myFloat> quantile(const vector<myFloat> &x, const vector<myFloat>& probs)
{
  myFloat eps = 100 * std::numeric_limits<myFloat>::epsilon();
  long np = probs.size();
  for (auto p : probs) {
    if (p < -eps || p > 1 + eps)
      throw std::range_error("quantile: 'probs' outside [0,1]");
  }
  vector<myFloat> qs(np);
  long n = x.size();
  if (n > 0 && np > 0) {
    vector<myFloat> x_sort = x;
    sort(x_sort.data(), x_sort.data()+n);
    for (int j=0; j<np; ++j) {
      myFloat index = (n - 1) * probs.at(j);
      int lo = static_cast<int>(floor(index));
      int hi = static_cast<int>(ceil(index));
      myFloat h = index - lo;
      qs.at(j) = (1-h) * x_sort.at(lo) + h * x_sort.at(hi);
    }
    return qs;
  } else {
    throw std::range_error("quantile: Both x and prob must be of length greater than 0");
  }
}

mutex kr_mutex;
mutex cout_mutex;

struct KappaResult {
  double alpha;
  double mad_rel_err;
  size_t m = 0;
  size_t m_ci = 0;
  vector<int> ns;
  void initialize(double alpha_in, const vector<int>& ns_in) {
    alpha = alpha_in;
    ns = ns_in;
    sum_dev.resize(ns_in.size(),0);
    sum_abs_dev.resize(ns_in.size(),0);
    conf_int.resize(ns_in.size());
    x_ci.resize(ns_in.size());
  }
  void update_dev(size_t m, size_t m_ci, const vector<myFloat>& dev,
                  const vector<myFloat>& abs_dev,
                  const vector<vector<double> >& x_in) {
    unique_lock<mutex> lock(kr_mutex);
    this->m += m;
    this->m_ci += m_ci;
    for (size_t i=0; i<ns.size(); ++i){
      sum_dev.at(i) += dev.at(i);
      sum_abs_dev.at(i) += abs_dev.at(i);
      x_ci.at(i).insert(x_ci.at(i).end(), x_in.at(i).begin(), x_in.at(i).end());
    }
  }
  void update_conf_int(double ci_level) {
    unique_lock<mutex> lock(kr_mutex);
    vector<double> p{ci_level/2,1-ci_level/2};
    for (size_t j=0; j<ns.size(); ++j) {
      vector<double> tmp = quantile(x_ci.at(j),p);
      conf_int.at(j) = tmp.at(1)-tmp.at(0);
      x_ci.at(j).clear();
    }
  }
  vector<myFloat> sum_dev;
  vector<myFloat> sum_abs_dev;
  vector<vector<double> > x_ci;

  double ci_rel_err;
  vector<double> conf_int;
  cpu_times elapsed_time;
};

ostream& operator<< (ostream& os, const KappaResult& k) {
  os << fixed << setw(7) << setprecision(2) << k.alpha
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

struct KappaResults {
  KappaResults(const vector<int>& ns, const vector<double> &alphas)
  : ns(ns), kr(alphas.size()) {
    for (size_t i=0; i<alphas.size(); ++i) {
      kr.at(i).initialize(alphas.at(i),ns);
    }
  }
  vector<int> ns;
  vector<KappaResult> kr;
};

ostream& operator<< (ostream& os, KappaResults& ks) {
  os << setw(7) << right << "alpha"
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
  return os;
}

template<typename Dist>
void calc_kappa(unsigned int thread_id, size_t m, size_t m_ci_limit,
                vector<int> ns, Dist dist,
                double ci_level,
                KappaResult* kp,
                bool verbose = false) {
  if (verbose) {
    unique_lock<mutex> lock(cout_mutex);
    cout << "Starting thread " << thread_id
         << ", with m = " << m
         << ", and m_ci_limit = " << m_ci_limit << endl;
  }
  mt19937 eng(thread_id);
  // Only thread 1 will update confidence intervals which require less data
  double mean = dist.mean();
  
  // Up to m_ci_limit we'll save the individual run results calculate
  // confidence intervals
  vector<vector<double> > x_ci(ns.size(), vector<double>(m_ci_limit));
  
  vector<myFloat> dev(ns.size(),0.);
  vector<myFloat> abs_dev(ns.size(),0.);

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

template<typename Dist>
void calculate_kappa(size_t m, vector<int> ns, Dist dist,
                     double ci_level,
                     KappaResult* kp,
                     bool verbose = false) {
  cpu_timer timer;
  
  size_t m_ci_limit = min(static_cast<size_t>(2000000/ci_level),m);
  /*
  calc_kappa<Dist> (1, m, ns, dist, ci_level, k, verbose);
   */
  // Subdivide m between eight threads giving at least
  // m_ci_limit to thread 1, which calculated the confidence intervals
  size_t m_unassigned = m;
  size_t m_assigned = m/8;
  size_t m_ci_assigned = m_ci_limit/8;
  thread thread1(calc_kappa<Dist>, 1, m_assigned, m_ci_assigned, ns, dist, ci_level, kp, verbose);
  m_unassigned -= m_assigned;
  m_ci_limit -= m_ci_assigned;
  m_assigned = m_unassigned/7;
  m_ci_assigned = m_ci_limit/7;
  thread thread2(calc_kappa<Dist>, 2, m_assigned, m_ci_assigned, ns, dist, ci_level, kp, verbose);
  m_unassigned -= m_assigned;
  m_ci_limit -= m_ci_assigned;
  m_assigned = m_unassigned/6;
  m_ci_assigned = m_ci_limit/6;
  thread thread3(calc_kappa<Dist>, 3, m_assigned, m_ci_assigned, ns, dist, ci_level, kp, verbose);
  m_unassigned -= m_assigned;
  m_ci_limit -= m_ci_assigned;
  m_assigned = m_unassigned/5;
  m_ci_assigned = m_ci_limit/5;
  thread thread4(calc_kappa<Dist>, 4, m_assigned, m_ci_assigned, ns, dist, ci_level, kp, verbose);
  m_unassigned -= m_assigned;
  m_ci_limit -= m_ci_assigned;
  m_assigned = m_unassigned/4;
  m_ci_assigned = m_ci_limit/4;
  thread thread5(calc_kappa<Dist>, 5, m_assigned, m_ci_assigned, ns, dist, ci_level, kp, verbose);
  m_unassigned -= m_assigned;
  m_ci_limit -= m_ci_assigned;
  m_assigned = m_unassigned/3;
  m_ci_assigned = m_ci_limit/3;
  thread thread6(calc_kappa<Dist>, 6, m_assigned, m_ci_assigned, ns, dist, ci_level, kp, verbose);
  m_unassigned -= m_assigned;
  m_ci_limit -= m_ci_assigned;
  m_assigned = m_unassigned/2;
  m_ci_assigned = m_ci_limit/2;
  thread thread7(calc_kappa<Dist>, 7, m_assigned, m_ci_assigned, ns, dist, ci_level, kp, verbose);
  m_unassigned -= m_assigned;
  m_ci_limit -= m_ci_assigned;
  m_assigned = m_unassigned;
  m_ci_assigned = m_ci_limit;
  thread thread8(calc_kappa<Dist>, 8, m_assigned, m_ci_assigned, ns, dist, ci_level, kp, verbose);
  thread1.join();
  thread2.join();
  thread3.join();
  thread4.join();
  thread5.join();
  thread6.join();
  thread7.join();
  thread8.join();
  kp->update_conf_int(ci_level);

  kp->mad_rel_err = rel_err(dist.mad(), static_cast<double>(kp->sum_abs_dev.at(0)/kp->m));
  kp->ci_rel_err = rel_err(dist.ci(ci_level), kp->conf_int.at(0));
  kp->elapsed_time = timer.elapsed();
}

void show_usage(path p) {
  cerr << "Usage: " << p.filename().string() << " m" << endl
  << "  where m is number of trials" << endl;
}


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

  vector<double> alphas;
  
  for (double alpha=1.25; alpha<=4; alpha+=.25)
    alphas.push_back(alpha);
  
  vector<int> ns{1, 2};
  KappaResults ks_pareto(ns, alphas);
  KappaResults ks_student(ns, alphas);
  
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
  
  for (size_t i=0; i<alphas.size(); ++i) {
    double alpha=alphas.at(i);
    size_t m_alpha = static_cast<size_t>(pow(.001,alpha/(1-alpha)));
    m_alpha = min(max(m_ci_limit, m_alpha), m);
    if (verbose)
      cout << "Starting Pareto Distribution with alpha = " << alphas.at(i) << endl;
    pareto_distribution<> pd(alpha);
    KappaResult kr;
    kr.initialize(alpha, ns);
    calculate_kappa(m_alpha, ns, pd, ci_level, &kr, verbose);
    ks_pareto.kr.at(i) = kr;
  }
  out << "Pareto Distribution" << endl << endl;
  out << ks_pareto;
  out.flush();
  
  for (size_t i=0; i<alphas.size(); ++i) {
    double alpha = alphas.at(i);
    size_t m_alpha = static_cast<size_t>(pow(.001,alpha/(1-alpha)));
    m_alpha = min(max(m_ci_limit, m_alpha), m);
    if (verbose)
      cout << "Starting Student T Distribution with alpha = " << alphas.at(i) << endl;
    student_t_distribution<> td(alpha);
    KappaResult kr;
    kr.initialize(alpha,ns);
    calculate_kappa(m_alpha, ns, td, ci_level, &kr, verbose);
    ks_student.kr.at(i) = kr;
  }
  out << "Student Distribution" << endl << endl;
  out << ks_student;
  return 0;
}
