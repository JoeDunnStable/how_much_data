/// \file adaptive_integration_impl.h
/// Implementation of the classes declared in adaptive_integration.h
/// Included in adaptive_integration.h when LIBRARY is defined.
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <limits>
#include <algorithm>
#include "gauss_kronrod.h"
#include <stdio.h>
#include <iomanip>

namespace adaptive_integration {
  
using std::endl;
using std::setw;
using std::setprecision;
using std::right;
using std::sort;
using std::make_heap;
using std::pop_heap;
using std::push_heap;
using std::max;
using std::min;
using namespace gauss_kronrod;
  

template<typename T>
class sort_data {
public:
    T a;
    int index;
    sort_data(T a, int index) : a(a), index(index){}
};

template<typename T>
bool sd_lt(const sort_data<T> &lhs,const sort_data<T>& rhs) {
    return lhs.a < rhs.a;
}

template<typename myFloat>
ostream& operator<< (ostream& os, const IntegrationController<myFloat>& ctl) {
  os << setw(15) << "noext = " << ctl.noext << endl
     << setw(15) << "epsabs = " << ctl.epsabs << endl
     << setw(15) << "epsrel = " << ctl.epsrel << endl
     << setw(15) << "limit = " << ctl.limit << endl
     << ctl.g_k << endl;
  return os;
    
}

template<typename myFloat>
void print_subs_summary(ostream& os, const vector<Subinterval<myFloat> >& subs, const int last, const vector<myFloat>& points) {
  if (last==0) return;
  if (subs.size()<last) {
    throw std::range_error("print_subs: last is greater than the length of subs.");
  }
  // Determine the geometric order of the subintervals
  
  vector<sort_data<myFloat> > srt_a;
  for (int i=0; i<last; i++) {
    sort_data<myFloat> elem(subs.at(i).a,i);
    srt_a.push_back(elem);
  }
  sort(srt_a.begin(), srt_a.end(),sd_lt<myFloat>);
  
  cout << setw(13) << right << "a"
  << setw(13) << right << "b"
  << setw(13) << right << "length"
  << setw(13) << right << "r"
  << setw(13) << right << "average"
  << setw(13) << right << "e"
  << setw(13) << right << "rabs"
  << setw(13) << right << "defabs"
  << setw(6) << right << "count" << endl;
  int j = 0;
  for (int i = 0; i < points.size()-1; i++) {
    myFloat a = points.at(i);
    myFloat b = points.at(i+1);
    
    myFloat r(0), e(0), rabs(0), defabs(0);
    int count(0);
    for (;j < last ; j++) {
      int k = srt_a.at(j).index;
      if ( subs.at(k).b > b ) break;
      r+= subs.at(k).r;
      e+= subs.at(k).e;
      rabs+= subs.at(k).rabs;
      defabs += subs.at(k).defabs;
      ++count;
    }
    cout << setw(13) << setprecision(5) << a
    << setw(13) << setprecision(5) << b
    << setw(13) << setprecision(5) << b-a
    << setw(13) << setprecision(5) << r
    << setw(13) << setprecision(5) << r/(b-a)
    << setw(13) << setprecision(5) << e
    << setw(13) << setprecision(5) << rabs
    << setw(13) << setprecision(5) << defabs
    << setw(6)  << count << endl;
  }
}

template<typename myFloat>
void print_subs(ostream& os, const vector<Subinterval<myFloat> >& subs, const int last, const vector<myFloat>& points) {
  if (last==0) return;
  if (subs.size()<last) {
    throw std::range_error("print_subs: last is greater than the length of subs.");
  }
  // Determine the geometric order of the subintervals
  
  vector<sort_data<myFloat> > srt_a;
  for (int i=0; i<last; i++) {
    sort_data<myFloat> elem(subs.at(i).a,i);
    srt_a.push_back(elem);
  }
  sort(srt_a.begin(), srt_a.end(),sd_lt<myFloat>);
  
  // Determine the error ranks for each subinterval
  
  vector<sort_data<myFloat> > srt_eord;
  for (int i=0; i<last; i++) {
    sort_data<myFloat> elem(subs.at(i).e,i);
    srt_eord.push_back(elem);
  }
  sort(srt_eord.begin(), srt_eord.end(),sd_lt<myFloat>);
  
  vector<sort_data<int> > srt_iord;
  for (int i=0; i<last; i++) {
    sort_data<int> elem(srt_eord.at(i).index,last-1-i);
    srt_iord.push_back(elem);
  }
  sort(srt_iord.begin(), srt_iord.end(),sd_lt<int>);
  
  os << " "
  << setw(13) << right << "a"
  << setw(13) << right << "b"
  << setw(13) << right << "length"
  << setw(13) << right << "r"
  << setw(13) << right << "average"
  << setw(13) << right << "e"
  << setw(13) << right << "rabs"
  << setw(13) << right << "defabs"
  << setw(5) << right << "rank"
  << setw(6) << right << "level" << endl;
  for (int i=0; i<last;i++){
    int j = srt_a[i].index;
    
    // Determine whether the left endpoint is in points.
    bool ispt = i==0;
    for (int ipt = 0; ipt<points.size(); ipt++) {
      ispt = ispt || subs.at(j).a==points.at(ipt);
      if (ispt) break;
    }
    if (ispt)
      cout << "*";
    else
      cout << " ";
    cout << setw(13) << setprecision(5) << subs.at(j).a
    << setw(13) << setprecision(5) << subs.at(j).b
    << setw(13) << setprecision(5) << subs.at(j).b-subs.at(j).a
    << setw(13) << setprecision(5) << subs.at(j).r
    << setw(13) << setprecision(5) << subs.at(j).r/(subs.at(j).b-subs.at(j).a)
    << setw(13) << setprecision(5) << subs.at(j).e
    << setw(13) << setprecision(5) << subs.at(j).rabs
    << setw(13) << setprecision(5) << subs.at(j).defabs
    << setw(5) << srt_iord.at(j).index+1
    << setw(6)  << subs.at(j).level << endl;
  }
}

/*
myFloat eps_mat[50][50];
*/

/// accelerate convergence by means of an epsilon table
///
/// the routine determines the limit of a given sequence of
/// approximations, by means of the epsilon algorithm of
/// p.wynn. an estimate of the absolute error is also given.
/// the condensed epsilon table is computed. only those
/// elements needed for the computation of the next diagonal
/// are preserved.
///
/// @author Joseph Dunn based on the original Fortran code of
/// @author piessens,robert,appl. math. & progr. div. - k.u.leuven
/// @author de doncker,elise,appl. math & progr. div. - k.u.leuven
///
template<typename myFloat>
class EpsilonTable {
public:
    int n;                         ///< index of new element in the first column
                                   ///< of the epsilon table
    vector<myFloat> epsilon_table;  ///< the elements of the two lower diagonals of
                                   ///< the triangular epsilon table. the elements are numbered
                                   ///< starting at the right-hand corner of the triangle.
    myFloat result;                 ///< the result of the extrapolation
    myFloat abserr;                 ///< estimate of absolue error of approximation based on
                                   ///< on result and the three previous results
    vector<myFloat> res3la;         ///< vector containing three previous results
    int nres;                      ///< number of calls to the routine.
    /// update the table with a new result 
    void update(myFloat new_result);
    EpsilonTable(int size) : n(0), epsilon_table(size), result(NAN), res3la(3), nres(0){}
};

template<typename myFloat>
void EpsilonTable<myFloat>::update(myFloat new_result) {
  //
  //           list of major variables
  //           -----------------------
  //
  //           e0     - the 4 elements on which the computation of a new
  //           e1       element in the epsilon table is based
  //           e2
  //           e3                 e0
  //                        e3    e1    new
  //                              e2
  //           newelm - number of elements to be computed in the new
  //                    diagonal
  //           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
  //           result - the element in the new diagonal with least value
  //                    of error
  //           limexp   is the maximum number of elements the epsilon
  //                    table can contain. if this number is reached, the upper
  //                    diagonal of the epsilon table is deleted.
  //
  //***first executable statement update
    
  myFloat oflow = std::numeric_limits<myFloat>::max();
  myFloat epmach = std::numeric_limits<myFloat>::epsilon();
  epsilon_table.at(n) = new_result;
  ++n;
  nres++;
  abserr = oflow;
  result = epsilon_table.at(n-1);
  if (n < 3) {
    abserr = std::max<myFloat>(abserr, 5 * epmach * fabs(result));
    return;
  }
  int limexp = static_cast<int>(epsilon_table.size())-2;
  epsilon_table.at(n+1) = epsilon_table.at(n-1);
  int newelm = (n - 1) / 2;
  epsilon_table.at(n-1) = oflow;
  int num = n;
  int k1 = n;
  for (int i=1; i<=newelm; i++) {
    int k2 = k1 - 1;
    int k3 = k1 - 2;
    myFloat res = epsilon_table.at(k1+1);
    myFloat e0 = epsilon_table.at(k3-1);
    myFloat e1 = epsilon_table.at(k2-1);
    myFloat e2 = res;
    myFloat e1abs = fabs(e1);
    myFloat delta2 = e2 - e1;
    myFloat err2 = fabs(delta2);
    myFloat tol2 = std::max<myFloat>(fabs(e2), e1abs) * epmach;
    myFloat delta3 = e1 - e0;
    myFloat err3 = fabs(delta3);
    myFloat tol3 = std::max<myFloat>(e1abs, fabs(e0)) * epmach;
    if (err2 <= tol2 && err3 <= tol3) {
      //
      //           if e0, e1 and e2 are equal to within machine
      //           accuracy, convergence is assumed.
      //           result = e2
      //           abserr = abs(e1-e0)+abs(e2-e1)
      //
      result = res;
      abserr = std::max<myFloat>(err2 + err3, 50 * epmach * fabs(result));
      return;
    }
    myFloat e3 = epsilon_table.at(k1-1);
    epsilon_table.at(k1-1) = e1;
    myFloat delta1 = e1 - e3;
    myFloat err1 = fabs(delta1);
    myFloat tol1 = std::max<myFloat>(e1abs, fabs(e3)) * epmach;
    //
    //           if two elements are very close to each other, omit
    //           a part of the table by adjusting the value of n
    //
    if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) {
      n = i + i - 1;
      break;
    }
    myFloat ss = 1 / delta1 + 1 / delta2 - 1 / delta3;
    myFloat epsinf = fabs(ss * e1);
    //
    //           test to detect irregular behaviour in the table, and
    //           eventually omit a part of the table adjusting the value
    //           of n.
    //
    if (epsinf <= 0.1e-03) {
      n = i + i - 1;
      break;
    }
    //
    //           compute a new element and eventually adjust
    //           the value of result.
    //
    res = e1 + 1 / ss;
    epsilon_table.at(k1-1) = res;
    k1 = k1 - 2;
    myFloat error = err2 + fabs(res - e2) + err3;
    if (error <= abserr) {
      abserr = error;
      result = res;
    }
  }
  //
  //           shift the table.
  //
  if (n == limexp) {
    n = 2 * (limexp / 2) - 1;
  }
  int ib = 1;
  if ((num / 2) * 2 == num) {
    ib = 2;
  }
  int ie = newelm + 1;
  for (int i=1; i<=ie; i++) {
    int ib2 = ib + 2;
    epsilon_table.at(ib-1) = epsilon_table.at(ib2-1);
    ib = ib2;
  }
  if (num != n) {
    int indx = num - n + 1;
    for (int i=1; i<=n; i++) {
      epsilon_table.at(i-1) = epsilon_table.at(indx-1);
      indx++;
    }
/*

    for (int i=0; i<num-2; i++) {
      for (int j=0; j<=num/2; j++)
        if (j>i+1)
          eps_mat[i][j]=0;
        else
          eps_mat[i][j]=eps_mat[i+num-n][j];
    }
*/
  }

  if (nres >= 4) {
    abserr = fabs(result - res3la.at(2)) + fabs(result - res3la.at(1)) + fabs(result - res3la.at(0));
    res3la.at(0) = res3la.at(1);
    res3la.at(1) = res3la.at(2);
    res3la.at(2) = result;
  } else{
    res3la.at(nres-1) = result;
    abserr = oflow;
  }
  //
  //           compute error estimate
  //
  abserr = std::max<myFloat>(abserr, 50 * epmach * fabs(result));
/*
  cout << "EpsilonTable::update: "
            << "result = " << setw(20) << setprecision(15) << result
            << ", abserr = " << setw(20) << setprecision(15) << abserr
            << endl;


  if (n==3) {
    for (int i=0; i<50; i++)
      for (int j=0; j<50; j++)
        eps_mat[i][j]=0.;
    eps_mat[0][1]=epsilon_table[0];
    eps_mat[0][0]=epsilon_table[1];
    eps_mat[1][0]=epsilon_table[2];
  } else {
    int row = (n-2);  //base zero
    int col = 0;  //base zero
    for (int i=n-1; i>=0; i-=2, row--, col++) {
      eps_mat[row][col]=epsilon_table[i];
    }
  }

  for (int row = 0; row <(n-1);row++){
    for (int col=0; col<(n+1)/2; col++){
      if (eps_mat[row][col]==0) break;
      cout << setw(11) << setprecision(5) << result-eps_mat[row][col];
    }
    cout << endl;
  }
*/
}

template<typename myFloat>
template<typename F>
void Subinterval<myFloat>::integrate(F& f,
                       const int n_gauss,
                       const std::vector<myFloat>& w_gauss,
                       const std::vector<myFloat>& x_kronrod,
                       const std::vector<myFloat>& w_kronrod,
                       std::vector<myFloat>& fv)
{
    //           list of major variables
    //           -----------------------
    //
    //           centr  - mid point of the interval
    //           hlgth  - half-length of the interval
    //           absc   - abscissa
    //           fval*  - function value
    //           resg   - result of the n_gauss point gauss formula
    //           resk   - result of the 2n_gauss-1 point kronrod formula
    //           reskh  - approximation to the mean value of f over (a,b),
    //                    i.e. to i/(b-a)
    //
    //           machine dependent constants
    //           ---------------------------
    //
    //           epmach is the largest relative spacing.
    //           uflow is the smallest positive magnitude.
    //
    //***first executable statement  integrate_myFloat
    //"
    
    myFloat epmach = std::numeric_limits<myFloat>::epsilon();
    myFloat uflow = std::numeric_limits<myFloat>::min();
    myFloat centr =  (a + b) / 2;
    myFloat hlgth = (b - a) / 2;
    myFloat dhlgth = fabs(hlgth);
    //
    //           compute the 2n_gauss+1-point kronrod approximation to
    //           the integral, and estimate the absolute error.
    //
    myFloat resg = 0.;
    myFloat resk = 0.;
    rabs = 0.;
    myFloat fval = 0.;

    for (int j=0; j<2*n_gauss+1; j++) {
        fv[j]=centr+hlgth * x_kronrod[j];
    }
    f(fv);
    for (int j=0; j<2*n_gauss+1; j++) {
        fval = fv[j];
        resk += w_kronrod[j] * fval;
        rabs += w_kronrod[j] * fabs(fval);
    }
    for (int j=0; j<n_gauss; j++) {
        fval = fv[2*j+1];
        resg += w_gauss[j]*fval;
    }
    myFloat reskh = resk / 2;
    defabs = 0;
    for (int j=0; j<2*n_gauss+1 ; j++) {
        defabs += w_kronrod[j] * fabs(fv[j] - reskh);
    }
    r = resk * hlgth;
    rabs = rabs * dhlgth;
    defabs = defabs * dhlgth;
    e = fabs((resk - resg) * hlgth);
    if (defabs != myFloat(0) && e != myFloat(0)) {
        e = defabs * min<myFloat>(static_cast<myFloat>(1), pow((200 * e / defabs), static_cast<myFloat>(1.5)));
    }
    if (rabs > uflow / (50 * epmach)) {
        e = std::max<myFloat>((epmach * 50) * rabs, e);
    }
}

template<typename myFloat>
template<typename F>
void IntegrationController<myFloat>::integrate(F& f,
                                 const vector<myFloat>& points,
                                 myFloat& result,
                                 myFloat& abserr,
                                 int& neval,
                                 TerminationCode& termination_code,
                                 int& last)
{
  //            list of major variables
  //            -----------------------
  //
  //           maxerr    - pointer to the interval with largest error
  //                       estimate
  //           errmax    - the error at maxerr
  //           erlast    - error on the interval currently subdivided
  //                       (before that subdivision has taken place)
  //           area      - sum of the integrals over the subintervals
  //           errsum    - sum of the errors over the subintervals
  //           errbnd    - requested accuracy max(epsabs,epsrel*
  //                       abs(result))
  //           *****1    - variable for the left subinterval
  //           *****2    - variable for the right subinterval
  //           last      - index for subdivision
  //           nres      - number of calls to the extrapolation routine
  //           erlarg    - sum of the errors over the intervals larger
  //                       than the smallest interval considered up to now
  //           keep_level_max    - logical variable denoting that the routine
  //                       is attempting to perform extrapolation. i.e.
  //                       before subdividing the smallest interval we
  //                       try to decrease the value of erlarg.
  //           noext     - logical variable denoting that extrapolation is
  //                       no longer allowed (true-value)
  //
  //            machine dependent constants
  //            ---------------------------
  //
  //           epmach is the largest relative spacing.
  //           uflow is the smallest positive magnitude.
  //           oflow is the largest positive magnitude.
  //
  //***first executable statement  IntegrationController::integrate_myFloat
/*
  for (int i =0; i<50; i++)
    for (int j=0; j<50; j++)
      eps_mat[i][j]=0;
*/
  myFloat epmach = std::numeric_limits<myFloat>::epsilon();
  myFloat uflow = std::numeric_limits<myFloat>::min();
  myFloat oflow = std::numeric_limits<myFloat>::max();
  myFloat min_length = oflow;
  //
  //            test on validity of parameters
  //            -----------------------------
  //
  termination_code = normal;
  neval = 0;
  last = 0;
  result = 0;
  abserr = 0;
  int npts = int(points.size());
  int nint = npts-1;
  if (npts < 2 || limit < nint) {
    cerr << "IntegrationController::integrate: Failed first test on input." << endl;
    if (npts < 2) cerr << "npts < 2" << endl;
    if (limit < nint) cerr << "limit < nint" << endl;
    termination_code = bad_input;
    result = NAN;
    return;
  }
  if (npts==2 && points.at(0)==points.at(1)) {
        result =0;
        return;
  }
  for (typename vector<myFloat>::const_iterator ppoint=points.begin(); ppoint < points.end()-1; ppoint++)
    if (*ppoint >= *(ppoint+1)) {
        cerr << "IntegrationController::integrate: points are either not distinct or not in ascending order" << endl;
        termination_code = bad_input;
        result = NAN;
        return;
    }

  //
  //            compute first integral and error approximations.
  //            ------------------------------------------------
  //
  myFloat resabs = 0;
  int i;
  for (i=1; i<=nint; i++) {
    subs.at(i-1).a = points.at(i-1);
    subs.at(i-1).b = points.at(i);
    subs.at(i-1).integrate(f, g_k.n_gauss, g_k.w_gauss, g_k.x_kronrod, g_k.w_kronrod, fv);
    abserr += subs.at(i-1).e;
    result += subs.at(i-1).r;
    if (subs.at(i-1).e == subs.at(i-1).defabs && subs.at(i-1).e != 0.0e+00) {
      subs.at(i-1).ndin = 1;
    } else {
      subs.at(i-1).ndin = 0;
    }
    resabs += subs.at(i-1).rabs;
    subs.at(i-1).level =0;
    min_length = std::min<myFloat>(min_length, fabs(subs.at(i-1).b-subs.at(i-1).a));
  }
  myFloat errsum = abserr;
  
  //
  //           test on accuracy.
  //
  last = nint;
  neval = (2*g_k.n_gauss+1) * nint;
  myFloat dres = fabs(result);
  myFloat errbnd = std::max<myFloat>(epsabs, epsrel * dres);
  if (abserr <= 100 * epmach * resabs && abserr > errbnd) {
    termination_code = roundoff_error;
  }
  SubintervalComp<myFloat> sub_comp(noext ? limit : 1);
  make_heap(subs.begin(),subs.begin()+nint, sub_comp);
  if (limit < npts) {
    termination_code = maximum_subdivisions;
  }
  if (termination_code != normal || abserr <= errbnd) {
    return;
  }
  //
  //           initialization
  //           --------------
  //
  myFloat correc = 0.;
  myFloat area = result;
  EpsilonTable<myFloat> eps_table(52);
  int ktmin = 0;
  bool keep_level_max = false;
  bool noext = this->noext;
  myFloat erlarg = errsum;
  myFloat ertest = errbnd;
  int iroff1 = 0;
  int iroff2 = 0;
  int iroff3 = 0;
  int ierro = 0;
  abserr = oflow;
  bool f_doesnt_change_sign =(dres >= (1 - 50 * epmach) * resabs);
  bool final_check{noext ? false : true},  need_sum{noext ? true : false};
  //
  //           main do-loop
  //           ------------
  //
  for (last=npts; last<=limit ; last++) {
    //           At this point last is the number of intervals after the biscection.
    //
    //           bisect the subinterval with the largest score
    //           estimate.
    //
    pop_heap(subs.begin(), subs.begin()+last-1, sub_comp);
    int maxerr = last - 1;
    if (!subs.at(maxerr-1).is_divisible()) {
        last--;
        break;
    }
    myFloat errmax = subs.at(maxerr-1).e;
    myFloat erlast = errmax;
    myFloat rold = subs.at(maxerr-1).r;
    int level_cur = subs.at(maxerr-1).level+1;
    myFloat a1 = subs.at(maxerr-1).a;
    myFloat b2 = subs.at(maxerr-1).b;
    myFloat b1 = (a1+b2)/2;
    myFloat a2 = b1;
    subs.at(maxerr-1).b = b1;
    subs.at(last-1).a = a2;
    subs.at(last-1).b = b2;
    subs.at(maxerr-1).integrate(f, g_k.n_gauss, g_k.w_gauss, g_k.x_kronrod, g_k.w_kronrod, fv);
    subs.at(last-1).integrate(f, g_k.n_gauss, g_k.w_gauss, g_k.x_kronrod, g_k.w_kronrod, fv);
    //
    //           improve previous approximations to integral
    //           and error and test for accuracy.
    
    neval += 2*(2*g_k.n_gauss+1);
    myFloat area12 = subs.at(maxerr-1).r + subs.at(last-1).r;
    myFloat erro12 = subs.at(maxerr-1).e + subs.at(last-1).e;
    errsum += erro12 - errmax;
    area += area12 - rold;
    if (subs.at(maxerr-1).defabs != subs.at(maxerr-1).e
                && subs.at(last-1).defabs != subs.at(last-1).e
                && erro12 > epsrel*fabs(area12)) {
      if ( fabs(rold - area12) <= static_cast<myFloat>(0.1e-04) * fabs(area12)
           && erro12 >= static_cast<myFloat>(0.99e+00) * errmax) {
        if (keep_level_max) {
          iroff2++;
        } else {
          iroff1++;
        }
      }
      if (last > 10 && erro12 > errmax) {
        iroff3++;
      }
    }
    subs.at(maxerr-1).ndin = 0;
    subs.at(maxerr-1).level=level_cur;
    subs.at(last-1).ndin = 0;
    subs.at(last-1).level=level_cur;
    push_heap(subs.begin(),subs.begin()+last-1, sub_comp);
    push_heap(subs.begin(),subs.begin()+last, sub_comp);
    errbnd = std::max<myFloat>(epsabs, epsrel * fabs(area));
    min_length = std::min<myFloat>(min_length,fabs(b2-a1)/2);
    //
    //           test for roundoff error and eventually set error flag.
    //
/*
    cout << "a1 = " << a1 << ", b2 = " << b2
              << ", level = " << level_cur <<endl
              << "rold = " << rold << ", rnew = " << area12 << endl
              << "errorold = " << errmax << ", errornew = " << erro12 << endl
              << "iroff1 = " << iroff1
              << ", iroff2 = " << iroff2
              << ", iroff3 = " << iroff3 << endl;
*/
    if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
      termination_code = roundoff_error;
    }
    if (iroff2 >= 5) {
      ierro = 3;
    }
    //
    //           set error flag in the case that the number of
    //           subintervals equals limit.
    //
    if (last == limit) {
      termination_code = maximum_subdivisions;
    }
    //
    //           set error flag in the case of bad integrand behaviour
    //           at a point of the integration range
    //
    if (max<myFloat>(fabs(a1), fabs(b2)) <= (1 + 100 * epmach) * (fabs(a2) + 1000 * uflow)) {
      termination_code = bad_integrand;
    }
    final_check= noext ? false : true;
    if (errsum <= errbnd) {
      // ***jump out of do-loop
 //     cout << "IntegrationController::integrate: Success errsum <= errbnd" << endl;
      need_sum = true;
      final_check = false;
      break;
    }
    if (termination_code != normal) {
      // ***jump out of do-loop
//      cout << "IntegrationController::integrate: Aborting with raw error code = " << termination_code << endl;
      break;
    }
    erlarg = erlarg - erlast;
      if (level_cur < sub_comp.level_max)
      erlarg += erro12;
    if (erlarg < std::max<myFloat>(epsabs,epsrel * fabs(result))) {
      if (noext) {
        continue;
      }
      //
      //           perform extrapolation.
      //
      eps_table.update(area);
      if (eps_table.n > 2) {
        myFloat reseps = eps_table.result;
        myFloat abseps = eps_table.abserr;
        ktmin++;
        if (ktmin > 10 && abserr < 0.1e-02 * errsum) {
          termination_code = extrapolation_roundoff_error;
        }
        if (abseps < abserr) {
//          cout << "Using reseps = " << reseps
//                    << ", vs result = " << result << endl;
          ktmin = 0;
          abserr = abseps;
          result = reseps;
          correc = erlarg;
          ertest = std::max<myFloat>(epsabs, epsrel * fabs(reseps));
          // ***jump out of do-loop
          if (abserr < ertest) {
 //           cout << "IntegrationController::integrate: Success.  abserr from dqelg < errtest." << endl;
            break;
          }
        }
        if (eps_table.n == 1) {
//          cout << "Performed extrapolation but eps_table.n == 1" << endl;
          noext = true;
          sub_comp.level_max = limit;
          make_heap(subs.begin(),subs.begin()+last, sub_comp);
        }
        if (termination_code >= extrapolation_roundoff_error) {
//          cout << "IntegrationController::integrate: Aborting because of roundoff error in extrapolation table" << endl;
          break;
        }
      }
      //
      //           prepare bisection of the smallest interval.
      //
      keep_level_max = false;
      erlarg = errsum;
      ++sub_comp.level_max;
      make_heap(subs.begin(),subs.begin()+last, sub_comp);
    } else {
      //
      //           keep working on intervals > level_max
      //
      keep_level_max = true;
    }
  } // main do loop
  //
  //           set the final result.
  //           ---------------------
  //
  if (final_check){
    bool test_for_divergence{true};
    if (abserr == oflow) {
      test_for_divergence=false;
      need_sum=true;
    } else if ((termination_code + ierro) == 0) {
      test_for_divergence=true;
    } else {
      if (ierro == 3) {
//        cout << "IntegrationController::integrate: ierro = 3 adding " << correc << " to abserr." << endl;
        abserr += correc;
      }
      if (termination_code == normal) {
        termination_code = roundoff_error;        // Roundoff error detected
      }
      if (result != 0.0e+00 && area != 0.0e+00) {
        if (abserr / fabs(result) > errsum / fabs(area)) {
          need_sum=true;
          test_for_divergence=false;
        } else {
          test_for_divergence=true;
        }
      } else if (abserr > errsum) {
         need_sum = true;
         test_for_divergence=false;
      } else if (area == 0.0e+00) {
        need_sum=false;
        test_for_divergence=false;
      } else {
        test_for_divergence=true;
      }
    }
    //
    //           test on divergence.
    //
    if (test_for_divergence){
      need_sum=false;
      if (f_doesnt_change_sign
          && max<myFloat>(fabs(result), fabs(area)) > resabs * 0.01) {
        if (   0.01 > (result / area)
            || (result / area) > 100
            || errsum > fabs(area)) {
          termination_code = divergent_integral;      // Divergent integral
        }
      }
    } // test for divergence
  } // final_check
  //
  //           compute global integral sum.
  //
  if (need_sum){
//    cout << "Calculating sum" << endl;
    result = 0.0e+00;
    errsum = 0;
    for (int k=1; k<=last; k++) {
      result += subs.at(k-1).r;
      errsum += subs.at(k-1).e;
    }
    abserr = errsum;
  }
}
  
  template<typename myFloat, typename F>
  myFloat Integral<myFloat, F>::operator() () {
    Fmt<myFloat> fmt;
    
    if (verbose==4){
      cout << endl
      << "Integral " << "with controller:" << endl << controller << endl;
    }
    controller.integrate(f, points,
                          result, abserr, neval, termination_code, last);
    
    if (verbose>=3){
      myFloat rsum=0, esum=0;
      for (int i=0; i<last; i++) {
        rsum += controller.subs.at(i).r;
        esum += controller.subs.at(i).e;
      }
      
      if (termination_code > 0)
        cout << msg() << ":" << endl;
      cout << "Integral from " << fmt << points.front()
      << " to " << fmt << points.back()
      << " = " << fmt << result
      << ", with absolute error = " << fmt << abserr << endl
      << "Number of evaluations = " << neval
      << ", subintervals = " << last << endl
      << "rsum = " << fmt << rsum << ", esum = " << fmt << esum << endl;
      print_subs_summary(cout, controller.subs, last, points);
    }
    if (verbose>=4){
      print_subs(cout, controller.subs, last, points);
    }
    return result;
  }
  
  
  
} // namespace adaptive_integration



