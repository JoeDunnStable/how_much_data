/// \file gauss_kronrod_impl.h
/// Implementation of routines to calculate nodes and weights
/// Included in file gauss_kronrod.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <Eigen/Eigenvalues>
#include <iomanip>
#include <vector>
#include <algorithm>

namespace gauss_kronrod {
using Eigen::SelfAdjointEigenSolver;

using std::endl;
using std::setw;
using std::setprecision;
using std::right;
using std::vector;

template<typename T>
class sort_data {
public:
  T a;
  int index;
  sort_data() : a(0), index(0){}
  sort_data(T a, int index) : a(a), index(index){}
};

template<typename T>
bool sd_lt(const sort_data<T> &lhs,const sort_data<T>& rhs) {
  return lhs.a < rhs.a;
}

// toms726 generates the nodes x and weights w for a quadrature rule
// given the recurrence factors a and b for the orthogonal polynomials
//
template<typename myFloat>
void toms726(const int n_gauss, const vector<myFloat>& a, const vector<myFloat>& b,
             vector<myFloat>& x, vector<myFloat>& w, const int verbose){
  Fmt<myFloat> fmt;
  using MatrixXF = Eigen::Matrix<myFloat,Eigen::Dynamic,Eigen::Dynamic>;
  using VectorXF = Eigen::Matrix<myFloat,Eigen::Dynamic,1>;
  int sumbgt0=0;
  for (int j=0; j<n_gauss; j++) sumbgt0+=(b.at(j)>0);
  if(sumbgt0!=n_gauss) {
    throw std::range_error("toms726: b is not all >0");
  }
  if (n_gauss == 1) {
    x.at(0)=a.at(0);
    w.at(0)=b.at(0);
    return;
  }
  MatrixXF J = MatrixXF::Zero(n_gauss,n_gauss);
  for (int k=0; k<n_gauss-1; k++){
    J(k,k)=a.at(k);
    J(k,k+1)=sqrt(b.at(k+1));
    J(k+1,k)=J(k,(k+1));
  }
  J(n_gauss-1,n_gauss-1)=a.at(n_gauss-1);
  SelfAdjointEigenSolver<MatrixXF> es;
  es.compute(J);
  VectorXF e_values = es.eigenvalues();   // These are as yet unsorted
  MatrixXF e_vectors = es.eigenvectors();
  vector<myFloat> e(n_gauss);
  for (int k=0; k<n_gauss; k++) {
    e.at(k) = b.at(0) * pow(e_vectors(0,k), 2);
  }
  vector<sort_data<myFloat> > srt_d(n_gauss);
  for (int i=0; i<n_gauss; i++) {
    sort_data<myFloat> elem(e_values(i),i);
    srt_d.at(i) = elem;
  }
  std::sort(srt_d.begin(), srt_d.end(),sd_lt<myFloat>);
  for (int i=0; i<n_gauss; i++){
    if (verbose)
      cout<< fmt  << srt_d.at(i).a
           << fmt << srt_d.at(i).index << endl;
    x.at(i)=srt_d.at(i).a;
    w.at(i)=e.at(srt_d.at(i).index);
  }
}

template<typename myFloat>
void r_jacobi01(const int n_gauss, const myFloat a, const myFloat b, vector<myFloat>& c,
                vector<myFloat>& d){
  if((n_gauss<=0)||(a<=-1)||(b<=-1)){
    throw std::range_error("r_jacobi01: parameter(s) out of range");
  }
  r_jacobi(n_gauss, a, b, c, d);
  for (int i=0; i<n_gauss; i++) {
    c.at(i)=(1+c.at(i))/2;
    if (i==0)
      d.at(i)=d.at(0)/pow(2,a+b+1);
    else
      d.at(i)=d.at(i)/4;
  }
}

template<typename myFloat>
void r_jacobi(const int n_gauss, const myFloat a, const myFloat b, vector<myFloat>& c,
              vector<myFloat>& d){
  if((n_gauss<=0)||(a<=-1)||(b<=-1)) {
    throw std::range_error("r_jacobi: parameter(s) out of range");
  }
  myFloat nu = (b-a)/(a+b+2);
  myFloat mu = pow(2,a+b+1)*tgamma(a+1)*tgamma(b+1)/tgamma(a+b+2);
  if (n_gauss==1) {
    c.at(0)=nu;
    d.at(0)=mu;
    return;
  }
  c.at(0)=nu;
  for (int n=1; n<n_gauss; n++){
    myFloat nab=2*n+a+b;
    c.at(n) = (pow(b,2)-pow(a,2))*1/(nab*(nab+2));
  }
  d.at(0) = mu;
  d.at(1) = 4*(a+1)*(b+1)/(pow(a+b+2,2)*(a+b+3));
  for (int n=2; n<n_gauss; n++) {
    myFloat nab=2*n+a+b;
    d.at(n)=4*(n+a)*(n+b)*n*(n+a+b)/(pow(nab,2)*(nab+1)*(nab-1));
  }
}

template<typename myFloat>
void r_kronrod(const int n_gauss, const vector<myFloat>& a0, const vector<myFloat>& b0,
               vector<myFloat>& a, vector<myFloat>& b) {
  if (a0.size()<ceil(3.*n_gauss/2.)+1) {
    throw std::range_error("r_kronrod: a0 is too short.");
  }
  for (int k=0; k<=ceil(3.*n_gauss/2.); k++) {
    a.at(k) = a0.at(k);
    b.at(k) = b0.at(k);
  }
  vector<myFloat> s((n_gauss/2)+2);
  vector<myFloat> t((n_gauss/2)+2);
  for (int k=0; k<(n_gauss/2)+2; k++)
    t.at(k)=s.at(k)=myFloat(0);
  t.at(1)=b[n_gauss+1];	// P.H.
  for (int m=0; m<=(n_gauss-2); m++) {
    myFloat cumsum = 0;
    vector<myFloat> tmp((n_gauss/2)+2);
    tmp.at(0)=s.at(0);
    for (int k=((m+1)/2); k>=0; k--){
      int l = m-k;
      cumsum=cumsum + (a.at(k+n_gauss+1) - a.at(l))*t.at(k+1) + b.at(k+n_gauss+1)*s.at(k)
               - b.at(l)*s.at(k+1);
      tmp.at(k+1)=cumsum;
    }
    for (int k=((m+1)/2)+2; k<((n_gauss/2)+2); k++)
      tmp.at(k)=s.at(k);
    for (int k=0; k<((n_gauss/2)+2); k++){
      s.at(k)=t.at(k);
      t.at(k)=tmp.at(k);
    }
  }
  for (int j=(n_gauss/2); j>=0; j--)
    s.at(j+1)=s.at(j);
  for (int m=n_gauss-1; m<=(2*n_gauss-3); m++){
    myFloat cumsum = 0;
    vector<myFloat> tmp((n_gauss/2)+2);
    for (int k=m+1-n_gauss; k<=((m-1)/2); k++){
      int l=m-k;
      int j=n_gauss-1-l;
      cumsum=cumsum+-(a.at(k+n_gauss+1)-a.at(l))*t.at(j+1)-b.at(k+n_gauss+1)*s.at(j+1)+b.at(l)*s.at(j+2);
      tmp.at(j+1)=cumsum;
    }
    int j{0};
    for (int k=m+1-n_gauss; k<=((m-1)/2); k++){
      j=n_gauss-1-m+k;
      s.at(j+1)=tmp.at(j+1);
    }
    int k=((m+1)/2);
    if ((m % 2)==0) {
      a.at(k+n_gauss+1)=a.at(k)+(s.at(j+1)-b.at(k+n_gauss+1)*s.at(j+2))/t.at(j+2);
    }
    else {
      b.at(k+n_gauss+1)=s.at(j+1)/s.at(j+2);
    }
    for (int j=0; j<((n_gauss/2)+2); j++) {
      myFloat swap = s.at(j);
      s.at(j)=t.at(j);
      t.at(j)=swap;
    }
  }
  a.at(2*n_gauss)=a.at(n_gauss-1)-b.at(2*n_gauss)*s.at(1)/t.at(1);
  return;
}
  
template<typename myFloat>
  Kronrod<myFloat>::Kronrod(const int n_gauss, const int verbose)
    : n_gauss(n_gauss), verbose(verbose) {
  x_gauss.resize(n_gauss);
  w_gauss.resize(n_gauss);
  x_kronrod.resize(2*n_gauss+1);
  w_kronrod.resize(2*n_gauss+1);
  int M = std::max(2*n_gauss,int(ceil(3.*n_gauss/2.)+1));
  vector<myFloat> a0(M);
  vector<myFloat> b0(M);
  r_jacobi01(M,myFloat(0),myFloat(0), a0, b0);
  toms726(n_gauss, a0, b0, x_gauss, w_gauss, verbose>4);
  vector<myFloat> a(2*n_gauss+1);
  vector<myFloat> b(2*n_gauss+1);
  r_kronrod(n_gauss, a0, b0, a, b);
  toms726(2*n_gauss+1, a, b, x_kronrod, w_kronrod, verbose > 4);
  for (int i = 0; i<n_gauss; i++) {
    x_gauss.at(i)=2*x_gauss.at(i)-1;
    w_gauss.at(i)=2*w_gauss.at(i);
  }
  
  for (int i = 0; i<2*n_gauss+1; i++) {
    x_kronrod.at(i) = 2*x_kronrod.at(i)-1;
    w_kronrod.at(i) = 2*w_kronrod.at(i);
  }
} // Kronrod constructor
  
/// Operator to print the kronrod nodes and wieghts
template<typename myFloat>
  ostream& operator<< (ostream& os, const Kronrod<myFloat>& k) {
    int d = std::numeric_limits<myFloat>::digits10;
    os << "Gauss nodes and weights" << endl << endl;
    os << setw(d+10) << right << "Node" << setw(d+10) << right << "Weight" << endl;
    for (int i=0; i<k.n_gauss; ++i) {
      os << setw(d+10) << setprecision(d) << k.x_gauss.at(i)
      << setw(d+10) << setprecision(d) << k.w_gauss.at(i) << endl;
    }
    os << endl << "Kronrod nodes and weights" << endl << endl;
    os << setw(d+10) << right << "Node" << setw(d+10) << right << "Weight" << endl;
    for (int i=0; i<2*k.n_gauss+1; ++i) {
      os << setw(d+10) << setprecision(d) << k.x_kronrod.at(i)
      << setw(d+10) << setprecision(d) << k.w_kronrod.at(i) << endl;
    }
    return os;
  }

  
} // namespace gauss_kronrod
