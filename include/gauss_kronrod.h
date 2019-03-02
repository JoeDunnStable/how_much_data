///
/// \file  gauss_kronrod.h
/// Routines calculating Gauss Kronrod nodes and weights
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef gauss_kronrod_h
#define gauss_kronrod_h
#include <vector>
#include <iostream>

namespace gauss_kronrod {
using std::vector;
using std::ostream;

/// Calculate integration nodes and weights given recursion coeffecients of orthogonal polynomials.
///  @see http://people.sc.fsu.edu/~jburkardt/f_src/toms726/toms726.html for the fortran version
///
template<typename myFloat>
void toms726(const int n_gauss,               ///< [in] the number of nodes
             const vector<myFloat>& a,  ///< [in] x-a is the coefficient of pj in the recursion
             const vector<myFloat>& b,  ///< [in] b is the coefficient of pj-1 in the recursion
             vector<myFloat>& x,        ///< [out] the nodes of the quadrature formula
             vector<myFloat>& w,        ///< [out] the weights attached to the nodes
             const int verbose          ///< [in] a flag indicating verbose output
);

/// Generate recurrence coefficients for monic Jacobi polynomials on [0,1].
///
///   r_jacobi01(n,a,b,c,d) generates the first n recurrence
///   coefficients for monic Jacobi polynomials on [0,1] with
///   parameters a and b. These are orthogonal on [0,1] relative
///   to the weight function w(t)=(1-t)^a t^b. The n alpha-
///   coefficients are stored in c, the n beta-
///   coefficients in d.
///
///   Translated to C++ by Joseph Dunn, Dec. 26, 2015
///
template<typename myFloat>
void r_jacobi01(const int n_gauss,         ///< [in] the number of recurrence coeficients to generate
                const myFloat a,     ///< [in] the exponent of (1-t) in the weight function
                const myFloat b,     ///<[in] the exponented of t in the weight function
                vector<myFloat>& c,  ///< [out] the alpha coefficients
                vector<myFloat>& d   ///< [out] the beta coefficients
);

/// Generate recurrence coefficients for monic Jacobi polynomials on [-1,1]
///
///   r_jacobi(n,a,b,c,d) generates the first n recurrence
///   coefficients for monic Jacobi polynomials with parameters
///   a and b. These are orthogonal on [-1,1] relative to the
///   weight function w(t)=(1-t)^a(1+t)^b. The n alpha-coefficients
///   are stored in c, the n beta-coefficients in d
///
///   Translated to C++ by Joseph Dunn, Dec. 26, 2015
///   Supplied by Dirk Laurie, 6-22-1998; edited by Walter
///   Gautschi, 4-4-2002.
///

template<typename myFloat>
void r_jacobi(const int n_gauss,        ///< [in] the number of recurrence coeficients to generate
              const myFloat a,    ///< [in] the exponent of (1-t) in the weight function
              const myFloat b,    ///<[in] the exponented of t in the weight function
              vector<myFloat>& c, ///< [out] the alpha coefficients
              vector<myFloat>& d  ///< [out] the beta coefficients
);

/// Generate the recurrence coeficients for the Jacobi-Kronrod matrix.
///
///   r_kronrod(n,a0,b0,a,b) produces the alpha- and beta-elements in
///   the Jacobi-Kronrod matrix of order 2n+1 for the weight
///   function (or measure) w. The input data for the weight
///   function w are the recurrence coefficients of the associated
///   orthogonal polynomials, which are stored in the vectors a0 & b0 of
///   dimension [ceiling(3*n/2)+1]. The alpha-elements are stored in
///   a, the beta-elements in b, of
///   dimension (2*n+1)
///
///   Supplied by Dirk Laurie, 6-22.1998
///   Edited by Pavel Holoborodko, n_gaussovember 7, 2011:
///
///   Translated to C++ by Joseph Dunn, Dec. 26 2015
///
///@see http://dip.sun.ac.za/~laurie/papers/kronrod/kronrod.ps
template<typename myFloat>
void r_kronrod(int n_gauss,                      ///< [in] the number of nodess in the Jacobi Gauss quadrature formula
               const vector<myFloat>& a0,  ///< [in] the alpha coeficients of the Jacobi Gauss recurrence
               const vector<myFloat>& b0,  ///< [in] the beta coefficients of the Jacobi Gauss recurrence
               vector<myFloat>& a,         ///< [out] the alpha coefficients of the Jacobi-Kronrod recurrence
               vector<myFloat>& b         ///< [out] the beta coefficients of the Jacobi-Kronrod recurrence
);
  
template<typename myFloat> class Kronrod;
/// Operator to print out class Kronrod
template<typename myFloat> ostream& operator<< (ostream& os, const Kronrod<myFloat>& k);

/// Class contain the x nodes and weights for gauss kronrod integration
/// For unit weight function
template<typename myFloat>
  class Kronrod {
  public:
    const int n_gauss;          ///< the number of nodes for Gauss integration
    const int verbose;          ///< a indicator for the level of trace output
    vector<myFloat> x_gauss;   ///< the gauss nodes
    vector<myFloat> w_gauss;   ///< the gauss weights
    vector<myFloat> x_kronrod; ///< the kronrod nodes
    vector<myFloat> w_kronrod; ///< the kronrod weights
    Kronrod(const int n_gauss, const int verbose=0);
    template<typename BigFloat>
    Kronrod(const Kronrod<BigFloat>& k_big)
    : n_gauss(k_big.n_gauss), verbose(k_big.verbose){
      x_gauss.resize(n_gauss);
      w_gauss.resize(n_gauss);
      x_kronrod.resize(2*n_gauss+1);
      w_kronrod.resize(2*n_gauss+1);
      for (int i=0; i<n_gauss; ++i) {
        x_gauss.at(i) = static_cast<myFloat>(k_big.x_gauss.at(i));
        reset_prec(x_gauss.at(i));
        w_gauss.at(i) = static_cast<myFloat>(k_big.w_gauss.at(i));
        reset_prec(w_gauss.at(i));
      }
      for (int i=0; i<2*n_gauss+1; ++i) {
        x_kronrod.at(i) = static_cast<myFloat>(k_big.x_kronrod.at(i));
        reset_prec(x_kronrod.at(i));
        w_kronrod.at(i) = static_cast<myFloat>(k_big.w_kronrod.at(i));
        reset_prec(w_kronrod.at(i));
      }
    };
    friend ostream& operator<< <> (ostream& os, const Kronrod<myFloat>& k);
  };
    
} // namespace gauss_kronrod

#include "gauss_kronrod_impl.h"

#endif // gauss_kronrod_h
