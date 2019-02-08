//
//  p_spline.h
//  spline_test
//
//  Created by Joseph Dunn on 2/2/19.
//  Copyright © 2019 Joseph Dunn. All rights reserved.
//

#ifndef p_spline_h
#define p_spline_h

#include <unsupported/Eigen/Splines>
using Eigen::Spline;
using Eigen::SplineFitting;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::DenseIndex;
using Eigen::ComputeThinU;
using Eigen::ComputeThinV;
#include <Eigen/SVD>

/// Fit a spline using linear least squares to fit parameterstemplate <typename SplineType>
template <typename SplineType, typename PointArrayType>
static SplineType SplineFitLeastSquare(const PointArrayType& pts,                    ///< the values at the fitting pts
                                       const typename SplineType::KnotVectorType& x, ///< the loaction for the values
                                       ///< Should be within the range of knot_parameters
                                       DenseIndex degree,                            ///< the degree of the polynomials
                                       const typename SplineType::KnotVectorType& knot_parameters) ///< knot for the spline
{
  typedef typename SplineType::KnotVectorType::Scalar Scalar;
  typedef typename SplineType::ControlPointVectorType ControlPointVectorType;
  
  typedef Matrix<Scalar,Dynamic,Dynamic> MatrixType;
  
  typename SplineType::KnotVectorType knots;
  KnotAveraging(knot_parameters, degree, knots);
  
  DenseIndex n = knot_parameters.size();
  DenseIndex m = pts.cols();
  
  MatrixType A = MatrixType::Zero(m,n);
  for (DenseIndex i=0; i<m; ++i)
  {
    const DenseIndex span = SplineType::Span(x[i], degree, knots);
    
    // The segment call should somehow be told the spline order at compile time.
    A.row(i).segment(span-degree, degree+1) = SplineType::BasisFunctions(x[i], degree, knots);
  }
  /*
   A(0,0) = 1.0;
   A(m-1,m-1) = 1.0;
   */
  
  // Here, we are creating a temporary due to an Eigen issue.
  ControlPointVectorType ctrls = A.bdcSvd(ComputeThinU | ComputeThinV).solve(MatrixType(pts.transpose())).transpose();
  
  return SplineType(knots, ctrls);
}

/// Fit a p-spline using linear least squares to fit parameterst
template <typename SplineType, typename PointArrayType>
static SplineType SplineFitPSpline(const PointArrayType& pts,                    ///< the values at the fitting pts
                                   const typename SplineType::KnotVectorType& x, ///< the location for the values
                                       ///< Should be within the range of knot_parameters
                                   DenseIndex degree,                            ///< the degree of the polynomials
                                   const typename SplineType::KnotVectorType& knot_parameters, ///< knot for the spline
                                   const typename SplineType::KnotVectorType::Scalar penality                       ///< penality on first differenct of ctrls
                                   )
{
  typedef typename SplineType::KnotVectorType::Scalar Scalar;
  typedef typename SplineType::ControlPointVectorType ControlPointVectorType;
  
  typedef Matrix<Scalar,Dynamic,Dynamic> MatrixType;
  
  typename SplineType::KnotVectorType knots;
  KnotAveraging(knot_parameters, degree, knots);
  
  DenseIndex n = knot_parameters.size();
  DenseIndex m = pts.cols();
  
  MatrixType A = MatrixType::Zero(m+n-1,n);
  for (DenseIndex i=0; i<m; ++i)
  {
    const DenseIndex span = SplineType::Span(x[i], degree, knots);
    
    // The segment call should somehow be told the spline order at compile time.
    A.row(i).segment(span-degree, degree+1) = SplineType::BasisFunctions(x[i], degree, knots);
  }
  A.block(m,0,n-1,n-1).diagonal().setConstant(sqrt(penality));
  A.block(m,0,n-1,n-1).diagonal(1).setConstant( -sqrt(penality));
  
  MatrixType rhs(m+n-1,pts.rows());
  rhs.block(0,0,m,pts.rows())=pts.transpose();
  rhs.block(m,0,n-1,pts.rows()).setZero();
  
  // Here, we are creating a temporary due to an Eigen issue.
  ControlPointVectorType ctrls = A.bdcSvd(ComputeThinU | ComputeThinV).solve(rhs).transpose();
  
  return SplineType(knots, ctrls);
}

#endif /* SplitFitLeastSquare_h */
