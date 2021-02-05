#ifndef __PESOPTEIGEN_H
#define __PESOPTEIGEN_H


#ifdef PESOPT_WITH_EIGEN

#include <includesEigen.h>

// modified versions of Eigen solvers
#include <EigenBiCGSTAB_New.h>
#include <EigenGMRES_New.h>
#include <EigenProjectedPCG_New.h>

#endif   // PESOPT_WITH_EIGEN



namespace pesopt{

// Returns minimum/maximium
template<class T> inline T Min ( const T a, const T b ) { return ( ( a < b ) ? a : b ); }
template<class T> inline T Max ( const T a, const T b ) { return ( ( a < b ) ? b : a ); }
template<class T> inline T Clamp ( const T Value, const T Min, const T Max ) { return ( pesopt::Max ( pesopt::Min ( Value, Max ), Min ) ); } // Returns Value clamped into [Min,Max].
template<class T> inline T Sqr (const T a) { return a * a; }
template<class T> inline T Cub (const T a) { return a * a * a; }
//! Signum function template, signum ( 0 ) = 0
// template<class T> inline T signum ( const T x ) {
//   if (x > 0.0) return 1;
//   else if (x == 0.0) return 0;
//   else if (x < 0.0) return -1;
//   return x; // NaN
// }
// int countDigitsOfNumber ( const int N ) {  return ( N != 0 ) ? static_cast<int>(floor(log10(static_cast<double>(std::abs(N)))))+1 : 1; }


template<typename RealType, typename VectorType>
void thresholdVector ( const VectorType & vec, VectorType & thresholdVec, RealType lowerBound, RealType upperBound, RealType threshold ){
    for( int i=0; i<vec.size(); ++i ){
      if( vec[i] < threshold ) thresholdVec[i] = lowerBound;
      else thresholdVec[i] = upperBound;
    }
}


template<typename RealType, typename VectorType>
void scaleVector ( const VectorType & vec, const RealType l, const RealType u,
                   VectorType & scaledVec, const RealType newL, const RealType newU ){
   const RealType a = (newU - newL ) / (u - l);
   const RealType b = ( u * newL - l * newU ) / (u - l);
   for( int i=0; i< scaledVec.size(); ++i ){
     scaledVec[i] = a * vec[i] + b;   
   }
}

// template<typename RealType, typename VectorType>
// void scaleVectorTo01 ( const VectorType & vec, VectorType & scaledVec){
//     RealType max = vec.maxCoeff();
// //     RealType min = vec.minCoeff();
// //     RealType range = max - min;
//     if( max > 1.e-16 ){ scaledVec = vec; scaledVec /= max; }
//     else scaledVec.setZero( );
// }


// template<typename RealType, typename VectorType>
// void scaleVectorToEps1 ( const VectorType & vec, VectorType & scaledVec, const RealType eps ){
//     RealType max = vec.maxCoeff();
//     if( max > 1.e-16 ){ scaledVec = vec; scaledVec /= max; }
//     else scaledVec.setZero( );
//     for( int i=0; i< scaledVec.size(); ++i ){
//      if( scaledVec[i] < 1.e-16 ) scaledVec[i] = eps;   
//     }
// }

//! Frobenius scalar product
template<typename RealType, typename Matrix>
RealType ddProd( const Matrix &MatA, const Matrix &MatB ){
    RealType out = 0.0;
    for( int i=0; i<MatA.rows(); ++i ) 
        for( int j=0; j<MatB.cols(); ++j )
            out += MatA(i,j) * MatB(i,j);
        
    return out;
}

} // namespace pesopt


#include <pesoptTensor.h>


#endif  //__PESOPTEIGEN_H
