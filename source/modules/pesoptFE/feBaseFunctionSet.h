#ifndef __FEBASEFUNCTIONSET_H
#define __FEBASEFUNCTIONSET_H

#include <initializer_list>
#include <pesopt_IO.h>


//! Base class
template <typename DataTypeContainer, class QuadRuleType >
class FEBaseFunctionSetBase  {
    
protected:
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
  QuadRuleType _quadRule;
  
public:
  FEBaseFunctionSetBase( ) {}

  int numQuadPoints( ) const { return QuadRuleType::numQuadPoints;}
  inline RealType getWeight ( int QuadPoint ) const { return _quadRule.getWeight ( QuadPoint );}
  inline const RealVecChart& getRefCoord ( int QuadPoint ) const { return _quadRule.getRefCoord ( QuadPoint );}
};



//! base class for cached version
template <typename DataTypeContainer, class QuadRuleType, int numBaseFuncs, typename Imp >
class FECachedBaseFunctionSetBase  {
    
protected:
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;

  /**** cache the values of the basis functions at the quadrature points ****/
  RealType     basisQuadValues   [numBaseFuncs][QuadRuleType::numQuadPoints];
  RealVecChart basisQuadGradients[numBaseFuncs][QuadRuleType::numQuadPoints];
  QuadRuleType _quadRule;
    
public:
  FECachedBaseFunctionSetBase( ) {}

  int numQuadPoints( ) const { return QuadRuleType::numQuadPoints;}
  inline RealType getWeight ( int QuadPoint ) const {return _quadRule.getWeight ( QuadPoint );}
  inline const RealVecChart& getRefCoord ( int QuadPoint ) const { return _quadRule.getRefCoord ( QuadPoint );}

  //! read the cached value of the basis function with number BaseFuncNum at the given quadrature point
  inline RealType evaluate ( int BaseFuncNum, int QuadPoint ) const { return basisQuadValues[BaseFuncNum][QuadPoint];}
  inline const RealVecChart& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const { return basisQuadGradients[BaseFuncNum][QuadPoint];}
  inline RealType evaluate ( int BaseFuncNum, const RealVecChart &RefCoord ) const { return asImp().evaluate ( BaseFuncNum, RefCoord );}
  inline void evaluateGradient ( int BaseFuncNum, const RealVecChart& RefCoord, RealVecChart& Gradient ) const {asImp().evaluateGradient ( BaseFuncNum, RefCoord, Gradient );}
  
protected:
      
  Imp &asImp() { return static_cast<Imp&> ( *this ); }
  const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

  void initializeQuadCache( ) {
    for ( int b = 0; b < numBaseFuncs; b++ )
      for ( int i = 0; i < QuadRuleType::numQuadPoints; i++ ) {
        basisQuadValues[b][i] = evaluate ( b,  _quadRule.getRefCoord ( i ) );
        evaluateGradient ( b, _quadRule.getRefCoord ( i ), basisQuadGradients[b][i] );
      }
  }

};





template <typename DataTypeContainer, class QuadRuleType, int numBaseFuncs, typename Imp >
class FECachedSecondOrderBaseFunctionSetBase  {
    
protected:
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
  typedef typename DataTypeContainer::DerivativeVectorValuedType DerivativeVectorValuedType;
    
   /**** cache the values of the basis functions at the quadrature points ****/
  RealType                            basisQuadValues   [numBaseFuncs][QuadRuleType::numQuadPoints];
  RealVecChart                        basisQuadGradients[numBaseFuncs][QuadRuleType::numQuadPoints];
  DerivativeVectorValuedType          basisQuadHessian  [numBaseFuncs][QuadRuleType::numQuadPoints];
  QuadRuleType _quadRule;
    
public:
  FECachedSecondOrderBaseFunctionSetBase ( ) {}

  int numQuadPoints( ) const { return QuadRuleType::numQuadPoints;}
  inline RealType getWeight ( int QuadPoint ) const {return _quadRule.getWeight ( QuadPoint );}
  inline const RealVecChart& getRefCoord ( int QuadPoint ) const { return _quadRule.getRefCoord ( QuadPoint );}

  //! read the cached value of the basis function with number BaseFuncNum at the given quadrature point
  inline RealType evaluate ( int BaseFuncNum, int QuadPoint ) const { return basisQuadValues[BaseFuncNum][QuadPoint];}
  inline const RealVecChart& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const { return basisQuadGradients[BaseFuncNum][QuadPoint];}
  inline const DerivativeVectorValuedType& evaluateHessian ( int BaseFuncNum, int QuadPoint ) const { return basisQuadHessian[BaseFuncNum][QuadPoint];
  }
  inline const RealType evaluateLaplace( int BaseFuncNum, int QuadPoint ) const {
      const int dimChartDomain = basisQuadHessian[BaseFuncNum][QuadPoint].rows();
      RealType lap = 0;
      for ( int i = 0; i<dimChartDomain; ++i )
          lap += basisQuadHessian[BaseFuncNum][QuadPoint](i, i);
      return lap;
  }
  
  //! functions to be implemented in the derived class
  inline RealType evaluate ( int BaseFuncNum, const RealVecChart &RefCoord ) const { return asImp().evaluate ( BaseFuncNum, RefCoord );}
  inline void evaluateGradient ( int BaseFuncNum, const RealVecChart& RefCoord, RealVecChart& Gradient ) const {asImp().evaluateGradient ( BaseFuncNum, RefCoord, Gradient );}
  inline void evaluateHessian ( int BaseFuncNum, const RealVecChart& RefCoord, DerivativeVectorValuedType& Hessian ) const {
    asImp().evaluateHessian( BaseFuncNum, RefCoord, Hessian );
  }

protected:
  Imp &asImp() { return static_cast<Imp&> ( *this ); }

  const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

  void initializeQuadCache( ) {
    for ( int b = 0; b < numBaseFuncs; b++ ) {
      for ( int i = 0; i < QuadRuleType::numQuadPoints; i++ ) {
        basisQuadValues[b][i] = evaluate ( b,  _quadRule.getRefCoord ( i ) );
        evaluateGradient ( b, _quadRule.getRefCoord ( i ), basisQuadGradients[b][i] );
        evaluateHessian ( b, _quadRule.getRefCoord ( i ), basisQuadHessian[b][i] );
      }
    }
  }

};



//! ========================================================================================
//! global basefunction set for gradient of scalar valued function 
//! \f$ u : \R^d \to \R, u(x) = <a,x> \f$
//! =======================================================================================

template <typename DataTypeContainer>
class GlobalAffineGradBaseFunctionSet2D {
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
    
  // b1(x) = (x,0)
  // b2(x) = (0,y)
  
  static void _b1  ( const RealVecChart &c, RealVecChart &dest ) { dest[0] = c[0]; dest[1] = 0.0; }
  static void _b2  ( const RealVecChart &c, RealVecChart &dest ) { dest[0] = 0.0;  dest[1] = c[1]; }
  
  static void _Grad_b1  ( const RealVecChart &c, RealVecChart &dest ) { dest.setZero(); dest(0) = 1.0; }
  static void _Grad_b2  ( const RealVecChart &c, RealVecChart &dest ) { dest.setZero(); dest(1) = 1.0; }

  typedef void ( *BASIS_FUNC_TYPE ) ( const RealVecChart &c, RealVecChart &dest );
  BASIS_FUNC_TYPE _basis[2];
  typedef void ( *BASIS_FUNC_TYPE_GRAD ) ( const RealVecChart &c, RealVecChart &dest );
  BASIS_FUNC_TYPE_GRAD _Grad_basis[2];

public:
  GlobalAffineGradBaseFunctionSet2D ( ) {
    _basis[0] = _b1; _basis[1] = _b2;
    _Grad_basis[0] = _Grad_b1; _Grad_basis[1] = _Grad_b2;
  }

  static const int numBaseFuncs = 2;
  
  void evaluate ( int BaseFuncNum, const RealVecChart &c, RealVecChart &dest ) const { _basis[BaseFuncNum] ( c, dest );}
  void evaluateGrad ( int BaseFuncNum, const RealVecChart &c, RealVecChart &dest ) const { _Grad_basis[BaseFuncNum] ( c, dest );}

};


template <typename DataTypeContainer>
class GlobalAffineGradBaseFunctionSet3D{

  typedef typename DataTypeContainer::RealType RealType;  
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
  
  static void _b1   ( const RealVecChart &c, RealVecChart & dest ) { dest[0] = c[0]; dest[1] = 0.0; dest[2] = 0.0; }
  static void _b2   ( const RealVecChart &c, RealVecChart & dest ) { dest[0] = 0.0; dest[1] = c[1]; dest[2] = 0.0; }
  static void _b3   ( const RealVecChart &c, RealVecChart & dest ) { dest[0] = 0.0; dest[1] = 0.0; dest[2] = c[2]; }

  static void _Grad_b1  ( const RealVecChart &c, RealVecChart &dest ) { dest.setZero(); dest(0) = 1.0; }
  static void _Grad_b2  ( const RealVecChart &c, RealVecChart &dest ) { dest.setZero(); dest(1) = 1.0; }
  static void _Grad_b3  ( const RealVecChart &c, RealVecChart &dest ) { dest.setZero(); dest(2) = 1.0; }

  typedef void ( *BASIS_FUNC_TYPE ) ( const RealVecChart &c, RealVecChart & dest );
  BASIS_FUNC_TYPE _basis[3];
  typedef void ( *BASIS_FUNC_TYPE_GRAD ) ( const RealVecChart &c, RealVecChart &dest );
  BASIS_FUNC_TYPE_GRAD _Grad_basis[6];

public:
  GlobalAffineGradBaseFunctionSet3D( ) {
    
    _basis[0] = _b1; _basis[1] = _b2; _basis[2] = _b3;
    _Grad_basis[0] = _Grad_b1; _Grad_basis[1] = _Grad_b2; _Grad_basis[2] = _Grad_b3;
  }

  static const int numBaseFuncs = 3;
  
  void evaluate ( int BaseFuncNum, const RealVecChart &c, RealVecChart &dest ) const {_basis[BaseFuncNum] ( c, dest );}
  void evaluateGrad ( int BaseFuncNum, const RealVecChart &c, RealVecChart &dest ) const { _Grad_basis[BaseFuncNum] ( c, dest );}

};



//! ========================================================================================
//! global affine basefunction set for functions 
//! \f$ u : R^d \to \R^{d x d}_{sym} \f$
//! =======================================================================================
template <typename DataTypeContainer>
class GlobalAffineSymGradBaseFunctionSet2D {
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
  typedef typename DataTypeContainer::DerivativeVectorValuedType DerivativeVectorValuedType;
    
  // b1(x) = (x,0)
  // b2(x) = (0,y)
  // b3(x) = (y,x)
  
  static void _b1  ( const RealVecChart &c, RealVecChart &dest ) { dest[0] = c[0]; dest[1] = 0.0; }
  static void _b2  ( const RealVecChart &c, RealVecChart &dest ) { dest[0] = 0.0;  dest[1] = c[1]; }
  static void _b3  ( const RealVecChart &c, RealVecChart &dest ) { dest[0] = c[1]; dest[1] = c[0]; }
  
  static void _symGrad_b1  ( const RealVecChart &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(0,0) = 1.0; }
  static void _symGrad_b2  ( const RealVecChart &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(1,1) = 1.0; }
  static void _symGrad_b3  ( const RealVecChart &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(0,1) = 1.0; dest(1,0) = 1.0; }
  
  static RealType _div_b1  ( const RealVecChart &c ) { return 1.; }
  static RealType _div_b2  ( const RealVecChart &c ) { return 1.; }
  static RealType _div_b3  ( const RealVecChart &c ) { return 0.; }
  
//   static RealType _dx_b1   ( const RealVecChart &c ) { return c[1] - 1.0; }
//   static RealType _dx_b2   ( const RealVecChart &c ) { return 1.0 - c[1]; }
//   static RealType _dx_b3   ( const RealVecChart &c ) { return -c[1]; }
// 
//   static RealType _dy_b1   ( const RealVecChart &c ) { return c[0] - 1.0; }
//   static RealType _dy_b2   ( const RealVecChart &c ) { return -c[0]; }
//   static RealType _dy_b3   ( const RealVecChart &c ) { return 1.0 - c[0]; }

  typedef void ( *BASIS_FUNC_TYPE ) ( const RealVecChart &c, RealVecChart &dest );
  BASIS_FUNC_TYPE _basis[3];
  typedef void ( *BASIS_FUNC_TYPE_SYMGRAD ) ( const RealVecChart &c, DerivativeVectorValuedType &dest );
  BASIS_FUNC_TYPE_SYMGRAD _symGrad_basis[3];
  typedef RealType ( *BASIS_FUNC_TYPE_DIV ) ( const RealVecChart &c );
  BASIS_FUNC_TYPE_DIV _div_basis[3];

public:
  GlobalAffineSymGradBaseFunctionSet2D ( ) {
    _basis[0] = _b1; _basis[1] = _b2; _basis[2] = _b3;
    _symGrad_basis[0] = _symGrad_b1; _symGrad_basis[1] = _symGrad_b2; _symGrad_basis[2] = _symGrad_b3;
    _div_basis[0] = _div_b1; _div_basis[1] = _div_b2; _div_basis[2] = _div_b3;
  }

  static const int numBaseFuncs = 3;
  
  void evaluate ( int BaseFuncNum, const RealVecChart &c, RealVecChart &dest ) const { _basis[BaseFuncNum] ( c, dest );}
  void evaluateSymGrad ( int BaseFuncNum, const RealVecChart &c, DerivativeVectorValuedType &dest ) const { 
      _symGrad_basis[BaseFuncNum] ( c, dest );
  }
  RealType evaluateDiv ( int BaseFuncNum, const RealVecChart &c ) const { return _div_basis[BaseFuncNum] (c);}

};


// base function set given by 
// (x,0,0), (0,y,0), (0,0,z), (y,x,0), (0,z,y), (z,0,x)
template <typename DataTypeContainer>
class GlobalAffineSymGradBaseFunctionSet3D{

  typedef typename DataTypeContainer::RealType                      RealType;  
  typedef typename DataTypeContainer::RealVecChart                  RealVecChart;
  typedef typename DataTypeContainer::DerivativeVectorValuedType    DerivativeVectorValuedType;
  
  static void _b1   ( const RealVecChart &c, RealVecChart & dest ) { dest[0] = c[0]; dest[1] = 0.0; dest[2] = 0.0; }
  static void _b2   ( const RealVecChart &c, RealVecChart & dest ) { dest[0] = 0.0; dest[1] = c[1]; dest[2] = 0.0; }
  static void _b3   ( const RealVecChart &c, RealVecChart & dest ) { dest[0] = 0.0; dest[1] = 0.0; dest[2] = c[2]; }
  static void _b4   ( const RealVecChart &c, RealVecChart & dest ) { dest[0] = c[1]; dest[1] = c[0]; dest[2] = 0.0; }
  static void _b6   ( const RealVecChart &c, RealVecChart & dest ) { dest[0] = c[2]; dest[1] = 0.0; dest[2] = c[0]; }
  static void _b5   ( const RealVecChart &c, RealVecChart & dest ) { dest[0] = 0.0; dest[1] = c[2]; dest[2] = c[1]; }
  
  static void _symGrad_b1  ( const RealVecChart &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(0,0) = 1.0; }
  static void _symGrad_b2  ( const RealVecChart &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(1,1) = 1.0; }
  static void _symGrad_b3  ( const RealVecChart &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(2,2) = 1.0; }
  static void _symGrad_b4  ( const RealVecChart &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(0,1) = 1.0; dest(1,0) = 1.0; }
  static void _symGrad_b6  ( const RealVecChart &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(0,2) = 1.0; dest(2,0) = 1.0; }
  static void _symGrad_b5  ( const RealVecChart &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(1,2) = 1.0; dest(2,1) = 1.0; }
  
  static RealType _div_b1  ( const RealVecChart &c ) { return 1.; }
  static RealType _div_b2  ( const RealVecChart &c ) { return 1.; }
  static RealType _div_b3  ( const RealVecChart &c ) { return 1.; }
  static RealType _div_b4  ( const RealVecChart &c ) { return 0.; }
  static RealType _div_b5  ( const RealVecChart &c ) { return 0.; }
  static RealType _div_b6  ( const RealVecChart &c ) { return 0.; }

  typedef void ( *BASIS_FUNC_TYPE ) ( const RealVecChart &c, RealVecChart & dest );
  BASIS_FUNC_TYPE _basis[6];
  typedef void ( *BASIS_FUNC_TYPE_SYMGRAD ) ( const RealVecChart &c, DerivativeVectorValuedType &dest );
  BASIS_FUNC_TYPE_SYMGRAD _symGrad_basis[6];
  typedef RealType ( *BASIS_FUNC_TYPE_DIV ) ( const RealVecChart &c );
  BASIS_FUNC_TYPE_DIV _div_basis[6];

public:
  GlobalAffineSymGradBaseFunctionSet3D( ) {
    _basis[0] = _b1; _basis[1] = _b2; _basis[2] = _b3; _basis[3] = _b4; _basis[4] = _b5; _basis[5] = _b6;
    _symGrad_basis[0] = _symGrad_b1; _symGrad_basis[1] = _symGrad_b2; _symGrad_basis[2] = _symGrad_b3; _symGrad_basis[3] = _symGrad_b4; _symGrad_basis[4] = _symGrad_b5; _symGrad_basis[5] = _symGrad_b6;
    _div_basis[0] = _div_b1; _div_basis[1] = _div_b2; _div_basis[2] = _div_b3; _div_basis[3] = _div_b4; _div_basis[4] = _div_b5; _div_basis[5] = _div_b6;
  }

  static const int numBaseFuncs = 6;
  
  void evaluate ( int BaseFuncNum, const RealVecChart &c, RealVecChart &dest ) const {_basis[BaseFuncNum] ( c, dest );}
  void evaluateSymGrad ( int BaseFuncNum, const RealVecChart &c, DerivativeVectorValuedType &dest ) const { _symGrad_basis[BaseFuncNum] ( c, dest );}
  RealType evaluateDiv ( int BaseFuncNum, const RealVecChart &c ) const { return _div_basis[BaseFuncNum] (c);}

};







#endif
