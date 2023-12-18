#ifndef __QUOCFESPLINEBASEFUNCTIONSETS_H
#define __QUOCFESPLINEBASEFUNCTIONSETS_H

#include <feBaseFunctionSet.h>


template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocSplineBaseFunctionSet1D 
: public FECachedSecondOrderBaseFunctionSetBase < DataTypeContainer, QuadType,  4,  QuocSplineBaseFunctionSet1D<DataTypeContainer, QuadType, ElementType> > {

public:
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
  typedef typename DataTypeContainer::DerivativeVectorValuedType DerivativeVectorValuedType;
  
  typedef FECachedSecondOrderBaseFunctionSetBase < DataTypeContainer, QuadType,  4,  QuocSplineBaseFunctionSet1D<DataTypeContainer, QuadType, ElementType> > BaseType;
  
  static const int numBaseFuncs = 4;
  static const unsigned short int numberOfDifferentBaseFunctionSets = 3;
  
private:
  
  //! basis functions for inner element
  static RealType _b1 ( const RealType RefCoord ) {  return 1./6. * pesopt::Cub( 1. - RefCoord ); }
  static RealType _b2 ( const RealType RefCoord ) { return 1./6. * ( pesopt::Cub( 2. - RefCoord ) - 4. * pesopt::Cub( 1. - RefCoord ) ); }
  static RealType _b3 ( const RealType RefCoord ) { return 1./6. * ( pesopt::Cub( 1. + RefCoord ) - 4. * pesopt::Cub( RefCoord ) ); }
  static RealType _b4 ( const RealType RefCoord ) { return 1./6. * pesopt::Cub( RefCoord ); }

  static RealType _d_b1 ( const RealType RefCoord ) { return -0.5 * pesopt::Sqr( 1. - RefCoord ); }
  static RealType _d_b2 ( const RealType RefCoord ) { return -0.5 * ( pesopt::Sqr( 2. - RefCoord ) - 4. * pesopt::Sqr( 1. - RefCoord ) ); }
  static RealType _d_b3 ( const RealType RefCoord ) { return 0.5 * ( pesopt::Sqr( 1. + RefCoord ) - 4. * pesopt::Sqr( RefCoord ) ); }
  static RealType _d_b4 ( const RealType RefCoord ) { return 0.5 * pesopt::Sqr( RefCoord ); }

  static RealType _dd_b1 ( const RealType RefCoord ) { return ( 1. - RefCoord ); }
  static RealType _dd_b2 ( const RealType RefCoord ) { return ( -2. + 3. * RefCoord ); }
  static RealType _dd_b3 ( const RealType RefCoord ) { return ( 1. - 3. * RefCoord ); }
  static RealType _dd_b4 ( const RealType RefCoord ) { return ( RefCoord ); }

  static RealType _ddd_b1 ( const RealType /*RefCoord*/ ) { return -1.; }
  static RealType _ddd_b2 ( const RealType /*RefCoord*/ ) { return 3.; }
  static RealType _ddd_b3 ( const RealType /*RefCoord*/ ) { return -3.; }
  static RealType _ddd_b4 ( const RealType /*RefCoord*/ ) { return 1.; }


  //! basis functions for element at left boundary
  static RealType _b1_left ( const RealType RefCoord ) { return _b2 ( RefCoord ) + 2. * _b1 ( RefCoord ); }
  static RealType _b2_left ( const RealType RefCoord ) { return _b3 ( RefCoord ) - _b1 ( RefCoord ); }
  static RealType _b3_left ( const RealType RefCoord ) { return _b4 ( RefCoord ); }
  static RealType _b4_left ( const RealType /*RefCoord*/ ) { return 0.; }

  static RealType _d_b1_left ( const RealType RefCoord ) { return _d_b2 ( RefCoord ) + 2. * _d_b1 ( RefCoord );}
  static RealType _d_b2_left ( const RealType RefCoord ) { return _d_b3 ( RefCoord ) - _d_b1 ( RefCoord ); }
  static RealType _d_b3_left ( const RealType RefCoord ) { return _d_b4 ( RefCoord ); }
  static RealType _d_b4_left ( const RealType /*RefCoord*/ ) { return 0.; }

  static RealType _dd_b1_left ( const RealType RefCoord ) { return _dd_b2 ( RefCoord ) + 2. * _dd_b1 ( RefCoord ); }
  static RealType _dd_b2_left ( const RealType RefCoord ) { return _dd_b3 ( RefCoord ) - _dd_b1 ( RefCoord ); }
  static RealType _dd_b3_left ( const RealType RefCoord ) { return _dd_b4 ( RefCoord ); }
  static RealType _dd_b4_left ( const RealType /*RefCoord*/ ) { return 0.; }

  static RealType _ddd_b1_left ( const RealType RefCoord ) { return _ddd_b2 ( RefCoord ) + 2. * _ddd_b1 ( RefCoord ); }
  static RealType _ddd_b2_left ( const RealType RefCoord ) { return _ddd_b3 ( RefCoord ) - _ddd_b1 ( RefCoord ); }
  static RealType _ddd_b3_left ( const RealType RefCoord ) { return _ddd_b4 ( RefCoord ); }
  static RealType _ddd_b4_left ( const RealType ) { return 0.; }


  //! basis functions for element at right boundary  
  static RealType _b1_right ( const RealType /*RefCoord*/ ) { return 0.; }
  static RealType _b2_right ( const RealType RefCoord ) { return _b1 ( RefCoord );}
  static RealType _b3_right ( const RealType RefCoord ) { return _b2 ( RefCoord ) - _b4 ( RefCoord ); }
  static RealType _b4_right ( const RealType RefCoord ) { return _b3 ( RefCoord ) + 2. * _b4 ( RefCoord ); }


  static RealType _d_b1_right ( const RealType /*RefCoord*/ ) { return 0.;}
  static RealType _d_b2_right ( const RealType RefCoord ) { return _d_b1 ( RefCoord ); }
  static RealType _d_b3_right ( const RealType RefCoord ) { return _d_b2 ( RefCoord ) - _d_b4 ( RefCoord ); }
  static RealType _d_b4_right ( const RealType RefCoord ) { return _d_b3 ( RefCoord ) + 2. * _d_b4 ( RefCoord ); }


  static RealType _dd_b1_right ( const RealType /*RefCoord*/ ) { return 0.; }
  static RealType _dd_b2_right ( const RealType RefCoord ) { return _dd_b1 ( RefCoord ); }
  static RealType _dd_b3_right ( const RealType RefCoord ) { return _dd_b2 ( RefCoord ) - _dd_b4 ( RefCoord ); }
  static RealType _dd_b4_right ( const RealType RefCoord ) { return _dd_b3 ( RefCoord ) + 2. * _dd_b4 ( RefCoord ); }

  static RealType _ddd_b1_right ( const RealType /*RefCoord*/ ) { return 0.; }
  static RealType _ddd_b2_right ( const RealType RefCoord ) { return _ddd_b1 ( RefCoord ); }
  static RealType _ddd_b3_right ( const RealType RefCoord ) { return _ddd_b2 ( RefCoord ) - _ddd_b4 ( RefCoord ); }
  static RealType _ddd_b4_right ( const RealType RefCoord ) { return _ddd_b3 ( RefCoord ) + 2. * _ddd_b4 ( RefCoord ); }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const RealType RefCoord );
  BASIS_FUNC_TYPE _thirdDeriv_basis[4];
  BASIS_FUNC_TYPE _secondDeriv_basis[4];
  BASIS_FUNC_TYPE _deriv_basis[4];
  BASIS_FUNC_TYPE _basis[4];  

  const RealType _hx;
  const unsigned short _splineElementType;
  
public:
  QuocSplineBaseFunctionSet1D ( const RealType hx, const unsigned short SplineElementType ) 
    : _hx ( hx ), _splineElementType ( SplineElementType ) {

    switch ( SplineElementType ){
    // element at left boundary
    case 0:{
      _basis[0] = _b1_left;
      _basis[1] = _b2_left;
      _basis[2] = _b3_left;
      _basis[3] = _b4_left;

      _deriv_basis[0] = _d_b1_left;
      _deriv_basis[1] = _d_b2_left;
      _deriv_basis[2] = _d_b3_left;
      _deriv_basis[3] = _d_b4_left;

      _secondDeriv_basis[0] = _dd_b1_left;
      _secondDeriv_basis[1] = _dd_b2_left;
      _secondDeriv_basis[2] = _dd_b3_left;
      _secondDeriv_basis[3] = _dd_b4_left;

      _thirdDeriv_basis[0] = _ddd_b1_left;
      _thirdDeriv_basis[1] = _ddd_b2_left;
      _thirdDeriv_basis[2] = _ddd_b3_left;
      _thirdDeriv_basis[3] = _ddd_b4_left;
      break;
    }
    // inner elements
    case 1:{
      _basis[0] = _b1;
      _basis[1] = _b2;
      _basis[2] = _b3;
      _basis[3] = _b4;

      _deriv_basis[0] = _d_b1;
      _deriv_basis[1] = _d_b2;
      _deriv_basis[2] = _d_b3;
      _deriv_basis[3] = _d_b4;

      _secondDeriv_basis[0] = _dd_b1;
      _secondDeriv_basis[1] = _dd_b2;
      _secondDeriv_basis[2] = _dd_b3;
      _secondDeriv_basis[3] = _dd_b4;

      _thirdDeriv_basis[0] = _ddd_b1;
      _thirdDeriv_basis[1] = _ddd_b2;
      _thirdDeriv_basis[2] = _ddd_b3;
      _thirdDeriv_basis[3] = _ddd_b4;
      break;
    }
    // element at right boundary
    case 2:{
      _basis[0] = _b1_right;
      _basis[1] = _b2_right;
      _basis[2] = _b3_right;
      _basis[3] = _b4_right;

      _deriv_basis[0] = _d_b1_right;
      _deriv_basis[1] = _d_b2_right;
      _deriv_basis[2] = _d_b3_right;
      _deriv_basis[3] = _d_b4_right;

      _secondDeriv_basis[0] = _dd_b1_right;
      _secondDeriv_basis[1] = _dd_b2_right;
      _secondDeriv_basis[2] = _dd_b3_right;
      _secondDeriv_basis[3] = _dd_b4_right;

      _thirdDeriv_basis[0] = _ddd_b1_right;
      _thirdDeriv_basis[1] = _ddd_b2_right;
      _thirdDeriv_basis[2] = _ddd_b3_right;
      _thirdDeriv_basis[3] = _ddd_b4_right;
      break;
    }
    default:
      cout << "SplineElementType = " << SplineElementType << endl;
      throw std::invalid_argument( pesopt::strprintf ( "Wrong SplineElementType in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      break;

    }

    this->initializeQuadCache( );
  }


  RealType evaluate ( int BaseFuncNum, const RealVecChart &c ) const { return _basis[BaseFuncNum] ( c[0] );}
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
      return BaseType::evaluate ( BaseFuncNum, QuadPoint );
  }
  void evaluateGradient ( int BaseFuncNum, const RealVecChart &RefCoord, RealVecChart &Gradient ) const {
      Gradient[0] = _deriv_basis[BaseFuncNum] ( RefCoord[0] ) / _hx;
  }
  inline const RealVecChart& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateGradient ( BaseFuncNum, QuadPoint );
  }
  void evaluateHessian ( int BaseFuncNum, const RealVecChart &RefCoord, DerivativeVectorValuedType &Hessian ) const {
      Hessian(0, 0) = _secondDeriv_basis[BaseFuncNum] ( RefCoord[0] ) / pesopt::Sqr( _hx );
  }
  inline const DerivativeVectorValuedType & evaluateHessian ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateHessian ( BaseFuncNum, QuadPoint );
  }

};















template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocSplineBaseFunctionSet2D 
: public FECachedSecondOrderBaseFunctionSetBase < DataTypeContainer, QuadType,  16,  QuocSplineBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType> > {
    
protected:
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
  typedef typename DataTypeContainer::DerivativeVectorValuedType DerivativeVectorValuedType;
    
public:
  typedef FECachedSecondOrderBaseFunctionSetBase < DataTypeContainer, QuadType,  16,  QuocSplineBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType> > BaseType;
  
  typedef Quoc1DDataTypeContainer DataTypeContainer1D;
  typedef typename DataTypeContainer1D::RealVecChart RealVecChart1D;
  typedef typename DataTypeContainer1D::DerivativeVectorValuedType DerivativeVectorValuedType1D;
  typedef typename QuadType::QuadRuleType1D QuadType1D;
  typedef QuocSplineBaseFunctionSet1D < DataTypeContainer1D, QuadType1D, QuocElement1D<DataTypeContainer1D> > BaseFunc1DType;
  
  static const unsigned short int numberOfDifferentBaseFunctionSets = 9;
  static const int numBaseFuncs = 16;
  
protected:
    
  const RealType _hx, _hy;
  const unsigned short _splineElementTypeX, _splineElementTypeY;
  const BaseFunc1DType _baseX, _baseY;

public:
  QuocSplineBaseFunctionSet2D ( const RealType hx, const RealType hy, 
                                const unsigned short SplineElementTypeX, 
                                const unsigned short SplineElementTypeY )
    : _hx ( hx ), _hy ( hy ), 
    _splineElementTypeX ( SplineElementTypeX ),
    _splineElementTypeY ( SplineElementTypeY ), 
    _baseX ( _hx, _splineElementTypeX ), 
    _baseY ( _hy, _splineElementTypeY ) {
        this->initializeQuadCache( );
    }


  RealType evaluate ( int BaseFuncNum, const RealVecChart &RefCoord ) const {
    RealVecChart1D RefCoordX,  RefCoordY;
    RefCoordX[0] = RefCoord[0];
    RefCoordY[0] = RefCoord[1];  
    return _baseX.evaluate ( getXBaseFuncNum ( BaseFuncNum ), RefCoordX ) * _baseY.evaluate ( getYBaseFuncNum ( BaseFuncNum ), RefCoordY );
  }
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
      return BaseType::evaluate ( BaseFuncNum, QuadPoint );
  }

  void evaluateGradient ( int BaseFuncNum, const RealVecChart &RefCoord, RealVecChart &Gradient ) const {
    const int xBaseFuncNum ( getXBaseFuncNum ( BaseFuncNum ) );
    const int yBaseFuncNum ( getYBaseFuncNum ( BaseFuncNum ) );
    RealVecChart1D RefCoordX,  RefCoordY;
    RefCoordX[0] = RefCoord[0];
    RefCoordY[0] = RefCoord[1];
    RealVecChart1D GradX; _baseX.evaluateGradient ( xBaseFuncNum, RefCoordX,  GradX );
    RealVecChart1D GradY; _baseY.evaluateGradient ( yBaseFuncNum, RefCoordY,  GradY );
    Gradient[0] = GradX[0] * _baseY.evaluate ( yBaseFuncNum, RefCoordY );
    Gradient[1] = _baseX.evaluate ( xBaseFuncNum, RefCoordX ) * GradY[0];
  }
  inline const RealVecChart& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateGradient ( BaseFuncNum, QuadPoint );
  }
  
  void evaluateHessian ( int BaseFuncNum, const RealVecChart &RefCoord, DerivativeVectorValuedType &Hessian ) const {
    const int xBaseFuncNum ( getXBaseFuncNum ( BaseFuncNum ) );
    const int yBaseFuncNum ( getYBaseFuncNum ( BaseFuncNum ) );
    RealVecChart1D RefCoordX,  RefCoordY;
    RefCoordX[0] = RefCoord[0];
    RefCoordY[0] = RefCoord[1];
    RealVecChart1D GradX; _baseX.evaluateGradient ( xBaseFuncNum, RefCoordX,  GradX );
    RealVecChart1D GradY; _baseY.evaluateGradient ( yBaseFuncNum, RefCoordY,  GradY );
    DerivativeVectorValuedType1D HessianX; _baseX.evaluateHessian ( xBaseFuncNum, RefCoordX,  HessianX );
    DerivativeVectorValuedType1D HessianY; _baseY.evaluateHessian ( yBaseFuncNum, RefCoordY,  HessianY );
    Hessian(0, 0) = HessianX(0, 0) * _baseY.evaluate ( yBaseFuncNum, RefCoordY );
    Hessian(1, 0) = Hessian(0, 1) = GradX[0] * GradY[0];
    Hessian(1, 1) = _baseX.evaluate ( xBaseFuncNum, RefCoordX ) * HessianY(0, 0);
  }
  inline const DerivativeVectorValuedType& evaluateHessian ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateHessian ( BaseFuncNum, QuadPoint );
  }
  

private:
  int getXBaseFuncNum ( const int BaseFuncNum ) const { return BaseFuncNum % 4; }
  int getYBaseFuncNum ( const int BaseFuncNum ) const { return BaseFuncNum / 4; }
};






template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocSplineBaseFunctionSet3D 
: public FECachedSecondOrderBaseFunctionSetBase < DataTypeContainer, QuadType,  64,  QuocSplineBaseFunctionSet3D<DataTypeContainer, QuadType, ElementType> > {
    
protected:
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
  typedef typename DataTypeContainer::DerivativeVectorValuedType DerivativeVectorValuedType;
    
public:
  typedef FECachedSecondOrderBaseFunctionSetBase < DataTypeContainer, QuadType,  64,  QuocSplineBaseFunctionSet3D<DataTypeContainer, QuadType, ElementType> > BaseType;

  typedef Quoc1DDataTypeContainer DataTypeContainer1D;
  typedef typename DataTypeContainer1D::RealVecChart RealVecChart1D;
  typedef typename DataTypeContainer1D::DerivativeVectorValuedType DerivativeVectorValuedType1D;
  typedef typename QuadType::QuadRuleType1D QuadType1D;
  typedef QuocSplineBaseFunctionSet1D < DataTypeContainer1D, QuadType1D, QuocElement1D<DataTypeContainer1D> > BaseFunc1DType;

  static const unsigned short int numberOfDifferentBaseFunctionSets = 27;
  static const int numBaseFuncs = 64;
  

private:
  const RealType _hx, _hy, _hz;
  const unsigned short _splineElementTypeX, _splineElementTypeY, _splineElementTypeZ;
  const BaseFunc1DType _baseX, _baseY, _baseZ;

public:
  QuocSplineBaseFunctionSet3D ( const RealType hx,  const RealType hy,  const RealType hz,
                              const unsigned short SplineElementTypeX, 
                              const unsigned short SplineElementTypeY,
                              const unsigned short SplineElementTypeZ )
    : _hx ( hx ), _hy ( hy ),  _hz ( hz ), 
      _splineElementTypeX ( SplineElementTypeX ), 
      _splineElementTypeY ( SplineElementTypeY ), 
      _splineElementTypeZ ( SplineElementTypeZ ),
      _baseX ( _hx, _splineElementTypeX ), 
      _baseY ( _hy, _splineElementTypeY ), 
      _baseZ ( _hz, _splineElementTypeZ ) {
    this->initializeQuadCache( );
  }

  
  RealType evaluate ( int BaseFuncNum, const RealVecChart &RefCoord ) const {
    RealVecChart1D RefCoordX,  RefCoordY,  RefCoordZ;
    RefCoordX[0] = RefCoord[0];
    RefCoordY[0] = RefCoord[1];
    RefCoordZ[0] = RefCoord[2];
    return _baseX.evaluate ( getXBaseFuncNum ( BaseFuncNum ), RefCoordX ) * _baseY.evaluate ( getYBaseFuncNum ( BaseFuncNum ), RefCoordY ) * _baseZ.evaluate ( getZBaseFuncNum ( BaseFuncNum ), RefCoordZ );
  }
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluate ( BaseFuncNum, QuadPoint );
  }

  void evaluateGradient ( int BaseFuncNum, const RealVecChart &RefCoord, RealVecChart &Gradient ) const {
    const int xBaseFuncNum ( getXBaseFuncNum ( BaseFuncNum ) );
    const int yBaseFuncNum ( getYBaseFuncNum ( BaseFuncNum ) );
    const int zBaseFuncNum ( getZBaseFuncNum ( BaseFuncNum ) );
    RealVecChart1D RefCoordX,  RefCoordY,  RefCoordZ;
    RefCoordX[0] = RefCoord[0];
    RefCoordY[0] = RefCoord[1];
    RefCoordZ[0] = RefCoord[2];
    RealVecChart1D GradX; _baseX.evaluateGradient ( xBaseFuncNum, RefCoordX,  GradX );
    RealVecChart1D GradY; _baseY.evaluateGradient ( yBaseFuncNum, RefCoordY,  GradY );
    RealVecChart1D GradZ; _baseZ.evaluateGradient ( zBaseFuncNum, RefCoordZ,  GradZ );
    Gradient[0] =  GradX(0) * _baseY.evaluate ( yBaseFuncNum, RefCoordY ) * _baseZ.evaluate ( zBaseFuncNum, RefCoordZ );
    Gradient[1] = _baseX.evaluate ( xBaseFuncNum, RefCoordX ) * GradY(0) * _baseZ.evaluate ( zBaseFuncNum, RefCoordZ );
    Gradient[2] = _baseX.evaluate ( xBaseFuncNum, RefCoordX ) * _baseY.evaluate ( yBaseFuncNum, RefCoordY ) * GradZ(0);
  }
  inline const RealVecChart& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateGradient ( BaseFuncNum, QuadPoint );
  }
  
  void evaluateHessian ( int BaseFuncNum, const RealVecChart &RefCoord, DerivativeVectorValuedType &Hessian ) const {
    const int xBaseFuncNum ( getXBaseFuncNum ( BaseFuncNum ) );
    const int yBaseFuncNum ( getYBaseFuncNum ( BaseFuncNum ) );
    const int zBaseFuncNum ( getZBaseFuncNum ( BaseFuncNum ) );
    RealVecChart1D RefCoordX,  RefCoordY,  RefCoordZ;
    RefCoordX[0] = RefCoord[0];
    RefCoordY[0] = RefCoord[1];
    RefCoordZ[0] = RefCoord[2];
    RealVecChart1D GradX; _baseX.evaluateGradient ( xBaseFuncNum, RefCoordX,  GradX );
    RealVecChart1D GradY; _baseY.evaluateGradient ( yBaseFuncNum, RefCoordY,  GradY );
    RealVecChart1D GradZ; _baseZ.evaluateGradient ( zBaseFuncNum, RefCoordZ,  GradZ );
    DerivativeVectorValuedType1D HessianX; _baseX.evaluateHessian ( xBaseFuncNum, RefCoordX,  HessianX );
    DerivativeVectorValuedType1D HessianY; _baseY.evaluateHessian ( yBaseFuncNum, RefCoordY,  HessianY );
    DerivativeVectorValuedType1D HessianZ; _baseZ.evaluateHessian ( zBaseFuncNum, RefCoordZ,  HessianZ );
    Hessian(0, 0) = HessianX(0, 0) * _baseY.evaluate ( yBaseFuncNum, RefCoordY ) * _baseZ.evaluate ( zBaseFuncNum, RefCoordZ );
    Hessian(1, 0) = Hessian(0, 1) = GradX(0) * GradY(0) * _baseZ.evaluate ( zBaseFuncNum, RefCoordZ );
    Hessian(1, 1) = _baseX.evaluate ( xBaseFuncNum, RefCoordX ) * HessianY(0, 0) * _baseZ.evaluate ( zBaseFuncNum, RefCoordZ );
    Hessian(2, 0) = Hessian(0, 2) = GradX(0) * _baseY.evaluate ( yBaseFuncNum, RefCoordY ) * GradZ(0);
    Hessian(2, 1) = Hessian(1, 2) = _baseX.evaluate ( xBaseFuncNum, RefCoordX ) * GradY(0) * GradZ(0);
    Hessian(2, 2) = _baseX.evaluate ( xBaseFuncNum, RefCoordX ) * _baseY.evaluate ( yBaseFuncNum, RefCoordY ) * HessianZ(0, 0);
  }
  inline const DerivativeVectorValuedType& evaluateHessian ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateHessian ( BaseFuncNum, QuadPoint );
  }



private:
  int getXBaseFuncNum ( const int BaseFuncNum ) const { return BaseFuncNum % 4; }
  int getYBaseFuncNum ( const int BaseFuncNum ) const { return ( BaseFuncNum % 16 ) / 4;}
  int getZBaseFuncNum ( const int BaseFuncNum ) const { return BaseFuncNum / 16; }
};



#endif //__QUOCSPLINEBASEFUNCTIONSETS_H
