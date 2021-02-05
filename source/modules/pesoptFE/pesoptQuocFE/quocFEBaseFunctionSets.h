#ifndef __QUOCFEBASEFUNCTIONSETS_H
#define __QUOCFEBASEFUNCTIONSETS_H

#include <feBaseFunctionSet.h>


template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocBaseFunctionSet1D 
: public FEBaseFunctionSetBase< DataTypeContainer, QuadType >  {

  typedef typename DataTypeContainer::RealType RealType;

  static RealType _b1   ( const RealType &c ) { return ( 1. - c ); }
  static RealType _b2   ( const RealType &c ) { return c[0]; }
  
  static RealType _dx_b1   ( const RealType &/*c*/ ) { return - 1.; }
  static RealType _dx_b2   ( const RealType &/*c*/ ) { return 1.; }
  

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const RealType &RefCoord );
  BASIS_FUNC_TYPE _basis[2];
  BASIS_FUNC_TYPE _deriv_basis[2];
  const RealType _hx;
  
public:
  QuocBaseFunctionSet1D( const RealType hx ) : _hx ( hx ){
    _basis[0] = _b1;
    _basis[1] = _b2;

    _deriv_basis[0] = _dx_b1;
    _deriv_basis[1] = _dx_b2;
  }

  static const int numBaseFuncs = 2;

  void setElement ( const ElementType &/*El*/ ) { }

  RealType evaluate ( int BaseFuncNum, const RealType &RefCoord ) const { return _basis[BaseFuncNum] ( RefCoord );}
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const { return evaluate ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );}
  
  void evaluateGradient ( int BaseFuncNum, const RealType &RefCoord, RealType &Gradient ) const {Gradient = _deriv_basis[BaseFuncNum] ( RefCoord ) / _hx;}
  inline RealType evaluateGradient ( int BaseFuncNum, int QuadPoint ) const { return RealType( _deriv_basis[BaseFuncNum] ( this->_quadRule.getRefCoord( QuadPoint ) ) / _hx );}
};



//! The basefunctionset for bilinear elements in 2d
template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocBaseFunctionSet2D : public FEBaseFunctionSetBase< DataTypeContainer, QuadType >{
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
    
  // b1 (0,0)
  // b2 (1,0)
  // b3 (0,1)
  // b4 (1,1)
  
  static RealType _b1   ( const RealVecChart &c ) { return ( 1.0 - c[0] ) * ( 1.0 - c[1] ); }
  static RealType _b2   ( const RealVecChart &c ) { return c[0] * ( 1.0 - c[1] ); }
  static RealType _b3   ( const RealVecChart &c ) { return ( 1.0 - c[0] ) * c[1]; }
  static RealType _b4   ( const RealVecChart &c ) { return c[0]*c[1]; }
  
  static RealType _dx_b1   ( const RealVecChart &c ) { return c[1] - 1.0; }
  static RealType _dx_b2   ( const RealVecChart &c ) { return 1.0 - c[1]; }
  static RealType _dx_b3   ( const RealVecChart &c ) { return -c[1]; }
  static RealType _dx_b4   ( const RealVecChart &c ) { return c[1]; }

  static RealType _dy_b1   ( const RealVecChart &c ) { return c[0] - 1.0; }
  static RealType _dy_b2   ( const RealVecChart &c ) { return -c[0]; }
  static RealType _dy_b3   ( const RealVecChart &c ) { return 1.0 - c[0]; }
  static RealType _dy_b4   ( const RealVecChart &c ) { return c[0]; }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const RealVecChart &c );
  BASIS_FUNC_TYPE _basis[4];
  BASIS_FUNC_TYPE _deriv_basis[2][4];
  const RealType _hx, _hy;

public:
  QuocBaseFunctionSet2D ( const RealType hx, const RealType hy ) : _hx ( hx ), _hy ( hy ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
  }

  static const int numBaseFuncs = 4;
  
  RealType evaluate ( int BaseFuncNum, const RealVecChart &c ) const { return _basis[BaseFuncNum] ( c );}
  inline const RealType evaluate( int BaseFuncNum, int QuadPoint ) const { return evaluate( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );}
  
  void evaluateGradient ( int BaseFuncNum, const RealVecChart &c, RealVecChart &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( c ) / _hx;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( c ) / _hy;
  }

  inline RealVecChart evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return RealVecChart( _deriv_basis[0][BaseFuncNum]( this->_quadRule.getRefCoord( QuadPoint ) ) / _hx, 
                       _deriv_basis[1][BaseFuncNum]( this->_quadRule.getRefCoord( QuadPoint ) ) / _hy );
  }

};


//! The basefunctionset for trilinear elements in 3d.
template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocBaseFunctionSet3D : public FEBaseFunctionSetBase< DataTypeContainer, QuadType >{

  typedef typename DataTypeContainer::RealType RealType;  
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
    
  static inline RealType _b1_1d ( RealType x ) { return 1. - x; }
  static inline RealType _b2_1d ( RealType x ) { return x; }

  static inline RealType _d_b1_1d ( RealType ) { return -1.; }
  static inline RealType _d_b2_1d ( RealType ) { return 1.; }

  // b1 (0,0,0)
  // b2 (1,0,0)
  // b3 (0,1,0)
  // b4 (1,1,0)
  // b5 (0,0,1)
  // b6 (1,0,1)
  // b7 (0,1,1)
  // b8 (1,1,1)
  
  static RealType _b1   ( const RealVecChart &c ) { return _b1_1d ( c[0] ) * _b1_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b2   ( const RealVecChart &c ) { return _b2_1d ( c[0] ) * _b1_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b3   ( const RealVecChart &c ) { return _b1_1d ( c[0] ) * _b2_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b4   ( const RealVecChart &c ) { return _b2_1d ( c[0] ) * _b2_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b5   ( const RealVecChart &c ) { return _b1_1d ( c[0] ) * _b1_1d ( c[1] ) * _b2_1d ( c[2] ); }
  static RealType _b6   ( const RealVecChart &c ) { return _b2_1d ( c[0] ) * _b1_1d ( c[1] ) * _b2_1d ( c[2] ); }
  static RealType _b7   ( const RealVecChart &c ) { return _b1_1d ( c[0] ) * _b2_1d ( c[1] ) * _b2_1d ( c[2] ); }
  static RealType _b8   ( const RealVecChart &c ) { return _b2_1d ( c[0] ) * _b2_1d ( c[1] ) * _b2_1d ( c[2] ); }

  static RealType _dx_b1 ( const RealVecChart &c ) { return _d_b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b2 ( const RealVecChart &c ) { return _d_b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b3 ( const RealVecChart &c ) { return _d_b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b4 ( const RealVecChart &c ) { return _d_b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b5 ( const RealVecChart &c ) { return _d_b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dx_b6 ( const RealVecChart &c ) { return _d_b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dx_b7 ( const RealVecChart &c ) { return _d_b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dx_b8 ( const RealVecChart &c ) { return _d_b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }

  static RealType _dy_b1 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b2 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b3 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b4 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b5 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dy_b6 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dy_b7 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dy_b8 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }

  static RealType _dz_b1 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b2 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b3 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b4 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b5 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b2_1d ( c[2] ); }
  static RealType _dz_b6 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b2_1d ( c[2] ); }
  static RealType _dz_b7 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b2_1d ( c[2] ); }
  static RealType _dz_b8 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b2_1d ( c[2] ); }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const RealVecChart &c );
  BASIS_FUNC_TYPE _deriv_basis[3][8];
  BASIS_FUNC_TYPE _basis[8];

  const RealType _hx, _hy, _hz;

public:
  QuocBaseFunctionSet3D( const RealType hx, const RealType hy, const RealType hz ) : _hx ( hx ), _hy ( hy ), _hz (hz) {
    
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;
    _basis[6] = _b7;
    _basis[7] = _b8;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;
    _deriv_basis[0][4] = _dx_b5;
    _deriv_basis[0][5] = _dx_b6;
    _deriv_basis[0][6] = _dx_b7;
    _deriv_basis[0][7] = _dx_b8;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
    _deriv_basis[1][4] = _dy_b5;
    _deriv_basis[1][5] = _dy_b6;
    _deriv_basis[1][6] = _dy_b7;
    _deriv_basis[1][7] = _dy_b8;

    _deriv_basis[2][0] = _dz_b1;
    _deriv_basis[2][1] = _dz_b2;
    _deriv_basis[2][2] = _dz_b3;
    _deriv_basis[2][3] = _dz_b4;
    _deriv_basis[2][4] = _dz_b5;
    _deriv_basis[2][5] = _dz_b6;
    _deriv_basis[2][6] = _dz_b7;
    _deriv_basis[2][7] = _dz_b8;

//     this->initializeQuadCache( );
  }

  static const int numBaseFuncs = 8;
  
  RealType evaluate ( int BaseFuncNum, const RealVecChart &c ) const {
    return _basis[BaseFuncNum] ( c );
  }
  
  inline const RealType evaluate( int BaseFuncNum, int QuadPoint ) const {
    return evaluate( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );
  }

//   inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
//     return BaseFunctionSetBase<RealType, RealVecChart, RealVecChart, 8, QuadRuleType, QuocBaseFunctionSet2D<RealType, qc::QC_3D, QuadRuleType> >::evaluate ( BaseFuncNum, QuadPoint );
//   }
//   
  void evaluateGradient ( int BaseFuncNum, const RealVecChart &c, RealVecChart &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( c ) / _hx;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( c ) / _hy;
    Gradient[2] = _deriv_basis[2][BaseFuncNum] ( c ) / _hz;
  }

  inline RealVecChart evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return RealVecChart( _deriv_basis[0][BaseFuncNum]( this->_quadRule.getRefCoord( QuadPoint ) ) / _hx, 
                       _deriv_basis[1][BaseFuncNum]( this->_quadRule.getRefCoord( QuadPoint ) ) / _hy,
                       _deriv_basis[2][BaseFuncNum]( this->_quadRule.getRefCoord( QuadPoint ) ) / _hz
                     );
  }
  
//   inline const RealVecChart& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
//     return BaseFunctionSetBase<RealType, RealVecChart, RealVecChart, 8, QuadRuleType, QuocBaseFunctionSet2D<RealType, qc::QC_3D, QuadRuleType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
//   }

};



//! ==================================================================================================================
//  ==================================================================================================================
//! ==================================================================================================================
//                           Cached Intefaces
//  ==================================================================================================================
//! ==================================================================================================================
//  ==================================================================================================================


template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocCachedBaseFunctionSet1D 
: public FECachedBaseFunctionSetBase< DataTypeContainer, QuadType, 2, QuocCachedBaseFunctionSet1D<DataTypeContainer,QuadType,ElementType> >  {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;

  static RealType _b1   ( const RealType &c ) { return ( 1. - c ); }
  static RealType _b2   ( const RealType &c ) { return c; }
  
  static RealType _dx_b1   ( const RealType &/*c*/ ) { return - 1.; }
  static RealType _dx_b2   ( const RealType &/*c*/ ) { return 1.; }
  

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const RealType &RefCoord );
  BASIS_FUNC_TYPE _basis[2];
  BASIS_FUNC_TYPE _deriv_basis[2];
  const RealType _hx;
  
public:
  QuocCachedBaseFunctionSet1D( const RealType hx ) : _hx ( hx ){
    _basis[0] = _b1;
    _basis[1] = _b2;

    _deriv_basis[0] = _dx_b1;
    _deriv_basis[1] = _dx_b2;
    
    this->initializeQuadCache( );
  }
  
  static const int numBaseFuncs = 2;

  void setElement ( const ElementType &/*El*/ ) { }
  
  RealType evaluate ( int BaseFuncNum, const RealVecChart &c ) const { return _basis[BaseFuncNum] ( c[0] );}
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return FECachedBaseFunctionSetBase<DataTypeContainer, QuadType, 2, QuocCachedBaseFunctionSet1D<DataTypeContainer, QuadType, ElementType> >::evaluate ( BaseFuncNum, QuadPoint );
  }
  void evaluateGradient ( int BaseFuncNum, const RealVecChart &RefCoord, RealVecChart &Gradient ) const {
      Gradient[0] = _deriv_basis[BaseFuncNum] ( RefCoord[0] ) / _hx;
  }
  inline const RealVecChart & evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return FECachedBaseFunctionSetBase<DataTypeContainer, QuadType, 2, QuocCachedBaseFunctionSet1D<DataTypeContainer, QuadType, ElementType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }
  
};


//! The basefunctionset for bilinear elements in 2d
template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocCachedBaseFunctionSet2D 
: public FECachedBaseFunctionSetBase< DataTypeContainer, QuadType, 4, QuocCachedBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType> >{
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
    
  static RealType _b1   ( const RealVecChart &c ) { return ( 1.0 - c[0] ) * ( 1.0 - c[1] ); }
  static RealType _b2   ( const RealVecChart &c ) { return c[0] * ( 1.0 - c[1] ); }
  static RealType _b3   ( const RealVecChart &c ) { return ( 1.0 - c[0] ) * c[1]; }
  static RealType _b4   ( const RealVecChart &c ) { return c[0]*c[1]; }
  
  static RealType _dx_b1   ( const RealVecChart &c ) { return c[1] - 1.0; }
  static RealType _dx_b2   ( const RealVecChart &c ) { return 1.0 - c[1]; }
  static RealType _dx_b3   ( const RealVecChart &c ) { return -c[1]; }
  static RealType _dx_b4   ( const RealVecChart &c ) { return c[1]; }

  static RealType _dy_b1   ( const RealVecChart &c ) { return c[0] - 1.0; }
  static RealType _dy_b2   ( const RealVecChart &c ) { return -c[0]; }
  static RealType _dy_b3   ( const RealVecChart &c ) { return 1.0 - c[0]; }
  static RealType _dy_b4   ( const RealVecChart &c ) { return c[0]; }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const RealVecChart &c );
  BASIS_FUNC_TYPE _basis[4];
  BASIS_FUNC_TYPE _deriv_basis[2][4];
  const RealType  _hx, _hy;

public:
  QuocCachedBaseFunctionSet2D ( const RealType hx, const RealType hy ) : _hx ( hx ), _hy ( hy ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
    
    this->initializeQuadCache( );
  }

  static const int numBaseFuncs = 4;
  
  RealType evaluate ( int BaseFuncNum, const RealVecChart &c ) const { return _basis[BaseFuncNum] ( c );}
  
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return FECachedBaseFunctionSetBase<DataTypeContainer, QuadType, 4, QuocCachedBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType> >::evaluate ( BaseFuncNum, QuadPoint );
  }
  
  void evaluateGradient ( int BaseFuncNum, const RealVecChart &c, RealVecChart &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( c ) / _hx;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( c ) / _hy;
  }
  
  inline const RealVecChart& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return FECachedBaseFunctionSetBase<DataTypeContainer, QuadType, 4, QuocCachedBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

};


//! The basefunctionset for trilinear elements in 3d.
template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocCachedBaseFunctionSet3D 
  : public FECachedBaseFunctionSetBase< DataTypeContainer, QuadType, 8, QuocCachedBaseFunctionSet3D<DataTypeContainer, QuadType, ElementType> >{

  typedef typename DataTypeContainer::RealType RealType;  
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
    
  static inline RealType _b1_1d ( RealType x ) { return 1. - x; }
  static inline RealType _b2_1d ( RealType x ) { return x; }

  static inline RealType _d_b1_1d ( RealType ) { return -1.; }
  static inline RealType _d_b2_1d ( RealType ) { return 1.; }

  // b1 (0,0,0)
  // b2 (1,0,0)
  // b3 (0,1,0)
  // b4 (1,1,0)
  // b5 (0,0,1)
  // b6 (1,0,1)
  // b7 (0,1,1)
  // b8 (1,1,1)
  
  static RealType _b1   ( const RealVecChart &c ) { return _b1_1d ( c[0] ) * _b1_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b2   ( const RealVecChart &c ) { return _b2_1d ( c[0] ) * _b1_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b3   ( const RealVecChart &c ) { return _b1_1d ( c[0] ) * _b2_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b4   ( const RealVecChart &c ) { return _b2_1d ( c[0] ) * _b2_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b5   ( const RealVecChart &c ) { return _b1_1d ( c[0] ) * _b1_1d ( c[1] ) * _b2_1d ( c[2] ); }
  static RealType _b6   ( const RealVecChart &c ) { return _b2_1d ( c[0] ) * _b1_1d ( c[1] ) * _b2_1d ( c[2] ); }
  static RealType _b7   ( const RealVecChart &c ) { return _b1_1d ( c[0] ) * _b2_1d ( c[1] ) * _b2_1d ( c[2] ); }
  static RealType _b8   ( const RealVecChart &c ) { return _b2_1d ( c[0] ) * _b2_1d ( c[1] ) * _b2_1d ( c[2] ); }

  static RealType _dx_b1 ( const RealVecChart &c ) { return _d_b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b2 ( const RealVecChart &c ) { return _d_b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b3 ( const RealVecChart &c ) { return _d_b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b4 ( const RealVecChart &c ) { return _d_b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b5 ( const RealVecChart &c ) { return _d_b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dx_b6 ( const RealVecChart &c ) { return _d_b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dx_b7 ( const RealVecChart &c ) { return _d_b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dx_b8 ( const RealVecChart &c ) { return _d_b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }

  static RealType _dy_b1 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b2 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b3 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b4 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b5 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dy_b6 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dy_b7 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dy_b8 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }

  static RealType _dz_b1 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b2 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b3 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b4 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b5 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b2_1d ( c[2] ); }
  static RealType _dz_b6 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b2_1d ( c[2] ); }
  static RealType _dz_b7 ( const RealVecChart &c ) { return _b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b2_1d ( c[2] ); }
  static RealType _dz_b8 ( const RealVecChart &c ) { return _b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b2_1d ( c[2] ); }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const RealVecChart &c );
  BASIS_FUNC_TYPE _deriv_basis[3][8];
  BASIS_FUNC_TYPE _basis[8];

  const RealType _hx, _hy, _hz;

public:
  QuocCachedBaseFunctionSet3D( const RealType hx, const RealType hy, const RealType hz ) : _hx ( hx ), _hy ( hy ), _hz (hz) {
    
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;
    _basis[6] = _b7;
    _basis[7] = _b8;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;
    _deriv_basis[0][4] = _dx_b5;
    _deriv_basis[0][5] = _dx_b6;
    _deriv_basis[0][6] = _dx_b7;
    _deriv_basis[0][7] = _dx_b8;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
    _deriv_basis[1][4] = _dy_b5;
    _deriv_basis[1][5] = _dy_b6;
    _deriv_basis[1][6] = _dy_b7;
    _deriv_basis[1][7] = _dy_b8;

    _deriv_basis[2][0] = _dz_b1;
    _deriv_basis[2][1] = _dz_b2;
    _deriv_basis[2][2] = _dz_b3;
    _deriv_basis[2][3] = _dz_b4;
    _deriv_basis[2][4] = _dz_b5;
    _deriv_basis[2][5] = _dz_b6;
    _deriv_basis[2][6] = _dz_b7;
    _deriv_basis[2][7] = _dz_b8;

    this->initializeQuadCache( );
  }

  static const int numBaseFuncs = 8;
  
  RealType evaluate ( int BaseFuncNum, const RealVecChart &c ) const {
    return _basis[BaseFuncNum] ( c );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return FECachedBaseFunctionSetBase<DataTypeContainer, QuadType, 8, QuocCachedBaseFunctionSet3D<DataTypeContainer, QuadType, ElementType> >::evaluate ( BaseFuncNum, QuadPoint );
  }
  
  void evaluateGradient ( int BaseFuncNum, const RealVecChart &c, RealVecChart &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( c ) / _hx;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( c ) / _hy;
    Gradient[2] = _deriv_basis[2][BaseFuncNum] ( c ) / _hz;
  }
  
  inline const RealVecChart& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return FECachedBaseFunctionSetBase<DataTypeContainer, QuadType, 8, QuocCachedBaseFunctionSet3D<DataTypeContainer, QuadType, ElementType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

};



#endif //__QUOCBASEFUNCTIONSETS_H
