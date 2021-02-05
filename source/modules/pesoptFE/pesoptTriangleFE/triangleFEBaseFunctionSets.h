#ifndef __TRIANGLEFEBASEFUNCTIONSETS_H
#define __TRIANGLEFEBASEFUNCTIONSETS_H

#include <pesopt_IO.h>




template <typename DataTypeContainer, typename QuadType, typename TriangleType>
class RefTriangMeshBaseFunctionSetP0 
: public FEBaseFunctionSetBase< DataTypeContainer, QuadType >  {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;

  const TriangleType *_triangle;
  
public:
  RefTriangMeshBaseFunctionSetP0( ) : _triangle ( NULL ) { }

  void setTriangle ( const TriangleType &T ) { _triangle = &T;}
  
  static const int numBaseFuncs = 1;
  
  RealType evaluate ( int BaseFuncNum, const RealVecChart &RefCoord ) const { return 1.;}
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const { return 1.;}
  
};




template <typename DataTypeContainer, typename QuadType, typename TriangleType, typename Imp>
class RefTriangMeshBaseFunctionSetPkBase 
: public FEBaseFunctionSetBase< DataTypeContainer, QuadType >  {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;

  const TriangleType *_triangle;

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
public:
  RefTriangMeshBaseFunctionSetPkBase( ) : _triangle ( NULL ) { }

  void setTriangle ( const TriangleType &T ) { _triangle = &T;}
  
  //****************************************************
  // on unit triangle
  //****************************************************
   
  //! base function, has to be provided in derived classes.
  RealType evaluateOnUnitTriang ( int BaseFuncNum, const RealVecChart &RefCoord ) const { 
     throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
     return this->asImp().evaluateOnUnitTriang ( BaseFuncNum, RefCoord );
  }

  inline const RealType evaluateOnUnitTriang ( int BaseFuncNum, int QuadPoint ) const { 
      return this->asImp().evaluateOnUnitTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );
  }
  
  //! base function, has to be provided in derived classes.
  void evaluateGradientOnUnitTriang ( int BaseFuncNum, const RealVecChart &RefCoord, RealVecChart &Gradient ) const {
     throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
     this->asImp().evaluateGradientOnUnitTriang ( BaseFuncNum, RefCoord, Gradient );
  }
  
  inline RealVecChart evaluateGradientOnUnitTriang ( int BaseFuncNum, int QuadPoint ) const {
     RealVecChart Gradient;
     this->asImp().evaluateGradientOnUnitTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), Gradient );
     return Gradient;
  }
  //*****************************************************
  
  
  //****************************************************
  // on reference triangle
  //****************************************************
  RealType evaluate ( int BaseFuncNum, const RealVecChart &RefCoord ) const { return this->asImp().evaluateOnUnitTriang(BaseFuncNum,RefCoord);}
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const { return this->evaluate ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );}
  
  // Gradient = Dx g^{-1} Grad_UnitTriang = Dx (Dx^T Dx)^-1 Grad_UnitTriang
  void evaluateGradient ( int BaseFuncNum, const RealVecChart &RefCoord, RealVecChart &Gradient ) const {
    RealVecChart GradUnitTriang; this->asImp().evaluateGradientOnUnitTriang( BaseFuncNum, RefCoord, GradUnitTriang );
    RealVecChart ginvDx = _triangle->getInverseMetric( ) * GradUnitTriang;  
    const RealVecChart dir0 = _triangle->getEdge(0,1);
    const RealVecChart dir1 = _triangle->getEdge(0,2); 
    Gradient = ginvDx[0] * dir0 + ginvDx[1] * dir1;
  }

  inline RealVecChart evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    RealVecChart Gradient;
    this->asImp().evaluateGradient ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), Gradient );
    return Gradient;
  }
  //*****************************************************
  
};







template <typename DataTypeContainer, typename QuadType, typename TriangleType>
class RefTriangMeshBaseFunctionSetP1 : 
public RefTriangMeshBaseFunctionSetPkBase< DataTypeContainer, QuadType, TriangleType, RefTriangMeshBaseFunctionSetP1<DataTypeContainer, QuadType, TriangleType> >  {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
    
  static RealType _b1   ( const RealVecChart &c ) { return 1. - c[0] - c[1]; }
  static RealType _b2   ( const RealVecChart &c ) { return c[0]; }
  static RealType _b3   ( const RealVecChart &c ) { return c[1]; }

  static RealType _d1_b1   ( const RealVecChart &/*c*/ ) { return - 1.; }
  static RealType _d1_b2   ( const RealVecChart &/*c*/ ) { return   1.; }
  static RealType _d1_b3   ( const RealVecChart &/*c*/ ) { return   0.; }

  static RealType _d2_b1   ( const RealVecChart &/*c*/ ) { return - 1.; }
  static RealType _d2_b2   ( const RealVecChart &/*c*/ ) { return   0.; }
  static RealType _d2_b3   ( const RealVecChart &/*c*/ ) { return   1.; }
  
  
  typedef RealType ( *BASIS_FUNC_TYPE ) ( const RealVecChart &RefCoord );
  BASIS_FUNC_TYPE _basis[3];
  BASIS_FUNC_TYPE _dbasis[2][3];

public:
  RefTriangMeshBaseFunctionSetP1(  ) :
  RefTriangMeshBaseFunctionSetPkBase< DataTypeContainer, QuadType, TriangleType, RefTriangMeshBaseFunctionSetP1<DataTypeContainer, QuadType, TriangleType> > () {
      _basis[0] = _b1; 
      _basis[1] = _b2; 
      _basis[2] = _b3;
      
      _dbasis[0][0] = _d1_b1;
      _dbasis[0][1] = _d1_b2;
      _dbasis[0][2] = _d1_b3;
    
      _dbasis[1][0] = _d2_b1;
      _dbasis[1][1] = _d2_b2;
      _dbasis[1][2] = _d2_b3;
  }
  
  static const int numBaseFuncs = 3;
  
  RealType evaluateOnUnitTriang ( int BaseFuncNum, const RealVecChart &RefCoord ) const { return _basis[BaseFuncNum] ( RefCoord );}
  
  void evaluateGradientOnUnitTriang ( int BaseFuncNum, const RealVecChart &RefCoord, RealVecChart &Gradient ) const {
    Gradient[0] = _dbasis[0][BaseFuncNum](RefCoord); Gradient[1] = _dbasis[1][BaseFuncNum](RefCoord);
  }
};





#endif
