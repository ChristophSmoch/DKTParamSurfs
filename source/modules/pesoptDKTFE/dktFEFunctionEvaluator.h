#ifndef __DKTFEFUNCTIONEVALUATOR_H
#define __DKTFEFUNCTIONEVALUATOR_H

#include <pesopt_IO.h>

template <typename ConfiguratorType, typename VectorType, typename NLTYPE >
struct DKTDiscreteFunctionLookup {

  void getLocalDof ( const ConfiguratorType & conf, const VectorType &dofs, const typename ConfiguratorType::ElementType &El,
                     const int bfNum,  const NLTYPE &bfValue, NLTYPE &aux ) {
      aux = bfValue;
      aux *= dofs[ conf.localToGlobal ( El, bfNum ) ];
  }

};


template <typename ConfiguratorType, typename VectorType = typename ConfiguratorType::VectorType >
class DKTFEScalarFunctionEvaluator {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::Matrix22  Matrix22;
  typedef typename ConfiguratorType::BaseFuncSetType BaseFuncSetType;

  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

  const ConfiguratorType & _conf;
  const VectorType & _dofs;

  DKTFEScalarFunctionEvaluator ( const ConfiguratorType & config, const VectorType & Dofs )
    : _conf ( config ),
      _dofs ( Dofs ) { }

  RealType evaluate ( const ElementType &El, const DomVecType& RefCoord ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    RealType w = 0.;
    RealType v, aux;
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateOnRefTriang ( b, RefCoord );
      DKTDiscreteFunctionLookup<ConfiguratorType, VectorType, RealType>().getLocalDof( _conf, _dofs, El, b, v, aux );
      w += aux;
    }
    return w;
  }

  RealType evaluateAtQuadPoint ( const ElementType &El, int QuadPoint ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    RealType w = 0.;
    RealType v, aux;
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateOnRefTriang ( b, QuadPoint );
      DKTDiscreteFunctionLookup<ConfiguratorType, VectorType, RealType>().getLocalDof ( _conf, _dofs, El, b, v, aux );
      w += aux;
    }
    return w;
  }

  void evaluateGradient ( const ElementType &El, const DomVecType& RefCoord, DomVecType& Grad ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    DomVecType v, aux;
    Grad.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      bfs.evaluateGradientOnRefTriang ( b, RefCoord, v );
      DKTDiscreteFunctionLookup<ConfiguratorType, VectorType, DomVecType>().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Grad += aux;
    }
  }

  void evaluateGradientAtQuadPoint ( const ElementType &El, int QuadPoint, DomVecType& Grad ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    DomVecType v, aux;
    Grad.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateGradientOnRefTriang ( b, QuadPoint );
      DKTDiscreteFunctionLookup<ConfiguratorType, VectorType, DomVecType >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Grad += aux;
    }
  }

  //! evaluation of approximative gradient and its derivative

  void evaluateApproxGradient ( const ElementType &El, const DomVecType& RefCoord, DomVecType& Grad ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    DomVecType v, aux;
    Grad.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      bfs.evaluateApproxGradientOnRefTriang ( b, RefCoord, v );
      DKTDiscreteFunctionLookup<ConfiguratorType, VectorType, DomVecType>().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Grad += aux;
    }
  }

  void evaluateApproxGradientAtQuadPoint ( const ElementType &El, int QuadPoint, DomVecType& Grad ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    DomVecType v, aux;
    Grad.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateApproxGradientOnRefTriang ( b, QuadPoint );
      DKTDiscreteFunctionLookup<ConfiguratorType, VectorType, DomVecType >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Grad += aux;
    }
  }


  //returns HessianVec = ( d_x theta_1 , d_y theta_2, d_y theta_1 + d_x theta_2  )
  void evaluateApproxHessianAsVec ( const ElementType &El, const DomVecType& RefCoord, TangentVecType & Hessian ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    TangentVecType v, aux;
    Hessian.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      bfs.evaluateApproxHessianAsVecOnRefTriang ( b, RefCoord, v );
      DKTDiscreteFunctionLookup<ConfiguratorType, VectorType, TangentVecType >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Hessian += aux;
    }
  }

  void evaluateApproxHessianAsVecAtQuadPoint ( const ElementType &El, int QuadPoint, TangentVecType & Hessian ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    TangentVecType v, aux;
    Hessian.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateApproxHessianAsVecOnRefTriang ( b, QuadPoint );
      DKTDiscreteFunctionLookup<ConfiguratorType, VectorType, TangentVecType >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Hessian += aux;
    }
  }

  //returns Hessian = ( d_x theta_1 ,  d_y theta_1 )
  //                    d_x theta_2 ,  d_y theta_2 )
  void evaluateApproxHessian ( const ElementType &El, const DomVecType& RefCoord, Matrix22 & Hessian ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    Matrix22 v, aux;
    Hessian.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      bfs.evaluateApproxHessianOnRefTriang ( b, RefCoord, v );
      DKTDiscreteFunctionLookup<ConfiguratorType, VectorType, Matrix22 >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Hessian += aux;
    }
  }

  void evaluateApproxHessianAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix22 & Hessian ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    Matrix22 v, aux;
    Hessian.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateApproxHessianOnRefTriang ( b, QuadPoint );
      DKTDiscreteFunctionLookup<ConfiguratorType, VectorType, Matrix22 >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Hessian += aux;
    }
  }

  //returns Hessian = ( d_x theta_1,                       1/2 (d_y theta_1 + d_x theta_2)
  //                    1/2( d_y theta_1 + d_x theta_2) ,  d_y theta_2,  )
  void evaluateApproxHessianSym ( const ElementType &El, const DomVecType& RefCoord, Matrix22 & Hessian ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    Matrix22 v, aux;
    Hessian.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      bfs.evaluateApproxHessianSymOnRefTriang ( b, RefCoord, v );
      DKTDiscreteFunctionLookup<ConfiguratorType, VectorType, Matrix22 >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Hessian += aux;
    }
  }
  void evaluateApproxHessianSymAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix22 & Hessian ) const {
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = _conf.getApproxGradientBaseFunctionSet ( El );
    Matrix22 v, aux;
    Hessian.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateApproxHessianSymOnRefTriang ( b, QuadPoint );
      DKTDiscreteFunctionLookup<ConfiguratorType, VectorType, Matrix22 >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Hessian += aux;
    }
  }

  const VectorType& getDofs (  ) const {return _dofs;}

};




template <typename ConfiguratorType>
class DKTFEVectorFunctionEvaluator {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::Matrix22  Matrix22;
  typedef typename ConfiguratorType::Matrix32  Matrix32;
  typedef typename ConfiguratorType::Matrix33  Matrix33;

  typedef typename ConfiguratorType::Tensor222Type  Tensor222Type;
  typedef typename ConfiguratorType::Tensor322Type  Tensor322Type;

  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType  DomVecType;

  typedef DKTFEScalarFunctionEvaluator<ConfiguratorType, Eigen::Ref<const VectorType> > DiscFuncType;

  const int _numComponents;
  mutable std::vector<Eigen::Ref<const VectorType> > _refs;
  mutable std::vector<DiscFuncType> _discrFuncs;

public:

  DKTFEVectorFunctionEvaluator ( const ConfiguratorType &Configurator,  const VectorType &Dofs, const int numComponents = 3 ) : _numComponents ( numComponents )
  {
     const int numGlobalDofs = Configurator.getNumGlobalDofs();
     _refs.reserve( _numComponents );
     _discrFuncs.reserve( _numComponents );

     _refs.push_back ( Dofs.segment( 0, numGlobalDofs) );
     _refs.push_back ( Dofs.segment( numGlobalDofs, numGlobalDofs) );
     _refs.push_back ( Dofs.segment( 2 * numGlobalDofs, numGlobalDofs) );

     _discrFuncs.push_back ( DiscFuncType ( Configurator, _refs[0] ) );
     _discrFuncs.push_back ( DiscFuncType ( Configurator, _refs[1] ) );
     _discrFuncs.push_back ( DiscFuncType ( Configurator, _refs[2] ) );
  }

  void evaluate( const ElementType &El, const DomVecType& RefCoord, Point3DType &Value ) const {
    for ( int c = 0; c < 3; ++c ) Value[c] = _discrFuncs[c].evaluate ( El, RefCoord );
  }

  void evaluateAtQuadPoint ( const ElementType &El, int QuadPoint, Point3DType &Value ) const {
    for ( int c = 0; c < 3; ++c ) Value[c] = _discrFuncs[c].evaluateAtQuadPoint ( El, QuadPoint );
  }

  void evaluateGradient ( const ElementType &El, const DomVecType& RefCoord, Matrix32 &Dx ) const {
    DomVecType v;
    for ( int c = 0; c < 3; c++ ) {
      _discrFuncs[c].evaluateGradient ( El, RefCoord, v );
      for( int i=0; i<v.size(); ++i ) Dx( c, i ) = v[i];
    }
  }

  void evaluateGradientAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix32 &Dx ) const {
    DomVecType v;
    for ( int c = 0; c < 3; ++c ) {
      _discrFuncs[c].evaluateGradientAtQuadPoint ( El, QuadPoint, v );
      for( int i=0; i<v.size(); ++i ) Dx( c, i ) = v[i];
    }
  }

// //   void evaluateGradientInTangentSpaceAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix32 &Dx ) const {
// //     Matrix32 Dx33; this->evaluateGradientAtQuadPoint( El, QuadPoint, Dx33 );
// //     //TODO
// //     Dx.col(0) = Dx33.col(0);
// //     Dx.col(1) = Dx33.col(1);
// //   }


  //! evaluation of approximative gradient and its derivative

  void evaluateApproxGradient ( const ElementType &El, const DomVecType& RefCoord, Matrix32 &Dx ) const {
    DomVecType v;
    for ( int c = 0; c < 3; c++ ) {
      _discrFuncs[c].evaluateApproxGradient ( El, RefCoord, v );
      for( int i=0; i<v.size(); ++i ) Dx( c, i ) = v[i];
    }
  }

  void evaluateApproxGradientAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix32 &Dx ) const {
    DomVecType v;
    for ( int c = 0; c < 3; c++ ) {
      _discrFuncs[c].evaluateApproxGradientAtQuadPoint ( El, QuadPoint, v );
      for( int i=0; i<v.size(); ++i ) Dx( c, i ) = v[i];
    }
  }

  // ddX = ( HVec_comp0      i.e. component corresponds to column
  //         HVec_comp1
  //         HVec_comp2 )
 void evaluateApproxHessianAsVec ( const ElementType &El, const DomVecType& RefCoord, Matrix33 & ddX ) const {
    TangentVecType hessian_c;
    for ( int c = 0; c < 3; ++c ) {
      _discrFuncs[c].evaluateApproxHessianAsVec ( El , RefCoord, hessian_c );
      for( int j=0; j<3; ++j )
          ddX( c, j ) = hessian_c[j];
    }
  }

  void evaluateApproxHessianAsVecAtQuadPoint ( const ElementType &El, int QuadPoint, Matrix33 & ddX ) const {
    TangentVecType hessian_c;
    for ( int c = 0; c < 3; ++c ) {
      _discrFuncs[c].evaluateApproxHessianAsVecAtQuadPoint ( El , QuadPoint, hessian_c );
      for( int j=0; j<3; ++j )
          ddX( c, j ) = hessian_c[j];
    }
  }

  void evaluateApproxHessianAtQuadPoint ( const ElementType &El, int QuadPoint, Tensor322Type & ddX ) const {
    Matrix22 hessian_c;
    for ( int c = 0; c < 3; ++c ) {
      _discrFuncs[c].evaluateApproxHessianAtQuadPoint ( El , QuadPoint, hessian_c );
      for( int j=0; j<2; ++j )
        for( int k=0; k<2; ++k )
          ddX.set ( c, j, k, hessian_c( j, k ) );
    }
  }


  void evaluateApproxHessianSym ( const ElementType &El, const DomVecType& RefCoord, Tensor322Type & ddX ) const {
    Matrix22 hessian_c;
    for ( int c = 0; c < 3; ++c ) {
      _discrFuncs[c].evaluateApproxHessianSym ( El , RefCoord, hessian_c );
      for( int j=0; j<2; ++j )
        for( int k=0; k<2; ++k )
          ddX.set ( c, j, k, hessian_c( j, k ) );
    }
  }
  void evaluateApproxHessian ( const ElementType &El, const DomVecType& RefCoord, Tensor322Type & ddX ) const {
    Matrix22 hessian_c;
    for ( int c = 0; c < 3; ++c ) {
      _discrFuncs[c].evaluateApproxHessian ( El , RefCoord, hessian_c );
      for( int j=0; j<2; ++j )
        for( int k=0; k<2; ++k )
          ddX.set ( c, j, k, hessian_c( j, k ) );
    }
  }
  void evaluateApproxHessianSymAtQuadPoint ( const ElementType &El, int QuadPoint, Tensor322Type & ddX ) const {
    Matrix22 hessian_c;
    for ( int c = 0; c < 3; ++c ) {
      _discrFuncs[c].evaluateApproxHessianSymAtQuadPoint ( El , QuadPoint, hessian_c );
      for( int j=0; j<2; ++j )
        for( int k=0; k<2; ++k )
          ddX.set ( c, j, k, hessian_c( j, k ) );
    }
  }


  const DiscFuncType& operator[] ( int i ) const { return _discrFuncs[i];}
  DiscFuncType& operator[] ( int i ) { return _discrFuncs[i];}

};






// template <typename ConfiguratorType, int dimDomain = ConfiguratorType::ElementType::DimDomainRefTriang>
// class PointWiseVectorFunctionEvaluatorShellFE {
//
// };




template <typename ConfiguratorType>
// class PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType,2> {
class PointWiseVectorFunctionEvaluatorShellFE {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::Matrix22  Matrix22;
  typedef typename ConfiguratorType::Matrix32  Matrix32;
  typedef typename ConfiguratorType::Matrix33  Matrix33;

  typedef typename ConfiguratorType::Tensor222Type  Tensor222Type;
  typedef typename ConfiguratorType::Tensor322Type  Tensor322Type;

  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType  DomVecType;


public:

  PointWiseVectorFunctionEvaluatorShellFE ( ) { }


  void evaluateFirstFundamentalForm ( const Matrix32 &Dx, Matrix22 &g ) const {
    g = Dx.transpose() * Dx;
  }

  RealType evaluateArea ( const Matrix32 &Dx ) const {
    Matrix22 g; g = Dx.transpose() * Dx;
    return std::sqrt( g.determinant() );
  }

  //TODO used for Dirichlet Energy : int g^1 Du \cdot Du
//   void evaluateFirstFundamentalFormInverseAsMatrix22 ( const Matrix22 &ginv, Matrix22 &ginvAsGradDeformType ) const {
//     ginvAsGradDeformType = ginv;
//   }



  // given chart X -> compute n \circ X  = d_1 X \cross d_2 X / |.|
  void evaluateNormal ( const Matrix32 &Dx, TangentVecType &normal ) const {
    TangentVecType unNormalizedNormal ( (Dx.col(0)).cross(Dx.col(1)) );
    normal = unNormalizedNormal.normalized();
  }

  void evaluateSecondFundamentalForm ( const Tensor322Type &ddX, const TangentVecType &normal, Matrix22& shapeTensor ) const {
      for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j ) {
          TangentVecType ddX_ij;
          ddX.getVector( ddX_ij, i, j );
          shapeTensor( i , j ) = normal.dot( ddX_ij );
        }
  }

  //! computes derivative of g: Deriv_ijk = \partial_k g_ij = \sum_comp  \partial_ik x^comp \partial_j x^comp + \partial_i x^comp \partial_jk x^comp
  // TODO choise of Hessian? approximation of first derivative?
  void evaluateApproxDerivativeOfFirstFundamentalForm ( const Matrix32 &dX, const Tensor322Type &ddX, Tensor222Type& Deriv ) const {
      for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j )
          for( int k=0; k<2; ++k ){
            RealType tmp = 0.0;
            for( int comp=0; comp<3; ++comp ){
              tmp += ddX( comp, i, k ) * dX( comp, j ) + dX( comp, i ) * ddX( comp, j, k );
            }
            Deriv.set ( i , j , k, tmp );
        }
  }

  //! Gamma_ijk = 1/2 ( \partial_j g_ki + \partial_i g_kj - \partial_k g_ij )
//   void evaluateApproxChristoffelSymbolsOfFirstKind ( const Matrix32 &dX, const Tensor322Type &ddX, Tensor222Type& ChristoffelSym ) const {
//       Tensor222Type Dg;
//       evaluateApproxDerivativeOfFirstFundamentalForm ( dX, ddX, Dg );
//       for ( int i = 0; i < 2; ++i )
//         for ( int j = 0; j < 2; ++j )
//           for( int k=0; k<2; ++k ){
//             ChristoffelSym.set ( i , j , k, 0.5 * ( Dg( k, i, j) + Dg(k, j, i) - Dg(i,j,k) )   );
//         }
//   }
  //TODO Gamma_ijk = partial_ij x cdot partial_k x
  void evaluateApproxChristoffelSymbolsOfFirstKind ( const Matrix32 &dX, const Tensor322Type &ddX, Tensor222Type& ChristoffelSym ) const {
      for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j )
          for( int k=0; k<2; ++k ){
              RealType aux = 0.;
              for( int a=0; a<3; ++a ){
                  aux += ddX( a, i, j ) * dX( a, k );
              }
              ChristoffelSym.set ( i , j , k, aux );
        }
  }

  //! Gamma_ij^m = \sum_k g^m,k Gamma_i,j,k = \sum_k 1/2 g^m,k ( \partial_j g_ki + \partial_i g_kj - \partial_k g_ij )
  void evaluateApproxChristoffelSymbolsOfSecondKind ( const Matrix22 &ginv, const Matrix32 &dX, const Tensor322Type &ddX, Tensor222Type &ChristoffelSym2 ) const {
    Tensor222Type ChristoffelSym1;
    evaluateApproxChristoffelSymbolsOfFirstKind ( dX, ddX, ChristoffelSym1 );
    for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j )
          for( int m=0; m<2; ++m ){
            RealType tmp = 0.0;
            for( int k=0; k<2; ++k ){
              tmp += 0.5 * ginv( m, k ) * ChristoffelSym1( i, j, k );
            }
            ChristoffelSym2.set ( i , j , m, tmp  );
        }
  }

  //! g^-1 Gamma
  void evaluateVectorForLaplacian ( const Matrix22 &ginv, const Matrix32 &dX, const Tensor322Type &ddX, DomVecType& vec ) const {
    vec.setZero();
    Tensor222Type ChristoffelSym2;
    evaluateApproxChristoffelSymbolsOfSecondKind ( ginv, dX, ddX, ChristoffelSym2 );
    for( int l=0; l<2; ++l )
      for ( int j=0; j<2; ++j )
        for( int k=0; k<2; ++k )
          vec[l] += ginv( j, k ) * ChristoffelSym2( j, k, l );
  }

  //! Laplace-Beltrami
  RealType evaluateLaplaceBeltrami ( const Matrix22 &gAinv, const DomVecType &dX, const Matrix22 &ddX, const DomVecType& vecForLaplaceA ) const {
      return pesopt::ddProd<RealType,Matrix22>( gAinv, ddX ) - vecForLaplaceA.dot( dX );
  }

  void evaluateLaplaceBeltrami ( const Matrix22 &gAinv, const Matrix32 &dX, const Tensor322Type &ddX, const DomVecType& vecForLaplaceA, Point3DType &laplace ) const {
      for( int comp = 0; comp <3; ++comp ) laplace[comp] = pesopt::ddProd<RealType,Matrix22>( gAinv, ddX[comp] ) - vecForLaplaceA.dot( dX.row(comp) );
  }

};



template <typename ConfiguratorType>
class SemiNonlinIsometryPointWiseVectorFunctionEvaluatorShellFE {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::Matrix22  Matrix22;
  typedef typename ConfiguratorType::Matrix32  Matrix32;
  typedef typename ConfiguratorType::Matrix33  Matrix33;

  typedef typename ConfiguratorType::Tensor222Type  Tensor222Type;
  typedef typename ConfiguratorType::Tensor322Type  Tensor322Type;

  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType  DomVecType;


public:

  SemiNonlinIsometryPointWiseVectorFunctionEvaluatorShellFE ( ) { }


  void evaluateA1 ( const Matrix32 &Dx, TangentVecType &a1 ) const {
    const TangentVecType a1unnormalized = Dx.col(0);
    a1 = a1unnormalized.normalized();
  }

  void evaluateA2Hat ( const Matrix32 &Dx, const TangentVecType &a1, TangentVecType &a2hat ) const {
    const TangentVecType partial2x = Dx.col(1);
    a2hat = partial2x - (a1.dot(partial2x)) * a1;
  }

  void evaluateA2 ( const TangentVecType &a2hat, TangentVecType &a2 ) const {
    a2 = a2hat.normalized();
  }

  void evaluateA0Tilde ( const Tensor322Type &ddX, const TangentVecType &normal, Matrix22& a0tilde ) const {
      for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j ) {
          TangentVecType ddX_ij;
          ddX.getVector( ddX_ij, i, j );
          a0tilde( i , j ) = normal.dot( ddX_ij );
        }
  }

  void evaluateA1Tilde ( const Tensor322Type &ddX, const Matrix32 &Dx, const TangentVecType &a1, const TangentVecType &a2, const TangentVecType &a2hat, Matrix22& a1tilde ) const {
      const TangentVecType partial1x = Dx.col(0);
      const TangentVecType partial2x = Dx.col(1);
      for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j ) {
          TangentVecType ddX_ij;
          ddX.getVector( ddX_ij, i, j );
          a1tilde( i , j ) = (a1.dot( ddX_ij ))/partial1x.norm() - (a2.dot( ddX_ij )/a2hat.norm())*(partial2x.dot( a1 )/partial1x.norm());
        }
  }

  void evaluateA2Tilde ( const Tensor322Type &ddX, const TangentVecType &a2, const TangentVecType &a2hat, Matrix22& a2tilde ) const {
      for ( int i = 0; i < 2; ++i )
        for ( int j = 0; j < 2; ++j ) {
          TangentVecType ddX_ij;
          ddX.getVector( ddX_ij, i, j );
          a2tilde( i , j ) = (a2.dot( ddX_ij ))/a2hat.norm();
        }
  }

};


// template <typename ConfiguratorType>
// class PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType,3> {
// public:
//   typedef typename ConfiguratorType::RealType RealType;
//   typedef typename ConfiguratorType::TangentVecType TangentVecType;
//   typedef typename ConfiguratorType::Point3DType Point3DType;
//   typedef typename ConfiguratorType::VectorType VectorType;
//   typedef typename ConfiguratorType::Matrix22  Matrix22;
//   typedef typename ConfiguratorType::Matrix32  Matrix32;
//   typedef typename ConfiguratorType::Matrix33  Matrix33;
//
//   typedef typename ConfiguratorType::Tensor222Type  Tensor222Type;
//   typedef typename ConfiguratorType::Tensor322Type  Tensor322Type;
//
//   typedef typename ConfiguratorType::ElementType ElementType;
//   typedef typename ElementType::DomVecType  DomVecType;
//   typedef typename ElementType::DomVecType DomVecType;
//   typedef typename ElementType::Matrix32 Matrix32;
//   typedef typename ElementType::Matrix22 Matrix22;
//
//
// public:
//
//   PointWiseVectorFunctionEvaluatorShellFE ( ) { }
//
//
// // TODO for Dx in R^3x3???????????
//   void evaluateFirstFundamentalForm ( const Matrix33 &Dx, Matrix22 &g ) const {
// //     g = Dx.transpose() * Dx;
//        g.setZero(); g(0,0) = 1.; g(1,1) = 1.;
//   }
//
//  //TODO used for Dirichlet Energy : int g^1 Du \cdot Du
//   void evaluateFirstFundamentalFormInverseAsMatrix22 ( const Matrix22 &/*ginv*/, Matrix22 &ginvAsGradDeformType ) const {
//     ginvAsGradDeformType.setZero(); ginvAsGradDeformType(0,0) = 1.; ginvAsGradDeformType(1,1) = 1.; ginvAsGradDeformType(2,2) = 1.;
//   }
//
//
//   // given chart X -> compute n \circ X  = d_1 X \cross d_2 X / |.|
// TODO
//   void evaluateNormal ( const Matrix33 &Dx, TangentVecType &normal ) const {
//     TangentVecType unNormalizedNormal ( (Dx.col(0)).cross(Dx.col(1)) );
//     normal = unNormalizedNormal.normalized();
//       normal = TangentVecType( 0., 0., 1. );
//   }
//
//
//   void evaluateSecondFundamentalForm ( const Tensor322Type &ddX, const TangentVecType &normal, Matrix22& shapeTensor ) const {
//       for ( int i = 0; i < 2; ++i )
//         for ( int j = 0; j < 2; ++j ) {
//           TangentVecType ddX_ij;
//           ddX.getVector( ddX_ij, i, j );
//           shapeTensor( i , j ) = normal.dot( ddX_ij );
//         }
//   }
/*

  //! computes derivative of g: Deriv_ijk = \partial_k g_ij = \sum_comp  \partial_ik x^comp \partial_j x^comp + \partial_i x^comp \partial_jk x^comp
  // TODO choise of Hessian? approximation of first derivative?
  void evaluateApproxDerivativeOfFirstFundamentalForm ( const Matrix33 &dX, const Tensor322Type &ddX, Tensor222Type& Deriv ) const {
//       for ( int i = 0; i < 2; ++i )
//         for ( int j = 0; j < 2; ++j )
//           for( int k=0; k<2; ++k ){
//             RealType tmp = 0.0;
//             for( int comp=0; comp<3; ++comp ){
//               tmp += ddX( comp, i, k ) * dX( comp, j ) + dX( comp, i ) * ddX( comp, j, k );
//             }
//             Deriv.set ( i , j , k, tmp );
//         }
//       cout << "TODO" << endl;
  }*/
 /*
  //! Gamma_ijk = 1/2 ( \partial_j g_ki + \partial_i g_kj - \partial_k g_ij )
  void evaluateApproxChristoffelSymbolsOfFirstKind ( const Matrix33 &dX, const Tensor322Type &ddX, Tensor222Type& ChristoffelSym ) const {
//       Tensor222Type Dg;
//       evaluateApproxDerivativeOfFirstFundamentalForm ( dX, ddX, Dg );
//       for ( int i = 0; i < 2; ++i )
//         for ( int j = 0; j < 2; ++j )
//           for( int k=0; k<2; ++k ){
//             ChristoffelSym.set ( i , j , k, 0.5 * ( Dg( k, i, j) + Dg(k, j, i) - Dg(i,j,k) )   );
//         }
//             cout << "TODO" << endl;
  }

  //! Gamma_ij^m = \sum_k g^m,k Gamma_i,j,k = \sum_k 1/2 g^m,k ( \partial_j g_ki + \partial_i g_kj - \partial_k g_ij )
  void evaluateApproxChristoffelSymbolsOfSecondKind ( const Matrix22 &ginv, const Matrix33 &dX, const Tensor322Type &ddX, Tensor222Type &ChristoffelSym2 ) const {
//     Tensor222Type ChristoffelSym1;
//     evaluateApproxChristoffelSymbolsOfFirstKind ( dX, ddX, ChristoffelSym1 );
//     for ( int i = 0; i < 2; ++i )
//         for ( int j = 0; j < 2; ++j )
//           for( int m=0; m<2; ++m ){
//             RealType tmp = 0.0;
//             for( int k=0; k<2; ++k ){
//               tmp += 0.5 * ginv( m, k ) * ChristoffelSym1( i, j, k );
//             }
//             ChristoffelSym2.set ( i , j , m, tmp  );
//         }
//       cout << "TODO" << endl;
  }

  //! g^-1 Gamma
  void evaluateVectorForLaplacian ( const Matrix22 &ginv, const Matrix33 &dX, const Tensor322Type &ddX, DomVecType& vec ) const {
//     vec.setZero();
//     Tensor222Type ChristoffelSym2;
//     evaluateApproxChristoffelSymbolsOfSecondKind ( ginv, dX, ddX, ChristoffelSym2 );
//     for( int l=0; l<2; ++l )
//       for ( int j=0; j<2; ++j )
//         for( int k=0; k<2; ++k )
//           vec[l] += ginv( j, k ) * ChristoffelSym2( j, k, l );
//             cout << "TODO" << endl;
  }

  //! Laplace-Beltrami
  void evaluateLaplaceBeltrami ( const Matrix22 &gAinv, const Matrix33 &dX, const Tensor322Type &ddX, const DomVecType& vecForLaplaceA, Point3DType &laplace ) const {
//       for( int comp = 0; comp <3; ++comp ) laplace[comp] = pesopt::ddProd<RealType,Matrix22>( ginv, ddX[comp] ) - vecForLaplace.dot( dX.row(comp) );
//             cout << "TODO" << endl;
  }

};


*/










//! Store all informations at QuadPoints
// HelperClass
template <typename NLType>
class DFDStorageStructure
{
  private:
      const int _numberOfElements;
      const int _numberOfQuadPointsPerElement;
      NLType **const _dataStorage;

  public:

      DFDStorageStructure ( const int &numberOfElements, const int &numberOfQuadPointsPerElement) :
      _numberOfElements {numberOfElements},
      _numberOfQuadPointsPerElement {numberOfQuadPointsPerElement},
      _dataStorage {new NLType* [numberOfElements]}{
          for (int elementIndex = 0; elementIndex < _numberOfElements; ++elementIndex)
              _dataStorage [elementIndex] = new NLType [_numberOfQuadPointsPerElement];
      }

    ~DFDStorageStructure (){
        for (int elementIndex = 0; elementIndex < _numberOfElements; ++elementIndex) delete [] _dataStorage [elementIndex];
        delete [] _dataStorage;
    }

    const NLType& at (const int &elementIndex, const int &localQuadPointIndex) const{ return _dataStorage [elementIndex][localQuadPointIndex];}
    NLType& at (const int &elementIndex, const int &localQuadPointIndex){ return _dataStorage [elementIndex][localQuadPointIndex];}

    const NLType* operator [] (const int &elementIndex) const{return _dataStorage [elementIndex];}
    NLType* operator [] (const int &elementIndex) { return _dataStorage [elementIndex];}

    void set ( const int &elementIndex, const int &localQuadPointIndex, const NLType value) { this->at (elementIndex, localQuadPointIndex) = value;}
};



// template <typename NLType>
// class DFDStorageStructureEdge
// {
//   private:
//       const int _numberOfEdges;
//       NLType *const _dataStorage;
//
//   public:
//
//       DFDStorageStructureEdge ( const int &numberOfEdges ) :
//       _numberOfEdges {numberOfEdges},
//       _dataStorage {new NLType* [numberOfEdges]}{
//       }
//
//     ~DFDStorageStructure (){
//         delete [] _dataStorage;
//     }
//
//     const NLType& at (const int &elementIndex, const int &localQuadPointIndex) const{ return _dataStorage [elementIndex][localQuadPointIndex];}
//     NLType& at (const int &elementIndex, const int &localQuadPointIndex){ return _dataStorage [elementIndex][localQuadPointIndex];}
//
//     const NLType* operator [] (const int &elementIndex) const{return _dataStorage [elementIndex];}
//     NLType* operator [] (const int &elementIndex) { return _dataStorage [elementIndex];}
//
//     void set ( const int &elementIndex, const int &localQuadPointIndex, const NLType value) { this->at (elementIndex, localQuadPointIndex) = value;}
// };




enum DiscreteFunctionCacheType {
  NoCache = 0,
  FirstOrder = 1,
  FirstAndSecondOrder = 2,
};





template <typename ConfiguratorType, DiscreteFunctionCacheType>
class DiscreteScalarFunctionStorage {};



template <typename ConfiguratorType>
class DiscreteScalarFunctionStorage<ConfiguratorType,NoCache>
{

protected:

    typedef typename ConfiguratorType::RealType                 RealType;
    typedef typename ConfiguratorType::ElementType ElementType;
    typedef typename ConfiguratorType::DomVecType  DomVecType;
    typedef typename ConfiguratorType::VectorType               VectorType;

    const ConfiguratorType &_conf;
    const DKTFEScalarFunctionEvaluator<ConfiguratorType> _discrFunc;
    const int _numberOfElements;
    const int _numberOfQuadPointsPerElement;

public :
    const VectorType& _dofs;

  public:
    DiscreteScalarFunctionStorage ( const ConfiguratorType &conf, const VectorType &dofs ) :
      _conf (conf),
      _discrFunc (conf, dofs ),
      _numberOfElements (conf.getInitializer ().getNumTriangs ()),
      _numberOfQuadPointsPerElement (conf.maxNumQuadPoints ()),
      _dofs( dofs ) {}
};



template <typename ConfiguratorType>
class DiscreteScalarFunctionStorage<ConfiguratorType,FirstOrder>
{
  protected:

    typedef typename ConfiguratorType::RealType                 RealType;
    typedef typename ConfiguratorType::ElementType              ElementType;
    typedef typename ConfiguratorType::DomVecType               DomVecType;
    typedef typename ConfiguratorType::VectorType               VectorType;

    const ConfiguratorType &_conf;
    const DKTFEScalarFunctionEvaluator<ConfiguratorType> _discrFunc;
    const int _numberOfElements;
    const int _numberOfQuadPointsPerElement;

protected :
    const VectorType& _dofs;
    //
    DFDStorageStructure <RealType>                        _function;
    DFDStorageStructure <DomVecType>                      _Gradient;

  public:
    DiscreteScalarFunctionStorage ( const ConfiguratorType &conf, const VectorType &dofs ) :
      _conf (conf),
      _discrFunc (conf, dofs ),
      _numberOfElements (conf.getInitializer ().getNumTriangs ()),
      _numberOfQuadPointsPerElement (conf.maxNumQuadPoints ()),
      _dofs( dofs ),
      //
      _function ( _numberOfElements, _numberOfQuadPointsPerElement ),
      _Gradient (_numberOfElements, _numberOfQuadPointsPerElement)
    {
        for ( int elementIdx = 0; elementIdx < _conf.getInitializer().getNumTriangs(); ++elementIdx){
            const ElementType& El ( _conf.getInitializer().getTriang( elementIdx ) );
            for ( int localQuadPointIndex = 0; localQuadPointIndex < _numberOfQuadPointsPerElement; ++localQuadPointIndex){
                //
                _function.at(elementIdx, localQuadPointIndex ) = _discrFunc.evaluateAtQuadPoint ( El, localQuadPointIndex  );
                _discrFunc.evaluateGradientAtQuadPoint (El, localQuadPointIndex, _Gradient.at (elementIdx, localQuadPointIndex));
            }
        }
    }

public :
    const VectorType& getDofs( ) const {return _dofs;};
    const RealType& getFunction( const int elementIdx, const int QuadPoint ) const { return _function[elementIdx][QuadPoint]; }
    const DomVecType& getGradient( const int elementIdx, const int QuadPoint ) const { return _Gradient[elementIdx][QuadPoint]; }

};



template <typename ConfiguratorType>
class DiscreteScalarFunctionStorage<ConfiguratorType,FirstAndSecondOrder> : public DiscreteScalarFunctionStorage<ConfiguratorType,FirstOrder>
{
 protected:

    typedef typename ConfiguratorType::RealType                 RealType;
    typedef typename ConfiguratorType::Matrix22                 Matrix22;
    typedef typename ConfiguratorType::Matrix32                 Matrix32;
    typedef typename ConfiguratorType::Matrix33                 Matrix33;
    typedef typename ConfiguratorType::TangentVecType           TangentVecType;
    typedef typename ConfiguratorType::Point3DType              Point3DType;
    typedef typename ConfiguratorType::VectorType               VectorType;
    typedef typename ConfiguratorType::ElementType ElementType;
    typedef typename ConfiguratorType::DomVecType  DomVecType;

    const ConfiguratorType &_conf;
    const DKTFEScalarFunctionEvaluator<ConfiguratorType> _discrFunc;
    const int _numberOfElements;
    const int _numberOfQuadPointsPerElement;

protected:

    //ApproxGradient
    DFDStorageStructure<DomVecType>            _ApproxGradient;
    DFDStorageStructure<TangentVecType>        _HessianAsVec;
    DFDStorageStructure<Matrix22>              _Hessian;
    //Mixed
//     DFDStorageStructure < TangentVecType >     _LaplaceBeltrami;


  public:
    DiscreteScalarFunctionStorage ( const ConfiguratorType &conf, const VectorType &dofs ) :
    DiscreteScalarFunctionStorage<ConfiguratorType,FirstOrder> ( conf, dofs ),
      _conf (conf),
      _discrFunc (conf, dofs ),
      _numberOfElements (conf.getInitializer ().getNumTriangs ()),
      _numberOfQuadPointsPerElement (conf.maxNumQuadPoints ()),

      //approx
      _ApproxGradient (_numberOfElements, _numberOfQuadPointsPerElement),
      _HessianAsVec ( _numberOfElements, _numberOfQuadPointsPerElement ),
      _Hessian ( _numberOfElements, _numberOfQuadPointsPerElement )
      //mixed
//       _LaplaceBeltrami ( _numberOfElements, _numberOfQuadPointsPerElement )
    {
        for ( int elementIdx = 0; elementIdx < _conf.getInitializer().getNumTriangs(); ++elementIdx){
            const ElementType& El ( _conf.getInitializer().getTriang( elementIdx ) );
            for ( int localQuadPointIndex = 0; localQuadPointIndex < _numberOfQuadPointsPerElement; ++localQuadPointIndex){
                //approx
                _discrFunc.evaluateApproxGradientAtQuadPoint (El, localQuadPointIndex, _ApproxGradient.at (elementIdx, localQuadPointIndex));
                _discrFunc.evaluateApproxHessianAsVecAtQuadPoint ( El, localQuadPointIndex, _HessianAsVec.at( elementIdx, localQuadPointIndex ) );
                _discrFunc.evaluateApproxHessianSymAtQuadPoint (El, localQuadPointIndex, _Hessian.at(elementIdx, localQuadPointIndex));
                //mixed
//                 _discrFuncs.evaluateLaplaceBeltrami ( _firstFFInv.at(elementIdx, localQuadPointIndex), _Gradient.at(elementIdx, localQuadPointIndex), _Hessian.at(elementIdx, localQuadPointIndex), _gInvChristoffel2.at(elementIdx, localQuadPointIndex), _LaplaceBeltrami.at(elementIdx, localQuadPointIndex ) );
            }
        }
    }

public:
    const DomVecType& getApproxGradient ( const int elementIdx, const int QuadPoint ) const {return _ApproxGradient[elementIdx][QuadPoint];}
    const TangentVecType& getHessianAsVec ( const int elementIdx, const int QuadPoint ) const {return _HessianAsVec[elementIdx][QuadPoint];}
    const Matrix22& getHessian ( const int elementIdx, const int QuadPoint ) const {return _Hessian[elementIdx][QuadPoint];}
};













template <typename ConfiguratorType, DiscreteFunctionCacheType>
class DiscreteVectorFunctionStorage {};



template <typename ConfiguratorType>
class DiscreteVectorFunctionStorage<ConfiguratorType,NoCache>{

protected:

    typedef typename ConfiguratorType::RealType                 RealType;
    typedef typename ConfiguratorType::Matrix22                 Matrix22;
    typedef typename ConfiguratorType::Matrix32                 Matrix32;
    typedef typename ConfiguratorType::Matrix33                 Matrix33;
    typedef typename ConfiguratorType::TangentVecType           TangentVecType;
    typedef typename ConfiguratorType::Point3DType              Point3DType;
    typedef typename ConfiguratorType::VectorType               VectorType;
    typedef typename ConfiguratorType::Tensor322Type            Tensor322Type;

    typedef typename ConfiguratorType::ElementType ElementType;
    typedef typename ConfiguratorType::DomVecType  DomVecType;

    const ConfiguratorType &_conf;
    const DKTFEVectorFunctionEvaluator <ConfiguratorType> _discrFuncs;
    PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _pointwiseEvaluator;
    const int _numberOfElements;
    const int _numberOfQuadPointsPerElement;

public :
    const VectorType& _dofs;

public:
    DiscreteVectorFunctionStorage ( const ConfiguratorType &conf, const VectorType &dofs, const int numComponents ) :
      _conf (conf),
      _discrFuncs (conf, dofs, numComponents),
      _numberOfElements (conf.getInitializer ().getNumTriangs ()),
      _numberOfQuadPointsPerElement (conf.maxNumQuadPoints ()),
      _dofs( dofs ) {}
};



template <typename ConfiguratorType>
class DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder>
{
  protected:

    typedef typename ConfiguratorType::RealType                 RealType;
    typedef typename ConfiguratorType::Matrix22                 Matrix22;
    typedef typename ConfiguratorType::Matrix32                 Matrix32;
    typedef typename ConfiguratorType::Matrix33                 Matrix33;
    typedef typename ConfiguratorType::TangentVecType           TangentVecType;
    typedef typename ConfiguratorType::Point3DType              Point3DType;
    typedef typename ConfiguratorType::VectorType               VectorType;
    typedef typename ConfiguratorType::Tensor322Type            Tensor322Type;

    typedef typename ConfiguratorType::ElementType ElementType;
    typedef typename ConfiguratorType::DomVecType  DomVecType;

    const ConfiguratorType &_conf;
    const DKTFEVectorFunctionEvaluator <ConfiguratorType> _discrFuncs;
    PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _pointwiseEvaluator;
    const int _numberOfElements;
    const int _numberOfQuadPointsPerElement;

protected :
    const VectorType& _dofs;
    //
    DFDStorageStructure <Point3DType>                        _coords;
    DFDStorageStructure <Matrix32>                           _Gradient;
    DFDStorageStructure <Matrix22>                           _firstFF;
    DFDStorageStructure <Matrix22>                           _firstFFInv;
    DFDStorageStructure <RealType>                           _detFirstFF;
    DFDStorageStructure <RealType>                           _area;
    DFDStorageStructure <TangentVecType>                     _normal;
    DFDStorageStructure <Matrix32>                           _GradientGInv;

  public:
    DiscreteVectorFunctionStorage ( const ConfiguratorType &conf, const VectorType &dofs, const int numComponents ) :
      _conf (conf),
      _discrFuncs (conf, dofs, numComponents),
      _numberOfElements (conf.getInitializer ().getNumTriangs ()),
      _numberOfQuadPointsPerElement (conf.maxNumQuadPoints ()),
      _dofs( dofs ),
      //
      _coords ( _numberOfElements, _numberOfQuadPointsPerElement ),
      _Gradient (_numberOfElements, _numberOfQuadPointsPerElement),
      _firstFF (_numberOfElements, _numberOfQuadPointsPerElement),
      _firstFFInv (_numberOfElements, _numberOfQuadPointsPerElement),
      _detFirstFF (_numberOfElements, _numberOfQuadPointsPerElement),
      _area ( _numberOfElements, _numberOfQuadPointsPerElement ),
      _normal (_numberOfElements, _numberOfQuadPointsPerElement),
      _GradientGInv( _numberOfElements, _numberOfQuadPointsPerElement )
    {
        for ( int elementIdx = 0; elementIdx < _conf.getInitializer().getNumTriangs(); ++elementIdx){
            const ElementType& El ( _conf.getInitializer().getTriang( elementIdx ) );
            for ( int q = 0; q < _numberOfQuadPointsPerElement; ++q){
                //
                _discrFuncs.evaluateAtQuadPoint ( El, q, _coords.at(elementIdx, q ) );
                _discrFuncs.evaluateGradientAtQuadPoint (El, q, _Gradient.at (elementIdx, q));
                _pointwiseEvaluator.evaluateFirstFundamentalForm ( _Gradient.at (elementIdx, q), _firstFF.at (elementIdx, q));
                _firstFFInv.at (elementIdx, q) = (_firstFF.at (elementIdx, q)).inverse();
                _detFirstFF.at (elementIdx, q) = _firstFF.at (elementIdx, q).determinant();
                _area.at (elementIdx, q) = std::sqrt (_detFirstFF.at (elementIdx, q));
                _pointwiseEvaluator.evaluateNormal ( _Gradient.at (elementIdx, q), _normal.at (elementIdx, q) );
                _GradientGInv.at ( elementIdx, q ) = _Gradient.at (elementIdx, q) * _firstFFInv.at (elementIdx, q);
            }
        }
    }

public :
    const VectorType& getDofs( ) const {return _dofs;};
    const Point3DType& getCoords( const int elementIdx, const int QuadPoint ) const { return _coords[elementIdx][QuadPoint]; }
    const Matrix32& getGradient( const int elementIdx, const int QuadPoint ) const { return _Gradient[elementIdx][QuadPoint]; }
    const Matrix22& getFirstFF( const int elementIdx, const int QuadPoint ) const { return _firstFF[elementIdx][QuadPoint]; }
    const Matrix22& getFirstFFInv( const int elementIdx, const int QuadPoint ) const { return _firstFFInv[elementIdx][QuadPoint]; }
    const RealType& getDetFirstFF( const int elementIdx, const int QuadPoint ) const { return _detFirstFF[elementIdx][QuadPoint]; }
    const RealType& getArea( const int elementIdx, const int QuadPoint ) const { return _area[elementIdx][QuadPoint]; }
    const TangentVecType& getNormal( const int elementIdx, const int QuadPoint ) const { return _normal[elementIdx][QuadPoint]; }
    const Matrix32& getGradientGInv( const int elementIdx, const int QuadPoint ) const { return _GradientGInv[elementIdx][QuadPoint]; }

};



template <typename ConfiguratorType>
class DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> : public DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder>
{
 protected:

    typedef typename ConfiguratorType::RealType                 RealType;
    typedef typename ConfiguratorType::Matrix22                 Matrix22;
    typedef typename ConfiguratorType::Matrix32                 Matrix32;
    typedef typename ConfiguratorType::Matrix33                 Matrix33;
    typedef typename ConfiguratorType::TangentVecType           TangentVecType;
    typedef typename ConfiguratorType::Point3DType              Point3DType;
    typedef typename ConfiguratorType::VectorType               VectorType;
    typedef typename ConfiguratorType::Tensor322Type            Tensor322Type;

    typedef typename ConfiguratorType::ElementType ElementType;
    typedef typename ConfiguratorType::DomVecType  DomVecType;

    const ConfiguratorType &_conf;
    const DKTFEVectorFunctionEvaluator <ConfiguratorType> _discrFuncs;
    PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _pointwiseEvaluator;
    SemiNonlinIsometryPointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _semiNonlinIsometryPointwiseEvaluator;
    const int _numberOfElements;
    const int _numberOfQuadPointsPerElement;

protected:

    //ApproxGradient
    DFDStorageStructure <Matrix32>            _ApproxGradient;
    DFDStorageStructure <Matrix22>            _ApproxFirstFF;
    DFDStorageStructure <Matrix22>            _ApproxFirstFFInv;
    DFDStorageStructure <RealType>            _ApproxDetFirstFF;
    DFDStorageStructure <RealType>            _ApproxArea;
    DFDStorageStructure <TangentVecType>      _ApproxNormal;
    DFDStorageStructure <Matrix33>            _HessianAsVec;
    DFDStorageStructure <Tensor322Type>       _Hessian;

    //Mixed
    DFDStorageStructure <Matrix22>            _secondFF;
    DFDStorageStructure <DomVecType>         _gInvChristoffel2;
//     DFDStorageStructure < TangentVecType >     _LaplaceBeltrami;

    // SemiNonlinIsometry
    DFDStorageStructure <TangentVecType>      _SemiNonlinIsometry_a1;
    DFDStorageStructure <TangentVecType>      _SemiNonlinIsometry_a2hat;
    DFDStorageStructure <TangentVecType>      _SemiNonlinIsometry_a2;
    DFDStorageStructure <Matrix22>      _SemiNonlinIsometry_a0tilde;
    DFDStorageStructure <Matrix22>      _SemiNonlinIsometry_a1tilde;
    DFDStorageStructure <Matrix22>      _SemiNonlinIsometry_a2tilde;

  public:
    DiscreteVectorFunctionStorage ( const ConfiguratorType &conf, const VectorType &dofs, const int numComponents ) :
    DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> ( conf, dofs, numComponents ),
      _conf (conf),
      _discrFuncs (conf, dofs, numComponents),
      _numberOfElements (conf.getInitializer ().getNumTriangs ()),
      _numberOfQuadPointsPerElement (conf.maxNumQuadPoints ()),

      //approx
      _ApproxGradient (_numberOfElements, _numberOfQuadPointsPerElement),
      _ApproxFirstFF (_numberOfElements, _numberOfQuadPointsPerElement),
      _ApproxFirstFFInv (_numberOfElements, _numberOfQuadPointsPerElement),
      _ApproxDetFirstFF (_numberOfElements, _numberOfQuadPointsPerElement),
      _ApproxArea ( _numberOfElements, _numberOfQuadPointsPerElement ),
      _ApproxNormal (_numberOfElements, _numberOfQuadPointsPerElement),
      _HessianAsVec ( _numberOfElements, _numberOfQuadPointsPerElement ),
      _Hessian ( _numberOfElements, _numberOfQuadPointsPerElement ),
      //mixed
      _secondFF( _numberOfElements, _numberOfQuadPointsPerElement ),
      _gInvChristoffel2 ( _numberOfElements, _numberOfQuadPointsPerElement ),
      //SemiNonlinIsometry
      _SemiNonlinIsometry_a1 ( _numberOfElements, _numberOfQuadPointsPerElement ),
      _SemiNonlinIsometry_a2hat ( _numberOfElements, _numberOfQuadPointsPerElement ),
      _SemiNonlinIsometry_a2 ( _numberOfElements, _numberOfQuadPointsPerElement ),
      _SemiNonlinIsometry_a0tilde ( _numberOfElements, _numberOfQuadPointsPerElement),
      _SemiNonlinIsometry_a1tilde ( _numberOfElements, _numberOfQuadPointsPerElement),
      _SemiNonlinIsometry_a2tilde ( _numberOfElements, _numberOfQuadPointsPerElement)

    {
        for ( int elementIdx = 0; elementIdx < _conf.getInitializer().getNumTriangs(); ++elementIdx){
            const ElementType& El ( _conf.getInitializer().getTriang( elementIdx ) );
            for ( int localQuadPointIndex = 0; localQuadPointIndex < _numberOfQuadPointsPerElement; ++localQuadPointIndex){
                //approx
                _discrFuncs.evaluateApproxGradientAtQuadPoint (El, localQuadPointIndex, _ApproxGradient.at (elementIdx, localQuadPointIndex));
                _pointwiseEvaluator.evaluateFirstFundamentalForm ( _ApproxGradient.at (elementIdx, localQuadPointIndex), _ApproxFirstFF.at (elementIdx, localQuadPointIndex));
                _ApproxFirstFFInv.at (elementIdx, localQuadPointIndex) = (_ApproxFirstFF.at (elementIdx, localQuadPointIndex)).inverse();
                _ApproxDetFirstFF.at (elementIdx, localQuadPointIndex) = _ApproxFirstFF.at (elementIdx, localQuadPointIndex).determinant ();
                _ApproxArea.at (elementIdx, localQuadPointIndex) = std::sqrt (_ApproxDetFirstFF.at (elementIdx, localQuadPointIndex));
                _pointwiseEvaluator.evaluateNormal ( _ApproxGradient.at (elementIdx, localQuadPointIndex), _ApproxNormal.at (elementIdx, localQuadPointIndex) );
                _discrFuncs.evaluateApproxHessianAsVecAtQuadPoint ( El, localQuadPointIndex, _HessianAsVec.at( elementIdx, localQuadPointIndex ) );
                _discrFuncs.evaluateApproxHessianSymAtQuadPoint (El, localQuadPointIndex, _Hessian.at(elementIdx, localQuadPointIndex));

                //mixed
                _pointwiseEvaluator.evaluateSecondFundamentalForm ( _Hessian.at (elementIdx, localQuadPointIndex), this->_normal.at (elementIdx, localQuadPointIndex), _secondFF.at (elementIdx, localQuadPointIndex));
                _pointwiseEvaluator.evaluateVectorForLaplacian( this->_firstFFInv.at(elementIdx, localQuadPointIndex), this->_Gradient.at(elementIdx, localQuadPointIndex), _Hessian.at(elementIdx, localQuadPointIndex), _gInvChristoffel2.at(elementIdx, localQuadPointIndex) );
//                 _discrFuncs.evaluateLaplaceBeltrami ( _firstFFInv.at(elementIdx, localQuadPointIndex), _Gradient.at(elementIdx, localQuadPointIndex), _Hessian.at(elementIdx, localQuadPointIndex), _gInvChristoffel2.at(elementIdx, localQuadPointIndex), _LaplaceBeltrami.at(elementIdx, localQuadPointIndex ) );


                // SemiNonlinIsometry
                _semiNonlinIsometryPointwiseEvaluator.evaluateA1(  this->_Gradient.at(elementIdx, localQuadPointIndex), _SemiNonlinIsometry_a1.at(elementIdx, localQuadPointIndex) );

                _semiNonlinIsometryPointwiseEvaluator.evaluateA2Hat(  this->_Gradient.at(elementIdx, localQuadPointIndex), _SemiNonlinIsometry_a1.at(elementIdx, localQuadPointIndex), _SemiNonlinIsometry_a2hat.at(elementIdx, localQuadPointIndex) );

                _semiNonlinIsometryPointwiseEvaluator.evaluateA2( _SemiNonlinIsometry_a2hat.at(elementIdx, localQuadPointIndex), _SemiNonlinIsometry_a2.at(elementIdx, localQuadPointIndex)  );

                _semiNonlinIsometryPointwiseEvaluator.evaluateA0Tilde( this->_Hessian.at(elementIdx, localQuadPointIndex), this->_normal.at(elementIdx, localQuadPointIndex), _SemiNonlinIsometry_a0tilde.at(elementIdx, localQuadPointIndex)  );

                _semiNonlinIsometryPointwiseEvaluator.evaluateA1Tilde( this->_Hessian.at(elementIdx, localQuadPointIndex), this->_Gradient.at(elementIdx, localQuadPointIndex), _SemiNonlinIsometry_a1.at(elementIdx, localQuadPointIndex), _SemiNonlinIsometry_a2.at(elementIdx, localQuadPointIndex), _SemiNonlinIsometry_a2hat.at(elementIdx, localQuadPointIndex), _SemiNonlinIsometry_a1tilde.at(elementIdx, localQuadPointIndex)  );

                _semiNonlinIsometryPointwiseEvaluator.evaluateA2Tilde( this->_Hessian.at(elementIdx, localQuadPointIndex), _SemiNonlinIsometry_a2.at(elementIdx, localQuadPointIndex), _SemiNonlinIsometry_a2hat.at(elementIdx, localQuadPointIndex), _SemiNonlinIsometry_a2tilde.at(elementIdx, localQuadPointIndex)  );


            }
        }
    }


public:
    const Matrix32& getApproxGradient ( const int elementIdx, const int QuadPoint ) const {return _ApproxGradient[elementIdx][QuadPoint];}
    const Matrix22& getApproxFirstFF ( const int elementIdx, const int QuadPoint ) const {return _ApproxFirstFF[elementIdx][QuadPoint];}
    const Matrix22& getApproxFirstFFInv ( const int elementIdx, const int QuadPoint ) const {return _ApproxFirstFFInv[elementIdx][QuadPoint];}
    const RealType& getApproxDetFirstFF ( const int elementIdx, const int QuadPoint ) const {return _ApproxDetFirstFF[elementIdx][QuadPoint];}
    const RealType& getApproxArea ( const int elementIdx, const int QuadPoint ) const {return _ApproxArea[elementIdx][QuadPoint];}
    const TangentVecType& getApproxNormal ( const int elementIdx, const int QuadPoint ) const {return _ApproxNormal[elementIdx][QuadPoint];}
    const Matrix33& getHessianAsVec ( const int elementIdx, const int QuadPoint ) const {return _HessianAsVec[elementIdx][QuadPoint];}
    const Tensor322Type& getHessian ( const int elementIdx, const int QuadPoint ) const {return _Hessian[elementIdx][QuadPoint];}

    const Matrix22& getSecondFF ( const int elementIdx, const int QuadPoint ) const {return _secondFF[elementIdx][QuadPoint];}
    const DomVecType& getGInvChristoffel2 ( const int elementIdx, const int QuadPoint ) const {return _gInvChristoffel2[elementIdx][QuadPoint];}

    const TangentVecType& getSemiNonlinIsometry_a1 ( const int elementIdx, const int QuadPoint ) const {return _SemiNonlinIsometry_a1[elementIdx][QuadPoint];}
    const TangentVecType& getSemiNonlinIsometry_a2hat ( const int elementIdx, const int QuadPoint ) const {return _SemiNonlinIsometry_a2hat[elementIdx][QuadPoint];}
    const TangentVecType& getSemiNonlinIsometry_a2 ( const int elementIdx, const int QuadPoint ) const {return _SemiNonlinIsometry_a2[elementIdx][QuadPoint];}
    const Matrix22& getSemiNonlinIsometry_a0tilde ( const int elementIdx, const int QuadPoint ) const {return _SemiNonlinIsometry_a0tilde[elementIdx][QuadPoint];}
    const Matrix22& getSemiNonlinIsometry_a1tilde ( const int elementIdx, const int QuadPoint ) const {return _SemiNonlinIsometry_a1tilde[elementIdx][QuadPoint];}
    const Matrix22& getSemiNonlinIsometry_a2tilde ( const int elementIdx, const int QuadPoint ) const {return _SemiNonlinIsometry_a2tilde[elementIdx][QuadPoint];}

};










template <typename ConfiguratorType, DiscreteFunctionCacheType>
class MixedDiscreteVectorFunctionStorage {};



// template <typename ConfiguratorType>
// class MixedDiscreteVectorFunctionStorage<ConfiguratorType,NoCache>{
//
// protected:
//
//     typedef typename ConfiguratorType::RealType                 RealType;
//     typedef typename ConfiguratorType::Matrix22                 Matrix22;
//     typedef typename ConfiguratorType::Matrix32                 Matrix32;
//     typedef typename ConfiguratorType::Matrix33                 Matrix33;
//     typedef typename ConfiguratorType::TangentVecType           TangentVecType;
//     typedef typename ConfiguratorType::Point3DType              Point3DType;
//     typedef typename ConfiguratorType::VectorType               VectorType;
//     typedef typename ConfiguratorType::Tensor322Type            Tensor322Type;
//
//     typedef typename ConfiguratorType::ElementType ElementType;
//     typedef typename ConfiguratorType::DomVecType  DomVecType;
//
//     const ConfiguratorType &_conf;
//     const DKTFEVectorFunctionEvaluator <ConfiguratorType> _discrFuncs;
//     PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _pointwiseEvaluator;
//     const int _numberOfElements;
//     const int _numberOfQuadPointsPerElement;
//
// public :
//     const VectorType& _dofs;
//
// public:
//     DiscreteVectorFunctionStorage ( const ConfiguratorType &conf, const VectorType &dofs, const int numComponents ) :
//       _conf (conf),
//       _discrFuncs (conf, dofs, numComponents),
//       _numberOfElements (conf.getInitializer ().getNumTriangs ()),
//       _numberOfQuadPointsPerElement (conf.maxNumQuadPoints ()),
//       _dofs( dofs ) {}
// };



// template <typename ConfiguratorType>
// class DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder>
// {
//   protected:
//
//     typedef typename ConfiguratorType::RealType                 RealType;
//     typedef typename ConfiguratorType::Matrix22                 Matrix22;
//     typedef typename ConfiguratorType::Matrix32                 Matrix32;
//     typedef typename ConfiguratorType::Matrix33                 Matrix33;
//     typedef typename ConfiguratorType::TangentVecType           TangentVecType;
//     typedef typename ConfiguratorType::Point3DType              Point3DType;
//     typedef typename ConfiguratorType::VectorType               VectorType;
//     typedef typename ConfiguratorType::Tensor322Type            Tensor322Type;
//
//     typedef typename ConfiguratorType::ElementType ElementType;
//     typedef typename ConfiguratorType::DomVecType  DomVecType;
//
//     const ConfiguratorType &_conf;
//     const DKTFEVectorFunctionEvaluator <ConfiguratorType> _discrFuncs;
//     PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _pointwiseEvaluator;
//     const int _numberOfElements;
//     const int _numberOfQuadPointsPerElement;
//
// protected :
//     const VectorType& _dofs;
//     //
//     DFDStorageStructure <Point3DType>                        _coords;
//     DFDStorageStructure <Matrix32>                           _Gradient;
//     DFDStorageStructure <Matrix22>                           _firstFF;
//     DFDStorageStructure <Matrix22>                           _firstFFInv;
//     DFDStorageStructure <RealType>                           _detFirstFF;
//     DFDStorageStructure <RealType>                           _area;
//     DFDStorageStructure <TangentVecType>                     _normal;
//     DFDStorageStructure <Matrix32>                           _GradientGInv;
//
//   public:
//     DiscreteVectorFunctionStorage ( const ConfiguratorType &conf, const VectorType &dofs, const int numComponents ) :
//       _conf (conf),
//       _discrFuncs (conf, dofs, numComponents),
//       _numberOfElements (conf.getInitializer ().getNumTriangs ()),
//       _numberOfQuadPointsPerElement (conf.maxNumQuadPoints ()),
//       _dofs( dofs ),
//       //
//       _coords ( _numberOfElements, _numberOfQuadPointsPerElement ),
//       _Gradient (_numberOfElements, _numberOfQuadPointsPerElement),
//       _firstFF (_numberOfElements, _numberOfQuadPointsPerElement),
//       _firstFFInv (_numberOfElements, _numberOfQuadPointsPerElement),
//       _detFirstFF (_numberOfElements, _numberOfQuadPointsPerElement),
//       _area ( _numberOfElements, _numberOfQuadPointsPerElement ),
//       _normal (_numberOfElements, _numberOfQuadPointsPerElement),
//       _GradientGInv( _numberOfElements, _numberOfQuadPointsPerElement )
//     {
//         for ( int elementIdx = 0; elementIdx < _conf.getInitializer().getNumTriangs(); ++elementIdx){
//             const ElementType& El ( _conf.getInitializer().getTriang( elementIdx ) );
//             for ( int q = 0; q < _numberOfQuadPointsPerElement; ++q){
//                 //
//                 _discrFuncs.evaluateAtQuadPoint ( El, q, _coords.at(elementIdx, q ) );
//                 _discrFuncs.evaluateGradientAtQuadPoint (El, q, _Gradient.at (elementIdx, q));
//                 _pointwiseEvaluator.evaluateFirstFundamentalForm ( _Gradient.at (elementIdx, q), _firstFF.at (elementIdx, q));
//                 _firstFFInv.at (elementIdx, q) = (_firstFF.at (elementIdx, q)).inverse();
//                 _detFirstFF.at (elementIdx, q) = _firstFF.at (elementIdx, q).determinant();
//                 _area.at (elementIdx, q) = std::sqrt (_detFirstFF.at (elementIdx, q));
//                 _pointwiseEvaluator.evaluateNormal ( _Gradient.at (elementIdx, q), _normal.at (elementIdx, q) );
//                 _GradientGInv.at ( elementIdx, q ) = _Gradient.at (elementIdx, q) * _firstFFInv.at (elementIdx, q);
//             }
//         }
//     }
//
// public :
//     const VectorType& getDofs( ) const {return _dofs;};
//     const Point3DType& getCoords( const int elementIdx, const int QuadPoint ) const { return _coords[elementIdx][QuadPoint]; }
//     const Matrix32& getGradient( const int elementIdx, const int QuadPoint ) const { return _Gradient[elementIdx][QuadPoint]; }
//     const Matrix22& getFirstFF( const int elementIdx, const int QuadPoint ) const { return _firstFF[elementIdx][QuadPoint]; }
//     const Matrix22& getFirstFFInv( const int elementIdx, const int QuadPoint ) const { return _firstFFInv[elementIdx][QuadPoint]; }
//     const RealType& getDetFirstFF( const int elementIdx, const int QuadPoint ) const { return _detFirstFF[elementIdx][QuadPoint]; }
//     const RealType& getArea( const int elementIdx, const int QuadPoint ) const { return _area[elementIdx][QuadPoint]; }
//     const TangentVecType& getNormal( const int elementIdx, const int QuadPoint ) const { return _normal[elementIdx][QuadPoint]; }
//     const Matrix32& getGradientGInv( const int elementIdx, const int QuadPoint ) const { return _GradientGInv[elementIdx][QuadPoint]; }
//
// };



template <typename ConfiguratorType>
class MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder>
{
 protected:

    typedef typename ConfiguratorType::RealType                 RealType;
    typedef typename ConfiguratorType::Matrix22                 Matrix22;
    typedef typename ConfiguratorType::Matrix32                 Matrix32;
    typedef typename ConfiguratorType::Matrix33                 Matrix33;
    typedef typename ConfiguratorType::TangentVecType           TangentVecType;
    typedef typename ConfiguratorType::Point3DType              Point3DType;
    typedef typename ConfiguratorType::VectorType               VectorType;
    typedef typename ConfiguratorType::Tensor322Type            Tensor322Type;

    typedef typename ConfiguratorType::ElementType ElementType;
    typedef typename ConfiguratorType::DomVecType  DomVecType;

    const ConfiguratorType &_conf;
    const DiscreteVectorFunctionStorage<ConfiguratorType, FirstAndSecondOrder> &_xAStorage,  &_xBStorage;
    PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _pointwiseEvaluator;
    SemiNonlinIsometryPointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _semiNonlinIsometryPointwiseEvaluator;
    const int _numberOfElements;
    const int _numberOfQuadPointsPerElement;

protected:

    //ApproxGradient
    DFDStorageStructure <Matrix22>            _gAInvSecondFFB;
    DFDStorageStructure <Tensor322Type>            _semiNonlinearityB;
    DFDStorageStructure <Tensor322Type>            _gAInvSemiNonlinearityB;

  public:
    MixedDiscreteVectorFunctionStorage ( const ConfiguratorType &conf,
                    const DiscreteVectorFunctionStorage<ConfiguratorType, FirstAndSecondOrder> &xAStorage,
                    const DiscreteVectorFunctionStorage<ConfiguratorType, FirstAndSecondOrder> &xBStorage,
                    const int numComponents ) :
      _conf (conf),
      _xAStorage ( xAStorage ),  _xBStorage ( xBStorage ),
      _numberOfElements (conf.getInitializer ().getNumTriangs ()),
      _numberOfQuadPointsPerElement (conf.maxNumQuadPoints ()),

      _semiNonlinearityB (_numberOfElements, _numberOfQuadPointsPerElement),
      _gAInvSemiNonlinearityB (_numberOfElements, _numberOfQuadPointsPerElement),
      _gAInvSecondFFB (_numberOfElements, _numberOfQuadPointsPerElement)

    {
        for ( int elementIdx = 0; elementIdx < _conf.getInitializer().getNumTriangs(); ++elementIdx){
            const ElementType& El ( _conf.getInitializer().getTriang( elementIdx ) );
            for ( int localQuadPointIndex = 0; localQuadPointIndex < _numberOfQuadPointsPerElement; ++localQuadPointIndex){
                //approx
                _gAInvSecondFFB.at (elementIdx, localQuadPointIndex) = _xAStorage.getFirstFFInv(elementIdx, localQuadPointIndex) * _xBStorage.getSecondFF(elementIdx, localQuadPointIndex );
                for ( int l = 0; l<3; ++l ) {
                         Matrix22 mat_temp; mat_temp.setZero();

                         Matrix22 gAinvD2xB = _xAStorage.getFirstFFInv( elementIdx, localQuadPointIndex ) * _xBStorage.getHessian( elementIdx, localQuadPointIndex )[l];
                         mat_temp += gAinvD2xB;

                         mat_temp -= _xBStorage.getNormal( elementIdx, localQuadPointIndex )[l] * _xAStorage.getFirstFFInv( elementIdx, localQuadPointIndex ) * _xAStorage.getSemiNonlinIsometry_a0tilde( elementIdx, localQuadPointIndex );

                         mat_temp -= _xBStorage.getGradient( elementIdx, localQuadPointIndex )(l, 0) * _xAStorage.getFirstFFInv( elementIdx, localQuadPointIndex ) * _xAStorage.getSemiNonlinIsometry_a1tilde( elementIdx, localQuadPointIndex );

                         mat_temp -= _xBStorage.getGradient( elementIdx, localQuadPointIndex )(l, 1) * _xAStorage.getFirstFFInv( elementIdx, localQuadPointIndex ) * _xAStorage.getSemiNonlinIsometry_a2tilde( elementIdx, localQuadPointIndex );

                         _semiNonlinearityB.at (elementIdx, localQuadPointIndex)[l] = mat_temp;
                         _gAInvSemiNonlinearityB.at (elementIdx, localQuadPointIndex)[l] = _xAStorage.getFirstFFInv( elementIdx, localQuadPointIndex ) * mat_temp;
                      }

            }
        }
    }

      // _semiNonlinearityB (_numberOfElements, _numberOfQuadPointsPerElement)

    // {
    //   for ( int elementIdx = 0; elementIdx < _conf.getInitializer().getNumTriangs(); ++elementIdx){
    //     const ElementType& El ( _conf.getInitializer().getTriang( elementIdx ) );
    //     for ( int localQuadPointIndex = 0; localQuadPointIndex < _numberOfQuadPointsPerElement; ++localQuadPointIndex){
    //
    //       for ( int l = 0; l<3; ++l ) {
    //          Matrix22 mat_temp; mat_temp.setZero();
    //
    //          Matrix22 gAinvD2xB = _xAStorage.getFirstFFInv( elementIdx, localQuadPointIndex ) * _xBStorage.getHessian( elementIdx, localQuadPointIndex )[l];
    //          mat_temp += gAinvD2xB;
    //
    //          mat_temp -= _xBStorage.getNormal( elementIdx, localQuadPointIndex )[l] * _xAStorage.getFirstFFInv( elementIdx, localQuadPointIndex ) * _xAStorage.getSemiNonlinIsometry_a0tilde( elementIdx, localQuadPointIndex );
    //
    //          mat_temp -= _xBStorage.getGradient( elementIdx, localQuadPointIndex )(l, 0) * _xAStorage.getFirstFFInv( elementIdx, localQuadPointIndex ) * _xAStorage.getSemiNonlinIsometry_a1tilde( elementIdx, localQuadPointIndex );
    //
    //          mat_temp -= _xBStorage.getGradient( elementIdx, localQuadPointIndex )(l, 1) * _xAStorage.getFirstFFInv( elementIdx, localQuadPointIndex ) * _xAStorage.getSemiNonlinIsometry_a2tilde( elementIdx, localQuadPointIndex );
    //
    //          _semiNonlinearityB.at (elementIdx, localQuadPointIndex)[l] = mat_temp;
    //       }
    //     }
    //   }
    //
    // }

public:
    const Matrix22& getGAInvSecondFFB ( const int elementIdx, const int QuadPoint ) const {
        return _gAInvSecondFFB[elementIdx][QuadPoint];
    }
    const Tensor322Type& getSemiNonlinearityB (const int elementIdx, const int QuadPoint) const{
        return _semiNonlinearityB[elementIdx][QuadPoint];
    }
    const Tensor322Type& getGAInvSemiNonlinearityB (const int elementIdx, const int QuadPoint) const{
        return _gAInvSemiNonlinearityB[elementIdx][QuadPoint];
    }

};






#endif
