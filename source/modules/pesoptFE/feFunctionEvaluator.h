#ifndef __FEFUNCTIONEVALUATOR_H
#define __FEFUNCTIONEVALUATOR_H

#include <initializer_list>
#include <pesopt_IO.h>
  
template <typename ConfiguratorType, typename VectorType, typename NLTYPE >
struct FEFunctionEvaluatorLookup {
  
  void getLocalDof ( const ConfiguratorType & conf, const VectorType &dofs, const typename ConfiguratorType::ElementType &El,  
                     const int bfNum,  const NLTYPE &bfValue, NLTYPE &aux ) {
      aux = bfValue;
      aux *= dofs[ conf.localToGlobal ( El, bfNum ) ];
  }

};
  

template <typename ConfiguratorType, typename VectorType = typename ConfiguratorType::VectorType >
class FEScalarFunctionEvaluator {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::RealVecChart  RealVecChart;
  typedef typename ConfiguratorType::DerivativeVectorValuedType  DerivativeVectorValuedType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::BaseFuncSetType BaseFuncSetType;

  const ConfiguratorType & _conf;
  const VectorType & _dofs;

  FEScalarFunctionEvaluator ( const ConfiguratorType & config, const VectorType & Dofs )
    : _conf ( config ),
      _dofs ( Dofs ) { }
    
  RealType evaluate ( const ElementType &El, const RealVecChart& RefCoord ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    RealType w = 0.;
    RealType v, aux;
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluate ( b, RefCoord );
      FEFunctionEvaluatorLookup<ConfiguratorType, VectorType, RealType>().getLocalDof( _conf, _dofs, El, b, v, aux );
      w += aux;
    }
    return w;
  }

  RealType evaluateAtQuadPoint ( const ElementType &El, int QuadPoint ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    RealType w = 0.;
    RealType v, aux;
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluate ( b, QuadPoint );
      FEFunctionEvaluatorLookup<ConfiguratorType, VectorType, RealType>().getLocalDof ( _conf, _dofs, El, b, v, aux );
      w += aux;
    }
    return w;
  }
  
  void evaluateGradient ( const ElementType &El, const RealVecChart& RefCoord, RealVecChart& Grad ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    RealVecChart v, aux;
    Grad.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      bfs.evaluateGradient ( b, RefCoord, v );
      FEFunctionEvaluatorLookup<ConfiguratorType, VectorType, RealVecChart>().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Grad += aux;
    }
  }
  
  void evaluateGradientAtQuadPoint ( const ElementType &El, int QuadPoint, RealVecChart& Grad ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    RealVecChart v, aux;
    Grad.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateGradient ( b, QuadPoint );
      FEFunctionEvaluatorLookup<ConfiguratorType, VectorType, RealVecChart >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Grad += aux;
    }
  }
  
  
  void evaluateHessian ( const ElementType &El, const RealVecChart& RefCoord, DerivativeVectorValuedType& Hessian ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    DerivativeVectorValuedType v, aux;
    Hessian.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      bfs.evaluateHessian ( b, RefCoord, v );
      FEFunctionEvaluatorLookup<ConfiguratorType, VectorType, DerivativeVectorValuedType>().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Hessian += aux;
    }
  }
  
  void evaluateHessianAtQuadPoint ( const ElementType &El, int QuadPoint, DerivativeVectorValuedType& Hessian ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    DerivativeVectorValuedType v, aux;
    Hessian.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateHessian ( b, QuadPoint );
      FEFunctionEvaluatorLookup<ConfiguratorType, VectorType, DerivativeVectorValuedType >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Hessian += aux;
    }
  }
  
  
  const VectorType& getDofs (  ) const {return _dofs;}

};


template <typename ConfiguratorType, typename VectorType = typename ConfiguratorType::VectorType >
class FEScalarFunctionEvaluatorFace {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::RealVecChart  RealVecChart;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::BaseFuncSetType BaseFuncSetType;

  const ConfiguratorType & _conf;
  const VectorType & _dofs;

  FEScalarFunctionEvaluatorFace ( const ConfiguratorType & config, const VectorType & Dofs , const VectorType &Dofs2)
    : _conf ( config ),
      _dofs ( Dofs ) { }
    
  RealType evaluate ( const ElementType &El, const RealVecChart& RefCoord ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    RealType w = 0.;
    RealType v, aux;
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluate ( b, RefCoord );
      FEFunctionEvaluatorLookup<ConfiguratorType, VectorType, RealType>().getLocalDof( _conf, _dofs, El, b, v, aux );
      w += aux;
    }
    return w;
  }

  RealType evaluateAtQuadPoint ( const ElementType &El, int QuadPoint ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    RealType w = 0.;
    RealType v, aux;
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluate ( b, QuadPoint );
      FEFunctionEvaluatorLookup<ConfiguratorType, VectorType, RealType>().getLocalDof ( _conf, _dofs, El, b, v, aux );
      w += aux;
    }
    return w;
  }
  
  void evaluateGradient ( const ElementType &El, const RealVecChart& RefCoord, RealVecChart& Grad ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    RealVecChart v, aux;
    Grad.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      bfs.evaluateGradient ( b, RefCoord, v );
      FEFunctionEvaluatorLookup<ConfiguratorType, VectorType, RealVecChart>().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Grad += aux;
    }
  }
  
  void evaluateGradientAtQuadPoint ( const ElementType &El, int QuadPoint, RealVecChart& Grad ) const {
    const BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    RealVecChart v, aux;
    Grad.setZero();
    for ( int b = 0; b < static_cast<int> ( _conf.getNumLocalDofs ( El ) ); ++b ) {
      v = bfs.evaluateGradient ( b, QuadPoint );
      FEFunctionEvaluatorLookup<ConfiguratorType, VectorType, RealVecChart >().getLocalDof ( _conf, _dofs, El, b, v, aux );
      Grad += aux;
    }
  }
  
  const VectorType& getDofs (  ) const {return _dofs;}

};




template <typename ConfiguratorType, int numComponents = ConfiguratorType::dimChartDomain>
class FEVectorFunctionEvaluator {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::RealVecChart  RealVecChart;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::Matrix22  Matrix22;
  typedef typename ConfiguratorType::Matrix32  Matrix32;
  typedef typename ConfiguratorType::Matrix33  Matrix33;
  typedef typename ConfiguratorType::DerivativeVectorValuedType        DerivativeVectorValuedType;
  
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef FEScalarFunctionEvaluator<ConfiguratorType, Eigen::Ref<const VectorType> > DiscFuncType;

  mutable std::vector<Eigen::Ref<const VectorType> > _refs;
  mutable std::vector<DiscFuncType> _discrFuncs;

public:
  
  FEVectorFunctionEvaluator ( const ConfiguratorType &Configurator,  const VectorType &Dofs )
  {
     const int numGlobalDofs = Configurator.getNumGlobalDofs();
     _refs.reserve( numComponents );
     _discrFuncs.reserve( numComponents );

     for( int c=0; c < numComponents; ++c )
         _refs.push_back ( Dofs.segment( c * numGlobalDofs, numGlobalDofs) );
     
     for( int c=0; c < numComponents; ++c )
     _discrFuncs.push_back ( DiscFuncType ( Configurator, _refs[c] ) );
 
  }

  void evaluate( const ElementType &El, const RealVecChart& RefCoord, RealVecChart &Value ) const {
    for ( int c = 0; c < numComponents; ++c )
      Value[c] = _discrFuncs[c].evaluate ( El, RefCoord );
  }
  
  void evaluateAtQuadPoint ( const ElementType &El, int QuadPoint, RealVecChart &Value ) const {
    for ( int c = 0; c < numComponents; ++c )
      Value[c] = _discrFuncs[c].evaluateAtQuadPoint ( El, QuadPoint );
  }

  void evaluateGradient ( const ElementType &El, const RealVecChart& RefCoord, DerivativeVectorValuedType &Dx ) const {
    RealVecChart v;
    for ( int c = 0; c < numComponents; c++ ) {
      _discrFuncs[c].evaluateGradient ( El, RefCoord, v );
      for( int d=0; d < numComponents; d++ ){
          Dx( c, d ) = v[d];
      }
    }
  }

  void evaluateGradientAtQuadPoint ( const ElementType &El, int QuadPoint, DerivativeVectorValuedType &Dx ) const {
    RealVecChart v;
    for ( int c = 0; c < numComponents; ++c ) {
      _discrFuncs[c].evaluateGradientAtQuadPoint ( El, QuadPoint, v );
      for( int d=0; d < numComponents; ++d ) Dx( c, d ) = v[d];
    }
  }
  
  // = (Dx + Dx^T)/2
  void evaluateSymmetrizedGradient( const ElementType &El, const RealVecChart& RefCoord, DerivativeVectorValuedType &DxSym ) const {
    DerivativeVectorValuedType Dx; this->evaluateGradient(El, RefCoord, Dx);
    DxSym = 0.5 * ( Dx + Dx.transpose() );
  }
  
  void evaluateSymmetrizedGradientAtQuadPoint ( const ElementType &El, int QuadPoint, DerivativeVectorValuedType &DxSym ) const {
    DerivativeVectorValuedType Dx; this->evaluateGradientAtQuadPoint(El,QuadPoint,Dx);
    DxSym = 0.5 * ( Dx + Dx.transpose() );
  }

  RealType evaluateDivergence ( const ElementType &El, const RealVecChart& RefCoord ) const {
    RealType divergence = 0.;
    RealVecChart v;
    for ( int c = 0; c < numComponents; ++c ) {
      _discrFuncs[c].evaluateGradient ( El, RefCoord, v );
      divergence += v[c];
    }
    return divergence;
  }
  
  RealType evaluateDivergenceAtQuadPoint ( const ElementType &El, int QuadPoint ) const {
    RealType divergence = 0.;
    RealVecChart v;
    for ( int c = 0; c < numComponents; ++c ) {
      _discrFuncs[c].evaluateGradientAtQuadPoint ( El, QuadPoint, v );
      divergence += v[c];
    }
    return divergence;
  }
  
  
  const DiscFuncType& operator[] ( int i ) const { return _discrFuncs[i];}
  DiscFuncType& operator[] ( int i ) { return _discrFuncs[i];}

};


























template <typename DataTypeContainer, const int dimChartDomain>
class FEAffineFunctionEvaluator { };



template <typename DataTypeContainer>
class FEAffineFunctionEvaluator<DataTypeContainer, 1> {
public:
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::PointType PointType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::DerivativeVectorValuedType  DerivativeVectorValuedType;

  const VectorType & _dofs;
  mutable DerivativeVectorValuedType _symGrad;
  mutable RealType _div;

  FEAffineFunctionEvaluator ( const VectorType & Dofs ) : _dofs ( Dofs ) { evaluateSymGradient(); evaluateDiv(); }
  
protected:
  void evaluateSymGradient ( ) const {
     _symGrad(0,0) = _dofs[0];   
  }
  
  void evaluateDiv ( ) const {
      _div = _dofs[0];
  }
  
public:
  DerivativeVectorValuedType & getSymGrad( ) const {return _symGrad;}
  RealType getDiv( ) const {return _div;}

  // computes u(x) = symgrad x for all nodes x
  template<typename MeshType>
  void getAffineDisplacementAtNodes( const MeshType &mesh, VectorType &dispAffine ) const {
      typename DataTypeContainer::DerivativeVectorValuedType symMat = this->getSymGrad();
      for( int nodeIdx = 0; nodeIdx < mesh.getNumVertices(); ++nodeIdx ){
        PointType coords = mesh.getVertex( nodeIdx );
        PointType offsetAffine = symMat * coords;
        for( int comp = 0; comp < coords.size(); ++comp ) 
            dispAffine[nodeIdx + comp * mesh.getNumVertices()] = offsetAffine[comp];
     }
  }
  
};

template <typename DataTypeContainer>
class FEAffineFunctionEvaluator<DataTypeContainer, 2> {
public:
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::PointType PointType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::DerivativeVectorValuedType  DerivativeVectorValuedType;

  const VectorType & _dofs;
  mutable DerivativeVectorValuedType _symGrad;
  mutable RealType _div;

  FEAffineFunctionEvaluator ( const VectorType & Dofs ) : _dofs ( Dofs ) { evaluateSymGradient(); evaluateDiv(); }
  
protected:
  void evaluateSymGradient ( ) const {
     _symGrad(0,0) = _dofs[0];   
     _symGrad(1,1) = _dofs[1]; 
     _symGrad(1,0) = _dofs[2]; 
     _symGrad(0,1) = _dofs[2]; 
  }
  
  void evaluateDiv ( ) const {
      _div = _dofs[0] + _dofs[1];
  }
  
public:
  DerivativeVectorValuedType & getSymGrad( ) const {return _symGrad;}
  RealType getDiv( ) const {return _div;}

  // computes u(x) = symgrad x for all nodes x
  template<typename MeshType>
  void getAffineDisplacementAtNodes( const MeshType &mesh, VectorType &dispAffine ) const {
      typename DataTypeContainer::DerivativeVectorValuedType symMat = this->getSymGrad();
      for( int nodeIdx = 0; nodeIdx < mesh.getNumVertices(); ++nodeIdx ){
        PointType coords = mesh.getVertex( nodeIdx );
        PointType offsetAffine = symMat * coords;
        for( int comp = 0; comp < coords.size(); ++comp ) 
            dispAffine[nodeIdx + comp * mesh.getNumVertices()] = offsetAffine[comp];
     }
  }
  
};


template <typename DataTypeContainer>
class FEAffineFunctionEvaluator<DataTypeContainer, 3> {
public:
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::PointType PointType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::DerivativeVectorValuedType  DerivativeVectorValuedType;

  const VectorType & _dofs;
  mutable DerivativeVectorValuedType _symGrad;
  mutable RealType _div;

  FEAffineFunctionEvaluator ( const VectorType & Dofs ) : 
   _dofs ( Dofs ) { 
       evaluateSymGradient(); 
       evaluateDiv(); 
    }
  
protected:
  // 0 3 5
  //   1 4
  //     2
  // instead of
  // 0 3 4
  //   1 5
  //     2
  void evaluateSymGradient ( ) const {
     _symGrad(0,0) = _dofs[0];   
     _symGrad(1,1) = _dofs[1]; 
     _symGrad(2,2) = _dofs[2]; 
     _symGrad(1,0) = _dofs[3]; 
     _symGrad(0,1) = _dofs[3]; 
//      _symGrad(2,0) = _dofs[4]; 
//      _symGrad(0,2) = _dofs[4]; 
//      _symGrad(1,2) = _dofs[5]; 
//      _symGrad(2,1) = _dofs[5]; 
     _symGrad(2,0) = _dofs[5]; 
     _symGrad(0,2) = _dofs[5]; 
     _symGrad(1,2) = _dofs[4]; 
     _symGrad(2,1) = _dofs[4]; 
  }
  
  void evaluateDiv ( ) const {
      _div = _dofs[0] + _dofs[1] + _dofs[2];
  }
  
public:
  DerivativeVectorValuedType & getSymGrad( ) const {return _symGrad;}
  RealType getDiv( ) const {return _div;}
  
  // computes u(x) = symgrad x for all nodes x
  template<typename MeshType>
  void getAffineDisplacementAtNodes( const MeshType &mesh, VectorType &dispAffine ) const {
      typename DataTypeContainer::DerivativeVectorValuedType symMat = this->getSymGrad();
      for( int nodeIdx = 0; nodeIdx < mesh.getNumVertices(); ++nodeIdx ){
        PointType coords = mesh.getVertex( nodeIdx );
        PointType offsetAffine = symMat * coords;
        for( int comp = 0; comp < coords.size(); ++comp ) 
            dispAffine[nodeIdx + comp * mesh.getNumVertices()] = offsetAffine[comp];
     }
  }
  

};



#endif
