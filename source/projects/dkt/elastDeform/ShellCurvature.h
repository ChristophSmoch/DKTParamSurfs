#ifndef __SHELLCURVATURE_H
#define __SHELLCURVATURE_H

#include <pesopt_DKT.h>

// computes int \sqrt(g) |det(g^-1 h )|
template<typename ConfiguratorType>
class GaussCurvatureL1 :
public RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xStorage;
  public:
    GaussCurvatureL1 ( const ConfiguratorType &conf,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xStorage ) :
      RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xStorage ( xStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gInv = _xStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& h = _xStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType GaussCurv = ( gInv * h ).determinant();
        return _xStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * std::abs( GaussCurv );
    }
};


// computes int \sqrt(g_A) |det(g_B^-1 h_B) - det(g_A^-1 h_A )|
template<typename ConfiguratorType>
class GaussCurvatureL1Diff :
public RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1Diff<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xBStorage;
  public:
    GaussCurvatureL1Diff ( const ConfiguratorType &conf,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xBStorage ) :
      RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1Diff<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xAStorage ( xAStorage ),
      _xBStorage ( xBStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gAInv = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& hA = _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType GaussCurvA = ( gAInv * hA ).determinant();
        const Matrix22& gBInv = _xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& hB = _xBStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType GaussCurvB = ( gBInv * hB ).determinant();
        return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * std::abs( GaussCurvB - GaussCurvA );
        // return _xBStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * std::abs( GaussCurvB );
    }
};


// computes int \sqrt(g) |tr(g^-1 h )|
template<typename ConfiguratorType>
class MeanCurvatureL1 :
public RefTriangleIntegrator<ConfiguratorType, MeanCurvatureL1<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xStorage;
  public:
    MeanCurvatureL1 ( const ConfiguratorType &conf,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xStorage ) :
      RefTriangleIntegrator<ConfiguratorType, MeanCurvatureL1<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xStorage ( xStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gInv = _xStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& h = _xStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType MeanCurv = ( gInv * h ).trace();
        return _xStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * std::abs( MeanCurv );
    }
};


// computes int \sqrt(g_A) |g_A^-1 (h_B -h_A )|^2
template<typename ConfiguratorType>
class RelativeShapeOperatorL2 :
public RefTriangleIntegrator<ConfiguratorType, RelativeShapeOperatorL2<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage, &_xBStorage;
  public:
    RelativeShapeOperatorL2 ( const ConfiguratorType &conf,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xBStorage ) :
      RefTriangleIntegrator<ConfiguratorType, RelativeShapeOperatorL2<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xAStorage ( xAStorage ),
      _xBStorage ( xBStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gAInv = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& hA = _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& hB = _xBStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType aux = ( gAInv * (hB - hA) ).squaredNorm();
        return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * aux;
    }
};



// computes int \sqrt(g_A) |D^2 (x_B - x_A)|^2
template<typename ConfiguratorType>
class SecondDerivativeEnergy :
public RefTriangleIntegrator<ConfiguratorType, SecondDerivativeEnergy<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage, &_xBStorage;
  public:
    SecondDerivativeEnergy ( const ConfiguratorType &conf,
                             const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                             const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xBStorage ) :
      RefTriangleIntegrator<ConfiguratorType, SecondDerivativeEnergy<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xAStorage ( xAStorage ),
      _xBStorage ( xBStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gAInv = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Tensor322Type& D2xA = _xAStorage.getHessian( El.getGlobalElementIdx(), QuadPoint );
        const Tensor322Type& D2xB = _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint );
        Tensor322Type D2u;
        for( int i=0; i<3; ++i )
            for( int j=0; j<3; ++j )
                for( int k = 0; k <2; ++k )
                    D2u.set(i,j,k, D2xB.get(i,j,k) - D2xA.get(i,j,k) );

        const RealType aux =  D2u.normSqr();
        return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * aux;
    }
};


#endif
