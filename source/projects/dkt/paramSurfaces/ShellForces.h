#ifndef __FORCESSHELLFE_H
#define __FORCESSHELLFE_H

#include <pesopt_DKT.h>

    
template< typename ConfiguratorType >
class createConstantLoad :
public RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createConstantLoad <ConfiguratorType> >
{      
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::TangentVecType TangentVecType;
  
    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
    const TangentVecType &_forceVec;

  public:
    createConstantLoad ( const ConfiguratorType &conf,
                         const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                         const TangentVecType & forceVec ) :
     RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createConstantLoad <ConfiguratorType> > (conf),
     _xAStorage ( xAStorage ),
     _forceVec( forceVec ) { }
      
    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, TangentVecType &NL) const {   
      NL = _forceVec;
      NL *= _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};


template< typename ConfiguratorType >
class createConstantLoadInCube :
public RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createConstantLoadInCube <ConfiguratorType> >
{      
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::DomVecType DomVecType;
    typedef typename ConfiguratorType::TangentVecType TangentVecType;
  
    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
    const TangentVecType &_forceVec;
    const DomVecType _xRange, _yRange, _zRange;

  public:
    createConstantLoadInCube ( const ConfiguratorType &conf,
                         const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                         const TangentVecType &forceVec,
                         const DomVecType &xRange, const DomVecType &yRange, const DomVecType &zRange ) :
     RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createConstantLoadInCube <ConfiguratorType> > (conf),
     _xAStorage ( xAStorage ),
     _forceVec ( forceVec ),
     _xRange ( xRange ), _yRange ( yRange ), _zRange ( zRange ) { }
     
    createConstantLoadInCube ( const ConfiguratorType &conf,
                         const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                         const TangentVecType &forceVec,
                         const TangentVecType &startPoint, const TangentVecType &endPoint ) :
     RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createConstantLoadInCube <ConfiguratorType> > (conf),
     _xAStorage ( xAStorage ),
     _forceVec ( forceVec ),
     _xRange ( startPoint[0], endPoint[0] ), _yRange ( startPoint[1], endPoint[1] ), _zRange ( startPoint[2], endPoint[2] ) { }
      
    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, TangentVecType &NL) const {   
      NL.setZero();
      TangentVecType coords ( _xAStorage.getCoords(El.getGlobalElementIdx(),QuadPoint) );
      if( (coords[0] >= _xRange[0]) && (coords[0] <= _xRange[1]) && (coords[1] >= _yRange[0]) && (coords[1] <= _yRange[1]) && (coords[2] >= _zRange[0]) && (coords[2] <= _zRange[1]) ){
        NL = _forceVec;
        NL *= _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
      }
    }
};



template< typename ConfiguratorType >
class createConstantLoadInBall :
public RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createConstantLoadInBall <ConfiguratorType> >
{      
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::DomVecType DomVecType;
    typedef typename ConfiguratorType::TangentVecType TangentVecType;
  
    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
    const TangentVecType &_forceVec;
    const RealType &_radius;
    const TangentVecType &_center;

  public:
    createConstantLoadInBall ( const ConfiguratorType &conf,
                         const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                         const TangentVecType &forceVec,
                         const RealType &radius, const TangentVecType &center ) :
     RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createConstantLoadInBall <ConfiguratorType> > (conf),
     _xAStorage ( xAStorage ),
     _forceVec ( forceVec ),
     _radius ( radius ), _center ( center ) { }
      
    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, TangentVecType &NL) const {   
      NL.setZero();
      TangentVecType coords ( _xAStorage.getCoords(El.getGlobalElementIdx(),QuadPoint) );
      const RealType distToCenter = (coords - _center).norm();
      if( distToCenter <= _radius ){
        NL = _forceVec;
        NL *= _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
      }
    }
};


template< typename ConfiguratorType >
class createNormalLoad :
public RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createNormalLoad <ConfiguratorType> >
{      
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::TangentVecType TangentVecType;
  
    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
    const RealType &_factor;

  public:
    createNormalLoad ( const ConfiguratorType &conf,
                         const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                         const RealType &factor ) :
     RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createNormalLoad <ConfiguratorType> > (conf),
     _xAStorage ( xAStorage ),
     _factor( _factor ) { }
      
    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, TangentVecType &NL) const {   
      NL = _xAStorage.getNormal(El.getGlobalElementIdx(),QuadPoint);
      NL *= _factor * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};


template< typename ConfiguratorType >
class createNormalLoadInCube :
public RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createNormalLoadInCube <ConfiguratorType> >
{      
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::DomVecType DomVecType;
    typedef typename ConfiguratorType::TangentVecType TangentVecType;
  
    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
    const RealType &_factor;
    const DomVecType _xRange, _yRange, _zRange;

  public:
    createNormalLoadInCube ( const ConfiguratorType &conf,
                         const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                         const RealType &factor,
                         const DomVecType &xRange, const DomVecType &yRange, const DomVecType &zRange ) :
     RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createNormalLoadInCube <ConfiguratorType> > (conf),
     _xAStorage ( xAStorage ),
     _factor ( factor ),
     _xRange ( xRange ), _yRange ( yRange ), _zRange ( zRange ) { }
     
    createNormalLoadInCube ( const ConfiguratorType &conf,
                         const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                         const RealType &factor,
                         const TangentVecType &startPoint, const TangentVecType &endPoint ) :
     RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createNormalLoadInCube <ConfiguratorType> > (conf),
     _xAStorage ( xAStorage ),
     _factor ( factor ),
     _xRange ( startPoint[0], endPoint[0] ), _yRange ( startPoint[1], endPoint[1] ), _zRange ( startPoint[2], endPoint[2] ) { }
      
    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, TangentVecType &NL) const {   
      NL.setZero();
      TangentVecType coords ( _xAStorage.getCoords(El.getGlobalElementIdx(),QuadPoint) );
      if( (coords[0] >= _xRange[0]) && (coords[0] <= _xRange[1]) && (coords[1] >= _yRange[0]) && (coords[1] <= _yRange[1]) && (coords[2] >= _zRange[0]) && (coords[2] <= _zRange[1]) ){
        NL = _xAStorage.getNormal(El.getGlobalElementIdx(),QuadPoint);
        NL *= _factor * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
      }
    }
};



template <typename ConfiguratorType>
void createForceShellFE ( const pesopt::BoostParser & parser,
                          const ConfiguratorType & conf,
                          const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                          typename ConfiguratorType::VectorType & rhs_force,
                          const typename ConfiguratorType::MaskType &boundaryMask,
                          const bool collapseBoundaryValues = true ) {
  
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::DomVecType DomVecType;
    typedef typename ConfiguratorType::TangentVecType TangentVecType;
    
    int numGlobalDofs = conf.getNumGlobalDofs();
    rhs_force.setZero();
    const int numLoads = parser.template get<int>( "Force.numLoads" );
    for( int load=1; load<=numLoads; ++load ){
        switch( parser.template get<int> ( pesopt::strprintf( "Force.Type%d", load ).c_str() ) ){
            //in cube
            case 1:{
                TangentVecType forceVec; parser.template getFixSizeVector<RealType,TangentVecType> ( pesopt::strprintf( "Force.Load%d", load ).c_str(), forceVec );
                TangentVecType startPoint; parser.template getFixSizeVector<RealType,TangentVecType> ( pesopt::strprintf( "Force.StartRange%d", load ).c_str(), startPoint );
                TangentVecType endPoint; parser.template getFixSizeVector<RealType,TangentVecType> ( pesopt::strprintf( "Force.EndRange%d", load ).c_str(), endPoint );
                createConstantLoadInCube<ConfiguratorType>( conf, xAStorage, forceVec, startPoint, endPoint ).assembleAdd( rhs_force );
            }break;
            
            //in ball
            case 2:{
             TangentVecType forceVec; parser.template getFixSizeVector<RealType,TangentVecType> ( pesopt::strprintf( "Force.Load%d", load ).c_str(), forceVec );
             RealType radius = parser.template get<RealType> ( pesopt::strprintf( "Force.RadiusBall%d", load ).c_str() );
             TangentVecType center; parser.template getFixSizeVector<RealType,TangentVecType> ( pesopt::strprintf( "Force.CenterBall%d", load ).c_str(), center );
             createConstantLoadInBall<ConfiguratorType>( conf, xAStorage, forceVec, radius, center ).assembleAdd( rhs_force );   
            }break;
            
            //in cube in normal direction
            case 11:{
             RealType factor = parser.template get<RealType> ( pesopt::strprintf( "Force.Factor%d", load ).c_str() );
             TangentVecType startPoint; parser.template getFixSizeVector<RealType,TangentVecType> ( pesopt::strprintf( "Force.StartRange%d", load ).c_str(), startPoint );
             TangentVecType endPoint; parser.template getFixSizeVector<RealType,TangentVecType> ( pesopt::strprintf( "Force.EndRange%d", load ).c_str(), endPoint );
             createNormalLoadInCube<ConfiguratorType>( conf, xAStorage, factor, startPoint, endPoint ).assembleAdd( rhs_force );   
            }break;

            default: 
                throw std::invalid_argument( pesopt::strprintf ( "Wrong NonlinearOptimizationType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
                break;
        };
    }

    //bc rhs_force 
    if( collapseBoundaryValues ){
        for( int i = 0; i < numGlobalDofs; ++i ){
          if ( boundaryMask[i] ){
            for( int comp=0; comp<3; ++comp )
              rhs_force[i + comp * numGlobalDofs] = 0.0;
          }
        } 
    }
}



//
// template< typename MatOpType >
// class createGravitationalLoad :
// public RefTriangleFENonlinVectorOpIntegrator <typename MatOpType::ConfiguratorType, createGravitationalLoad <MatOpType> >
// {      
//   
//   protected:
//     typedef typename MatOpType::ConfiguratorType ConfiguratorType;
//     typedef typename MatOpType::ConfiguratorTypePf ConfiguratorTypePf;
//     typedef typename ConfiguratorType::RealType RealType;
//     typedef typename ConfiguratorType::TangentVecType TangentVecType;
//     typedef typename ConfiguratorType::VectorType VectorType;
//   
//     const MatOpType &_matOpConf;
//     const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
//     DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
//     const RealType _fz;
// 
//   public:
//     createGravitationalLoad ( const MatOpType &matOpConf,
//                               const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
//                               const VectorType &pf,
//                               const RealType force_z ) :
//      RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, createGravitationalLoad <MatOpType> > (matOpConf._conf),
//      _matOpConf ( matOpConf ),
//      _xAStorage ( xAStorage ),
//      _pf ( matOpConf._confpf, pf ),
//      _fz ( force_z ) { }
//       
//     void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, TangentVecType &NL) const {   
//       NL[0] = 0.; NL[1] = 0.; NL[2] = _fz;
//       NL *= _matOpConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ) ) * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
//     }
// };

#endif
