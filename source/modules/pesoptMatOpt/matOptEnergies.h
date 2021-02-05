#ifndef __MATOPTENERGIES_H
#define __MATOPTENERGIES_H

#include <pesopt_IO.h>
#include <energyDefines.h>

#include <pesopt_fe.h>

#include "matOptDefines.h" //for approxCharFct, doubleWell


    
//==========================================================================================================================
//! computes \f$ volume(m) = \int_\Omega \chi(m) \, dx \f$
//==========================================================================================================================
template <typename ConfiguratorType>
class PfOp_Volume :
public FEIntegrator < ConfiguratorType, PfOp_Volume<ConfiguratorType> >{
  
protected :
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  const ConfiguratorType &_conf;
  PhaseFieldFunctions<RealType> _pfFcts;
  mutable FEScalarFunctionEvaluator<ConfiguratorType> *_PfPtr;
  
  public:
    PfOp_Volume ( const ConfiguratorType &conf ) : 
     FEIntegrator<ConfiguratorType, PfOp_Volume<ConfiguratorType> > ( conf ),
     _conf ( conf ), 
     _PfPtr ( NULL ) {}
     
  ~PfOp_Volume() {
     delete _PfPtr;
  };
  
  void apply( const VectorType & v, RealType& energy ) const {   
      energy = 0.;
      delete _PfPtr;
     _PfPtr = new FEScalarFunctionEvaluator<ConfiguratorType> ( _conf, v );
      FEIntegrator< ConfiguratorType, PfOp_Volume<ConfiguratorType>  >::assembleAdd( energy );
  }

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
     return _pfFcts.approxCharFct_vol( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

//==========================================================================================================================
//! computes \f$  D Volume(v) (\hat v) = \int_\Omega \chi'(v) (\hat v) \, dx \f$
//==========================================================================================================================
template <typename ConfiguratorType>
class PfOp_VolumeDerivative 
: public FEIntegratorVec< ConfiguratorType, PfOp_VolumeDerivative<ConfiguratorType> > {
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    
    const ConfiguratorType & _conf;
    PhaseFieldFunctions<RealType> _pfFcts;
    mutable FEScalarFunctionEvaluator<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_VolumeDerivative ( const ConfiguratorType & conf ) 
  : FEIntegratorVec< ConfiguratorType, PfOp_VolumeDerivative<ConfiguratorType>  > ( conf ),
    _conf ( conf ), _PfPtr ( NULL ) {}
     
  ~PfOp_VolumeDerivative() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, VectorType& Deriv ) const {   
      Deriv.setZero();
      delete _PfPtr;
     _PfPtr = new FEScalarFunctionEvaluator<ConfiguratorType> ( _conf, v );
      FEIntegratorVec< ConfiguratorType, PfOp_VolumeDerivative<ConfiguratorType> >::assembleAdd( Deriv );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      return _pfFcts.approxCharFct_vol_Derivative ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
  }
};


//! \f$ D^2 Volume(v) (hat v) (hat w) = \int \chi''(v) (\hat v)(\hat w) \f$
template <typename ConfiguratorType>
class PfOp_VolumeSecondDerivative 
: public FEWeightedMassIntegrator< ConfiguratorType, PfOp_VolumeSecondDerivative<ConfiguratorType> > {
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const ConfiguratorType &_conf;
    PhaseFieldFunctions<RealType> _pfFcts;
    mutable FEScalarFunctionEvaluator<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_VolumeSecondDerivative ( const ConfiguratorType &conf ) 
  : FEWeightedMassIntegrator< ConfiguratorType, PfOp_VolumeSecondDerivative<ConfiguratorType>  > ( conf ),
    _conf ( conf ), 
    _PfPtr ( NULL ) {}
     
  ~PfOp_VolumeSecondDerivative() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, SparseMatrixType& Hessian ) const {   
      Hessian.setZero();
      delete _PfPtr;
     _PfPtr = new FEScalarFunctionEvaluator<ConfiguratorType> ( _conf, v );
      FEWeightedMassIntegrator< ConfiguratorType, PfOp_VolumeSecondDerivative<ConfiguratorType> >::assemble( Hessian );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      return _pfFcts.approxCharFct_vol_SecondDerivative ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
  }
};


//! \f$ Volume(v) = \int_\Omega \chi(v) \, d x \f$ and derivatives
template <typename ConfiguratorType>
class VolumeOp 
: public pesopt::NonlinearEnergyOp< typename ConfiguratorType::DTContainer  >
{
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  
protected:
    const ConfiguratorType &_conf;
    const FEBoundaryConditionHandler<ConfiguratorType> &_bdryHandlerMaterial;
    const int _numDofs;

public:
  VolumeOp ( const ConfiguratorType &conf, 
             const FEBoundaryConditionHandler<ConfiguratorType> &bdryHandlerMaterial, 
             const int numDofs )  :
    _conf ( conf ), _bdryHandlerMaterial ( bdryHandlerMaterial ), _numDofs ( numDofs ) { } 
      
    void evaluateEnergy( const VectorType & v, RealType& energy ) const {
         energy = 0;
        VectorType vCollabsedAndExtended( v );
        _bdryHandlerMaterial.collabseVector( vCollabsedAndExtended );
        _bdryHandlerMaterial.extendVector( vCollabsedAndExtended );
         PfOp_Volume<ConfiguratorType> ( _conf ).apply( vCollabsedAndExtended, energy );
    }
    
    void evaluateJacobian( const VectorType & v, VectorType& Deriv ) const {
        VectorType vCollabsedAndExtended( v );
        _bdryHandlerMaterial.collabseVector( vCollabsedAndExtended );
        _bdryHandlerMaterial.extendVector( vCollabsedAndExtended );
         PfOp_VolumeDerivative<ConfiguratorType> ( _conf ).apply( vCollabsedAndExtended, Deriv );
         _bdryHandlerMaterial.collabseVectorAdditive( Deriv );
    }
    
    void evaluateHessian( const VectorType & /*v*/, SparseMatrixType& /*Hessian*/ ) const {
         throw std::invalid_argument( pesopt::strprintf ( "Hessian is so far not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
    void evaluateTripletListHessian( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const {
        throw std::invalid_argument( pesopt::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
    void evaluateTripletListHessianSym( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const {
        throw std::invalid_argument( pesopt::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
    const int getNumDofs ( ) const { return _numDofs; }
      
    //TODO getEnergyInfo
};



template <typename ConfiguratorType>
class VolumeConstraint 
: public pesopt::NonlinearConstraintOps< typename ConfiguratorType::DTContainer  >
{
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::DTContainer DataTypeContainer;
  typedef pesopt::BoostParser ParameterParserType;
  
  protected:
        const ConfiguratorType &_conf;
        const FEBoundaryConditionHandler<ConfiguratorType> &_bdryHandlerMaterial;
        const RealType _lowerBound, _upperBound;
        mutable RealType _lastVolume;

  public:
  VolumeConstraint ( const ConfiguratorType &conf,
                               const FEBoundaryConditionHandler<ConfiguratorType> &bdryHandlerMaterial,
                               const RealType lowerBound, const RealType upperBound ) 
  : pesopt::NonlinearConstraintOps< typename ConfiguratorType::DTContainer  > ( 1 ), 
    _conf ( conf ),
    _bdryHandlerMaterial ( bdryHandlerMaterial ),
    _lowerBound ( lowerBound ), _upperBound ( upperBound ) { } 
      
    void evaluateEnergy( const int numConstraint, const VectorType & v, RealType& energy ) const override{
         energy = 0;
        VectorType vCollabsedAndExtended( v );
        _bdryHandlerMaterial.collabseVector( vCollabsedAndExtended );
        _bdryHandlerMaterial.extendVector( vCollabsedAndExtended );
         PfOp_Volume<ConfiguratorType> ( _conf ).apply( vCollabsedAndExtended, _lastVolume );
         energy = _lastVolume;
    }
    
    void evaluateJacobian( const int numConstraint, const VectorType & v, VectorType& Deriv ) const override{
        VectorType vCollabsedAndExtended( v );
        _bdryHandlerMaterial.collabseVector( vCollabsedAndExtended );
        _bdryHandlerMaterial.extendVector( vCollabsedAndExtended );
         PfOp_VolumeDerivative<ConfiguratorType> ( _conf ).apply( vCollabsedAndExtended, Deriv );
         _bdryHandlerMaterial.collabseVectorAdditive( Deriv );
    }
    
    void evaluateHessian( const int numConstraint, const VectorType & /*v*/, SparseMatrixType& /*Hessian*/ ) const override{
         throw std::invalid_argument( pesopt::strprintf ( "Hessian is so far not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }

    const RealType getLowerBound( const int numConstraint) const override { return _lowerBound; }
    const RealType getUpperBound( const int numConstraint) const override { return _upperBound; }
    
    void getConstraintInfo( const VectorType & v, ParameterParserType & energyInfo ) const override{
        RealType energy=0.;  this->evaluateEnergy( 0, v, energy );
        energyInfo.set ( "VolumeConstraint.volume", _lastVolume );
   }
      
};


 
//==========================================================================================================================
//! \f$ barycenter_i(m) = (1/int_domain \chi(m) ) \int_domain chi(m) x_i  -  c_i \f$
//==========================================================================================================================
template <typename ConfiguratorType>
class PfOp_Barycenter :
public FEIntegrator < ConfiguratorType, PfOp_Barycenter<ConfiguratorType> >{
  
protected :
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::PointType PointType;
  typedef typename ConfiguratorType::QuadRuleType QuadRuleType;
  
  const ConfiguratorType & _conf;
  PhaseFieldFunctions<RealType> _pfFcts;
  const int _direction; 
  const RealType _c;
  mutable FEScalarFunctionEvaluator<ConfiguratorType> *_PfPtr;
  QuadRuleType _quadRule;
  
  public:
    PfOp_Barycenter ( const ConfiguratorType & conf, const int direction, const RealType c ) : 
     FEIntegrator<ConfiguratorType, PfOp_Barycenter<ConfiguratorType> > ( conf ),
     _conf ( conf ), _direction ( direction ), _c ( c ),
     _PfPtr ( NULL ) {}
     
  ~PfOp_Barycenter() {
     delete _PfPtr;
  };
  
  void apply( const VectorType & v, RealType& energy ) const {   
      energy = 0.;
      delete _PfPtr;
     _PfPtr = new FEScalarFunctionEvaluator<ConfiguratorType> ( _conf, v );
      FEIntegrator< ConfiguratorType, PfOp_Barycenter<ConfiguratorType>  >::assembleAdd( energy );
  }

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
     PointType GlobalCoord; El.getGlobalCoord ( _quadRule.getRefCoord( QuadPoint ), GlobalCoord );
     return _pfFcts.approxCharFct_vol( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) * ( GlobalCoord[_direction] - _c);
  }
};


template <typename ConfiguratorType>
class PfOp_BarycenterDerivative 
: public FEIntegratorVec< ConfiguratorType, PfOp_BarycenterDerivative<ConfiguratorType> > {
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::PointType PointType;
    typedef typename ConfiguratorType::QuadRuleType QuadRuleType;
    
    const ConfiguratorType & _conf;
    PhaseFieldFunctions<RealType> _pfFcts;
    const int _direction; 
    const RealType _c;
    mutable FEScalarFunctionEvaluator<ConfiguratorType> *_PfPtr;
    QuadRuleType _quadRule;
    
  public:
  PfOp_BarycenterDerivative ( const ConfiguratorType & conf, const int direction, const RealType c  ) 
  : FEIntegratorVec< ConfiguratorType, PfOp_BarycenterDerivative<ConfiguratorType>  > ( conf ),
    _conf ( conf ), _direction ( direction ), _c ( c ), _PfPtr ( NULL ) {}
     
  ~PfOp_BarycenterDerivative() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, VectorType& Deriv ) const {   
      Deriv.setZero();
      delete _PfPtr;
     _PfPtr = new FEScalarFunctionEvaluator<ConfiguratorType> ( _conf, v );
      FEIntegratorVec< ConfiguratorType, PfOp_BarycenterDerivative<ConfiguratorType> >::assembleAdd( Deriv );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      PointType GlobalCoord; El.getGlobalCoord ( _quadRule.getRefCoord( QuadPoint ), GlobalCoord );
      return _pfFcts.approxCharFct_vol_Derivative ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) * ( GlobalCoord[_direction] - _c);
  }
};


// D^2 Area(v) (hat v) (dhat v) = \int chi''(v) (hat v)(dhat v)
template <typename ConfiguratorType>
class PfOp_BarycenterSecondDerivative 
: public FEWeightedMassIntegrator< ConfiguratorType, PfOp_BarycenterSecondDerivative<ConfiguratorType> > {
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    typedef typename ConfiguratorType::PointType PointType;
    
    const ConfiguratorType & _conf;
    PhaseFieldFunctions<RealType> _pfFcts;
    const int _direction; 
    const RealType _c;
    mutable FEScalarFunctionEvaluator<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_BarycenterSecondDerivative ( const ConfiguratorType & conf, const int direction, const RealType c  ) 
  : FEWeightedMassIntegrator< ConfiguratorType, PfOp_BarycenterSecondDerivative<ConfiguratorType>  > ( conf ),
    _conf ( conf ), _direction ( direction ), _c ( c ),  _PfPtr ( NULL ) {}
     
  ~PfOp_BarycenterSecondDerivative() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, SparseMatrixType& Hessian ) const {   
      Hessian.setZero();
      delete _PfPtr;
     _PfPtr = new FEScalarFunctionEvaluator<ConfiguratorType> ( _conf, v );
      FEWeightedMassIntegrator< ConfiguratorType, PfOp_BarycenterSecondDerivative<ConfiguratorType> >::assemble( Hessian );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      PointType GlobalCoord; _conf.getGlobalCoords ( El, QuadPoint, GlobalCoord );
      return _pfFcts.approxCharFct_vol_SecondDerivative ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) * ( GlobalCoord[_direction] - _c);
  }
};


template <typename ConfiguratorType>
class BaryCenterConstraint 
: public pesopt::NonlinearConstraintOps< typename ConfiguratorType::DTContainer  >
{
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::PointType PointType;
  typedef typename ConfiguratorType::DTContainer DataTypeContainer;
  typedef pesopt::BoostParser ParameterParserType;
  
  protected:
        const ConfiguratorType &_conf;
        const FEBoundaryConditionHandler<ConfiguratorType> & _bdryHandlerMaterial;
        const PointType &_c;
        const PointType &_lowerBounds, &_upperBounds;

  public:
  BaryCenterConstraint ( const ConfiguratorType &conf, 
                                   const FEBoundaryConditionHandler<ConfiguratorType> &bdryHandlerMaterial, 
                                   const PointType &c,
                                   const PointType &lowerBounds, const PointType &upperBounds
                                 ) 
  : pesopt::NonlinearConstraintOps< typename ConfiguratorType::DTContainer  > ( c.size() ),
  _conf ( conf ), _bdryHandlerMaterial ( bdryHandlerMaterial ), _c ( c ), _lowerBounds ( lowerBounds ), _upperBounds ( upperBounds )  { }
      
    void evaluateEnergy( const int numConstraint, const VectorType & v, RealType& energy ) const {
         energy = 0;
        VectorType vCollabsedAndExtended( v );
        _bdryHandlerMaterial.collabseVector( vCollabsedAndExtended );
        _bdryHandlerMaterial.extendVector( vCollabsedAndExtended );
         PfOp_Barycenter<ConfiguratorType> ( _conf, numConstraint, _c[numConstraint] ).apply( vCollabsedAndExtended, energy );
    }
    
    void evaluateJacobian( const int numConstraint, const VectorType & v, VectorType& Deriv ) const {
        VectorType vCollabsedAndExtended( v );
        _bdryHandlerMaterial.collabseVector( vCollabsedAndExtended );
        _bdryHandlerMaterial.extendVector( vCollabsedAndExtended );
         PfOp_BarycenterDerivative<ConfiguratorType> ( _conf, numConstraint, _c[numConstraint] ).apply( vCollabsedAndExtended, Deriv );
         _bdryHandlerMaterial.collabseVectorAdditive( Deriv );
    }
    
    void evaluateHessian( const int numConstraint, const VectorType & /*v*/, SparseMatrixType& /*Hessian*/ ) const {
         throw std::invalid_argument( pesopt::strprintf ( "Hessian is so far not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
    const PointType& getBarycenterPoint( ) const {return _c;};
      
    const RealType getLowerBound( const int numConstraint) const override { return _lowerBounds[numConstraint]; }
    const RealType getUpperBound( const int numConstraint) const override { return _upperBounds[numConstraint]; }
    
    void getConstraintInfo ( const VectorType &Arg, ParameterParserType & energyInfo) const {
       PointType barycenterDiff;
       for( int i=0; i<this->getNumConstraints(); ++i ) this->evaluateEnergy( i, Arg, barycenterDiff[i] );
       PointType barycenter ( this->getBarycenterPoint() + barycenterDiff );
       energyInfo.template setFixSizeVector<PointType>( "BarycenterConstraint.barycenter", barycenter, false ); 
    }
};





//======================================================================================================================================
//================================= ModicaMortola and Derivative ============================================================================
//======================================================================================================================================
//! class to compute perimeter_eps(v) 
// = 1/2 \int_S eps |nabla v|^2 + 1/eps _factor (v^2-1)^2 d x
// Dirichlet Part
template <typename ConfiguratorType>
class PfOp_ModicaMortolaPart1 
: public FEIntegrator<ConfiguratorType, PfOp_ModicaMortolaPart1<ConfiguratorType> >{

protected:  
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::RealVecChart RealVecChart;
    
    const ConfiguratorType & _conf;
    const RealType _eps;
    mutable FEScalarFunctionEvaluator<ConfiguratorType> *_PfPtr;

public:
  PfOp_ModicaMortolaPart1 ( const ConfiguratorType & conf, const RealType eps_area ) : 
     FEIntegrator<ConfiguratorType, PfOp_ModicaMortolaPart1<ConfiguratorType> >( conf ), 
     _conf ( conf ), 
     _eps ( eps_area ),
     _PfPtr ( NULL ) { }
     
  ~PfOp_ModicaMortolaPart1() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, RealType& energy ) const {   
      energy = 0.;
      delete _PfPtr;
     _PfPtr = new FEScalarFunctionEvaluator<ConfiguratorType> ( _conf, v );
      FEIntegrator< ConfiguratorType, PfOp_ModicaMortolaPart1<ConfiguratorType> >::assembleAdd( energy );
  }
  
  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const {
      RealVecChart Dv;
      _PfPtr->evaluateGradientAtQuadPoint (El, QuadPoint, Dv );
      RealType Part1 = ( Dv ).dot( Dv );                
      return 0.5 * _eps * Part1;
  }
};


template <typename ConfiguratorType>
class PfOp_ModicaMortolaPart2 
: public FEIntegrator<ConfiguratorType, PfOp_ModicaMortolaPart2<ConfiguratorType> >{

protected:  
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::RealVecChart RealVecChart;
    
    const ConfiguratorType & _conf;
    PhaseFieldFunctions<RealType> _pfFcts;
    const RealType _eps;
    mutable FEScalarFunctionEvaluator<ConfiguratorType> *_PfPtr;

public:
  PfOp_ModicaMortolaPart2 ( const ConfiguratorType & conf, const RealType eps_area ) : 
     FEIntegrator<ConfiguratorType, PfOp_ModicaMortolaPart2<ConfiguratorType> >( conf ), 
     _conf ( conf ),
     _eps ( eps_area ),
     _PfPtr ( NULL ) { }
     
  ~PfOp_ModicaMortolaPart2() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, RealType& energy ) const {   
      energy = 0.;
      delete _PfPtr;
     _PfPtr = new FEScalarFunctionEvaluator<ConfiguratorType> ( _conf, v );
      FEIntegrator< ConfiguratorType, PfOp_ModicaMortolaPart2<ConfiguratorType>  >::assembleAdd( energy );
  }
  
  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const {
      return 0.5 * ( _pfFcts.doubleWell( _PfPtr->evaluateAtQuadPoint ( El, QuadPoint ) ) / _eps );
  }
};


template <typename ConfiguratorType>
class PfOp_ModicaMortola {

  protected:
    
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    
    const ConfiguratorType & _conf;
    const RealType _eps_area;
  
public:
  PfOp_ModicaMortola( const ConfiguratorType & conf, const RealType eps_area ) : 
       _conf ( conf ), _eps_area( eps_area ) { }
  
  void apply( const VectorType& v, RealType& energy ) const {
    RealType temp = 0;
    PfOp_ModicaMortolaPart1<ConfiguratorType> ( _conf, _eps_area ).apply( v, temp );
    energy = temp;
    PfOp_ModicaMortolaPart2<ConfiguratorType> ( _conf, _eps_area ).apply( v, temp );
    energy += temp;
  }
  
};



//! class to compute \int \sqrt(eps nabla m nabla b_i
template<typename ConfiguratorType>
class PfOp_ModicaMortolaDerivativePart1 :
public FEDiffOpIntegratorVec <ConfiguratorType, PfOp_ModicaMortolaDerivativePart1 <ConfiguratorType> >
{      
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::RealVecChart RealVecChart;
  
    const ConfiguratorType & _conf;
    const RealType _eps;
    mutable FEScalarFunctionEvaluator<ConfiguratorType> *_PfPtr;

  public:
    PfOp_ModicaMortolaDerivativePart1 ( const ConfiguratorType & conf, const RealType eps_area ) :
     FEDiffOpIntegratorVec<ConfiguratorType, PfOp_ModicaMortolaDerivativePart1<ConfiguratorType> > ( conf ),
      _conf ( conf ), _eps ( eps_area ),   _PfPtr ( NULL ) { }
      
     ~PfOp_ModicaMortolaDerivativePart1() {
           delete _PfPtr;
     };
  
   void apply( const VectorType & v, VectorType& Deriv ) const {   
      Deriv.setZero();
      delete _PfPtr;
     _PfPtr = new FEScalarFunctionEvaluator<ConfiguratorType> ( _conf, v );
      FEDiffOpIntegratorVec <ConfiguratorType, PfOp_ModicaMortolaDerivativePart1 <ConfiguratorType> > ::assembleAdd( Deriv );
  }
      
    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, RealVecChart &NL) const {     
      _PfPtr->evaluateGradientAtQuadPoint( El, QuadPoint, NL );
      NL *= _eps; 
      
    }
}; 



//! doubleWell part
template <typename ConfiguratorType>
class PfOp_ModicaMortolaDerivativePart2 : 
public FEIntegratorVec< ConfiguratorType, PfOp_ModicaMortolaDerivativePart2<ConfiguratorType> > {

  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    const ConfiguratorType & _conf;
    PhaseFieldFunctions<RealType> _pfFcts;
    const RealType _eps;
    mutable FEScalarFunctionEvaluator<ConfiguratorType> *_PfPtr;

  public:
  PfOp_ModicaMortolaDerivativePart2 ( const ConfiguratorType & conf, const RealType eps_area ) 
  : FEIntegratorVec< ConfiguratorType, PfOp_ModicaMortolaDerivativePart2<ConfiguratorType>  > ( conf ),
    _conf ( conf ), _eps ( eps_area ), _PfPtr ( NULL ) { }
      
     ~PfOp_ModicaMortolaDerivativePart2() {
           delete _PfPtr;
     };
  
   void apply( const VectorType & v, VectorType& Deriv ) const {   
      Deriv.setZero();
      delete _PfPtr;
     _PfPtr = new FEScalarFunctionEvaluator<ConfiguratorType> ( _conf, v );
      FEIntegratorVec <ConfiguratorType, PfOp_ModicaMortolaDerivativePart2 <ConfiguratorType> > ::assembleAdd( Deriv );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El,  int QuadPoint ) const {   
      return 0.5 * _pfFcts.doubleWellDerivative( _PfPtr->evaluateAtQuadPoint ( El, QuadPoint ) ) / _eps;
  }
};


//! ModicaMortola Derivative
template <typename ConfiguratorType>
class PfOp_ModicaMortolaDerivative {

  protected:
    
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    
    const ConfiguratorType & _conf;
    const RealType _eps_area;
  
public:
  PfOp_ModicaMortolaDerivative( const ConfiguratorType & conf, const RealType eps_area ) : 
       _conf ( conf ),  _eps_area( eps_area ) { }
  
  void apply( const VectorType& v, VectorType& Deriv ) const {
    VectorType temp ( Deriv.size() );
    PfOp_ModicaMortolaDerivativePart1<ConfiguratorType> ( _conf, _eps_area ).apply( v, temp );
    Deriv = temp;
    PfOp_ModicaMortolaDerivativePart2<ConfiguratorType> ( _conf, _eps_area ).apply( v, temp );
    Deriv += temp;
  }
  
};





/*
// Dirichlet-Part
template <typename ConfiguratorType>
class PfOp_ModicaMortolaSecondDerivativePart1 
: public UnitTriangleFELinAsymMatrixWeightedStiffIntegrator< ConfiguratorType, PfOp_ModicaMortolaSecondDerivativePart1<ConfiguratorType> > {
  
  protected:
    
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::Matrix22 Matrix22;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const ConfiguratorType & _conf;
    
    const RealType _eps_area;
    mutable FEScalarFunctionEvaluator<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_ModicaMortolaSecondDerivativePart1 ( const ConfiguratorType & conf, const RealType eps_area ) 
  : UnitTriangleFELinAsymMatrixWeightedStiffIntegrator< ConfiguratorType, PfOp_ModicaMortolaSecondDerivativePart1<ConfiguratorType>  > ( conf ),
    _conf ( conf ), 
    _eps_area ( eps_area ),
    _PfPtr ( NULL ) {}
     
  ~PfOp_ModicaMortolaSecondDerivativePart1() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, SparseMatrixType& Hessian ) const {   
      Hessian.setZero();
      delete _PfPtr;
     _PfPtr = new FEScalarFunctionEvaluator<ConfiguratorType> ( _conf, v );
      UnitTriangleFELinAsymMatrixWeightedStiffIntegrator< ConfiguratorType, PfOp_ModicaMortolaSecondDerivativePart1<ConfiguratorType> >::assemble( Hessian );
  }
  
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, Matrix22 &Matrix ) const {
    Matrix *= * _eps_area;
  }
  
};

// Double-Well-Part
template <typename ConfiguratorType>
class PfOp_ModicaMortolaSecondDerivativePart2 : 
public FEWeightedMassIntegrator< ConfiguratorType, PfOp_ModicaMortolaSecondDerivativePart2<ConfiguratorType> > {
// ,public pesopt::Op< ConfiguratorType::VectorType, ConfiguratorType::SparseMatrixType > {
  
  protected:
    
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const ConfiguratorType & _conf;
    PhaseFieldFunctions<RealType> _pfFcts;
    const RealType _eps_area;
    mutable FEScalarFunctionEvaluator<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_ModicaMortolaSecondDerivativePart2 ( const ConfiguratorType & conf,  const RealType eps_area ) 
  : FEWeightedMassIntegrator< ConfiguratorType, PfOp_ModicaMortolaSecondDerivativePart2<ConfiguratorType>  > ( conf ),
    _conf ( conf ), 
    
    _eps_area ( eps_area ),
    _PfPtr ( NULL ) {}
     
  ~PfOp_ModicaMortolaSecondDerivativePart2() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, SparseMatrixType& Hessian ) const {   
      Hessian.setZero();
      delete _PfPtr;
     _PfPtr = new FEScalarFunctionEvaluator<ConfiguratorType> ( _conf, v );
      FEWeightedMassIntegrator< ConfiguratorType, PfOp_ModicaMortolaSecondDerivativePart2<ConfiguratorType> >::assemble( Hessian );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      return 0.5 * _pfFcts.doubleWellSecondDerivative ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) / _eps_area;
  }
};

//! ModicaMortola Second Derivative
template <typename ConfiguratorType>
class PfOp_ModicaMortolaSecondDerivative {

  protected:
    
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const ConfiguratorType & _conf;
    
    const RealType _eps_area;
  
public:
  PfOp_ModicaMortolaSecondDerivative( const ConfiguratorType & conf,
                                             
                                             const RealType eps_area ) : 
       _conf ( conf ),         
       _eps_area( eps_area ) {}
  
  void apply( const VectorType& v, SparseMatrixType& Hessian ) const {
    SparseMatrixType temp ( Hessian.rows(), Hessian.cols() );
    PfOp_ModicaMortolaSecondDerivativePart1<ConfiguratorType> ( _conf,_eps_area ).apply( v, temp );
    Hessian = temp;
    PfOp_ModicaMortolaSecondDerivativePart2<ConfiguratorType> ( _conf, _eps_area ).apply( v, temp );
    Hessian += temp;
  }
  
};*/



#endif
