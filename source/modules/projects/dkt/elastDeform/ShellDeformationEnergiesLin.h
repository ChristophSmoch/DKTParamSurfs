#ifndef __LINEARDEFORMATIONENERGIESWITHMATERIAL_H
#define __LINEARDEFORMATIONENERGIESWITHMATERIAL_H

#include "ShellDeformationEnergiesInterfaces.h"


// ************************************************************************************
// file contains
// Interfaces: 
//             LinElastEnergy_EvaluationHelper
//             LinElastEnergyOpInterface 
// For Plates:
//             DirichletEnergy
//             DirichletEnergyWithApproxGrad 
//             LaplaceEnergy
// General Shells :
//             Dirichlet Energy
//             Laplace Energy 
//             KirchhoffLove Energy              
//*************************************************************************************




//*************************************************************************************
//*************************************************************************************
//************************************************************************************
//                                  Plates:             
//*************************************************************************************
//*************************************************************************************
//*************************************************************************************




//-----------------------------------------------------------------
//-----------------------------------------------------------------
//-----------------------------------------------------------------
/**  Dirichlet energy on plate
  *   E[u] = 1/2 \int_M chi(m) |nabla_M(u)|^2
  **/
//-----------------------------------------------------------------
//-----------------------------------------------------------------


template<typename MatOptConfType>
class Plate_DirichletEnergy_EvaluationHelper : public LinElastEnergy_EvaluationHelper<MatOptConfType,Plate_DirichletEnergy_EvaluationHelper<MatOptConfType>> {

protected: 

  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  typedef typename ConfiguratorType::ElementType ElementType;
  
public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstOrder;
  
protected:
  
  const MatOptConfType &_matOptConf;
  const ConfiguratorType &_conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  mutable RealType _factorMembraneEnergy, _factorBendingEnergy;
  
public:
    Plate_DirichletEnergy_EvaluationHelper ( const MatOptConfType &matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage ) :
    LinElastEnergy_EvaluationHelper<MatOptConfType,Plate_DirichletEnergy_EvaluationHelper<MatOptConfType>> ( matOpConf ),
    _matOptConf ( matOpConf ), _conf ( matOpConf._conf ),
    _xAStorage ( xAStorage ),     
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
    _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ), _factorBendingEnergy( 0. ) {}
   
   RealType getFactorMembraneEnergy (  ) const { return _factorMembraneEnergy; } 
   RealType getFactorBendingEnergy (  ) const { return _factorBendingEnergy; }
   
   
  RealType evaluateMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const Matrix32 &dX_u, const Matrix32 &dX_v, const RealType pf ) const {
      RealType aux = pesopt::ddProd<RealType, Matrix32> ( dX_u, dX_v );
      return _factorMembraneEnergy * _matOptConf.approxCharFct_material_Derivative( pf, _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus()) * aux;
  }
   
  RealType evaluateMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const VectorType &u, const VectorType &v, const RealType pf ) const {
      Matrix32 dX_u, dX_v;
      DKTFEVectorFunctionEvaluator<ConfiguratorType> uDFD ( _conf, u, 3 ), vDFD( _conf, v, 3 );
      uDFD.evaluateGradientAtQuadPoint( El, QuadPoint, dX_u );
      vDFD.evaluateGradientAtQuadPoint( El, QuadPoint, dX_v );
      return this->evaluateMixedSecondDerivativeAtQuadPoint( El, QuadPoint, dX_u, dX_v, pf );
  }

  RealType evaluateMembraneStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> uDFD ( _conf, displacement, 3 );
      RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluate( El, RefCoord ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      Matrix32 dX_u; uDFD.evaluateGradient( El, RefCoord, dX_u );
      RealType aux = pesopt::ddProd<RealType,Matrix32> ( dX_u, dX_u );
      return _factorMembraneEnergy * aux * chi;
  }
  
  RealType evaluateMembraneStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> uDFD ( _matOptConf._conf, displacement, 3 );
      RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      Matrix32 dX_u; uDFD.evaluateGradientAtQuadPoint( El, QuadPoint, dX_u );
      RealType aux = pesopt::ddProd<RealType,Matrix32> ( dX_u, dX_u );
      return _factorMembraneEnergy * aux * chi;
  }
      
  RealType evaluateBendingStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{ return 0.;}
  RealType evaluateBendingStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{ return 0.;}

};






template <typename MatOptConf>
class Plate_DirichletEnergy : 
public RefTriangleIntegrator < typename MatOptConf::ConfiguratorType, Plate_DirichletEnergy<MatOptConf> > {
  
protected :
  typedef typename MatOptConf::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConf::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
  
public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstOrder;
protected:
  
  const MatOptConf & _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  DKTFEVectorFunctionEvaluator<ConfiguratorType> _dispDFD;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  mutable RealType _factorMembraneEnergy;
  
  public:
    Plate_DirichletEnergy ( const MatOptConf & matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, const VectorType &disp, const VectorType &pf ) : 
     RefTriangleIntegrator<ConfiguratorType, Plate_DirichletEnergy<MatOptConf> > ( matOpConf._conf ),
     _matOptConf ( matOpConf ), _xAStorage ( xAStorage ),  _dispDFD ( matOpConf._conf, disp, 3 ), _pf( matOpConf._confpf, pf ),
     _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ) {}
     
  ~Plate_DirichletEnergy() { };

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
    RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
     Matrix32 GradDisp; _dispDFD.evaluateGradientAtQuadPoint( El, QuadPoint, GradDisp );
     RealType aux = pesopt::ddProd<RealType,Matrix32> ( GradDisp, GradDisp );
     return 0.5 * aux * _factorMembraneEnergy * chi;  // * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) (=1);
  }
};





template<typename MatOptConf>
class Plate_DirichletEnergy_SubHessian : 
public RefTriangleFELinAsymMatrixWeightedStiffIntegrator<typename MatOptConf::ConfiguratorType, Plate_DirichletEnergy_SubHessian <MatOptConf> > {
protected:
  typedef typename MatOptConf::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConf::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstOrder;
protected:
  const MatOptConf & _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  mutable RealType _factorMembraneEnergy, _factor;
  
public:
    Plate_DirichletEnergy_SubHessian(const MatOptConf &matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, const VectorType &pf) : 
    RefTriangleFELinAsymMatrixWeightedStiffIntegrator< ConfiguratorType, Plate_DirichletEnergy_SubHessian<MatOptConf>>(matOpConf._conf),
     _matOptConf ( matOpConf ), _xAStorage ( xAStorage ), _pf( matOpConf._confpf, pf ),
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
    _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ) {}

    inline void getCoeffMatrix(const typename ConfiguratorType::ElementType &El,  const int &QuadPoint, Matrix22 &matrix) const {
        RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
        DomVecType tmp;
        matrix.setZero(); 
        for( int i=0; i<tmp.size(); ++i ) matrix(i,i) = 1.;
        matrix *= _factorMembraneEnergy * chi;
    }
};




template <typename MatOptConfType>
class Plate_DirichletEnergyOp : public LinElastEnergyOpInterface<MatOptConfType, Plate_DirichletEnergyOp<MatOptConfType> > {
     
 protected: 
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType; 
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::MaskType MaskType;
  
 public:
  typedef Plate_DirichletEnergy_EvaluationHelper<MatOptConfType> EvaluationHelper;
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstOrder;
 protected:
  const MatOptConfType &_matOptConf; 
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  const VectorType & _pf;
     
   public:
    Plate_DirichletEnergyOp ( const MatOptConfType &matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, const VectorType &pf ) : 
    LinElastEnergyOpInterface<MatOptConfType, Plate_DirichletEnergyOp<MatOptConfType> > ( matOpConf ),
    _matOptConf ( matOpConf ), _xAStorage ( xAStorage ), _pf ( pf ) {}
  
  void evaluateElasticEnergies( const VectorType &disp, RealType &membraneEnergy, RealType &bendingEnergy, RealType &totalEnergy ) const {
      Plate_DirichletEnergy<MatOptConfType> energyOp( _matOptConf, _xAStorage, disp, _pf );
      membraneEnergy = 0.; energyOp.assembleAdd( membraneEnergy );
      bendingEnergy = 0.;
      totalEnergy = membraneEnergy + bendingEnergy;
  }
  
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    std::vector<TripletType>  tripletListPart;
    Plate_DirichletEnergy_SubHessian<MatOptConfType> ( _matOptConf, _xAStorage, _pf ).assembleTripletListDirichlet( tripletListPart, boundaryMask, Factor );

    tripletListMasked.reserve( 3 * ( tripletListPart.size() ) );
    for( int k=0; k<3; ++k ){
        const int offset = k * _matOptConf._conf.getNumGlobalDofs();
        for( int iter=0; iter<tripletListPart.size(); ++iter )
            tripletListMasked.push_back( TripletType( tripletListPart[iter].row() + offset, tripletListPart[iter].col() + offset, tripletListPart[iter].value() ) );
    }
  }

};





//-----------------------------------------------------------------
//-----------------------------------------------------------------
//-----------------------------------------------------------------
/**  Dirichlet energy on plate with approx gradient
  *   E[u] = 1/2 \int_M chi(m) | theta(u) |^2
  **/
//-----------------------------------------------------------------
//-----------------------------------------------------------------


template<typename MatOptConfType>
class Plate_DirichletEnergyApproxGrad_EvaluationHelper : public LinElastEnergy_EvaluationHelper<MatOptConfType,Plate_DirichletEnergyApproxGrad_EvaluationHelper<MatOptConfType>> {

protected: 

  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
  const MatOptConfType &_matOptConf;
  const ConfiguratorType &_conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  mutable RealType _factorMembraneEnergy, _factorBendingEnergy;
  
public:
    Plate_DirichletEnergyApproxGrad_EvaluationHelper ( const MatOptConfType &matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage ) :
    LinElastEnergy_EvaluationHelper<MatOptConfType,Plate_DirichletEnergyApproxGrad_EvaluationHelper<MatOptConfType>> ( matOpConf ),
    _matOptConf ( matOpConf ), _conf ( matOpConf._conf ),
    _xAStorage ( xAStorage ),     
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
    _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ), _factorBendingEnergy( 0. ) {}
   
   RealType getFactorMembraneEnergy (  ) const { return _factorMembraneEnergy; } 
   RealType getFactorBendingEnergy (  ) const { return _factorBendingEnergy; }
   
  RealType evaluateMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const Matrix32 &dX_u, const Matrix32 &dX_v, const RealType pf ) const {
      RealType aux = pesopt::ddProd<RealType, Matrix32> ( dX_u, dX_v );
      return _factorMembraneEnergy * _matOptConf.approxCharFct_material_Derivative( pf, _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus()) * aux;
  }
   
  RealType evaluateMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const VectorType &u, const VectorType &v, const RealType pf ) const {
      Matrix32 dX_u, dX_v;
      DKTFEVectorFunctionEvaluator<ConfiguratorType> uDFD ( _conf, u, 3 ), vDFD( _conf, v, 3 );
      uDFD.evaluateApproxGradientAtQuadPoint( El, QuadPoint, dX_u );
      vDFD.evaluateApproxGradientAtQuadPoint( El, QuadPoint, dX_v );
      return this->evaluateMixedSecondDerivativeAtQuadPoint( El, QuadPoint, dX_u, dX_v, pf );
  }

  RealType evaluateMembraneStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> uDFD ( _conf, displacement, 3 );
      RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluate( El, RefCoord ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      Matrix32 dX_u; uDFD.evaluateApproxGradient( El, RefCoord, dX_u );
      RealType aux = pesopt::ddProd<RealType,Matrix32> ( dX_u, dX_u );
      return _factorMembraneEnergy * aux * chi;
  }
  
  RealType evaluateMembraneStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> uDFD ( _matOptConf._conf, displacement, 3 );
      RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      Matrix32 dX_u; uDFD.evaluateApproxGradientAtQuadPoint( El, QuadPoint, dX_u );
      RealType aux = pesopt::ddProd<RealType,Matrix32> ( dX_u, dX_u );
      return _factorMembraneEnergy * aux * chi;
  }
      
  RealType evaluateBendingStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{ return 0.;}
  RealType evaluateBendingStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{ return 0.;}
  
protected:   

};


template <typename MatOptConf>
class Plate_DirichletEnergyApproxGrad : 
public RefTriangleIntegrator < typename MatOptConf::ConfiguratorType, Plate_DirichletEnergyApproxGrad<MatOptConf> > {
  
protected :
  typedef typename MatOptConf::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConf::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
  const MatOptConf & _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  DKTFEVectorFunctionEvaluator<ConfiguratorType> _dispDFD;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  mutable RealType _factorMembraneEnergy;
  
  public:
    Plate_DirichletEnergyApproxGrad ( const MatOptConf & matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, const VectorType &disp, const VectorType &pf ) : 
     RefTriangleIntegrator<ConfiguratorType, Plate_DirichletEnergyApproxGrad<MatOptConf> > ( matOpConf._conf ),
     _matOptConf ( matOpConf ), _xAStorage ( xAStorage ),  _dispDFD ( matOpConf._conf, disp, 3 ), _pf( matOpConf._confpf, pf ),
     _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ){}
     
  ~Plate_DirichletEnergyApproxGrad() { };

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
    RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
    Matrix32 GradDisp; _dispDFD.evaluateApproxGradientAtQuadPoint( El, QuadPoint, GradDisp );
    RealType aux = pesopt::ddProd<RealType,Matrix32> ( GradDisp, GradDisp );
    return 0.5 * aux * _factorMembraneEnergy * chi;
  }
};





template<typename MatOptConf>
class Plate_DirichletEnergyApproxGrad_Hessian : 
public RefTriangleFELinAsymMatrixWeightedApproxStiffIntegrator<typename MatOptConf::ConfiguratorType, Plate_DirichletEnergyApproxGrad_Hessian <MatOptConf> > {
protected:
  typedef typename MatOptConf::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConf::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
  
public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
  const MatOptConf & _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
   
  const Material<RealType> & _HardMaterial;
  const Material<RealType> & _SoftMaterial;
  mutable RealType _factorMembraneEnergy;
  
public:
    Plate_DirichletEnergyApproxGrad_Hessian(const MatOptConf &matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, const VectorType &pf) : 
    RefTriangleFELinAsymMatrixWeightedApproxStiffIntegrator< ConfiguratorType, Plate_DirichletEnergyApproxGrad_Hessian<MatOptConf>>(matOpConf._conf),
     _matOptConf ( matOpConf ), _xAStorage ( xAStorage ), _pf( matOpConf._confpf, pf ),
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
    _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy )  {}

    inline void getCoeffMatrix(const typename ConfiguratorType::ElementType &El,  const int &QuadPoint, Matrix22 &matrix) const {
        RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
        DomVecType tmp;
        matrix.setZero();  for( int i=0; i<tmp.size(); ++i ) matrix(i,i) = 1.;
        matrix *= _factorMembraneEnergy * chi;
    }
};




template <typename MatOptConfType>
class Plate_DirichletEnergyApproxGradOp : public LinElastEnergyOpInterface<MatOptConfType, Plate_DirichletEnergyApproxGradOp<MatOptConfType> > {
     
 protected: 
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType; 
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::MaskType MaskType;
  
 public:
  typedef Plate_DirichletEnergyApproxGrad_EvaluationHelper<MatOptConfType> EvaluationHelper;
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
  const MatOptConfType &_matOptConf; 
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  const VectorType & _pf;
     
   public:
    Plate_DirichletEnergyApproxGradOp ( const MatOptConfType &matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, const VectorType &pf ) : 
    LinElastEnergyOpInterface<MatOptConfType, Plate_DirichletEnergyApproxGradOp<MatOptConfType> > ( matOpConf ),
    _matOptConf ( matOpConf ), _xAStorage ( xAStorage ), _pf ( pf ) {}
  
  void evaluateElasticEnergies( const VectorType &disp, RealType &membraneEnergy, RealType &bendingEnergy, RealType &totalEnergy ) const {
      Plate_DirichletEnergyApproxGrad<MatOptConfType> energyOp( _matOptConf, _xAStorage, disp, _pf );
      membraneEnergy = 0.; energyOp.assembleAdd( membraneEnergy );
      bendingEnergy = 0.;
      totalEnergy = membraneEnergy + bendingEnergy;
  }
  
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    std::vector<TripletType>  tripletListPart;
    Plate_DirichletEnergyApproxGrad_Hessian<MatOptConfType> ( _matOptConf, _xAStorage, _pf ).assembleTripletListDirichlet( tripletListPart, boundaryMask, Factor );

    tripletListMasked.reserve( 3 * ( tripletListPart.size() ) );
    for( int k=0; k<3; ++k ){
        const int offset = k * _matOptConf._conf.getNumGlobalDofs();
        for( int iter=0; iter<tripletListPart.size(); ++iter )
            tripletListMasked.push_back( TripletType( tripletListPart[iter].row() + offset, tripletListPart[iter].col() + offset, tripletListPart[iter].value() ) );
    }
  }

};





//-----------------------------------------------------------------
//-----------------------------------------------------------------
//-----------------------------------------------------------------
//   Plate Laplace energy and derivatives
// E = 1/24 \int_Plate |D^2 u|^2 = 1/24 \int_Plate |D^2_11 u|^2 + |D^2_22 u|^ + 0.5 |D^2_12 u + D^2_21 u|^2
//-----------------------------------------------------------------
//-----------------------------------------------------------------


template<typename MatOptConfType>
class Plate_LaplaceEnergy_EvaluationHelper : public LinElastEnergy_EvaluationHelper<MatOptConfType,Plate_LaplaceEnergy_EvaluationHelper<MatOptConfType>> {

protected: 

  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType; 
  typedef typename ConfiguratorType::LocalVectorType LocalVectorType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::ElementType ElementType; 
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
  const MatOptConfType &_matOptConf;
  const ConfiguratorType &_conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  mutable RealType _factorMembraneEnergy, _factorBendingEnergy;
  
public:
    Plate_LaplaceEnergy_EvaluationHelper ( const MatOptConfType &matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage ) :
    LinElastEnergy_EvaluationHelper<MatOptConfType,Plate_LaplaceEnergy_EvaluationHelper<MatOptConfType>> ( matOpConf ),
    _matOptConf ( matOpConf ), _conf ( matOpConf._conf ),
    _xAStorage ( xAStorage ),     
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
    _factorMembraneEnergy( 0. ), _factorBendingEnergy ( matOpConf._materialInfo._factorBendingEnergy ) {}
  
  RealType evaluateMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const Matrix33 &ddXAsVec_u, const Matrix33 &ddXAsVec_v, const RealType pf ) const {
      RealType aux = 0.;
      for( int comp = 0; comp < 3; ++comp ){
            Point3DType hessianComp_u; hessianComp_u = ddXAsVec_u.row(comp);
            Point3DType hessianComp_v; hessianComp_v = ddXAsVec_v.row(comp);
            aux += hessianComp_u(0) * hessianComp_v(0)  + hessianComp_u(1) * hessianComp_v(1) + 0.5 * hessianComp_u(2) * hessianComp_v(2);
      }
      return (1./12.) * _factorBendingEnergy * _matOptConf.approxCharFct_material_Derivative( pf, _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus()) * aux;
  }
    
  RealType evaluateMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const VectorType &u, const VectorType &v, const RealType pf ) const {
      Matrix33 ddX_u, ddX_v;
      DKTFEVectorFunctionEvaluator<ConfiguratorType> uDFD ( _conf, u, 3 ), vDFD( _conf, v, 3 );
      uDFD.evaluateApproxHessianAsVecAtQuadPoint( El, QuadPoint, ddX_u );
      vDFD.evaluateApproxHessianAsVecAtQuadPoint( El, QuadPoint, ddX_v );
      return this->evaluateMixedSecondDerivativeAtQuadPoint( El, QuadPoint, ddX_u, ddX_v, pf );
  }
  
  RealType getFactorMembraneEnergy (  ) const { return _factorMembraneEnergy; }  
  RealType getFactorBendingEnergy (  ) const { return _factorBendingEnergy; }
  
  RealType evaluateMembraneStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{ return 0.; }
  RealType evaluateMembraneStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{ return 0.;}
  
  RealType evaluateBendingStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( _matOptConf._conf, displacement, 3 );
      Matrix33 ddXAsVec_u; dispDFD.evaluateApproxHessianAsVec( El, RefCoord, ddXAsVec_u );
      RealType aux = 0.;
      for( int comp = 0; comp < 3; ++comp ){
            Point3DType hessianComp_u; hessianComp_u = ddXAsVec_u.row(comp);
            aux += pesopt::Sqr( hessianComp_u(0) ) + pesopt::Sqr( hessianComp_u(1) ) + 0.5 * pesopt::Sqr( hessianComp_u(2) );
      }
      RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluate( El, RefCoord ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return _factorBendingEnergy * (1./12.) * 0.5 * aux * chi;
  }
  
   RealType evaluateBendingStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( _matOptConf._conf, displacement, 3 );
      Matrix33 ddXAsVec_u; dispDFD.evaluateApproxHessianAsVecAtQuadPoint( El, QuadPoint, ddXAsVec_u );
      RealType aux = 0.;
      for( int comp = 0; comp < 3; ++comp ){
            Point3DType hessianComp_u; hessianComp_u = ddXAsVec_u.row(comp);
            aux += pesopt::Sqr( hessianComp_u(0) ) + pesopt::Sqr( hessianComp_u(1) ) + 0.5 * pesopt::Sqr( hessianComp_u(2) );
      }
      RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return _factorBendingEnergy * (1./12.) * 0.5 * aux * chi;
  }
  
  
};
    


template<typename MatOptConfType>
class Plate_LaplaceEnergy :
public RefTriangleIntegrator<typename MatOptConfType::ConfiguratorType, Plate_LaplaceEnergy<MatOptConfType> >,
public Plate_LaplaceEnergy_EvaluationHelper<MatOptConfType>
{
  protected:
       typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
       typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
  
       typedef Plate_LaplaceEnergy_EvaluationHelper<MatOptConfType> EvaluationHelper;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const MatOptConfType &_matOptConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _disp;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> &_HardMaterial, &_SoftMaterial;
       mutable RealType _factorBendingEnergy;
  public:
    Plate_LaplaceEnergy ( const MatOptConfType &matOpConf, 
                    const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                    const VectorType &displacement,
                    const VectorType &pf ) :
      RefTriangleIntegrator<ConfiguratorType, Plate_LaplaceEnergy<MatOptConfType>> (matOpConf._conf),
      Plate_LaplaceEnergy_EvaluationHelper<MatOptConfType>( matOpConf, xAStorage ), 
      _matOptConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorBendingEnergy ( matOpConf._materialInfo._factorBendingEnergy ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        Matrix33 HessianAsVec; _disp.evaluateApproxHessianAsVecAtQuadPoint( El, QuadPoint, HessianAsVec );
        RealType aux = 0.;
        for( int comp = 0; comp < 3; ++comp ){
            Point3DType hessianComp; hessianComp = HessianAsVec.row(comp);
            aux += pesopt::Sqr( hessianComp(0) ) + pesopt::Sqr( hessianComp(1) ) + 0.5 * pesopt::Sqr( hessianComp(2) );
        }
        RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
        return (1./12.) * 0.5 * aux * chi * _factorBendingEnergy;//* _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) (= 1);
    }
};


template<typename MatOptConfType>
class Plate_LaplaceEnergy_SubHessian :
public RefTriangleFELinMatrixDiff2AsVecIntegrator <typename MatOptConfType::ConfiguratorType, Plate_LaplaceEnergy_SubHessian <MatOptConfType> > {
   protected:
       typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
       typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;

public:
     static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
     const MatOptConfType &_matOptConf;   
     const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
     DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
     const Material<RealType> &_HardMaterial, &_SoftMaterial;
     mutable RealType _factorBendingEnergy;
 
   public:
     Plate_LaplaceEnergy_SubHessian ( const MatOptConfType &matOpConf, 
                                      const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, 
                                      const VectorType &pf ) :
      RefTriangleFELinMatrixDiff2AsVecIntegrator<ConfiguratorType, Plate_LaplaceEnergy_SubHessian<MatOptConfType>> (matOpConf._conf),
      _matOptConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorBendingEnergy ( matOpConf._materialInfo._factorBendingEnergy ) {}
 
     inline void getCoeffMatrix ( const ElementType &El, int QuadPoint, Matrix33 &Matrix ) const{
       const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
       Matrix.setZero(); Matrix(0,0) = 1.; Matrix(1,1) = 1.; Matrix(2,2) = 0.5;
       Matrix *= (1./12.) * _factorBendingEnergy * chi;
     }
 };

 
template <typename MatOptConfType>
class Plate_LaplaceEnergyOp : public LinElastEnergyOpInterface<MatOptConfType, Plate_LaplaceEnergyOp<MatOptConfType> >{
     
 protected: 
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType; 
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::MaskType MaskType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  
 public:
  typedef Plate_LaplaceEnergy_EvaluationHelper<MatOptConfType> EvaluationHelper;
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
  const MatOptConfType &_matOptConf; 
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  const VectorType & _pf;
     
 public:
  Plate_LaplaceEnergyOp ( const MatOptConfType &matOpConf, 
                          const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, 
                          const VectorType &pf ) :
    LinElastEnergyOpInterface<MatOptConfType, Plate_LaplaceEnergyOp<MatOptConfType> > ( matOpConf ),
    _matOptConf ( matOpConf ), _xAStorage ( xAStorage ), _pf ( pf ) {}
  
  void evaluateElasticEnergies( const VectorType &disp, RealType &membraneEnergy, RealType &bendingEnergy, RealType &totalEnergy ) const {
      Plate_LaplaceEnergy<MatOptConfType> energyOp( _matOptConf, _xAStorage, disp, _pf );
      membraneEnergy = 0.;
      bendingEnergy = 0.; energyOp.assembleAdd( bendingEnergy );
      totalEnergy = membraneEnergy + bendingEnergy;
  }
  
//   void evaluateGradient( const VectorType &disp, VectorType &Deriv ) const {
//     Deriv.setZero();
//     Plate_LaplaceEnergyGradient_Part1<MatOptConfType> ( _matOptConf, _xA, disp, _pf ).assembleAdd( Deriv );
//     Plate_LaplaceEnergyGradient_Part2<MatOptConfType> ( _matOptConf, _xA, disp, _pf ).assembleAdd( Deriv );
//   }
  
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    std::vector<TripletType> tripletListPart1;
    Plate_LaplaceEnergy_SubHessian<MatOptConfType> ( _matOptConf, _xAStorage, _pf ).assembleTripletListDirichlet( tripletListPart1, boundaryMask, Factor );

    tripletListMasked.reserve( 3 * tripletListPart1.size() );
    for( int comp=0; comp<3; ++comp ){
        const int offset = comp * _matOptConf._conf.getNumGlobalDofs();
        for( int iter=0; iter<tripletListPart1.size(); ++iter )
            tripletListMasked.push_back( TripletType( tripletListPart1[iter].row() + offset, tripletListPart1[iter].col() + offset, tripletListPart1[iter].value() ) );
    }
  }
     
};




template<typename MatOptConfType>
class Plate_LaplaceEnergyMixedSecondDerivative
: public RefTriangleFENonlinOpIntegrator< typename MatOptConfType::ConfiguratorTypePf, Plate_LaplaceEnergyMixedSecondDerivative<MatOptConfType> >,
  public Plate_LaplaceEnergy_EvaluationHelper<MatOptConfType> {

protected: 
  
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename MatOptConfType::RealType RealType;
  typedef typename ConfiguratorTypePf::VectorType VectorType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
//   typedef Plate_LaplaceEnergy_EvaluationHelper<MatOptConfType> EvaluationHelper;
  
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage;
  const VectorType &_disp;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const VectorType &_adjoint;
  
public:
  Plate_LaplaceEnergyMixedSecondDerivative( const MatOptConfType &matOptConf, 
                                            const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                                            const VectorType &disp,
                                            const VectorType &pf,
                                            const VectorType &adjoint  ) : 
     RefTriangleFENonlinOpIntegrator<ConfiguratorTypePf, Plate_LaplaceEnergyMixedSecondDerivative<MatOptConfType> > ( matOptConf._confpf ),
     Plate_LaplaceEnergy_EvaluationHelper<MatOptConfType> ( matOptConf, xAStorage ),
     _xAStorage(xAStorage),
     _disp ( disp ),
     _pf( matOptConf._confpf, pf ),
     _adjoint( adjoint ) {}
      

  RealType getNonlinearity ( const typename ConfiguratorTypePf::ElementType &El, int QuadPoint ) const {
    RealType NL = this->evaluateMixedSecondDerivativeAtQuadPoint( El, QuadPoint, _disp, _adjoint, _pf.evaluateAtQuadPoint( El, QuadPoint ) );
    NL *= _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    return NL;
  }

};












//*************************************************************************************
//*************************************************************************************
//************************************************************************************
//                                  General Shells:             
//*************************************************************************************
//*************************************************************************************
//*************************************************************************************





//-----------------------------------------------------------------
//-----------------------------------------------------------------
//-----------------------------------------------------------------
/**  Dirichlet energy
  *   E[u] = 1/2 \int_M chi(m) |nabla_M(u)|^2
  *        = 1/2 \int_chart sqrt(det(gA)) chi(m) |(nabla_M u ) \circ xA )|^2 
  *        = 1/2 \int_chart sqrt(det(gA)) chi(m) \sum_comp gA^{-1} Du_comp \cdot Du_comp
  **/
//-----------------------------------------------------------------
//-----------------------------------------------------------------


template<typename MatOptConfType>
class DirichletEnergy_EvaluationHelper : public LinElastEnergy_EvaluationHelper<MatOptConfType,DirichletEnergy_EvaluationHelper<MatOptConfType>> {

 protected: 
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::ElementType ElementType; 
  typedef typename ConfiguratorType::VectorType VectorType;
 public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstOrder;
 protected:
  const MatOptConfType &_matOptConf;
  const ConfiguratorType &_conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  mutable RealType _factorMembraneEnergy, _factorBendingEnergy;
  
 public:
    DirichletEnergy_EvaluationHelper ( const MatOptConfType &matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage ) :
    LinElastEnergy_EvaluationHelper<MatOptConfType,DirichletEnergy_EvaluationHelper<MatOptConfType>> ( matOpConf ),
    _matOptConf ( matOpConf ), _conf ( matOpConf._conf ),
    _xAStorage ( xAStorage ),     
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
    _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ), _factorBendingEnergy( 0. ) {}
   
   RealType getFactorMembraneEnergy (  ) const { return _factorMembraneEnergy; } 
   RealType getFactorBendingEnergy (  ) const { return _factorBendingEnergy; }
   
   
  RealType evaluateMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const VectorType &u, const VectorType &v, const RealType pf ) const {

      Matrix32 dX_u, dX_v;
      DKTFEVectorFunctionEvaluator<ConfiguratorType> uDFD ( _conf, u, 3 ), vDFD( _conf, v, 3 );
      uDFD.evaluateGradientAtQuadPoint( El, QuadPoint, dX_u );
      vDFD.evaluateGradientAtQuadPoint( El, QuadPoint, dX_v );
      
      Matrix33 grad_u, grad_v;
      Matrix32 DxAGAinv = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint) * _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
      this->evaluateGradient( DxAGAinv, dX_u, grad_u );
      this->evaluateGradient( DxAGAinv, dX_v, grad_v );
      
      RealType Dchi = _matOptConf.approxCharFct_material_Derivative( pf, _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus());
      
      return _factorMembraneEnergy * Dchi * pesopt::ddProd<RealType,Matrix33>( grad_u, grad_v );
  }

  RealType evaluateMembraneStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( _matOptConf._conf, displacement, 3 );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> xADFD ( _matOptConf._conf, xA, 3 );
      PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _pointwiseEvaluator;
      Matrix32 DxA; xADFD.evaluateGradient(El,RefCoord,DxA);
      Matrix22 gA; _pointwiseEvaluator.evaluateFirstFundamentalForm( DxA, gA );
      Matrix32 DxAGAinv = DxA * (gA.inverse());
      Matrix32 dX; dispDFD.evaluateGradient( El, RefCoord, dX );
      Matrix33 grad_u;
      this->evaluateGradient ( DxAGAinv, dX, grad_u );
      RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluate( El, RefCoord ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return _factorMembraneEnergy * pesopt::ddProd<RealType,Matrix33>(grad_u,grad_u) * chi;
  }
  
  RealType evaluateMembraneStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( _matOptConf._conf, displacement, 3 );
      Matrix32 dX; dispDFD.evaluateGradientAtQuadPoint( El, QuadPoint, dX );
      Matrix32 DxAGAinv = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint) * _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
      Matrix33 grad_u;
      this->evaluateGradient ( DxAGAinv, dX, grad_u );
      RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return _factorMembraneEnergy * pesopt::ddProd<RealType,Matrix33>(grad_u,grad_u) * chi;
  }
      
  RealType evaluateBendingStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{ return 0.;}
  RealType evaluateBendingStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{ return 0.;}
  
protected:   
  
  void evaluateGradient ( const Matrix32 &DxAGAinv, const Matrix32 &dX, Matrix33 &gradient ) const {
      gradient = DxAGAinv * dX.transpose();
  }

};


// Dirichlet Energy
template<typename MatOptConfType>
class DirichletEnergy :
public RefTriangleIntegrator<typename MatOptConfType::ConfiguratorType, DirichletEnergy<MatOptConfType> >,
public DirichletEnergy_EvaluationHelper<MatOptConfType>
{
  protected:
       typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
       typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::VectorType VectorType;
       typedef DirichletEnergy_EvaluationHelper<MatOptConfType> EvaluationHelper;
       
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstOrder;
 protected:
       const MatOptConfType &_matOptConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _disp;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> &_HardMaterial, &_SoftMaterial;
       mutable RealType _factorMembraneEnergy;
  public:
    DirichletEnergy ( const MatOptConfType &matOpConf, 
                    const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                    const VectorType &displacement,
                    const VectorType &pf ) :
      RefTriangleIntegrator<ConfiguratorType, DirichletEnergy<MatOptConfType>> (matOpConf._conf),
      DirichletEnergy_EvaluationHelper<MatOptConfType>( matOpConf, xAStorage ), 
      _matOptConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
      Matrix32 dX; _disp.evaluateGradientAtQuadPoint( El, QuadPoint, dX );
      Matrix32 DxAGAinv = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint) * _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
      
      Matrix33 gradient;
      this->evaluateGradient ( DxAGAinv, dX, gradient );
      RealType aux = pesopt::ddProd<RealType,Matrix33>(gradient,gradient);
      
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return 0.5 * aux * chi * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * _factorMembraneEnergy;
    }
};



//! \brief Assembles \f$ \sum_{EL} \int_T \sqrt{\det g} g^{-1} D v_i \cdot D v_j \f$
template<typename MatOptConfType>
class DirichletEnergy_Hessian : public RefTriangleFELinAsymMatrixWeightedStiffIntegrator<typename MatOptConfType::ConfiguratorType, DirichletEnergy_Hessian <MatOptConfType> >
{
 protected:
    typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
    typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::TripletType TripletType;
    typedef typename ConfiguratorType::MaskType MaskType;
    typedef typename ConfiguratorType::RealType   RealType;
    typedef typename ConfiguratorType::Matrix22   Matrix22;
    typedef typename ConfiguratorType::Matrix32   Matrix32;
 public:
    static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstOrder;
 protected:
    const MatOptConfType &_matOptConf;
    const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
    const DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
    const Material<RealType> &_HardMaterial, &_SoftMaterial;
    mutable RealType _factorMembraneEnergy;

 public:
    DirichletEnergy_Hessian(const MatOptConfType &matOpConf,  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> & xAStorage, const VectorType &pf ) : 
    RefTriangleFELinAsymMatrixWeightedStiffIntegrator< ConfiguratorType, DirichletEnergy_Hessian<MatOptConfType>>(matOpConf._conf),
    _matOptConf( matOpConf ),
    _xAStorage ( xAStorage ), _pf( matOpConf._confpf, pf),
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
    _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ) {}

    inline void getCoeffMatrix(const typename ConfiguratorType::ElementType &El,  const int &QuadPoint, Matrix22 &matrix) const {
        RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
        matrix  = _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
        matrix *= _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * chi * _factorMembraneEnergy;
    }
};


template <typename MatOptConfType>
class DirichletEnergyOp : public LinElastEnergyOpInterface<MatOptConfType, DirichletEnergyOp<MatOptConfType> > {
     
 protected: 
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType; 
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::MaskType MaskType;
 public:
  typedef DirichletEnergy_EvaluationHelper<MatOptConfType> EvaluationHelper;
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstOrder;
 protected:
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage; 
  const VectorType &_pf; 
 public:
  DirichletEnergyOp ( const MatOptConfType &matOpConf, 
                      const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, 
                      const VectorType &pf ) : 
    LinElastEnergyOpInterface<MatOptConfType, DirichletEnergyOp<MatOptConfType> > ( matOpConf ),
    _xAStorage( xAStorage ), _pf( pf ) {}
  
  void evaluateElasticEnergies( const VectorType &disp, RealType &membraneEnergy, RealType &bendingEnergy, RealType &totalEnergy ) const {
      DirichletEnergy<MatOptConfType> energyOp( this->_matOptConf, this->_xAStorage, disp, this->_pf );
      membraneEnergy = 0.; energyOp.assembleAdd( membraneEnergy );
      bendingEnergy = 0.;
      totalEnergy = membraneEnergy + bendingEnergy;
  }
  
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    std::vector<TripletType>  tripletListPart;
    DirichletEnergy_Hessian<MatOptConfType> ( this->_matOptConf, this->_xAStorage, this->_pf ).assembleTripletListDirichlet( tripletListPart, boundaryMask, Factor );

    tripletListMasked.reserve( 3 * ( tripletListPart.size() ) );
    for( int k=0; k<3; ++k ){
        const int offset = k * this->_matOptConf._conf.getNumGlobalDofs();
        for( int iter=0; iter<tripletListPart.size(); ++iter )
            tripletListMasked.push_back( TripletType( tripletListPart[iter].row() + offset, tripletListPart[iter].col() + offset, tripletListPart[iter].value() ) );
    }
  }
     
};








//-----------------------------------------------------------------
//-----------------------------------------------------------------
//-----------------------------------------------------------------
/**  Dirichlet energy with approx gradient
  *   E[u] = 1/2 \int_M chi(m) | theta(u) |^2
  *        = 1/2 \int_chart sqrt(det(gA)) chi(m) |(theta(u)) \circ xA )|^2 
  *        = 1/2 \int_chart sqrt(det(gA)) chi(m) \sum_comp gA^{-1} theta(u)_comp \cdot theta(u)_comp
  **/
//-----------------------------------------------------------------
//-----------------------------------------------------------------


template<typename MatOptConfType>
class DirichletEnergyApproxGradient_EvaluationHelper : public LinElastEnergy_EvaluationHelper<MatOptConfType,DirichletEnergyApproxGradient_EvaluationHelper<MatOptConfType>> {

protected: 

  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::ElementType ElementType; 
  typedef typename ConfiguratorType::VectorType VectorType;
public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
  const MatOptConfType &_matOptConf;
  const ConfiguratorType &_conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  mutable RealType _factorMembraneEnergy, _factorBendingEnergy;
  
public:
    DirichletEnergyApproxGradient_EvaluationHelper ( const MatOptConfType &matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage ) :
    LinElastEnergy_EvaluationHelper<MatOptConfType,DirichletEnergyApproxGradient_EvaluationHelper<MatOptConfType>> ( matOpConf ),
    _matOptConf ( matOpConf ), _conf ( matOpConf._conf ),
    _xAStorage ( xAStorage ),     
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
    _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ), _factorBendingEnergy( 0. ) {}
   

   RealType getFactorMembraneEnergy (  ) const { return _factorMembraneEnergy; } 
   RealType getFactorBendingEnergy (  ) const { return _factorBendingEnergy; }
   
   
  RealType evaluateMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const VectorType &u, const VectorType &v, const RealType pf ) const {

      Matrix32 dX_u, dX_v;
      DKTFEVectorFunctionEvaluator<ConfiguratorType> uDFD ( _conf, u, 3 ), vDFD( _conf, v, 3 );
      uDFD.evaluateApproxGradientAtQuadPoint( El, QuadPoint, dX_u );
      vDFD.evaluateApproxGradientAtQuadPoint( El, QuadPoint, dX_v );
      
      Matrix33 grad_u, grad_v;
      Matrix32 DxAGAinv = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint) * _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
      this->evaluateGradient( DxAGAinv, dX_u, grad_u );
      this->evaluateGradient( DxAGAinv, dX_v, grad_v );
      
      RealType Dchi = _matOptConf.approxCharFct_material_Derivative( pf, _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus());
      
      return _factorMembraneEnergy * Dchi * pesopt::ddProd<RealType,Matrix33>( grad_u, grad_v );
  }

  
  RealType evaluateMembraneStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( _matOptConf._conf, displacement, 3 );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> xADFD ( _matOptConf._conf, xA, 3 );
      PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _pointwiseEvaluator;
      Matrix32 DxA; xADFD.evaluateGradient(El,RefCoord,DxA);
      Matrix22 gA; _pointwiseEvaluatorevaluateFirstFundamentalForm( DxA, gA );
      Matrix32 DxAGAinv = DxA * (gA.inverse());
      Matrix32 dX; dispDFD.evaluateGradient( El, RefCoord, dX );
      Matrix33 grad_u;
      this->evaluateGradient ( DxAGAinv, dX, grad_u );
      RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluate( El, RefCoord ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return _factorMembraneEnergy * pesopt::ddProd<RealType,Matrix33>(grad_u,grad_u) * chi;
  }
  
  RealType evaluateMembraneStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( _matOptConf._conf, displacement, 3 );
      Matrix32 dX; dispDFD.evaluateApproxGradientAtQuadPoint( El, QuadPoint, dX );
      Matrix32 DxAGAinv = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint) * _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
      Matrix33 grad_u;
      this->evaluateGradient ( DxAGAinv, dX, grad_u );
      RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() ); 
      return _factorMembraneEnergy * pesopt::ddProd<RealType,Matrix33>(grad_u,grad_u) * chi;
  }
      
  RealType evaluateBendingStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{ return 0.;}
  RealType evaluateBendingStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{ return 0.;}
  
protected:   
  
  void evaluateGradient ( const Matrix32 &DxAGAinv, const Matrix32 &dX, Matrix33 &gradient ) const {
      gradient = DxAGAinv * dX.transpose();
  }

};


// Dirichlet Energy
template<typename MatOptConfType>
class DirichletEnergyApproxGradient :
public RefTriangleIntegrator<typename MatOptConfType::ConfiguratorType, DirichletEnergyApproxGradient<MatOptConfType> >,
public DirichletEnergyApproxGradient_EvaluationHelper<MatOptConfType>
{
  protected:

       typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
       typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::VectorType VectorType;
  
       typedef DirichletEnergyApproxGradient_EvaluationHelper<MatOptConfType> EvaluationHelper;
       public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
       const MatOptConfType &_matOptConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _disp;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> &_HardMaterial, &_SoftMaterial;
       mutable RealType _factorMembraneEnergy;
  public:
    DirichletEnergyApproxGradient ( const MatOptConfType &matOpConf, 
                    const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                    const VectorType &displacement,
                    const VectorType &pf ) :
      RefTriangleIntegrator<ConfiguratorType, DirichletEnergyApproxGradient<MatOptConfType>> (matOpConf._conf),
      DirichletEnergyApproxGradient_EvaluationHelper<MatOptConfType>( matOpConf, xAStorage ), 
      _matOptConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
     
      Matrix32 dX; _disp.evaluateApproxGradientAtQuadPoint( El, QuadPoint, dX );
      Matrix32 DxAGAinv = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint) * _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
      
      Matrix33 gradient;
      this->evaluateGradient ( DxAGAinv, dX, gradient );
      RealType aux = pesopt::ddProd<RealType,Matrix33>(gradient,gradient);
      
      RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return 0.5 * aux * chi * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * _factorMembraneEnergy;
    }
};



//! \brief Assembles \f$ \sum_{EL} \int_T \sqrt{\det g} g^{-1} D v_i \cdot D v_j \f$
template<typename MatOptConfType>
class DirichletEnergyApproxGradient_Hessian : public RefTriangleFELinAsymMatrixWeightedApproxStiffIntegrator<typename MatOptConfType::ConfiguratorType, DirichletEnergyApproxGradient_Hessian <MatOptConfType> >
{
 protected:
    typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
    typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::TripletType TripletType;
    typedef typename ConfiguratorType::MaskType MaskType;
    typedef typename ConfiguratorType::RealType   RealType;
    typedef typename ConfiguratorType::Matrix22   Matrix22;
    typedef typename ConfiguratorType::Matrix32   Matrix32;
 public:
    static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
    const MatOptConfType &_matOptConf;
    const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
    const DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
    const Material<RealType> &_HardMaterial, &_SoftMaterial;
    mutable RealType _factorMembraneEnergy;
public:
    DirichletEnergyApproxGradient_Hessian(const MatOptConfType &matOpConf,  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> & xAStorage, const VectorType &pf ) : 
    RefTriangleFELinAsymMatrixWeightedApproxStiffIntegrator< ConfiguratorType, DirichletEnergyApproxGradient_Hessian<MatOptConfType>>(matOpConf._conf),
    _matOptConf( matOpConf ),
    _xAStorage ( xAStorage ), _pf( matOpConf._confpf, pf),
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
    _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ) {}

    inline void getCoeffMatrix(const typename ConfiguratorType::ElementType &El,  const int &QuadPoint, Matrix22 &matrix) const {
        RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
        matrix  = _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
        matrix *= _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * chi * _factorMembraneEnergy;
    }
};


template <typename MatOptConfType>
class DirichletEnergyApproxGradientOp : public LinElastEnergyOpInterface<MatOptConfType, DirichletEnergyApproxGradientOp<MatOptConfType> >{
     
 protected: 
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType; 
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::MaskType MaskType;
  
 public:
  typedef DirichletEnergyApproxGradient_EvaluationHelper<MatOptConfType> EvaluationHelper;
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
     
  DirichletEnergyApproxGradientOp ( const MatOptConfType &matOpConf, 
                        const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, 
                        const VectorType &pf ) : 
   LinElastEnergyOpInterface<MatOptConfType, DirichletEnergyApproxGradientOp<MatOptConfType> > ( matOpConf ) {}
  
  void evaluateElasticEnergies( const VectorType &disp, RealType &membraneEnergy, RealType &bendingEnergy, RealType &totalEnergy ) const {
      DirichletEnergyApproxGradient<MatOptConfType> energyOp( this->_matOptConf, this->_xAStorage, disp, this->_pf );
      membraneEnergy = 0.; energyOp.assembleAdd( membraneEnergy );
      bendingEnergy = 0.;
      totalEnergy = membraneEnergy + bendingEnergy;
  }
  
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    std::vector<TripletType>  tripletListPart;
    DirichletEnergyApproxGradient_Hessian<MatOptConfType> ( this->_matOptConf, this->_xAStorage, this->_pf ).assembleTripletListDirichlet( tripletListPart, boundaryMask, Factor );

    tripletListMasked.reserve( 3 * ( tripletListPart.size() ) );
    for( int k=0; k<3; ++k ){
        const int offset = k * this->_matOptConf._conf.getNumGlobalDofs();
        for( int iter=0; iter<tripletListPart.size(); ++iter )
            tripletListMasked.push_back( TripletType( tripletListPart[iter].row() + offset, tripletListPart[iter].col() + offset, tripletListPart[iter].value() ) );
    }
  }
     
};















//-----------------------------------------------------------------
//-----------------------------------------------------------------
//-----------------------------------------------------------------
//   Laplace energy and derivatives
// E = 1/24 \int_M |lap_M (xB - xA)|^2 = 1/24 \int_chart sqrt(det(gA)) |(lap_M (xB - xA) ) \circ xA )|^2
//-----------------------------------------------------------------
//-----------------------------------------------------------------


template<typename MatOptConfType>
class LaplaceEnergy_EvaluationHelper : public LinElastEnergy_EvaluationHelper<MatOptConfType,LaplaceEnergy_EvaluationHelper<MatOptConfType>> {

protected: 
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType; 
  typedef typename ConfiguratorType::LocalVectorType LocalVectorType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::ElementType ElementType; 
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
  const MatOptConfType &_matOptConf;
  const ConfiguratorType &_conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  mutable RealType _factorMembraneEnergy, _factorBendingEnergy;
  
public:
    LaplaceEnergy_EvaluationHelper ( const MatOptConfType &matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage ) :
    LinElastEnergy_EvaluationHelper<MatOptConfType,LaplaceEnergy_EvaluationHelper<MatOptConfType>> ( matOpConf ),
    _matOptConf ( matOpConf ), _conf ( matOpConf._conf ),
    _xAStorage ( xAStorage ),     
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
    _factorMembraneEnergy( 0. ), _factorBendingEnergy ( matOpConf._materialInfo._factorBendingEnergy ) {}
  
  //This is for example used for the derivative of the Energy wrt. both the displacment and the material, s.t. the integrand has to be computed in this way
  RealType evaluateMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const Matrix32 &dX_u, const Tensor322Type &ddX_u, const Matrix32 &dX_v, const Tensor322Type &ddX_v, const RealType pf ) const {
      TangentVecType laplace_u, laplace_v;
      this->evaluateLaplaceBeltrami( _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint), dX_u, ddX_u, _xAStorage.getGInvChristoffel2(El.getGlobalElementIdx(),QuadPoint), laplace_u );
      this->evaluateLaplaceBeltrami( _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint), dX_v, ddX_v, _xAStorage.getGInvChristoffel2(El.getGlobalElementIdx(),QuadPoint), laplace_v );
      return (1./12.) * _factorBendingEnergy * _matOptConf.approxCharFct_material_Derivative( pf, _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus()) * laplace_u.dot(laplace_v);
  }
    
  RealType evaluateMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const VectorType &u, const VectorType &v, const RealType pf ) const {

      Matrix32 dX_u, dX_v; Tensor322Type ddX_u, ddX_v;
      
      DKTFEVectorFunctionEvaluator<ConfiguratorType> uDFD ( _conf, u, 3 ), vDFD( _conf, v, 3 );
      uDFD.evaluateGradientAtQuadPoint( El, QuadPoint, dX_u );
      uDFD.evaluateApproxHessianSymAtQuadPoint( El, QuadPoint, ddX_u );
      vDFD.evaluateGradientAtQuadPoint( El, QuadPoint, dX_v );
      vDFD.evaluateApproxHessianSymAtQuadPoint( El, QuadPoint, ddX_v );
      
      return this->evaluateMixedSecondDerivativeAtQuadPoint( El, QuadPoint, dX_u, ddX_u, dX_v, ddX_v, pf );
  }
  
  
  RealType getFactorMembraneEnergy (  ) const { return _factorMembraneEnergy; }  
  RealType getFactorBendingEnergy (  ) const { return _factorBendingEnergy; }
  
  RealType evaluateMembraneStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{ return 0.; }
  RealType evaluateMembraneStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{ return 0.;}
  
  RealType evaluateBendingStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( _matOptConf._conf, displacement, 3 );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> xADFD ( _matOptConf._conf, xA, 3 );
      PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _pointwiseEvaluator;
      Matrix32 dXA; xADFD.evaluateGradient( El, RefCoord, dXA ); 
      Matrix22 gA, gAinv; _pointwiseEvaluator.evaluateFirstFundamentalForm( dXA, gA ); gAinv = gA.inverse();
      Tensor322Type ddXA; xADFD.evaluateApproxHessianSym( El, RefCoord, ddXA );
      DomVecType vecForLaplaceA; _pointwiseEvaluator.evaluateVectorForLaplacian ( gAinv, dXA, ddXA, vecForLaplaceA );
      Matrix32 dX; dispDFD.evaluateGradient( El, RefCoord, dX );
      Tensor322Type ddX; dispDFD.evaluateApproxHessianSym( El, RefCoord, ddX );
      TangentVecType laplaceX;
      this->evaluateLaplaceBeltrami ( gAinv, dX, ddX, vecForLaplaceA, laplaceX );
      RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluate( El, RefCoord ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return _factorBendingEnergy * (1./12.) * laplaceX.squaredNorm() * chi;
  }
  
   RealType evaluateBendingStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( _matOptConf._conf, displacement, 3 );
      Matrix32 dX; dispDFD.evaluateGradientAtQuadPoint( El, QuadPoint, dX );
      Tensor322Type ddX; dispDFD.evaluateApproxHessianSymAtQuadPoint( El, QuadPoint, ddX );
      TangentVecType laplaceX;
      this->evaluateLaplaceBeltrami ( _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint), dX, ddX, _xAStorage.getGInvChristoffel2(El.getGlobalElementIdx(),QuadPoint), laplaceX );
      RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return _factorBendingEnergy * (1./12.) * laplaceX.squaredNorm() * chi;
  }
  
protected:   
  
  void evaluateLaplaceBeltrami ( const Matrix22 &gAinv, const Matrix32 &dX, const Tensor322Type &ddX, const DomVecType& vecForLaplaceA, TangentVecType &laplace ) const {
      for( int comp = 0; comp <3; ++comp ) laplace[comp] = pesopt::ddProd<RealType,Matrix22>( gAinv, ddX[comp] ) - vecForLaplaceA.dot( dX.row(comp) );
  }

};



// Laplace-Beltrami(u) = g^jk partial_jk (u \circ x) - g^jk Gamma_{jk}^l \partial_l (u \circ x )
template<typename MatOptConfType>
class LaplaceEnergy :
public RefTriangleIntegrator<typename MatOptConfType::ConfiguratorType, LaplaceEnergy<MatOptConfType> >,
public LaplaceEnergy_EvaluationHelper<MatOptConfType>
{
  protected:
       typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
       typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
       typedef LaplaceEnergy_EvaluationHelper<MatOptConfType> EvaluationHelper;
public:
       static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
       const MatOptConfType &_matOptConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _disp;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> & _HardMaterial, & _SoftMaterial;
       mutable RealType _factorBendingEnergy;
  public:
    LaplaceEnergy ( const MatOptConfType &matOpConf, 
                    const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                    const VectorType &displacement,
                    const VectorType &pf ) :
      RefTriangleIntegrator<ConfiguratorType, LaplaceEnergy<MatOptConfType>> (matOpConf._conf),
      LaplaceEnergy_EvaluationHelper<MatOptConfType>( matOpConf, xAStorage ), 
      _matOptConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorBendingEnergy ( matOpConf._materialInfo._factorBendingEnergy ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
     
      Matrix32 GradDisp; _disp.evaluateGradientAtQuadPoint( El, QuadPoint, GradDisp );
      
      Tensor322Type HessianDisp; _disp.evaluateApproxHessianSymAtQuadPoint( El, QuadPoint, HessianDisp );
    
      Point3DType laplace;
      this->evaluateLaplaceBeltrami ( _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint), GradDisp, HessianDisp, _xAStorage.getGInvChristoffel2(El.getGlobalElementIdx(),QuadPoint), laplace );
      RealType aux = laplace.squaredNorm();
      
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return (1./12.) * 0.5 * aux * chi * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * _factorBendingEnergy;
    }
};




template<typename MatOptConfType>
class LaplaceEnergyGradient_Part1 :
public RefTriangleMVDiffOpIntegrator<typename MatOptConfType::ConfiguratorType, LaplaceEnergyGradient_Part1<MatOptConfType> >,
public LaplaceEnergy_EvaluationHelper<MatOptConfType>
{
  protected:
       typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
       typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
       typedef LaplaceEnergy_EvaluationHelper<MatOptConfType> EvaluationHelper;
public:
       static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
       const MatOptConfType &_matOptConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _disp;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> &_HardMaterial, &_SoftMaterial;
       mutable RealType _factorBendingEnergy;
  public:
    LaplaceEnergyGradient_Part1 ( const MatOptConfType &matOpConf, 
                    const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                    const VectorType &displacement,
                    const VectorType &pf ) :
      RefTriangleMVDiffOpIntegrator<ConfiguratorType, LaplaceEnergyGradient_Part1<MatOptConfType>> (matOpConf._conf),
      LaplaceEnergy_EvaluationHelper<MatOptConfType>( matOpConf, xAStorage ), 
      _matOptConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorBendingEnergy ( matOpConf._materialInfo._factorBendingEnergy ) {}
      
    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, Matrix32 &NL) const {     
      
      Matrix32 GradDisp; _disp.evaluateGradientAtQuadPoint( El, QuadPoint, GradDisp );
      Tensor322Type HessianDisp; _disp.evaluateApproxHessianSymAtQuadPoint( El, QuadPoint, HessianDisp );
        
      for( int comp = 0; comp <3; ++comp ){
        RealType aux1 = pesopt::ddProd<RealType, Matrix22>( _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint), HessianDisp[comp] );
        RealType aux2 = (_xAStorage.getGInvChristoffel2(El.getGlobalElementIdx(),QuadPoint)).dot( GradDisp.row(comp) );
        NL.row(comp) = (aux2 - aux1) * ( _xAStorage.getGInvChristoffel2(El.getGlobalElementIdx(),QuadPoint) );
      }
      
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      
      NL *= (1./12.) * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * chi * _factorBendingEnergy;
    }
}; 



template<typename MatOptConfType>
class LaplaceEnergyGradient_Part2 :
public RefTriangleMVDiff2OpIntegrator<typename MatOptConfType::ConfiguratorType, LaplaceEnergyGradient_Part2<MatOptConfType> >,
public LaplaceEnergy_EvaluationHelper<MatOptConfType>
{ 
  protected:
       typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
       typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
       typedef LaplaceEnergy_EvaluationHelper<MatOptConfType> EvaluationHelper;
 public:
       static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const MatOptConfType &_matOptConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _disp;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> &_HardMaterial, &_SoftMaterial;
       mutable RealType _factorBendingEnergy;
  public:
    LaplaceEnergyGradient_Part2 ( const MatOptConfType &matOpConf, 
                    const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                    const VectorType &displacement,
                    const VectorType &pf ) :
      RefTriangleMVDiff2OpIntegrator<ConfiguratorType, LaplaceEnergyGradient_Part2<MatOptConfType>> (matOpConf._conf),
      LaplaceEnergy_EvaluationHelper<MatOptConfType>( matOpConf, xAStorage ), 
      _matOptConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorBendingEnergy ( matOpConf._materialInfo._factorBendingEnergy ) {}

    void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint, Tensor322Type &NL) const {     
      
      Matrix32 GradDisp; _disp.evaluateGradientAtQuadPoint( El, QuadPoint, GradDisp );
      Tensor322Type HessianDisp; _disp.evaluateApproxHessianSymAtQuadPoint( El, QuadPoint, HessianDisp );

      for( int comp = 0; comp <3; ++comp ){
        RealType aux1 = pesopt::ddProd<RealType,Matrix22>(_xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint), HessianDisp[comp] );
        RealType aux2 = (_xAStorage.getGInvChristoffel2(El.getGlobalElementIdx(),QuadPoint)).dot( GradDisp.row(comp) );
        NL[comp] = (aux1 - aux2) *_xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
      }
      
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      
      NL *= (1./12.) * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * chi * _factorBendingEnergy; 
    }
}; 


// int sqrt(det_g) g^-1 : D^2(u_i \circ x)   g^-1 : D^2(u_j \circ x)
template<typename MatOptConfType>
class LaplaceEnergySubHessian_Part1 :
public RefTriangleFELinMatrixSeperatedDiff2Integrator <typename MatOptConfType::ConfiguratorType, LaplaceEnergySubHessian_Part1 <MatOptConfType>> 
 {
   protected:
       typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
       typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
     const MatOptConfType &_matOptConf;   
     const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
     DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
     const Material<RealType> &_HardMaterial, &_SoftMaterial;
     mutable RealType _factorBendingEnergy;
 
   public:
     LaplaceEnergySubHessian_Part1 ( const MatOptConfType &matOpConf, 
                                     const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, 
                                     const VectorType &pf ) :
      RefTriangleFELinMatrixSeperatedDiff2Integrator <ConfiguratorType, LaplaceEnergySubHessian_Part1<MatOptConfType>> (matOpConf._conf),
      _matOptConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorBendingEnergy ( matOpConf._materialInfo._factorBendingEnergy ) {}
 
     inline void getCoeffMatrices ( const ElementType &El, int QuadPoint, Matrix22 &MatrixA, Matrix22 &MatrixB ) const{
       const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
       MatrixA = _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
       MatrixA *= (1./12.) * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * _factorBendingEnergy * chi;
       MatrixB = _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
     }
 };

// int sqrt(det_g) sum_jkl g^jk Gamma_jk^l partial_l (u_j circ x)   g^-1 : D^2(u_i \circ x)
template<typename MatOptConfType>
class LaplaceEnergySubHessian_Part2 :
public RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator <typename MatOptConfType::ConfiguratorType, LaplaceEnergySubHessian_Part2 <MatOptConfType> > 
{
 protected:
       typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
       typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::DomVecType DomVecType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
     static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
     const MatOptConfType &_matOptConf;    
     const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
     DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
     const Material<RealType> &_HardMaterial, &_SoftMaterial;
     mutable RealType _factorBendingEnergy;
 
   public:
     LaplaceEnergySubHessian_Part2 (  const MatOptConfType &matOpConf,
                                      const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, 
                                      const VectorType &pf ) :
      RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator <ConfiguratorType, LaplaceEnergySubHessian_Part2<MatOptConfType>> (matOpConf._conf),
      _matOptConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorBendingEnergy ( matOpConf._materialInfo._factorBendingEnergy ) {}
 
     inline void getCoeffMatrixAndVector ( const ElementType &El, int QuadPoint,
                                           DomVecType &Vector, Matrix22 &Matrix ) const{
       const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
                                               
       Matrix = _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
       Matrix *= -1. * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * _factorBendingEnergy * chi;
       
       Vector = (1./12.) * _xAStorage.getGInvChristoffel2(El.getGlobalElementIdx(),QuadPoint);

     }
 };
 
 
 // int sqrt(det_g)  sum_jkl g^jk Gamma_jk^l partial_l (u_i circ x)  sum_jkl g^jk Gamma_jk^l partial_l (u_j circ x)
template<typename MatOptConfType>
class LaplaceEnergySubHessian_Part3 :
public RefTriangleFELinSeperatedVectorWeightedStiffIntegrator <typename MatOptConfType::ConfiguratorType, LaplaceEnergySubHessian_Part3 <MatOptConfType>> 
{
   protected:
       typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
       typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::DomVecType DomVecType;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
     static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
     const MatOptConfType &_matOptConf;  
     const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
     DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
     const Material<RealType> &_HardMaterial, &_SoftMaterial;
     mutable RealType _factorBendingEnergy;
 
   public:
     LaplaceEnergySubHessian_Part3 (  const MatOptConfType &matOpConf,
                                      const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, 
                                      const VectorType &pf ) :
      RefTriangleFELinSeperatedVectorWeightedStiffIntegrator <ConfiguratorType, LaplaceEnergySubHessian_Part3<MatOptConfType>> (matOpConf._conf),
      _matOptConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorBendingEnergy ( matOpConf._materialInfo._factorBendingEnergy ) {}
 
     inline void getCoeffVectors ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                                   DomVecType &vec1, DomVecType &vec2 ) const{
       const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
                                       
       vec1 = _xAStorage.getGInvChristoffel2(El.getGlobalElementIdx(),QuadPoint);
       vec1 *= _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * _factorBendingEnergy * chi;
      
       vec2 = (1./12.) * _xAStorage.getGInvChristoffel2(El.getGlobalElementIdx(),QuadPoint);
     }
};


template <typename MatOptConfType>
class LaplaceEnergyOp : public LinElastEnergyOpInterface<MatOptConfType, LaplaceEnergyOp<MatOptConfType> > {
     
 protected: 
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType; 
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::MaskType MaskType;
  
 public:
  typedef LaplaceEnergy_EvaluationHelper<MatOptConfType> EvaluationHelper;
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
     
 protected:
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage; 
  const VectorType &_pf; 
 public:
  
  LaplaceEnergyOp ( const MatOptConfType &matOpConf, 
                      const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, 
                      const VectorType &pf ) : 
    LinElastEnergyOpInterface<MatOptConfType, LaplaceEnergyOp<MatOptConfType> > ( matOpConf ),
    _xAStorage( xAStorage ), _pf( pf ){}
  
  void evaluateElasticEnergies( const VectorType &disp, RealType &membraneEnergy, RealType &bendingEnergy, RealType &totalEnergy ) const {
      LaplaceEnergy<MatOptConfType> energyOp( this->_matOptConf, this->_xAStorage, disp, this->_pf );
      membraneEnergy = 0.;
      bendingEnergy = 0.; energyOp.assembleAdd( bendingEnergy );
      totalEnergy = membraneEnergy + bendingEnergy;
  }
  
  void evaluateGradient( const VectorType &disp, VectorType &Deriv ) const {
    Deriv.setZero();
    LaplaceEnergyGradient_Part1<MatOptConfType> ( this->_matOptConf, this->_xAStorage, disp, this->_pf ).assembleAdd( Deriv );
    LaplaceEnergyGradient_Part2<MatOptConfType> ( this->_matOptConf, this->_xAStorage, disp, this->_pf ).assembleAdd( Deriv );
  }
  
  
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    std::vector<TripletType>  tripletListPart1, tripletListPart2, tripletListPart3;
    LaplaceEnergySubHessian_Part1<MatOptConfType> ( this->_matOptConf, this->_xAStorage, this->_pf ).assembleTripletListDirichlet( tripletListPart1, boundaryMask, Factor );
    LaplaceEnergySubHessian_Part2<MatOptConfType> ( this->_matOptConf, this->_xAStorage, this->_pf ).assembleTripletListDirichlet( tripletListPart2, boundaryMask, Factor );
    LaplaceEnergySubHessian_Part3<MatOptConfType> ( this->_matOptConf, this->_xAStorage, this->_pf ).assembleTripletListDirichlet( tripletListPart3, boundaryMask, Factor );

    tripletListMasked.reserve( 3 * ( tripletListPart1.size() + 2 * tripletListPart2.size() + tripletListPart3.size() ) );
    for( int k=0; k<3; ++k ){
        const int offset = k * this->_matOptConf._conf.getNumGlobalDofs();
        for( int iter=0; iter<tripletListPart1.size(); ++iter )
            tripletListMasked.push_back( TripletType( tripletListPart1[iter].row() + offset, tripletListPart1[iter].col() + offset, tripletListPart1[iter].value() ) );
        for( int iter=0; iter<tripletListPart2.size(); ++iter ){
              tripletListMasked.push_back( TripletType( tripletListPart2[iter].row() + offset, tripletListPart2[iter].col() + offset, tripletListPart2[iter].value() ) );
              tripletListMasked.push_back( TripletType( tripletListPart2[iter].col() + offset, tripletListPart2[iter].row() + offset, tripletListPart2[iter].value() ) );
        }
        for( int iter=0; iter<tripletListPart3.size(); ++iter )
            tripletListMasked.push_back( TripletType( tripletListPart3[iter].row() + offset, tripletListPart3[iter].col() + offset, tripletListPart3[iter].value() ) );
    }
  }
     
};














//-----------------------------------------------------------------
//-----------------------------------------------------------------
//-----------------------------------------------------------------
//  KirchhoffLove
// E = 1/2 \int_M |g_B - g_A|^2 + |h_B - h_A|^2 with linearized g,h
//-----------------------------------------------------------------
//-----------------------------------------------------------------
    
template<typename MatOptConfType>
class KirchhoffLoveEnergy_EvaluationHelper : public LinElastEnergy_EvaluationHelper<MatOptConfType,KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType>> {

protected: 

  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType; 
  typedef typename ConfiguratorType::LocalVectorType LocalVectorType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::ElementType ElementType; 
  typedef typename ConfiguratorType::VectorType VectorType;
public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
protected:
  const MatOptConfType &_matOptConf;
  const ConfiguratorType &_conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  const Material<RealType> & _HardMaterial, & _SoftMaterial;
  mutable RealType _factorMembraneEnergy, _factorBendingEnergy;
  mutable RealType _membCoefHard, _bendCoefHard, _membCoefSoft, _bendCoefSoft;
  
public:
    KirchhoffLoveEnergy_EvaluationHelper ( const MatOptConfType &matOpConf, const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage ) :
    LinElastEnergy_EvaluationHelper<MatOptConfType,KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType>> ( matOpConf ),
    _matOptConf( matOpConf ), _conf ( matOpConf._conf ),
    _xAStorage ( xAStorage ),
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
    _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ), _factorBendingEnergy ( matOpConf._materialInfo._factorBendingEnergy ),
    _membCoefHard ( computeMembraneCoef( _HardMaterial.getElastModulus(), _HardMaterial.getPoissonRatio() ) ), 
    _bendCoefHard ( computeBendingCoef(  _HardMaterial.getElastModulus(), _HardMaterial.getPoissonRatio() ) ), 
    _membCoefSoft ( computeMembraneCoef( _SoftMaterial.getElastModulus(), _SoftMaterial.getPoissonRatio() ) ),
    _bendCoefSoft ( computeBendingCoef(  _SoftMaterial.getElastModulus(), _SoftMaterial.getPoissonRatio() ) ) { }  
  
  //The factor is possible, because this can either used only for the integrand or for another interface which need a factor e.g. from the chain rule
  void assembleAddLocalMatrix( const ElementType& El, const int QuadPoint, 
                               const Matrix33 & materialMatrixMembrane, const Matrix33 & materialMatrixBending, 
                               LocalMatrixType (&localMatrix)[3][3],
                               const RealType factor) const{
        
    const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet(El); 
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &approxbfs = _conf.getApproxGradientBaseFunctionSet(El); 
    const int numLocDofs = _conf.getNumLocalDofs(El); 

    //Cache membrane and bending matrix for all local dofs
    Matrix33 membMatrix [ConfiguratorType::maxNumLocalDofs], bendMatrix [ConfiguratorType::maxNumLocalDofs];
    for (int dof = 0; dof < numLocDofs; ++dof) {
        this->computeMembraneMatrix( bfs.evaluateGradientOnRefTriang( dof, QuadPoint ), _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint), membMatrix[dof] ); 
        this->computeBendingMatrix ( bfs.evaluateGradientOnRefTriang( dof, QuadPoint ), approxbfs.evaluateApproxHessianAsVecOnRefTriang( dof, QuadPoint ), _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint), _xAStorage.getHessianAsVec(El.getGlobalElementIdx(),QuadPoint), _xAStorage.getNormal(El.getGlobalElementIdx(),QuadPoint), _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint), bendMatrix[dof]);
    }   
    
    Matrix33 totEnergyMat, tempMembMat, tempBendMat, membEnergyMat, bendEnergyMat; 
    for (int locIdx1 = 0; locIdx1 < numLocDofs; ++locIdx1) {

        tempMembMat = membMatrix[locIdx1].transpose() * materialMatrixMembrane;
        tempBendMat = bendMatrix[locIdx1].transpose() * materialMatrixBending; 

        for (int locIdx2 = 0; locIdx2 < numLocDofs; ++locIdx2) {

          membEnergyMat = tempMembMat * membMatrix[locIdx2];
          bendEnergyMat = tempBendMat * bendMatrix[locIdx2];
          totEnergyMat  = factor * ( membEnergyMat + bendEnergyMat );

          // Add contribution of locIdx1 and locIdx2 to localMatrix. 
          for (int argComp = 0; argComp < 3; ++argComp)
            for (int destComp = 0; destComp < 3; ++destComp)
              localMatrix[argComp][destComp]( locIdx1, locIdx2 ) += totEnergyMat(argComp,destComp); 
          
        } // End iteration of locIdx2.
      } // End iteration of locIdx1.
  }
  
  
  void assembleAddLocalMatrix( const ElementType& El, const DomVecType &RefCoord, 
                               const Matrix32 &dXA, const RealType sqrtDetgAinv, const Matrix33 &ddXA, const TangentVecType &normalA,
                               const Matrix33 & materialMatrixMembrane, const Matrix33 & materialMatrixBending, 
                               LocalMatrixType (&localMatrix)[3][3],
                               const RealType factor) const{
        
    const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet(El); 
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &approxbfs = _conf.getApproxGradientBaseFunctionSet(El); 
    const int numLocDofs = _conf.getNumLocalDofs(El); 

    //Cache membrane and bending matrix for all local dofs
    Matrix33 membMatrix [ConfiguratorType::maxNumLocalDofs], bendMatrix [ConfiguratorType::maxNumLocalDofs];
    for (int dof = 0; dof < numLocDofs; ++dof) {
        DomVecType gradientBf; bfs.evaluateGradientOnRefTriang( dof, RefCoord, gradientBf );
        TangentVecType hessianBf; approxbfs.evaluateApproxHessianAsVecOnRefTriang( dof, RefCoord, hessianBf );
        this->computeMembraneMatrix( gradientBf, dXA, membMatrix[dof] ); 
        this->computeBendingMatrix ( gradientBf, hessianBf, dXA, ddXA, normalA, sqrtDetgAinv, bendMatrix[dof]);
    }   
    
    Matrix33 totEnergyMat, tempMembMat, tempBendMat, membEnergyMat, bendEnergyMat; 
    for (int locIdx1 = 0; locIdx1 < numLocDofs; ++locIdx1) {

        tempMembMat = membMatrix[locIdx1].transpose() * materialMatrixMembrane;
        tempBendMat = bendMatrix[locIdx1].transpose() * materialMatrixBending; 

        for (int locIdx2 = 0; locIdx2 < numLocDofs; ++locIdx2) {

          membEnergyMat = tempMembMat * membMatrix[locIdx2];
          bendEnergyMat = tempBendMat * bendMatrix[locIdx2];
          totEnergyMat  = factor * ( membEnergyMat + bendEnergyMat );

          // Add contribution of locIdx1 and locIdx2 to localMatrix. 
          for (int argComp = 0; argComp < 3; ++argComp)
            for (int destComp = 0; destComp < 3; ++destComp)
              localMatrix[argComp][destComp]( locIdx1, locIdx2 ) += totEnergyMat(argComp,destComp); 
          
        } // End iteration of locIdx2.
      } // End iteration of locIdx1.
  }
  
  
  //===========================================================================================================================
  RealType applyLocalMatrix ( const ElementType &El, const int QuadPoint, 
                              const Matrix33& materialMatrixMembrane, const Matrix33& materialMatrixBending, 
                              const VectorType &u, const VectorType & v, 
                              const RealType factor ) const {
    
    LocalMatrixType localMatrix[3][3];
    for (int argComp = 0; argComp < 3; ++argComp)
        for (int destComp = 0; destComp < 3; ++destComp)
           localMatrix[argComp][destComp].setZero();  
    
    this->assembleAddLocalMatrix( El, QuadPoint, materialMatrixMembrane, materialMatrixBending, localMatrix, factor );
    
    // apply vectors
    const int numGlobalDofs = _conf.getNumGlobalDofs();
    const int numLocDofs = _conf.getNumLocalDofs(El); 
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( int i = 0; i < numLocDofs; ++i ) globalDofs[ i ] = _conf.localToGlobal ( El, i );
    RealType aux = 0.;
    for (int i = 0; i < numLocDofs; ++i)
        for (int j = 0; j < numLocDofs; ++j) 
          for (int argComp = 0; argComp < 3; ++argComp)
            for (int destComp = 0; destComp < 3; ++destComp)
              aux += u[globalDofs[j] + numGlobalDofs * argComp] * v[globalDofs[i] + numGlobalDofs * destComp ] * localMatrix[argComp][destComp]( i, j ); 
    
    return aux;
  }
  
  
  RealType applyLocalMatrix ( const ElementType &El, const DomVecType &RefCoord, 
                              const Matrix32 &dXA, const RealType sqrtDetgAinv, const Matrix33 &ddXA, const TangentVecType &normalA,
                              const Matrix33& materialMatrixMembrane, const Matrix33& materialMatrixBending, 
                              const VectorType &u, const VectorType & v, 
                              const RealType factor ) const {
    LocalMatrixType localMatrix[3][3];
    for (int argComp = 0; argComp < 3; ++argComp)
        for (int destComp = 0; destComp < 3; ++destComp)
           localMatrix[argComp][destComp].setZero();  
        
    this->assembleAddLocalMatrix( El, RefCoord, dXA, sqrtDetgAinv, ddXA, normalA,  materialMatrixMembrane, materialMatrixBending, localMatrix, factor );
    
    // apply vectors
    const int numGlobalDofs = _conf.getNumGlobalDofs();
    const int numLocDofs = _conf.getNumLocalDofs(El); 
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( int i = 0; i < numLocDofs; ++i ) globalDofs[ i ] = _conf.localToGlobal ( El, i );
    
    RealType aux = 0.;
    for (int i = 0; i < numLocDofs; ++i)
        for (int j = 0; j < numLocDofs; ++j) 
          for (int argComp = 0; argComp < 3; ++argComp)
            for (int destComp = 0; destComp < 3; ++destComp)
              aux += u[globalDofs[j] + numGlobalDofs * argComp] * v[globalDofs[i] + numGlobalDofs * destComp ] * localMatrix[argComp][destComp]( i, j ); 
    
    return aux;
  }

  RealType getFactorMembraneEnergy (  ) const { return _factorMembraneEnergy; }  
  RealType getFactorBendingEnergy (  ) const { return _factorBendingEnergy; }

  
  RealType evaluateMembraneStress( const ElementType &El, const DomVecType &RefCoord, 
                                   const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( _matOptConf._conf, displacement, 3 );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> xADFD ( _matOptConf._conf, xA, 3 );
      PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _pointwiseEvaluator;
      
      const RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluate( El, RefCoord ), 1., 0. ); //E and nu are in material matrix
      Matrix33 materialMatrixMembrane, materialMatrixBending; materialMatrixBending.setZero();  
      
      Matrix32 dXA; xADFD.evaluateGradient( El, RefCoord, dXA );
      Matrix22 gA, gAinv; _pointwiseEvaluator.evaluateFirstFundamentalForm( dXA, gA ); gAinv = gA.inverse();
      RealType sqrtDetgAinv = sqrt( gAinv.determinant() );
      Matrix33 ddXA; xADFD.evaluateApproxHessianAsVec( El, RefCoord, ddXA );
      TangentVecType normalA; _pointwiseEvaluator.evaluateNormal( dXA, normalA );
      Matrix33 materialMatrixHard, materialMatrixSoft;
      this->computeMaterialMatrix( gAinv, _HardMaterial.getPoissonRatio(), materialMatrixHard ); 
      this->computeMaterialMatrix( gAinv, _SoftMaterial.getPoissonRatio(), materialMatrixSoft );
      materialMatrixMembrane = chi * _membCoefHard * materialMatrixHard + (1. - chi) * _membCoefSoft * materialMatrixSoft;

      RealType aux = this->applyLocalMatrix ( El, RefCoord, dXA, sqrtDetgAinv, ddXA, normalA, materialMatrixMembrane, materialMatrixBending, displacement, displacement, 1.0 );
      return 0.5 * aux * sqrt( gA.determinant() );
  }
  
  RealType evaluateMembraneStressAtQuadPoint( const ElementType &El, const int QuadPoint,
                                              const VectorType &material, const VectorType &displacement ) const{
        DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      
        const RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluateAtQuadPoint( El, QuadPoint ), 1., 0. ); //E and nu are in material matrix
        Matrix33 materialMatrixMembrane, materialMatrixBending; 
        this->computeMaterialMatricesWithPf( El, QuadPoint, chi, materialMatrixMembrane, materialMatrixBending );
        materialMatrixBending.setZero();
        RealType aux = this->applyLocalMatrix ( El, QuadPoint, materialMatrixMembrane, materialMatrixBending, displacement, displacement, 1.0 );  
        return 0.5 * aux * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
  }
  
  RealType evaluateBendingStress( const ElementType &El, const DomVecType &RefCoord, 
                                  const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{
      DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
      DKTFEVectorFunctionEvaluator<ConfiguratorType> xADFD ( _matOptConf._conf, xA, 3 );
      PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> _pointwiseEvaluator;
      
      const RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluate( El, RefCoord ), 1., 0. ); //E and nu are in material matrix
      Matrix33 materialMatrixMembrane, materialMatrixBending; materialMatrixMembrane.setZero();  
      Matrix32 dXA; xADFD.evaluateGradient( El, RefCoord, dXA );
      Matrix22 gA, gAinv; _pointwiseEvaluator.evaluateFirstFundamentalForm( dXA, gA ); gAinv = gA.inverse();
      RealType sqrtDetgAinv = sqrt( gAinv.determinant() );
      Matrix33 ddXA; xADFD.evaluateApproxHessianAsVec( El, RefCoord, ddXA );
      TangentVecType normalA; _pointwiseEvaluator.evaluateNormal( dXA, normalA );
      Matrix33 materialMatrixHard, materialMatrixSoft;
      this->computeMaterialMatrix( gAinv, _HardMaterial.getPoissonRatio(), materialMatrixHard ); 
      this->computeMaterialMatrix( gAinv, _SoftMaterial.getPoissonRatio(), materialMatrixSoft );
      materialMatrixBending = chi * _bendCoefHard * materialMatrixHard + (1. - chi) * _bendCoefSoft * materialMatrixSoft;

      RealType aux = this->applyLocalMatrix ( El, RefCoord, dXA, sqrtDetgAinv, ddXA, normalA, materialMatrixMembrane, materialMatrixBending, displacement, displacement, 1.0 );
      return 0.5 * aux * sqrt( gA.determinant() );
  }
  
  
  RealType evaluateBendingStressAtQuadPoint( const ElementType &El, const int QuadPoint, 
                                             const VectorType &material, const VectorType &displacement ) const{
        DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> materialDFD ( _matOptConf._confpf, material );
        DKTFEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( _matOptConf._conf, displacement, 3 );
        const RealType chi = _matOptConf.approxCharFct_material ( materialDFD.evaluateAtQuadPoint( El, QuadPoint ), 1., 0. ); //E and nu are in material matrix
        Matrix33 materialMatrixMembrane, materialMatrixBending; 
        this->computeMaterialMatricesWithPf( El, QuadPoint, chi, materialMatrixMembrane, materialMatrixBending );
        materialMatrixMembrane.setZero();
        RealType aux = this->applyLocalMatrix ( El, QuadPoint, materialMatrixMembrane, materialMatrixBending, displacement, displacement, 1.0 );
        return 0.5 * aux * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
  }
      
      
  //for example used for D^2_{disp,material} E, where the integrand has to be computed in this way
   RealType evaluateMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const VectorType &u, const VectorType & v, const RealType pf ) const {
      Matrix33 materialMatrixMembrane, materialMatrixBending;
      computeMaterialMatricesWithPf_Derivative ( El, QuadPoint, materialMatrixMembrane, materialMatrixBending );
      return applyLocalMatrix ( El, QuadPoint, materialMatrixMembrane, materialMatrixBending, u, v, 1.0) * _matOptConf.approxCharFct_material_Derivative( pf, 1., 0. );
  }
  
  RealType evaluateBendingMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const VectorType &u, const VectorType & v, const RealType pf ) const {
      Matrix33 materialMatrixMembrane, materialMatrixBending;
      computeMaterialMatricesWithPf_Derivative ( El, QuadPoint, materialMatrixMembrane, materialMatrixBending );
      materialMatrixMembrane.setZero();
      return applyLocalMatrix ( El, QuadPoint, materialMatrixMembrane, materialMatrixBending, u, v, 1.0) * _matOptConf.approxCharFct_material_Derivative( pf, 1., 0. );
  }
  
  //===========================================================================================================================
  void assembleLocalVector_Gradient ( const ElementType &El, const int QuadPoint, 
                                      const Matrix33& materialMatrixMembrane, const Matrix33& materialMatrixBending, 
                                      const VectorType &u,  
                                      LocalVectorType (&localVector)[3],
                                      const RealType factor ) const {
    
    for (int destComp = 0; destComp < 3; ++destComp)
        localVector[destComp].setZero();  
                                          
    LocalMatrixType localMatrix[3][3];
    for (int argComp = 0; argComp < 3; ++argComp)
        for (int destComp = 0; destComp < 3; ++destComp)
           localMatrix[argComp][destComp].setZero();  
    
    this->assembleAddLocalMatrix( El, QuadPoint, materialMatrixMembrane, materialMatrixBending, localMatrix, factor );
    
    // apply vectors
    const int numGlobalDofs = _conf.getNumGlobalDofs();
    const int numLocDofs = _conf.getNumLocalDofs(El); 
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( int i = 0; i < numLocDofs; ++i )  globalDofs[ i ] = _conf.localToGlobal ( El, i );
    
    for (int i = 0; i < numLocDofs; ++i)
        for (int j = 0; j < numLocDofs; ++j) 
          for (int argComp = 0; argComp < 3; ++argComp)
            for (int destComp = 0; destComp < 3; ++destComp)
              localVector[destComp](i) += u[globalDofs[j] + numGlobalDofs * argComp] * localMatrix[argComp][destComp]( i, j ); 
  }


 //===========================================================================================================================
protected:   

  void computeMaterialMatrix( const Matrix22 &gInv, const RealType poissonRatio, Matrix33 &materialMatrix) const {
    for( int comp=0; comp<2; ++comp ){
      materialMatrix(comp,comp) = pesopt::Sqr(gInv(comp, comp));
      materialMatrix(comp,2) = gInv(comp,comp) * gInv(0,1);
    }
    materialMatrix(0,1) =            poissonRatio * gInv(0,0) * gInv(1,1)  + (1. - poissonRatio)       * pesopt::Sqr(gInv(0,1));
    materialMatrix(2,2) = 0.5 * (1.-poissonRatio) * gInv(0,0) * gInv(1,1)  + 0.5 * (1. + poissonRatio) * pesopt::Sqr(gInv(0,1));
    // use symmetry 
    for( int k=1; k<3; ++k )
      for( int l=0; l<k; ++l )
        materialMatrix(k,l) = materialMatrix(l,k);
  }
  
  // computes chi C_h + (1-chi) C_s
  void computeMaterialMatricesWithPf ( const ElementType& El, const int QuadPoint, const RealType chi, Matrix33& materialMatrixMembrane, Matrix33& materialMatrixBending  ) const {
      Matrix33 materialMatrixHard, materialMatrixSoft;
      this->computeMaterialMatrix( _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint), _HardMaterial.getPoissonRatio(), materialMatrixHard ); 
      this->computeMaterialMatrix( _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint), _SoftMaterial.getPoissonRatio(), materialMatrixSoft ); 
      materialMatrixMembrane = chi * _membCoefHard * materialMatrixHard + (1. - chi) * _membCoefSoft * materialMatrixSoft;
      materialMatrixBending  = chi * _bendCoefHard * materialMatrixHard + (1. - chi) * _bendCoefSoft * materialMatrixSoft;
  }; 
      
  // computes  C_h - C_s
  void computeMaterialMatricesWithPf_Derivative ( const ElementType& El, const int QuadPoint, Matrix33& materialMatrixMembrane, Matrix33& materialMatrixBending  ) const {
      Matrix33 materialMatrixHard, materialMatrixSoft;
      this->computeMaterialMatrix( _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint), _HardMaterial.getPoissonRatio(), materialMatrixHard ); 
      this->computeMaterialMatrix( _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint), _SoftMaterial.getPoissonRatio(), materialMatrixSoft ); 
      materialMatrixMembrane  = _membCoefHard * materialMatrixHard - _membCoefSoft * materialMatrixSoft;
      materialMatrixBending   = _bendCoefHard * materialMatrixHard - _bendCoefSoft * materialMatrixSoft;
  };

  void computeMembraneMatrix( const DomVecType &gradientBf, const Matrix32 &dxA,  Matrix33 &membMatrix) const {
    for( int i=0; i<2; ++i )
        for( int j=0; j<3; ++j )
            membMatrix(i,j) = dxA(j,i) * gradientBf[i];
        
    for( int j=0; j<3; ++j )
        membMatrix(2,j) = dxA(j,1) * gradientBf[0] + dxA(j,0) * gradientBf[1];
  }

  void computeBendingMatrix( const DomVecType &gradientBf, const TangentVecType &hessianBf,
                             const Matrix32 &dxA, const Matrix33 &ddxA, const TangentVecType &normalA, const RealType &area, 
                             Matrix33 &bendMatrix) const {
    
    TangentVecType dxA_Col[2]; TangentVecType ddxA_Col[3];
    dxA_Col[0] = TangentVecType( dxA.col(0) ); dxA_Col[1] = TangentVecType( dxA.col(1) );
    ddxA_Col[0] = TangentVecType( ddxA.col(0) ); ddxA_Col[1] = TangentVecType( ddxA.col(1) ); ddxA_Col[2] = TangentVecType( ddxA.col(2) );
                                 
    TangentVecType helper( gradientBf[0] * (dxA_Col[1]).cross(normalA) + gradientBf[1] * normalA.cross(dxA_Col[0]) );
    
    for( int comp=0; comp<3; ++comp ){
        TangentVecType b (normalA);
        b *= hessianBf(comp); //i.e. for comp=2: b*=H(0,1)+H(1,0)
        for( int l=0; l<3; ++l ){
            RealType factor = 1. / area;
            b += factor * (  gradientBf[0] * (ddxA_Col[l]).cross( dxA_Col[1] ) + gradientBf[1] * dxA_Col[0].cross( ddxA_Col[l] ) + normalA.dot(ddxA_Col[l]) * helper );
        }
        for( int col=0; col<3; ++col )
            bendMatrix(comp,col) = b[col];
    }
  }

  RealType computeMembraneCoef( const RealType E, const RealType nu ) const { return _factorMembraneEnergy * E / (1. - pesopt::Sqr(nu)); }
  RealType computeBendingCoef( const RealType E, const RealType nu ) const {return _factorBendingEnergy * E / ( 12 * (1. - pesopt::Sqr(nu)) ); }

};


template<typename MatOptConfType>
class KirchhoffLoveMembraneEnergy :
public RefTriangleIntegrator<typename MatOptConfType::ConfiguratorType, KirchhoffLoveMembraneEnergy<MatOptConfType> >,
public KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType>
{
  protected:
       typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
       typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
       typedef KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType> EvaluationHelper;
 public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const MatOptConfType &_matOptConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
       const VectorType &_disp;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _dispDFD;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> & _HardMaterial;
       const Material<RealType> & _SoftMaterial;
       mutable RealType _factorMembraneEnergy;
  public:
    KirchhoffLoveMembraneEnergy ( const MatOptConfType &matOpConf, 
                    const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                    const VectorType &displacement,
                    const VectorType &pf ) :
      RefTriangleIntegrator<ConfiguratorType, KirchhoffLoveMembraneEnergy<MatOptConfType>> (matOpConf._conf),
      KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType>( matOpConf, xAStorage ), 
      _matOptConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( displacement ), _dispDFD( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ),
      _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorMembraneEnergy ( matOpConf._materialInfo._factorMembraneEnergy ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), 1., 0. ); //E and nu are in material matrix
        Matrix33 materialMatrixMembrane, materialMatrixBending; 
        this->computeMaterialMatricesWithPf( El, QuadPoint, chi, materialMatrixMembrane, materialMatrixBending );
        materialMatrixBending.setZero();
        RealType aux = this->applyLocalMatrix ( El, QuadPoint, materialMatrixMembrane, materialMatrixBending, _disp, _disp, 1.0 );
        return 0.5 * aux * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};

    

template<typename MatOptConfType>
class KirchhoffLoveBendingEnergy :
public RefTriangleIntegrator<typename MatOptConfType::ConfiguratorType, KirchhoffLoveBendingEnergy<MatOptConfType> >,
public KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType>
{
  protected:
       typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
       typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
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
       const MatOptConfType &_matOptConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
       const VectorType &_disp;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _dispDFD;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> & _HardMaterial;
       const Material<RealType> & _SoftMaterial;
       mutable RealType _factorBendingEnergy;
  public:
    KirchhoffLoveBendingEnergy ( const MatOptConfType &matOpConf, 
                    const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                    const VectorType &displacement,
                    const VectorType &pf ) :
      RefTriangleIntegrator<ConfiguratorType, KirchhoffLoveBendingEnergy<MatOptConfType>> (matOpConf._conf),
      KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType>( matOpConf, xAStorage ), 
      _matOptConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( displacement ), _dispDFD( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ),
      _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _factorBendingEnergy ( matOpConf._materialInfo._factorBendingEnergy ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), 1., 0. ); //E and nu are in material matrix
        Matrix33 materialMatrixMembrane, materialMatrixBending; 
        this->computeMaterialMatricesWithPf( El, QuadPoint, chi, materialMatrixMembrane, materialMatrixBending );
        materialMatrixMembrane.setZero();
        RealType aux = this->applyLocalMatrix ( El, QuadPoint, materialMatrixMembrane, materialMatrixBending, _disp, _disp, 1.0 );
        return 0.5 * aux * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};




//! D^2_(m,u) E(m, u(m)) (hat m)(P)
// see also ComplianceEnergy_DerivativePartAdjoint in ComplianceEnergies.h
template<typename MatOptConfType>
class KirchhoffLoveBendingEnergyMixedSecondDerivative
: public RefTriangleFENonlinOpIntegrator< typename MatOptConfType::ConfiguratorTypePf, KirchhoffLoveBendingEnergyMixedSecondDerivative<MatOptConfType> >,
 public KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType> {

protected: 
  
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename MatOptConfType::RealType RealType;
  typedef typename ConfiguratorTypePf::VectorType VectorType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType> EvaluationHelper;
  
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage;
  const VectorType &_disp;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const VectorType &_adjoint;
  
public:
  KirchhoffLoveBendingEnergyMixedSecondDerivative( const MatOptConfType &matOptConf, 
                                                   const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                                                   const VectorType &disp,
                                                   const VectorType &pf,
                                                   const VectorType &adjoint  ) : 
     RefTriangleFENonlinOpIntegrator<ConfiguratorTypePf, KirchhoffLoveBendingEnergyMixedSecondDerivative<MatOptConfType> > ( matOptConf._confpf ),
     KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType> ( matOptConf, xAStorage ),
     _xAStorage(xAStorage),
     _disp ( disp ),
     _pf( matOptConf._confpf, pf ),
     _adjoint( adjoint ) {}
      

  RealType getNonlinearity ( const typename ConfiguratorTypePf::ElementType &El, int QuadPoint ) const {
    RealType NL = this->evaluateBendingMixedSecondDerivativeAtQuadPoint( El, QuadPoint, _disp, _adjoint, _pf.evaluateAtQuadPoint( El, QuadPoint ) );
    NL *= _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    return NL;
  }

};


template<typename MatOptConfType>
class KirchhoffLoveBendingEnergyHessian : 
   public BlockMatrixValuedIntegratorBase< typename MatOptConfType::ConfiguratorType, KirchhoffLoveBendingEnergyHessian<MatOptConfType> >, 
   public KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType> 
{

 protected: 

  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType; 
  typedef typename ConfiguratorType::Matrix33  Matrix33;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::ElementType ElementType; 
  typedef typename ConfiguratorType::VectorType VectorType;
  
 public:
  typedef KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType>  EvaluationHelper;
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;

 protected:
  const MatOptConfType &_matOptConf; 
  const ConfiguratorType &_conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  const VectorType &_pf;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pfDFD;

 public:
    KirchhoffLoveBendingEnergyHessian ( const MatOptConfType &matOpConf, 
                      const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, 
                      const VectorType &pf )
    : BlockMatrixValuedIntegratorBase< ConfiguratorType, KirchhoffLoveBendingEnergyHessian<MatOptConfType> > ( matOpConf._conf ), 
      KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType> ( matOpConf, xAStorage ),
    _matOptConf ( matOpConf ),
    _conf ( matOpConf._conf ),
    _xAStorage ( xAStorage ),
    _pf( pf ),
    _pfDFD ( matOpConf._confpf, pf ) {}

    void prepareLocalMatrix( const ElementType &El, LocalMatrixType (&localMatrix)[3][3]) const {

        for (int argComp = 0; argComp < 3; ++argComp)
          for (int destComp = 0; destComp < 3; ++destComp)
            localMatrix[argComp][destComp].setZero();
            
        const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet(El);
            
        for (int QuadPoint = 0; QuadPoint < _conf.maxNumQuadPoints(); ++QuadPoint) {
            const RealType chi = _matOptConf.approxCharFct_material ( _pfDFD.evaluateAtQuadPoint( El, QuadPoint ), 1., 0. ); //E and nu are in material matrix
            Matrix33 materialMatrixMembrane, materialMatrixBending; 
            this->computeMaterialMatricesWithPf( El, QuadPoint, chi, materialMatrixMembrane, materialMatrixBending );
            materialMatrixMembrane.setZero();
            this->assembleAddLocalMatrix( El, QuadPoint, materialMatrixMembrane, materialMatrixBending, localMatrix, _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * bfs.getWeight(QuadPoint) );
        }
    }
};


template<typename MatOptConfType>
class KirchhoffLoveEnergyHessian : 
   public BlockMatrixValuedIntegratorBase< typename MatOptConfType::ConfiguratorType, KirchhoffLoveEnergyHessian<MatOptConfType> >, 
   public KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType> 
{

 protected: 

  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType; 
  typedef typename ConfiguratorType::Matrix33  Matrix33;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::ElementType ElementType; 
  typedef typename ConfiguratorType::VectorType VectorType;
  
 public:
  typedef KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType>  EvaluationHelper;
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;

 protected:
  const MatOptConfType &_matOptConf; 
  const ConfiguratorType &_conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
  const VectorType &_pf;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pfDFD;

 public:
    KirchhoffLoveEnergyHessian ( const MatOptConfType &matOpConf, 
                      const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, 
                      const VectorType &pf )
    : BlockMatrixValuedIntegratorBase< ConfiguratorType, KirchhoffLoveEnergyHessian<MatOptConfType> > ( matOpConf._conf ), 
      KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType> ( matOpConf, xAStorage ),
    _matOptConf ( matOpConf ),
    _conf ( matOpConf._conf ),
    _xAStorage ( xAStorage ),
    _pf( pf ),
    _pfDFD ( matOpConf._confpf, pf ) {}

    void prepareLocalMatrix( const ElementType &El, LocalMatrixType (&localMatrix)[3][3]) const {

        for (int argComp = 0; argComp < 3; ++argComp)
          for (int destComp = 0; destComp < 3; ++destComp)
            localMatrix[argComp][destComp].setZero();
            
        const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet(El);
            
        for (int QuadPoint = 0; QuadPoint < _conf.maxNumQuadPoints(); ++QuadPoint) {
            const RealType chi = _matOptConf.approxCharFct_material ( _pfDFD.evaluateAtQuadPoint( El, QuadPoint ), 1., 0. ); //E and nu are in material matrix
            Matrix33 materialMatrixMembrane, materialMatrixBending; 
            this->computeMaterialMatricesWithPf( El, QuadPoint, chi, materialMatrixMembrane, materialMatrixBending );
            this->assembleAddLocalMatrix( El, QuadPoint, materialMatrixMembrane, materialMatrixBending, localMatrix, _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * bfs.getWeight(QuadPoint) );
        }
    }
};



template<typename MatOptConfType>
class KirchhoffLoveEnergyOp : public LinElastEnergyOpInterface<MatOptConfType, KirchhoffLoveEnergyOp<MatOptConfType> > {

 protected: 

  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType; 
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::MaskType MaskType;
 public:
  typedef KirchhoffLoveEnergy_EvaluationHelper<MatOptConfType>  EvaluationHelper;
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
  
  const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage; 
  const VectorType &_pf; 
  
  
  KirchhoffLoveEnergyOp ( const MatOptConfType &matOpConf, 
                          const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage, 
                          const VectorType &pf )
    : LinElastEnergyOpInterface<MatOptConfType, KirchhoffLoveEnergyOp<MatOptConfType> > ( matOpConf ),
    _xAStorage( xAStorage ), _pf( pf ) {}

   void evaluateElasticEnergies( const VectorType &disp, RealType &membraneEnergy, RealType &bendingEnergy, RealType &totalEnergy  ) const {
      KirchhoffLoveMembraneEnergy<MatOptConfType> membraneEnergyOp( this->_matOptConf, this->_xAStorage, disp, this->_pf );
      membraneEnergy = 0.; membraneEnergyOp.assembleAdd( membraneEnergy );
      KirchhoffLoveBendingEnergy<MatOptConfType> bendingEnergyOp( this->_matOptConf, this->_xAStorage, disp, this->_pf );
      bendingEnergy = 0.; bendingEnergyOp.assembleAdd( bendingEnergy );
      totalEnergy = membraneEnergy + bendingEnergy;
  }
  
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    KirchhoffLoveEnergyHessian<MatOptConfType> ( this->_matOptConf, this->_xAStorage, this->_pf ).assembleTripletListDirichlet( tripletListMasked, boundaryMask, Factor );   
  }

};



#endif // __LINEARDEFORMATIONENERGIESWITHMATERIAL_H
