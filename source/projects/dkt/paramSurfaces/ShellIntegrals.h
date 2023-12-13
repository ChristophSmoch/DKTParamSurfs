#ifndef __SHELLINTEGRALS_H
#define __SHELLINTEGRALS_H

#include <pesopt_DKT.h>
# include <linearSystemSolver.h>

    
//! Area(v) = \int_S 1 dH^2 = \sum_El \int_T  sqrt( det( Dx^T Dx ) ) \circ x dL^2 , where T is reference triangle
template <typename ConfiguratorType>
class AreaOnManifold :
public RefTriangleIntegrator < ConfiguratorType, AreaOnManifold<ConfiguratorType> >{
  
protected :
  typedef typename ConfiguratorType::RealType RealType;
  
  const ConfiguratorType & _conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
  
  public:
    AreaOnManifold ( const ConfiguratorType & conf,
                     const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage  ) : 
     RefTriangleIntegrator<ConfiguratorType, AreaOnManifold<ConfiguratorType> > ( conf ),
     _conf ( conf ), _xAStorage ( xAStorage ) {}

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
     return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
  }
};
    
    
    
// int_M chi(m) u_i u_j = int_chart sqrt(det_g) chi(m) u_i u_j
template<typename MatOpType>
class MassEnergy_Hessian : public RefTriangleFELinWeightedMassIntegrator<typename MatOpType::ConfiguratorType, MassEnergy_Hessian <MatOpType> >
{
protected:
    typedef typename MatOpType::ConfiguratorType ConfiguratorType;
    typedef typename MatOpType::ConfiguratorTypePf ConfiguratorTypePf;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::TripletType TripletType;
    typedef typename ConfiguratorType::MaskType MaskType;
    typedef typename ConfiguratorType::RealType   RealType;
    typedef typename ConfiguratorType::Matrix22   Matrix22;
    typedef typename ConfiguratorType::Matrix32   Matrix32;
    
    const MatOpType &_matOptConf;
    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
    const DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
     
    const Material<RealType> & _HardMaterial;
    const Material<RealType> & _SoftMaterial;
     
public:
    MassEnergy_Hessian(const MatOpType &matOpConf,  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> & xAStorage, const VectorType &pf ) : 
    RefTriangleFELinWeightedMassIntegrator< ConfiguratorType, MassEnergy_Hessian<MatOpType>>(matOpConf._conf),
    _matOptConf( matOpConf ),
    _xAStorage ( xAStorage ), _pf( matOpConf._confpf, pf),
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ) {}

    inline RealType getNonlinearity(const typename ConfiguratorType::ElementType &El,  const int &QuadPoint ) const {
        RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
        return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * chi;
    }
}; 

    

// int_Surface |f|^p = int_chart sqrt det_g |f|^p
template<typename ConfiguratorTypeChart, typename ConfiguratorTypeFct>
class ScalarFctLpNormToP :
public RefTriangleIntegrator<ConfiguratorTypeFct, ScalarFctLpNormToP<ConfiguratorTypeChart,ConfiguratorTypeFct> >
{
  protected:

       typedef typename ConfiguratorTypeFct::RealType RealType; 
       typedef typename ConfiguratorTypeFct::ElementType ElementType; 
       typedef typename ConfiguratorTypeFct::VectorType VectorType;
       
       const ConfiguratorTypeFct &_conf;   
       const DiscreteVectorFunctionStorage<ConfiguratorTypeChart,FirstOrder> &_xAStorage;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypeFct> _f;
       const RealType _p;
  public:
    ScalarFctLpNormToP ( const ConfiguratorTypeFct &conf, 
                const DiscreteVectorFunctionStorage<ConfiguratorTypeChart,FirstOrder> &xAStorage,
                const VectorType &f,
                const RealType p = 2 ) :
      RefTriangleIntegrator<ConfiguratorTypeFct, ScalarFctLpNormToP<ConfiguratorTypeChart,ConfiguratorTypeFct>> (conf),
      _conf ( conf ),
      _xAStorage ( xAStorage ),
      _f( _conf, f ),
      _p(p) {}

    RealType evaluateIntegrand ( const typename ConfiguratorTypeFct::ElementType &El, int QuadPoint ) const{
      RealType F = _f.evaluateAtQuadPoint( El, QuadPoint );
      RealType aux = std::pow( std::abs( F ), _p );
      return aux * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};


// int_Surface |f|^p = int_chart sqrt det_g |f|^p
template<typename ConfiguratorTypeChart, typename ConfiguratorTypeFct>
class LpNormToP :
public RefTriangleIntegrator<ConfiguratorTypeFct, LpNormToP<ConfiguratorTypeChart,ConfiguratorTypeFct> >
{
  protected:

       typedef typename ConfiguratorTypeFct::RealType RealType; 
       typedef typename ConfiguratorTypeFct::Point3DType Point3DType;
       typedef typename ConfiguratorTypeFct::ElementType ElementType; 
       typedef typename ConfiguratorTypeFct::VectorType VectorType;
       
       const ConfiguratorTypeFct &_conf;   
       const DiscreteVectorFunctionStorage<ConfiguratorTypeChart,FirstOrder> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorTypeFct> _f;
       const RealType _p;
  public:
    LpNormToP ( const ConfiguratorTypeFct &conf, 
                const DiscreteVectorFunctionStorage<ConfiguratorTypeChart,FirstOrder> &xAStorage,
                const VectorType &f,
                const RealType p = 2 ) :
      RefTriangleIntegrator<ConfiguratorTypeFct, LpNormToP<ConfiguratorTypeChart,ConfiguratorTypeFct>> (conf),
      _conf ( conf ),
      _xAStorage ( xAStorage ),
      _f( _conf, f, 3 ),
      _p(p) {}

    RealType evaluateIntegrand ( const typename ConfiguratorTypeFct::ElementType &El, int QuadPoint ) const{
      Point3DType F; _f.evaluateAtQuadPoint( El, QuadPoint, F );
      RealType aux = 0.;
      for( int comp=0; comp<3; ++comp) aux += std::pow( std::abs(F[comp]), _p );
      return aux * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};


// int chi(m) |f|^p
template<typename MatOpType>
class PfWeightedLpNormToP :
public RefTriangleIntegrator<typename MatOpType::ConfiguratorType, PfWeightedLpNormToP<MatOpType> >
{
  protected:

       typedef typename MatOpType::ConfiguratorType ConfiguratorType;
       typedef typename MatOpType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
       
       const MatOpType &_matOpConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _disp;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> & _HardMaterial, _SoftMaterial;
       const RealType _p;
  public:
    PfWeightedLpNormToP ( const MatOpType &matOpConf, 
                const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                const VectorType &displacement,
                const VectorType &pf,
                const RealType p = 2 ) :
      RefTriangleIntegrator<ConfiguratorType, PfWeightedLpNormToP<MatOpType>> (matOpConf._conf),
      _matOpConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _p(p) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
      Point3DType Disp; _disp.evaluateAtQuadPoint( El, QuadPoint, Disp );
      RealType aux = 0.;
      for( int comp=0; comp<3; ++comp) aux += std::pow( std::abs(Disp[comp]), _p );
      RealType chi = _matOpConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return aux * chi * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};


// \| f \|_L^\infty
template<typename MatOpType>
class LInftyNorm :
public RefTriangleIntegrator<typename MatOpType::ConfiguratorType, LInftyNorm<MatOpType> >
{
  protected:

       typedef typename MatOpType::ConfiguratorType ConfiguratorType;
       typedef typename MatOpType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
       
       const MatOpType &_matOpConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _disp;
//        DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
//        const Material<RealType> & _HardMaterial, _SoftMaterial;
       mutable RealType _LInftyNorm;
  public:
    LInftyNorm ( const MatOpType &matOpConf, 
                const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                const VectorType &displacement
//                 const VectorType &pf 
               ) :
      RefTriangleIntegrator<ConfiguratorType, LInftyNorm<MatOpType>> (matOpConf._conf),
      _matOpConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( matOpConf._conf, displacement, 3 ),
//       _pf ( matOpConf._confpf, pf ),
//       _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _LInftyNorm ( 0. ) {}

  void assembleAdd( RealType& norm ) const {   
      RealType energyDummy = 0.;
      RefTriangleIntegrator< ConfiguratorType, LInftyNorm<MatOpType> >::assembleAdd( energyDummy );
      norm = _LInftyNorm;
  }
      
  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
      Point3DType Disp; _disp.evaluateAtQuadPoint( El, QuadPoint, Disp );
      RealType dispLInfty = 0.;
      for( int comp=0; comp<3; ++comp){
          RealType abs = std::abs(Disp[comp]);
          if( abs > dispLInfty ) dispLInfty = abs;
      }
//       RealType chi = _matOpConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
//       RealType aux = dispLInfty * chi * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
      RealType aux = dispLInfty;
      if( aux > _LInftyNorm ) _LInftyNorm = aux;
      return aux;
    }
};
    
    
    
// int_M chi(m) |Df|^p
template<typename MatOpType>
class LpNormOfGradToP :
public RefTriangleIntegrator<typename MatOpType::ConfiguratorType, LpNormOfGradToP<MatOpType> >
{
  protected:

       typedef typename MatOpType::ConfiguratorType ConfiguratorType;
       typedef typename MatOpType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix32 Matrix32;
       typedef typename ConfiguratorType::Matrix33 Matrix33;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
       
       const MatOpType &_matOpConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _disp;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> & _HardMaterial, _SoftMaterial;
       const RealType _p;
  public:
    LpNormOfGradToP ( const MatOpType &matOpConf, 
                const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                const VectorType &displacement,
                const VectorType &pf,
                const RealType p = 2 ) :
      RefTriangleIntegrator<ConfiguratorType, LpNormOfGradToP<MatOpType>> (matOpConf._conf),
      _matOpConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _p(p) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
      Matrix32 dX; _disp.evaluateGradientAtQuadPoint( El, QuadPoint, dX );
      Matrix32 DxAGAinv = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint) * _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
      Matrix33 grad = DxAGAinv * dX.transpose();
      RealType aux = 0.;
      for( int comp=0; comp<3; ++comp)
          for(int i=0; i<3; ++i )
              aux += std::pow( std::abs(grad(comp,i)), _p );
      RealType chi = _matOpConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return aux * chi * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};


// int_M chi(m) |D_h f|^p
template<typename MatOpType>
class LpNormOfApproxGradToP :
public RefTriangleIntegrator<typename MatOpType::ConfiguratorType, LpNormOfApproxGradToP<MatOpType> >
{
  protected:

       typedef typename MatOpType::ConfiguratorType ConfiguratorType;
       typedef typename MatOpType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix32 Matrix32;
       typedef typename ConfiguratorType::Matrix33 Matrix33;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
       
       const MatOpType &_matOpConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _disp;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> & _HardMaterial, _SoftMaterial;
       const RealType _p;
  public:
    LpNormOfApproxGradToP ( const MatOpType &matOpConf, 
                const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                const VectorType &displacement,
                const VectorType &pf,
                const RealType p = 2 ) :
      RefTriangleIntegrator<ConfiguratorType, LpNormOfApproxGradToP<MatOpType>> (matOpConf._conf),
      _matOpConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _p(p) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
      Matrix32 dX; _disp.evaluateApproxGradientAtQuadPoint( El, QuadPoint, dX );
      Matrix32 DxAGAinv = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint) * _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
      Matrix33 grad = DxAGAinv * dX.transpose();
      
      RealType aux = 0.;
      for( int comp=0; comp<3; ++comp) 
          for( int i=0; i<3; ++i )
              aux += std::pow( std::abs(grad(comp,i)), _p );
      RealType chi = _matOpConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return aux * chi * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};







// int_M chi(m) |D_h f - D f|^p
template<typename MatOpType>
class LpNormOfDiffGradVsApproxGradToP :
public RefTriangleIntegrator<typename MatOpType::ConfiguratorType, LpNormOfDiffGradVsApproxGradToP<MatOpType> >
{
  protected:

       typedef typename MatOpType::ConfiguratorType ConfiguratorType;
       typedef typename MatOpType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix32 Matrix32;
       typedef typename ConfiguratorType::Matrix33 Matrix33;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
       
       const MatOpType &_matOpConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _disp;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> & _HardMaterial, _SoftMaterial;
       const RealType _p;
  public:
    LpNormOfDiffGradVsApproxGradToP ( const MatOpType &matOpConf, 
                const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                const VectorType &displacement,
                const VectorType &pf,
                const RealType p = 2 ) :
      RefTriangleIntegrator<ConfiguratorType, LpNormOfDiffGradVsApproxGradToP<MatOpType>> (matOpConf._conf),
      _matOpConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _p(p) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
      Matrix32 DxAGAinv = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint) * _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
     
      Matrix32 dX; _disp.evaluateGradientAtQuadPoint( El, QuadPoint, dX );
      Matrix33 grad = DxAGAinv * dX.transpose();
      
      Matrix32 dXApprox; _disp.evaluateApproxGradientAtQuadPoint( El, QuadPoint, dXApprox );
      Matrix33 gradApprox = DxAGAinv * dXApprox.transpose();
      
      RealType aux = 0.;
      for( int comp=0; comp<3; ++comp) 
          for( int i=0; i<3; ++i )
              aux += std::pow( std::abs(grad(comp,i) - gradApprox(comp,i)), _p );
      RealType chi = _matOpConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return aux * chi * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};







// int_M chi(m) |D D_h f|^p
template<typename MatOpType>
class LpNormOfGradOfApproxGradToP :
public RefTriangleIntegrator<typename MatOpType::ConfiguratorType, LpNormOfGradOfApproxGradToP<MatOpType> >
{
  protected:

       typedef typename MatOpType::ConfiguratorType ConfiguratorType;
       typedef typename MatOpType::ConfiguratorTypePf ConfiguratorTypePf;
       typedef typename ConfiguratorType::RealType RealType; 
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix33 Matrix33;
       typedef typename ConfiguratorType::ElementType ElementType; 
       typedef typename ConfiguratorType::VectorType VectorType;
       
       const MatOpType &_matOpConf;   
       const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
       DKTFEVectorFunctionEvaluator<ConfiguratorType> _disp;
       DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
       const Material<RealType> & _HardMaterial, _SoftMaterial;
       const RealType _p;
  public:
    LpNormOfGradOfApproxGradToP ( const MatOpType &matOpConf, 
                const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                const VectorType &displacement,
                const VectorType &pf,
                const RealType p = 2 ) :
      RefTriangleIntegrator<ConfiguratorType, LpNormOfGradOfApproxGradToP<MatOpType>> (matOpConf._conf),
      _matOpConf ( matOpConf ),
      _xAStorage ( xAStorage ),
      _disp( matOpConf._conf, displacement, 3 ),
      _pf ( matOpConf._confpf, pf ),
      _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ),
      _p(p) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
      Matrix33 GradDisp; _disp.evaluateApproxHessianAsVecAtQuadPoint( El, QuadPoint, GradDisp );
      RealType aux = 0.;
      for( int comp=0; comp<3; ++comp) 
          for( int i=0; i<3; ++i )
              aux += std::pow( std::abs(GradDisp(comp,i)), _p );
      RealType chi = _matOpConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      return aux * chi * _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};
    
    
    
    
//! \brief Assembles \f$ \sum_{EL} \int_T \sqrt{\det g} <v_i, v_j> \f$
template<typename ConfiguratorType>
class ScalarValuedMassOp : public RefTriangleFELinWeightedMassIntegrator<ConfiguratorType, ScalarValuedMassOp<ConfiguratorType> >
{
protected:
    typedef typename ConfiguratorType::RealType   RealType;
    typedef typename ConfiguratorType::Matrix22   Matrix22;
    typedef typename ConfiguratorType::Matrix32   Matrix32;
    
    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;

public:
    ScalarValuedMassOp (const ConfiguratorType &conf, const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> & xAStorage ) : 
    RefTriangleFELinWeightedMassIntegrator< ConfiguratorType, ScalarValuedMassOp<ConfiguratorType>>(conf),
    _xAStorage ( xAStorage ) {}
    
    inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
       return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
    
};

//! \brief Assembles \f$ \sum_{EL} \int_T \sqrt{\det g} < g^{-1} Dv_i, Dv_j > \f$
template<typename ConfiguratorType>
class ScalarValuedStiffOp : public RefTriangleFELinAsymMatrixWeightedStiffIntegrator<ConfiguratorType, ScalarValuedStiffOp<ConfiguratorType> >
{
protected:
    typedef typename ConfiguratorType::RealType   RealType;
    typedef typename ConfiguratorType::Matrix22   Matrix22;
    typedef typename ConfiguratorType::Matrix32   Matrix32;
    
    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;

public:
    ScalarValuedStiffOp (const ConfiguratorType &conf, const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> & xAStorage ) : 
    RefTriangleFELinAsymMatrixWeightedStiffIntegrator< ConfiguratorType, ScalarValuedStiffOp<ConfiguratorType>>(conf),
    _xAStorage ( xAStorage ) {}
    
    inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, Matrix22 & NL ) const {
       NL = _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint);
       NL *= _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
    
};
    

//! \brief Assembles \f$ \sum_{EL} \int_T \sqrt{\det g} <v_i, v_j> \f$
template<typename ConfiguratorType>
class VectorValuedMassOp : public RefTriangleFELinWeightedBlockMassIntegrator<ConfiguratorType, VectorValuedMassOp <ConfiguratorType> >
{
protected:
    typedef typename ConfiguratorType::RealType   RealType;
    typedef typename ConfiguratorType::Matrix22   Matrix22;
    typedef typename ConfiguratorType::Matrix32   Matrix32;
    
    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;

public:
    VectorValuedMassOp(const ConfiguratorType &conf,  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> & xAStorage ) : 
    RefTriangleFELinWeightedBlockMassIntegrator< ConfiguratorType, VectorValuedMassOp<ConfiguratorType>>(conf),
    _xAStorage ( xAStorage ) {}
    
    inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
       return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
    
};



//! \brief Assembles \f$ \sum_{EL} \int_T \sqrt{\det g} F u_i \f$
// for F wrt ConfiguratorTypeForce, u_i base functions wrt ConfiguratorType
template< typename ConfiguratorType, typename ConfiguratorTypeForce >
class ScalarValuedForceOp_MixedConfigurators :
 public RefTriangleFENonlinOpIntegrator<ConfiguratorType, ScalarValuedForceOp_MixedConfigurators<ConfiguratorType,ConfiguratorTypeForce> >
{
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::TangentVecType TangentVecType;

    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
    const DKTFEScalarFunctionEvaluator<ConfiguratorTypeForce> _forceDFD;
    
public:
    ScalarValuedForceOp_MixedConfigurators( const ConfiguratorType &conf, const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage, 
                                            const ConfiguratorTypeForce &confForce, const VectorType &force )
    : RefTriangleFENonlinOpIntegrator <ConfiguratorType, ScalarValuedForceOp_MixedConfigurators<ConfiguratorType,ConfiguratorTypeForce> > (conf),
      _xAStorage(xAStorage),
      _forceDFD( confForce, force ) {}

    RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {   
      RealType NL = _forceDFD.evaluateAtQuadPoint( El, QuadPoint );
      NL *= _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
      return NL;
    }
};


//! \brief Assembles \f$ \sum_{EL} \int_T \sqrt{\det g} g^-1 grad(F) \cdot grad(u_i) \f$
// for F wrt ConfiguratorTypeForce, u_i base functions wrt ConfiguratorType
template< typename ConfiguratorType, typename ConfiguratorTypeForce >
class ScalarValuedGradForceOp_MixedConfigurators :
 public RefTriangleFENonlinDiffOpIntegrator<ConfiguratorType, ScalarValuedGradForceOp_MixedConfigurators<ConfiguratorType,ConfiguratorTypeForce> >
{
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::TangentVecType TangentVecType;
    typedef typename ConfiguratorType::DomVecType DomVecType;

    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
    const DKTFEScalarFunctionEvaluator<ConfiguratorTypeForce> _forceDFD;
    
public:
    ScalarValuedGradForceOp_MixedConfigurators( const ConfiguratorType &conf, 
                                                const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage, 
                                                const ConfiguratorTypeForce &confForce, const VectorType &force )
    : RefTriangleFENonlinDiffOpIntegrator <ConfiguratorType, ScalarValuedGradForceOp_MixedConfigurators<ConfiguratorType,ConfiguratorTypeForce> > (conf),
      _xAStorage(xAStorage),
      _forceDFD( confForce, force ) {}

    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, DomVecType &NL ) const {   
      DomVecType dF; _forceDFD.evaluateGradientAtQuadPoint( El, QuadPoint, dF );
      NL = _xAStorage.getFirstFFInv(El.getGlobalElementIdx(),QuadPoint) * dF;
      NL *= _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};


//! Assembly operator for constant RHS: \int <f, v_i> = \int sqrt(g) <f \circ x_A, v_i \circ x_A>
template< typename ConfiguratorType, typename ConfiguratorTypeForce >
class VectorValuedForceOp_MixedConfigurators :
    public RefTriangleFENonlinVectorOpIntegrator<ConfiguratorType, VectorValuedForceOp_MixedConfigurators<ConfiguratorType,ConfiguratorTypeForce> >
{
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::TangentVecType TangentVecType;

    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
    const DKTFEVectorFunctionEvaluator<ConfiguratorTypeForce> _forceDFD;
    
public:
    VectorValuedForceOp_MixedConfigurators(const ConfiguratorType &conf, const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage, 
                                   const ConfiguratorTypeForce &confForce, const VectorType &force)
    : RefTriangleFENonlinVectorOpIntegrator <ConfiguratorType, VectorValuedForceOp_MixedConfigurators<ConfiguratorType,ConfiguratorTypeForce> > (conf),
      _xAStorage(xAStorage),
      _forceDFD( confForce, force, 3 ) {}

    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, TangentVecType &NL) const {   
      NL.setZero();
      _forceDFD.evaluateAtQuadPoint( El, QuadPoint, NL );
      NL *= _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint);
    }
};






//aproximate P1 function by DKT function (scalar valued) wrt L2Norm
template<typename ConfiguratorTypeDKT,typename ConfiguratorTypeP1>
void ApproximateP1ByDKTFunction_wrtL2Norm( const ConfiguratorTypeP1 &confP1, 
                                 const typename ConfiguratorTypeP1::VectorType &P1Fct,
                                 const ConfiguratorTypeDKT &confDKT,
                                 const DiscreteVectorFunctionStorage<ConfiguratorTypeDKT,FirstOrder> &xAStorage, 
                                 typename ConfiguratorTypeDKT::VectorType &DKTFct,
                                 const typename ConfiguratorTypeDKT::MaskType &boundaryMask,
                                 bool useBoundaryMask ){
    
       typedef typename ConfiguratorTypeDKT::RealType RealType;
       typedef typename ConfiguratorTypeDKT::VectorType VectorType;
       typedef typename ConfiguratorTypeDKT::SparseMatrixType SparseMatrixType;
       typedef typename ConfiguratorTypeDKT::DTContainer DataTypeContainer;
    
       SparseMatrixType MassMatrix( DKTFct.size(), DKTFct.size() );
       ScalarValuedMassOp<ConfiguratorTypeDKT> MassOp( confDKT, xAStorage );
       if( useBoundaryMask ) MassOp.assembleDirichlet ( MassMatrix, boundaryMask );
       else MassOp.assemble( MassMatrix );
       MassMatrix.makeCompressed();  
       
       VectorType RHS( DKTFct.size() ); RHS.setZero();
       ScalarValuedForceOp_MixedConfigurators<ConfiguratorTypeDKT, ConfiguratorTypeP1> RHSOp( confDKT, xAStorage, confP1, P1Fct );   
       if( useBoundaryMask ) RHSOp.assembleDirichlet( RHS, boundaryMask );
       else RHSOp.assembleAdd( RHS );
       
    
       IterativeLinearSystemSolver<DataTypeContainer, EigenConjugateGradient_NEW > iterativeLinearSystemSolver ( 1.e-14, 1., "massop" );
       iterativeLinearSystemSolver.solve( MassMatrix, DKTFct, RHS );
}


//aproximate P1 function by DKT function (scalar valued) wrt H1-Norm
template<typename ConfiguratorTypeDKT,typename ConfiguratorTypeP1>
void ApproximateP1ByDKTFunction_wrtH1Norm( const ConfiguratorTypeP1 &confP1, 
                                           const typename ConfiguratorTypeP1::VectorType &P1Fct,
                                           const ConfiguratorTypeDKT &confDKT,
                                           const DiscreteVectorFunctionStorage<ConfiguratorTypeDKT,FirstOrder> &xAStorage, 
                                           typename ConfiguratorTypeDKT::VectorType &DKTFct,
                                           const typename ConfiguratorTypeDKT::MaskType &boundaryMask,
                                           bool useBoundaryMask ){
    
       typedef typename ConfiguratorTypeDKT::RealType RealType;
       typedef typename ConfiguratorTypeDKT::VectorType VectorType;
       typedef typename ConfiguratorTypeDKT::SparseMatrixType SparseMatrixType;
       typedef typename ConfiguratorTypeDKT::DTContainer DataTypeContainer;
    
       SparseMatrixType MassMatrix( DKTFct.size(), DKTFct.size() );
       ScalarValuedMassOp<ConfiguratorTypeDKT> MassOp( confDKT, xAStorage );
       if( useBoundaryMask ) MassOp.assembleDirichlet ( MassMatrix, boundaryMask );
       else MassOp.assemble( MassMatrix );
       
       SparseMatrixType StiffMatrix( DKTFct.size(), DKTFct.size() );
       ScalarValuedStiffOp<ConfiguratorTypeDKT> StiffOp( confDKT, xAStorage );
       if( useBoundaryMask ) StiffOp.assembleDirichlet ( StiffMatrix, boundaryMask );
       else MassOp.assemble( StiffMatrix );
       
       SparseMatrixType SystemMatrix( DKTFct.size(), DKTFct.size() );
       SystemMatrix = MassMatrix + StiffMatrix;
       SystemMatrix.makeCompressed();
       
       VectorType RHS_Mass( DKTFct.size() ); RHS_Mass.setZero();
       ScalarValuedForceOp_MixedConfigurators<ConfiguratorTypeDKT, ConfiguratorTypeP1> RHSMassOp( confDKT, xAStorage, confP1, P1Fct );   
       if( useBoundaryMask ) RHSMassOp.assembleDirichlet( RHS_Mass, boundaryMask );
       else RHSMassOp.assembleAdd( RHS_Mass );
       
       VectorType RHS_Stiff( DKTFct.size() ); RHS_Stiff.setZero();
       ScalarValuedGradForceOp_MixedConfigurators<ConfiguratorTypeDKT, ConfiguratorTypeP1> RHSStiffOp( confDKT, xAStorage, confP1, P1Fct );   
       if( useBoundaryMask ) RHSStiffOp.assembleDirichlet( RHS_Stiff, boundaryMask );
       else RHSStiffOp.assembleAdd( RHS_Stiff );
       
       VectorType RHS( DKTFct.size() ); RHS = RHS_Mass + RHS_Stiff;
       
       IterativeLinearSystemSolver<DataTypeContainer,  EigenConjugateGradient_NEW > iterativeLinearSystemSolver ( 1.e-14, 1., "ApproximateP1ByDKTFunction_wrtH1Norm" );
       iterativeLinearSystemSolver.solve( SystemMatrix, DKTFct, RHS );
}


#endif
