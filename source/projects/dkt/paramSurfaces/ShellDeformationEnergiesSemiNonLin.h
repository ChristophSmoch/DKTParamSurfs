#ifndef __SEMINONLINEARDEFORMATIONENERGIESWITHMATERIAL_H
#define __SEMINONLINEARDEFORMATIONENERGIESWITHMATERIAL_H

# include <energyDefines.h>

#include "ShellDeformationEnergiesInterfaces.h"
#include "ShellDeformationEnergiesNonLin.h"



//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
//                                 FEM based energies
//*************************************************************************************
//*************************************************************************************
//*************************************************************************************


//-----------------------------------------------------------------
//-----------------------------------------------------------------
//              Bending (Relative shape Operator ) energy and derivatives
//              E_bend(x_B) = \int_chart sqrt( gA ) |S_rel|^2
//              where S_rel = gAinv (SB - SA)
//-----------------------------------------------------------------
//-----------------------------------------------------------------

template <typename MatOptConfType>
class SemiNonlinearBendingEnergy
: public RefTriangleIntegrator<typename MatOptConfType::ConfiguratorType, SemiNonlinearBendingEnergy<MatOptConfType> >
{
 protected :
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::InitType MeshType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const MatOptConfType& _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage, &_xBStorage;
  const MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xABStorage;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  const RealType _factorBendingEnergy;
  const RealType _thicknessHard, _thicknessSoft;

  public:
    SemiNonlinearBendingEnergy ( const MatOptConfType &matOptConf,
                             const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                             const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                             const MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xABStorage,
                             const VectorType &pf,
                             const RealType factorBendingEnergy  ) :
     RefTriangleIntegrator <ConfiguratorType, SemiNonlinearBendingEnergy<MatOptConfType> > ( matOptConf._conf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage), _xBStorage ( xBStorage ), _xABStorage (xABStorage),
     _pf( matOptConf._confpf, pf ),
    _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
      const RealType materialFactor = 1./24. * chi * _factorBendingEnergy * deltaSqr;
      RealType temp = 0.0;
      Matrix22 tempMat;
      for( int l = 0; l < 3; ++l){
        tempMat = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[l] * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) ;
        temp += pesopt::ddProd<RealType,Matrix22> (tempMat,  _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[l] );
        // cout << "temp=" << temp << endl;
        // for( int i = 1; i < 2; ++i)
        //   for( int j = 1; j < 2; ++j){
        //     temp += tempMat(i,j);
        //   }
      }

      return materialFactor * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint ) * ( temp - 2 * pesopt::ddProd<RealType,Matrix22>( _xBStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) , _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) ) ) ;
      // return materialFactor * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint ) * (  (- 2) * pesopt::ddProd<RealType,Matrix22>( _xBStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) , _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) ) ) ;
      // return materialFactor * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint ) * ( temp  ) ;
    }

};



template<typename MatOptConfType>
class SemiNonlinearBendingEnergyGradient_Part1 :
public RefTriangleMVDiff2OpIntegrator<typename MatOptConfType::ConfiguratorType, SemiNonlinearBendingEnergyGradient_Part1<MatOptConfType> >
{
 protected :
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::InitType MeshType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const MatOptConfType& _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage, &_xBStorage;
  const MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xABStorage;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  const RealType _factorBendingEnergy;
  const RealType _thicknessHard, _thicknessSoft;

 public:
    SemiNonlinearBendingEnergyGradient_Part1 ( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xABStorage,
                              const VectorType &pf,
                              const RealType factorBendingEnergy  ) :
     RefTriangleMVDiff2OpIntegrator<ConfiguratorType, SemiNonlinearBendingEnergyGradient_Part1<MatOptConfType> > ( matOptConf._conf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _xABStorage ( xABStorage ),
     _pf( matOptConf._confpf, pf ),
     _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
     _factorBendingEnergy ( factorBendingEnergy ),
     _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ) {}

     void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, Tensor322Type &NL) const {

        // material factors
        const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
        const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
        const RealType materialFactor = 1./24. * chi * _factorBendingEnergy * deltaSqr;
        //Zweite Ableitungstest
        for ( int l = 0; l<3; ++l){
          NL[l] = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint) * _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint)[l] * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint);
          // NL[l].setZero();
          NL[l] -= _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint)[l] * (_xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint) * _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint));
        }
        // NL.setZero();
        NL *= materialFactor * 2.0 * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );
   }
};


template<typename MatOptConfType>
class SemiNonlinearBendingEnergyGradient_Part2 :
public RefTriangleMVDiffOpIntegrator< typename MatOptConfType::ConfiguratorType, SemiNonlinearBendingEnergyGradient_Part2<MatOptConfType> >
{
 protected :
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::InitType MeshType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const MatOptConfType& _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage, &_xBStorage;
  const MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xABStorage;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  const RealType _factorBendingEnergy;
  const RealType _thicknessHard, _thicknessSoft;

  public:
    SemiNonlinearBendingEnergyGradient_Part2 ( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xABStorage,
                              const VectorType &pf,
                              const RealType factorBendingEnergy  ) :
     RefTriangleMVDiffOpIntegrator<ConfiguratorType, SemiNonlinearBendingEnergyGradient_Part2<MatOptConfType> > ( matOptConf._conf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _xABStorage ( xABStorage ),
     _pf( matOptConf._confpf, pf ),
     _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
     _factorBendingEnergy ( factorBendingEnergy ),
     _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ) {}

     void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, Matrix32 &NL) const{


       // material factors
       const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
       const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
       const RealType materialFactor = 1./24. * chi * _factorBendingEnergy * deltaSqr;
       //
       //
       // //Zweite Ableitungstest
       // DomVecType gAInvGuA0tildeDu; gAInvGuA0tildeDu.setZero();
       //
       // for( int m=0; m<2; ++m){
       //   Tensor322Type a0tildedmxB;
       //   for ( int l = 0; l<3; ++l ) {
       //     a0tildedmxB [l] = _xAStorage.getSemiNonlinIsometry_a0tilde(El.getGlobalElementIdx(), QuadPoint) * _xBStorage.getGradient(El.getGlobalElementIdx(), QuadPoint)(l,m);
       //   }
       //   gAInvGuA0tildeDu(m) = _xABStorage.getGAInvSemiNonlinearityB(El.getGlobalElementIdx(), QuadPoint).ddprod(a0tildedmxB);
       // }
       DomVecType temp;
       temp.setZero();
       Matrix22 gAInvhAgAInv = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint) * _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint);
       for ( int i = 0; i<2; ++i){
         for ( int m = 0; m<3; ++m ) {
           temp(i) +=  pesopt::ddProd<RealType,Matrix22> ( gAInvhAgAInv, _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint)[m] ) * _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint)(m,i) ;
         }
       }
      //  derivativetest
       NL = _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint) * (_xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint) * temp).transpose();

       NL *= materialFactor * 2.0 * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );
      //  NL.setZero();

}
};



template<typename MatOptConfType>
class SemiNonlinearBendingEnergyGradient {
 protected :
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  const MatOptConfType& _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage, &_xBStorage;
  const MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xABStorage;
  const VectorType &_pf;
  const RealType _factorBendingEnergy;

 public:
    SemiNonlinearBendingEnergyGradient ( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xABStorage,
                              const VectorType &pf,
                              const RealType factorBendingEnergy  ) :
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _xABStorage ( xABStorage ),
     _pf( pf ),
    _factorBendingEnergy ( factorBendingEnergy ) {}

  void assembleAdd( VectorType &Deriv ) const {
    SemiNonlinearBendingEnergyGradient_Part1<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _xABStorage, _pf, _factorBendingEnergy ).assembleAdd( Deriv );
    SemiNonlinearBendingEnergyGradient_Part2<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _xABStorage, _pf, _factorBendingEnergy ).assembleAdd( Deriv );
  }
};




template<typename MatOptConfType>
class SemiNonlinearBendingEnergySubHessian_PartDiff1 :
public RefTriangleFELinAsymMatrixWeightedStiffIntegrator<typename MatOptConfType::ConfiguratorType, SemiNonlinearBendingEnergySubHessian_PartDiff1<MatOptConfType> >
{
 protected :
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::InitType MeshType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const MatOptConfType& _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage, &_xBStorage;
  const MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xABStorage;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  const RealType _factorBendingEnergy;
  const RealType _thicknessHard, _thicknessSoft;
  const int _k, _l;

 public:
  SemiNonlinearBendingEnergySubHessian_PartDiff1 ( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xABStorage,
                              const VectorType &pf,
                              const RealType factorBendingEnergy,
                              const int k, const int l  ) :
     RefTriangleFELinAsymMatrixWeightedStiffIntegrator<ConfiguratorType, SemiNonlinearBendingEnergySubHessian_PartDiff1<MatOptConfType> > ( matOptConf._conf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _xABStorage ( xABStorage),
     _pf( matOptConf._confpf, pf ),
    _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ),
    _k(k), _l(l) {}
    // _k(l), _l(k) {}

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, const int &QuadPoint, Matrix22 &Matrix ) const{

     const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
     const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
     const RealType materialFactor = 1./24. * chi * _factorBendingEnergy * deltaSqr;

     Matrix22 gAInvhAgAInv = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );

     Matrix22 Matrix_Part1 = _xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
     Matrix_Part1 *= pesopt::ddProd<RealType,Matrix22>( gAInvhAgAInv, _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[_l] ) * _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )( _k );

     Matrix22 Matrix_Part2 = _xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
    //  ghg = gAInvhAgAInv
     DomVecType ghgHessBGradB;
     ghgHessBGradB.setZero();
     for(int p=0; p<2; ++p){
       for (int m=0; m<3; ++m){
         ghgHessBGradB(p) += pesopt::ddProd<RealType,Matrix22>( gAInvhAgAInv, _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[m] ) * _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(m,p);
       }
     }
     Matrix_Part2 *= (-1.) * _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )( _k ) * ghgHessBGradB.transpose() * _xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint ).row( _l).transpose();

     Matrix22 Matrix_Part3;
     Matrix_Part3 = (_xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * ghgHessBGradB ) * (_xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint ).row( _k).transpose()).transpose();
     Matrix_Part3 *= (-1.) * _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )( _l );

     Matrix22 Matrix_Part4;
     Matrix_Part4 = (_xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint ).row( _l).transpose()) * (_xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * ghgHessBGradB ).transpose() ;
     Matrix_Part4 *= (-1.) * _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )( _k );

      //! sum up
      Matrix.setZero();
      //Zweite Ableitungstest
      Matrix += Matrix_Part1;
      Matrix += Matrix_Part2;
      Matrix += Matrix_Part3;
      Matrix += Matrix_Part4;


      // cout << "_k = " << _k << endl;
      // cout << "_l = " << _l << endl;
      // cout << "norm(Matrix - Matrix^T)^2 = " << (Matrix - Matrix.transpose()).squaredNorm() << endl;
      // cout << "norm(Matrix_Part221 - Matrix_Part221^T)^2 = " << (Matrix_Part221 - Matrix_Part221.transpose()).squaredNorm() << endl;
      // cout << "norm(Matrix_Part21 - Matrix_Part222^T)^2 = " << (Matrix_Part21 - Matrix_Part222.transpose()).squaredNorm() << endl;
      // cout << "norm(Matrix_Part233 - Matrix_Part32^T)^2 = " << (Matrix_Part233 - Matrix_Part32.transpose()).squaredNorm() << endl;
      // Matrix = Matrix.transpose();


      Matrix *= 2.0 * materialFactor * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );
      // Matrix.setZero();
    }
 };


// Part 2.1 = Part 2.2
template<typename MatOptConfType>
class SemiNonlinearBendingEnergySubHessian_PartDiffMixed1 :
public RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator<typename MatOptConfType::ConfiguratorType, SemiNonlinearBendingEnergySubHessian_PartDiffMixed1 <MatOptConfType> >
{
 protected :
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::InitType MeshType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const MatOptConfType& _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage, &_xBStorage;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  const RealType _factorBendingEnergy;
  const RealType _thicknessHard, _thicknessSoft;
  const int _k, _l;

 public:
  SemiNonlinearBendingEnergySubHessian_PartDiffMixed1 ( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const VectorType &pf,
                              const RealType factorBendingEnergy,
                              const int k, const int l  ) :
     RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator<ConfiguratorType, SemiNonlinearBendingEnergySubHessian_PartDiffMixed1<MatOptConfType> > ( matOptConf._conf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _pf( matOptConf._confpf, pf ),
    _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ),
    _k(k), _l(l) {}

     inline void getCoeffMatrixAndVector ( const typename ConfiguratorType::ElementType &El,
                                           const int &QuadPoint,
                                           DomVecType &Vector, Matrix22 &Matrix ) const{
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
      const RealType materialFactor = 1./24. * chi * _factorBendingEnergy * deltaSqr;


       Matrix = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
       Vector = _xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint) * _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint).row( _k ).transpose() * _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint)( _l );

       Vector *= 2.0 * materialFactor * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );
      //  Matrix.setZero();
      //  Vector.setZero();
     }
 };



// // Part 1.3 = Part 1.2^T
template<typename MatOptConfType>
class SemiNonlinearBendingEnergySubHessian_PartDiffMixed2 :
public RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator_Tensor<typename MatOptConfType::ConfiguratorType, SemiNonlinearBendingEnergySubHessian_PartDiffMixed2 <MatOptConfType>>
{
 protected :
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::InitType MeshType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::Tensor222Type Tensor222Type;
  typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const MatOptConfType& _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage, &_xBStorage;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  const RealType _factorBendingEnergy;
  const RealType _thicknessHard, _thicknessSoft;
  const int _k, _l;

 public:
  SemiNonlinearBendingEnergySubHessian_PartDiffMixed2 ( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const VectorType &pf,
                              const RealType factorBendingEnergy,
                              const int k, const int l  ) :
     RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator_Tensor<ConfiguratorType, SemiNonlinearBendingEnergySubHessian_PartDiffMixed2<MatOptConfType> > ( matOptConf._conf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _pf( matOptConf._confpf, pf ),
    _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ),
    _k(k), _l(l) {}

     inline void getCoeffTensor ( const typename ConfiguratorType::ElementType &El,
                                  const int QuadPoint,
                                  Tensor222Type &Tensor ) const {
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
      const RealType materialFactor = 1./24. * chi * _factorBendingEnergy * deltaSqr;

      //Zweite Ableitungstest
        Matrix22 gAInv2A1Tilde;
        gAInv2A1Tilde = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) *  _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a1tilde( El.getGlobalElementIdx(), QuadPoint );
        Matrix22 gAInv2A2Tilde;
        gAInv2A2Tilde = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) *  _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a2tilde( El.getGlobalElementIdx(), QuadPoint );
        for( int p=0; p <2 ; ++p )
          for( int q=0; q<2 ; ++q ){
            if (_l == _k){
              Tensor.set(p, q, 0, gAInv2A1Tilde(p,q));
              Tensor.set(p, q, 1, gAInv2A2Tilde(p,q));

              // Tensor.set(0, p, q, gAInv2A1Tilde(p,q));
              // Tensor.set(1, p, q, gAInv2A2Tilde(p,q));
            }
            else {
              Tensor.set( p, q, 0, 0.0);
              Tensor.set( p, q, 1, 0.0);

              // Tensor.set( 0, p, q, 0.0);
              // Tensor.set( 1, p, q, 0.0);
            }
          }
        Tensor.setZero();
        // Tensor *= -1. * 2. * materialFactor * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );


     }
 };




// Part 1.1
template<typename MatOptConfType>
class SemiNonlinearBendingEnergySubHessian_PartDiff2 :
public RefTriangleFELinMatrixMatrixDiff2Integrator<typename MatOptConfType::ConfiguratorType, SemiNonlinearBendingEnergySubHessian_PartDiff2<MatOptConfType> >
{
 protected :
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::InitType MeshType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const MatOptConfType& _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage, &_xBStorage;
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  const RealType _factorBendingEnergy;
  const RealType _thicknessHard, _thicknessSoft;
  const int _k, _l;

 public:
  SemiNonlinearBendingEnergySubHessian_PartDiff2 ( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const VectorType &pf,
                              const RealType factorBendingEnergy,
                              const int k, const int l  ) :
     RefTriangleFELinMatrixMatrixDiff2Integrator<ConfiguratorType, SemiNonlinearBendingEnergySubHessian_PartDiff2<MatOptConfType> > ( matOptConf._conf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _pf( matOptConf._confpf, pf ),
    _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ),
    _k(k), _l(l) {}

   inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, const int QuadPoint,  Matrix22 &Matrix1, Matrix22 &Matrix2 ) const{
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
      const RealType materialFactor = 1./24. * chi * _factorBendingEnergy * deltaSqr;

      // Zweite Ableitungstest
      if( _l == _k){
        Matrix1 = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        Matrix2 = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        //
        // for (int i = 0; i<2; ++i){
        //   for (int j = 0; j<2; ++j){
        //       Matrix(i,j) = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint )(i,j) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint )(i,j);
        //   }
        //   Matrix(i,2) =  _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint )(i,i) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint )(0,1);
        //   Matrix(2,i) = Matrix(i,2);
        // }
        // Matrix(2,2) =(_xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint )(0,0) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint )(1,1) + _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint )(0,1) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint )(0,1) );
        // Matrix = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
      }
      else{
        Matrix1.setZero();
        Matrix2.setZero();
      }
      // Matrix.setZero();
      Matrix1 *=  2. * materialFactor * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );

   }
 };



template<typename MatOptConfType>
class SemiNonlinearBendingEnergySubHessian {

 protected :
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const MatOptConfType& _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage, &_xBStorage;
  const MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xABStorage;
  const VectorType &_pf;
  const RealType _factorBendingEnergy;
  const int _k, _l;

 public:
  SemiNonlinearBendingEnergySubHessian( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xABStorage,
                              const VectorType &pf,
                              const RealType factorBendingEnergy,
                              const int k, const int l  ) :
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _xABStorage ( xABStorage),
     _pf( pf ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _k(k), _l(l) {}

     void assembleTripletList ( std::vector<TripletType> &tripletList ) const{

       //Part 1.1
       std::vector<TripletType> tripletList_Part1;
       SemiNonlinearBendingEnergySubHessian_PartDiff2<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _pf, _factorBendingEnergy, _k, _l).assembleTripletList( tripletList_Part1, 1. );

       //Parts 2.1 + 2.2 + 2.3.2 + 2.3.3 + 2.4 +  3.2 + 3.3
       std::vector<TripletType> tripletList_Part2;
       SemiNonlinearBendingEnergySubHessian_PartDiff1<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _xABStorage, _pf, _factorBendingEnergy, _k, _l).assembleTripletList( tripletList_Part2, 1. );

       //Zweite Ableitungstest
       //Part 1.3 = 3.1^T
       std::vector<TripletType> tripletList_Part3;
       SemiNonlinearBendingEnergySubHessian_PartDiffMixed2<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _pf, _factorBendingEnergy, _k, _l).assembleTripletList( tripletList_Part3, 2. );

       //Part 1.2 = 2.3.1
       std::vector<TripletType> tripletList_Part4;
       SemiNonlinearBendingEnergySubHessian_PartDiffMixed1<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _pf, _factorBendingEnergy, _k, _l).assembleTripletList( tripletList_Part4, 2. );

       tripletList.reserve( tripletList_Part1.size() + tripletList_Part2.size() + tripletList_Part3.size() + tripletList_Part4.size() );
       for( int i=0; i<tripletList_Part1.size(); ++i ) tripletList.push_back( tripletList_Part1[i] );
       for( int i=0; i<tripletList_Part2.size(); ++i ) tripletList.push_back( tripletList_Part2[i] );
       for( int i=0; i<tripletList_Part3.size(); ++i ) tripletList.push_back( tripletList_Part3[i] );
       for( int i=0; i<tripletList_Part4.size(); ++i ) tripletList.push_back( tripletList_Part4[i] );

     }

 };





template <typename MatOptConfType>
class SemiNonlinearBendingEnergyOp {

protected :
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::InitType MeshType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const MatOptConfType& _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage;
  const VectorType& _xA;
  const VectorType& _pf;
  const RealType _factorBendingEnergy;

public:

  SemiNonlinearBendingEnergyOp ( const MatOptConfType& matOptConf,
                             const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                             const VectorType& pf,
                             const RealType factorBendingEnergy )
  : _matOptConf( matOptConf ), _xAStorage(xAStorage), _xA(xAStorage.getDofs()), _pf( pf ),
    _factorBendingEnergy ( factorBendingEnergy ) {}

public:

  void evaluateEnergy( const VectorType& Displacement, RealType& Dest ) const {
      Dest = 0.;
      VectorType xB ( Displacement.size() ); xB = _xA + Displacement;
      DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xBStorage ( _matOptConf._conf, xB, 3 );
      MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xABStorage ( _matOptConf._conf, _xAStorage, xBStorage, 3 );
      SemiNonlinearBendingEnergy<MatOptConfType> ( _matOptConf, _xAStorage, xBStorage, xABStorage, _pf, _factorBendingEnergy ).assembleAdd( Dest );
  }

  void evaluateStressOnElements( const VectorType &Displacement, VectorType &bendingStressVec ) const{
      bendingStressVec.setZero();
      VectorType xB ( Displacement.size() ); xB = _xA + Displacement;
      DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xBStorage ( _matOptConf._conf, xB, 3 );
      MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xABStorage ( _matOptConf._conf, _xAStorage, xBStorage, 3 );
      SemiNonlinearBendingEnergy<MatOptConfType> ( _matOptConf, _xAStorage, xBStorage, xABStorage, _pf, _factorBendingEnergy ).assembleOnElements( bendingStressVec );
  }

  void evaluateGradient(const VectorType& Displacement, VectorType& Dest) const {
      Dest.setZero();
      VectorType xB ( Displacement.size() ); xB = _xA + Displacement;
      DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xBStorage ( _matOptConf._conf, xB, 3 );

      MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xABStorage ( _matOptConf._conf, _xAStorage, xBStorage, 3 );

      SemiNonlinearBendingEnergyGradient<MatOptConfType> ( _matOptConf, _xAStorage, xBStorage, xABStorage, _pf, _factorBendingEnergy ).assembleAdd( Dest );
  }


  void evaluateHessian( const VectorType& Displacement, SparseMatrixType& Hessian ) const {
    Hessian.setZero();
    std::vector<TripletType> tripletList;

    this->assembleTripletListHessian( Displacement, tripletList );
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
    Hessian.makeCompressed();
  }

  //
  void assembleTripletListHessian( const VectorType& Displacement, std::vector<TripletType>& tripletList, const RealType Fac = 1. ) const  {
      VectorType xB ( Displacement.size() ); xB = _xA + Displacement;
      DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xBStorage ( _matOptConf._conf, xB, 3 );
      MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xABStorage ( _matOptConf._conf, _xAStorage, xBStorage, 3 );
      const int numGlobalDofs = _matOptConf._conf.getNumGlobalDofs();
      for( int k=0; k<3; ++k )
          for( int l=0; l<3; ++l ){
              std::vector<TripletType> tripletList_kl;
              SemiNonlinearBendingEnergySubHessian<MatOptConfType> ( _matOptConf, _xAStorage, xBStorage, xABStorage, _pf, _factorBendingEnergy, k, l ).assembleTripletList( tripletList_kl );
              for( int i=0; i<tripletList_kl.size(); ++i ){
                tripletList.push_back( TripletType( tripletList_kl[i].row() + k * numGlobalDofs, tripletList_kl[i].col() + l * numGlobalDofs, 0.5 * tripletList_kl[i].value() ) );
                tripletList.push_back( TripletType( tripletList_kl[i].col() + l * numGlobalDofs, tripletList_kl[i].row() + k * numGlobalDofs, 0.5 * tripletList_kl[i].value() ) );
              }
          }
  }

};





//!==========================================================================================================
//! Nonlinear Membrane + Bending
//!==========================================================================================================


template<typename MatOptConfType>
class SemiNonlinearMembraneBendingEnergyOp
: public NonLinElastEnergyOpInterface<MatOptConfType, SemiNonlinearMembraneBendingEnergyOp<MatOptConfType> > {

 protected:
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::MaskType MaskType;
public:
  static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;

protected :
    const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
    const VectorType &_xA;
    const VectorType &_pf;
    const MaskType &_mask;
    const VectorType &_force;
    const RealType _factorMembraneEnergy, _factorBendingEnergy;

public:

  SemiNonlinearMembraneBendingEnergyOp ( const MatOptConfType &matOptConf,
                                     const MaskType &mask,
                                     const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                                     const VectorType &pf,
                                     const VectorType &force  )
    : NonLinElastEnergyOpInterface<MatOptConfType, SemiNonlinearMembraneBendingEnergyOp<MatOptConfType> > ( matOptConf ),
    _xAStorage ( xAStorage ), _xA( xAStorage.getDofs() ), _pf( pf ),  _mask( mask ),  _force( force ),
    _factorMembraneEnergy ( matOptConf._materialInfo._factorMembraneEnergy ), _factorBendingEnergy ( matOptConf._materialInfo._factorBendingEnergy ) {}

  void evaluateElasticEnergies( const VectorType &disp, RealType &membraneEnergy, RealType &bendingEnergy, RealType &potEnergy, RealType &totalEnergy  ) const {
      NonlinearMembraneEnergyOp<MatOptConfType> membraneEnergyOp( this->_matOptConf, this->_xAStorage, this->_pf, _factorMembraneEnergy );
      membraneEnergy = 0.; membraneEnergyOp.evaluateEnergy( disp, membraneEnergy );
      SemiNonlinearBendingEnergyOp<MatOptConfType> bendingEnergyOp( this->_matOptConf, this->_xAStorage, this->_pf, _factorBendingEnergy );
      bendingEnergy = 0.; bendingEnergyOp.evaluateEnergy( disp, bendingEnergy );
      potEnergy = _force.dot(disp);
      totalEnergy = membraneEnergy + bendingEnergy - potEnergy;
  }

  void evaluateJacobian( const VectorType &disp, VectorType &Deriv ) const override{
      NonlinearMembraneEnergyOp<MatOptConfType> membraneEnergyOp( this->_matOptConf, this->_xAStorage, this->_pf, _factorMembraneEnergy );
      VectorType membraneEnergyGrad ( Deriv.size() ); membraneEnergyOp.evaluateGradient( disp, membraneEnergyGrad );
      SemiNonlinearBendingEnergyOp<MatOptConfType> bendingEnergyOp( this->_matOptConf, this->_xAStorage, this->_pf, _factorBendingEnergy );
      VectorType bendingEnergyGrad ( Deriv.size() ); bendingEnergyOp.evaluateGradient( disp, bendingEnergyGrad );
      Deriv = membraneEnergyGrad + bendingEnergyGrad - _force;
      for( int i = 0; i < _mask.size(); ++i ){
        if ( _mask[i] ){
            for( int comp=0; comp<3; ++comp ) Deriv[i + comp * _mask.size()] = 0.0;}
      }
  }

//   void evaluateMixedSecondDerivative(const VectorType& disp, const VectorType& adjointSol, VectorType& Deriv) const {
//       NonlinearMembraneEnergyOp<MatOptConfType> membraneEnergyOp( this->_matOptConf, this->_xAStorage, this->_pf, _factorMembraneEnergy );
//       VectorType membraneEnergyGrad ( Deriv.size() ); membraneEnergyOp.evaluateMixedSecondDerivative( disp, adjointSol, membraneEnergyGrad );
//       SemiNonlinearBendingEnergyOp<MatOptConfType> bendingEnergyOp( this->_matOptConf, this->_xAStorage, this->_pf, _factorBendingEnergy );
//       VectorType bendingEnergyGrad ( Deriv.size() ); bendingEnergyOp.evaluateMixedSecondDerivative( disp, adjointSol, bendingEnergyGrad );
//       Deriv = membraneEnergyGrad + bendingEnergyGrad;
//   }

  void evaluateStressOnElements( const VectorType &disp, VectorType &membraneStressVec, VectorType &bendingStressVec, VectorType &totalStressVec ) const{
      membraneStressVec.setZero(); bendingStressVec.setZero(); totalStressVec.setZero();
      NonlinearMembraneEnergyOp<MatOptConfType> membraneEnergyOp( this->_matOptConf, this->_xAStorage, this->_pf, _factorMembraneEnergy );
      membraneEnergyOp.evaluateStressOnElements( disp, membraneStressVec );
      SemiNonlinearBendingEnergyOp<MatOptConfType> bendingEnergyOp( this->_matOptConf, this->_xAStorage, this->_pf, _factorBendingEnergy );
      bendingEnergyOp.evaluateStressOnElements( disp, bendingStressVec );
      totalStressVec = membraneStressVec + bendingStressVec;
    }

  void assembleTripletListHessian ( const VectorType &disp,
                                    std::vector<TripletType> & tripletListMasked,
                                    const RealType Factor = 1.0 ) const {

   std::vector<TripletType> tripletListMembrane, tripletListBending;
   NonlinearMembraneEnergyOp<MatOptConfType>( this->_matOptConf, this->_xAStorage, this->_pf, _factorMembraneEnergy ).assembleTripletListHessian( disp, tripletListMembrane, Factor );
   SemiNonlinearBendingEnergyOp<MatOptConfType>( this->_matOptConf, this->_xAStorage, this->_pf, _factorBendingEnergy ).assembleTripletListHessian( disp, tripletListBending, Factor );

   tripletListMasked.reserve( tripletListMembrane.size() + tripletListBending.size() );

   const int numDofs = _mask.size();
   for( unsigned iter=0; iter < tripletListMembrane.size(); ++iter ){
      if( (_mask[tripletListMembrane[iter].row() % numDofs ]) || (_mask[tripletListMembrane[iter].col() % numDofs ] ) ){
       //Boundary node!
      } else {
        tripletListMasked.push_back( tripletListMembrane[iter] );
      }
   }

   for( unsigned iter=0; iter < tripletListBending.size(); ++iter ){
      if( (_mask[tripletListBending[iter].row() % numDofs] ) || (_mask[tripletListBending[iter].col() % numDofs]) ){
       //Boundary node!
      } else {
        tripletListMasked.push_back( tripletListBending[iter] );
      }
   }

   for ( int i = 0; i < numDofs; ++i ){
      if ( _mask[i] ){ for ( int Comp = 0; Comp < 3; ++Comp ) tripletListMasked.push_back( TripletType( i + Comp * numDofs, i + Comp * numDofs, 1.0 ) );}
   }

  }

};



#endif
