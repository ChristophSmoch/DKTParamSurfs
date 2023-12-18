#ifndef __SEMINONLINEARDEFORMATIONENERGIESWITHMATERIAL_H
#define __SEMINONLINEARDEFORMATIONENERGIESWITHMATERIAL_H

# include <energyDefines.h>

#include "../elastDeform/ShellDeformationEnergiesInterfaces.h"
#include "../elastDeform/ShellDeformationEnergiesNonLin.h"



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
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  const RealType _factorBendingEnergy;
  const RealType _thicknessHard, _thicknessSoft;

  public:
    SemiNonlinearBendingEnergy ( const MatOptConfType &matOptConf,
                             const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                             const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                             const VectorType &pf,
                             const RealType factorBendingEnergy  ) :
     RefTriangleIntegrator <ConfiguratorType, SemiNonlinearBendingEnergy<MatOptConfType> > ( matOptConf._conf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage), _xBStorage ( xBStorage ),
     _pf( matOptConf._confpf, pf ),
    _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
      // material factors
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
      const RealType materialFactor = 1./24. * chi * _factorBendingEnergy * deltaSqr;


      RealType aux = 0.;


      for ( int l = 0; l<3; ++l ) {
         Matrix22 mat_temp; mat_temp.setZero();

         Matrix22 gAinvD2xB = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[l];
         mat_temp += gAinvD2xB;

         // mat_temp -= _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )[l] * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getSemiNonlinIsometry_a0tilde( El.getGlobalElementIdx(), QuadPoint );
         //
         // mat_temp -= _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, 0) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getSemiNonlinIsometry_a1tilde( El.getGlobalElementIdx(), QuadPoint );
         //
         // mat_temp -= _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, 1) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getSemiNonlinIsometry_a2tilde( El.getGlobalElementIdx(), QuadPoint );

         // mat_temp += _xBStorage.getFirstFF( El.getGlobalElementIdx(), QuadPoint );
         // mat_temp -= _xAStorage.getFirstFF( El.getGlobalElementIdx(), QuadPoint );


         mat_temp -= _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )[l] * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a0tilde( El.getGlobalElementIdx(), QuadPoint );

         mat_temp -= _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, 0) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a1tilde( El.getGlobalElementIdx(), QuadPoint );

         mat_temp -= _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, 1) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a2tilde( El.getGlobalElementIdx(), QuadPoint );

         aux +=  mat_temp.squaredNorm( );
      }

      return materialFactor * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint ) * aux;
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
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  const RealType _factorBendingEnergy;
  const RealType _thicknessHard, _thicknessSoft;

 public:
    SemiNonlinearBendingEnergyGradient_Part1 ( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const VectorType &pf,
                              const RealType factorBendingEnergy  ) :
     RefTriangleMVDiff2OpIntegrator<ConfiguratorType, SemiNonlinearBendingEnergyGradient_Part1<MatOptConfType> > ( matOptConf._conf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _pf( matOptConf._confpf, pf ),
    _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ) {}

   void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, Tensor322Type &NL) const {

      // material factors
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
      const RealType materialFactor = 1./24. * chi * _factorBendingEnergy * deltaSqr;


      for ( int l = 0; l<3; ++l ) {
         Matrix22 mat_temp; mat_temp.setZero();

         mat_temp +=  _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[l];

         mat_temp -= _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )[l] * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a0tilde( El.getGlobalElementIdx(), QuadPoint );

         mat_temp -= _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, 0) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a1tilde( El.getGlobalElementIdx(), QuadPoint );

         mat_temp -= _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, 1) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a2tilde( El.getGlobalElementIdx(), QuadPoint );

         for( int i=0; i<2; ++i )
           for( int j=0; j<2; ++j ) {
             NL.set( l, i, j, mat_temp(i,j));
           }
      }
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

      cout << "gradient part 2: getNonlinearity " << endl;

      // TEST MixedDiscreteVectorFunctionStorage
      // cout <<  "gAinvSecondFFB = " <<  _xABStorage.getGAInvSecondFFB(El.getGlobalElementIdx(),  QuadPoint ) << endl;
      cout << "test wieso geht das nicht..." << endl;
      cout << "G[u] = " << _xABStorage.getSemiNonlinearityB(El.getGlobalElementIdx(), QuadPoint ) << endl;

      // material factors
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
      const RealType materialFactor = 1./24. * chi * _factorBendingEnergy * deltaSqr;

      cout << "test 1" << endl;
      VectorType gAInvGuA0tildeDu; gAInvGuA0tildeDu.setZero();
      for( int m=0; m<2; ++m){

        RealType aux = 0.0;
        for ( int l = 0; l<3; ++l ) {
           Matrix22 mat_temp; mat_temp.setZero();

           mat_temp +=  _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[l];

           mat_temp -= _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )[l] * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a0tilde( El.getGlobalElementIdx(), QuadPoint );

           mat_temp -= _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, 0) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a1tilde( El.getGlobalElementIdx(), QuadPoint );

           mat_temp -= _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, 1) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a2tilde( El.getGlobalElementIdx(), QuadPoint );

           RealType w;
           w = pesopt::ddProd<RealType,Matrix22> ( mat_temp, _xAStorage.getSemiNonlinIsometry_a0tilde( El.getGlobalElementIdx(), QuadPoint ) ) ;
           w *=  _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, m);

           aux += w;

           }
        gAInvGuA0tildeDu[m] = aux;


      }
      NL = _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint ) * (_xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * gAInvGuA0tildeDu).transpose();


      Matrix32 NL1;

      cout << "test 2" << endl;
      for ( int l = 0; l<3; ++l ) {
         Matrix22 mat_temp; mat_temp.setZero();

         mat_temp +=  _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[l];

         mat_temp -= _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )[l] * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a0tilde( El.getGlobalElementIdx(), QuadPoint );

         mat_temp -= _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, 0) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a1tilde( El.getGlobalElementIdx(), QuadPoint );

         mat_temp -= _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, 1) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a2tilde( El.getGlobalElementIdx(), QuadPoint );

         RealType gAGuA1tilde;
         gAGuA1tilde = pesopt::ddProd<RealType,Matrix22> ( mat_temp, _xAStorage.getSemiNonlinIsometry_a1tilde( El.getGlobalElementIdx(), QuadPoint ));

         NL1(l,0) = gAGuA1tilde;
    }
    cout << "test 3" << endl;
      for ( int l = 0; l<3; ++l ) {
         Matrix22 mat_temp; mat_temp.setZero();

         mat_temp +=  _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[l];

         mat_temp -= _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )[l] * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a0tilde( El.getGlobalElementIdx(), QuadPoint );

         mat_temp -= _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, 0) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a1tilde( El.getGlobalElementIdx(), QuadPoint );

         mat_temp -= _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint )(l, 1) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getSemiNonlinIsometry_a2tilde( El.getGlobalElementIdx(), QuadPoint );

         RealType gAGuA2tilde;
         gAGuA2tilde = pesopt::ddProd<RealType,Matrix22> ( mat_temp, _xAStorage.getSemiNonlinIsometry_a2tilde( El.getGlobalElementIdx(), QuadPoint ));

         NL1(l,1) =  gAGuA2tilde;
    }
    NL -= NL1;
    NL *= materialFactor * 2.0 * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );

    cout << "finished gradient part 2: getNonlinearity " << endl;

///////////////////////////////////

    //   Matrix22 gAInvGAInvHBMinHA = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * ( _xBStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) - _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) );
    //   TangentVecType v;
    //   for (int comp = 0; comp < 3; ++comp){
    //     v(comp) = pesopt::ddProd<RealType,Matrix22> ( gAInvGAInvHBMinHA, _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[comp] );
    //   }
    //   Matrix33 normalTensorV = _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint ) * v.transpose();
    //   NL = normalTensorV * _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint );
    //
    //   NL *= materialFactor * 2.0 * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );
    // }
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
    SemiNonlinearBendingEnergyGradient_Part1<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _pf, _factorBendingEnergy ).assembleAdd( Deriv );
    SemiNonlinearBendingEnergyGradient_Part2<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _xABStorage, _pf, _factorBendingEnergy ).assembleAdd( Deriv );
  }
};






template<typename MatOptConfType>
class SemiNonlinearBendingEnergyMixedSecondDerivative_Part1
: public RefTriangleFENonlinOpIntegrator<typename MatOptConfType::ConfiguratorTypePf, SemiNonlinearBendingEnergyMixedSecondDerivative_Part1<MatOptConfType> >
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
  DKTFEVectorFunctionEvaluator<ConfiguratorType> _adjoint;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  const RealType _factorBendingEnergy;
  const RealType _thicknessHard, _thicknessSoft;

  public:
    SemiNonlinearBendingEnergyMixedSecondDerivative_Part1 ( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const VectorType &pf,
                              const VectorType &adjoint,
                              const RealType factorBendingEnergy  ) :
     RefTriangleFENonlinOpIntegrator<ConfiguratorTypePf, SemiNonlinearBendingEnergyMixedSecondDerivative_Part1<MatOptConfType> > ( matOptConf._confpf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _pf( matOptConf._confpf, pf ),
     _adjoint( matOptConf._conf, adjoint, 3 ),
    _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ){}

  RealType getNonlinearity ( const typename ConfiguratorTypePf::ElementType &El, int QuadPoint ) const {
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      const RealType Dchi = _matOptConf.approxCharFct_material_Derivative ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
      const RealType DdeltaSqr = _matOptConf.approxCharFct_thicknessSqr_Derivative ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );

      Matrix22 mat;
      mat = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) *  ( _xBStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) -_xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) );
      Tensor322Type NL;
      for( int i=0; i<3; ++i )
        for( int j=0; j<2; ++j )
          for( int k=0; k<2; ++k )
            NL.set( i, j, k, _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )[i] * mat(j,k) );
      NL *= ( deltaSqr * Dchi + DdeltaSqr * chi ) * _factorBendingEnergy * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );

      Tensor322Type HessianAdjoint; _adjoint.evaluateApproxHessianSymAtQuadPoint( El, QuadPoint, HessianAdjoint );
      return 1./12. * NL.ddprod(HessianAdjoint);
  }

};



template<typename MatOptConfType>
class SemiNonlinearBendingEnergyMixedSecondDerivative_Part2
: public RefTriangleFENonlinOpIntegrator<typename MatOptConfType::ConfiguratorTypePf, SemiNonlinearBendingEnergyMixedSecondDerivative_Part2<MatOptConfType> >
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
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> _pf;
  DKTFEVectorFunctionEvaluator<ConfiguratorType> _adjoint;
  const Material<RealType> &_HardMaterial, &_SoftMaterial;
  const RealType _factorBendingEnergy;
  const RealType _thicknessHard, _thicknessSoft;

  public:
    SemiNonlinearBendingEnergyMixedSecondDerivative_Part2 ( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const VectorType &pf,
                              const VectorType &adjoint,
                              const RealType factorBendingEnergy  ) :
     RefTriangleFENonlinOpIntegrator<ConfiguratorTypePf, SemiNonlinearBendingEnergyMixedSecondDerivative_Part2<MatOptConfType> > ( matOptConf._confpf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _pf( matOptConf._confpf, pf ),
     _adjoint( matOptConf._conf, adjoint, 3 ),
    _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ){}

  RealType getNonlinearity ( const typename ConfiguratorTypePf::ElementType &El, int QuadPoint ) const {
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      const RealType Dchi = _matOptConf.approxCharFct_material_Derivative ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
      const RealType DdeltaSqr = _matOptConf.approxCharFct_thicknessSqr_Derivative ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );

      Matrix22 gAInvGAInvHBMinHA = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * ( _xBStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) - _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) );
      TangentVecType v;
      for (int comp = 0; comp < 3; ++comp){
        v(comp) = pesopt::ddProd<RealType,Matrix22> ( gAInvGAInvHBMinHA, _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[comp] );
      }
      //normalTensorV.makeTensorProduct ( _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint ) , v);
      Matrix33 normalTensorV = _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint ) * v.transpose();
      Matrix32 NL = normalTensorV * _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint );
      NL *= _xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
      NL *= -1. * ( deltaSqr * Dchi + DdeltaSqr * chi ) * _factorBendingEnergy * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );

      Matrix32 GradientAdjoint; _adjoint.evaluateGradientAtQuadPoint( El, QuadPoint, GradientAdjoint );
      return 1./12. * pesopt::ddProd<RealType, Matrix32 > ( NL, GradientAdjoint );
  }

};





template<typename MatOptConfType>
class SemiNonlinearBendingEnergyMixedSecondDerivative {
 protected :
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfType::ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  const MatOptConfType& _matOptConf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &_xAStorage, &_xBStorage;
  const VectorType &_pf;
  const VectorType &_adjoint;
  const RealType _factorBendingEnergy;

 public:
    SemiNonlinearBendingEnergyMixedSecondDerivative ( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const VectorType &pf,
                              const VectorType &adjoint,
                              const RealType factorBendingEnergy  ) :
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _pf( pf ), _adjoint( adjoint ),
    _factorBendingEnergy ( factorBendingEnergy ) {}

  void assembleAdd( VectorType &Deriv ) const {
    SemiNonlinearBendingEnergyMixedSecondDerivative_Part1<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _pf, _adjoint, _factorBendingEnergy ).assembleAdd( Deriv );
    SemiNonlinearBendingEnergyMixedSecondDerivative_Part2<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _pf, _adjoint, _factorBendingEnergy ).assembleAdd( Deriv );
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
  SemiNonlinearBendingEnergySubHessian_PartDiff1 ( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const VectorType &pf,
                              const RealType factorBendingEnergy,
                              const int k, const int l  ) :
     RefTriangleFELinAsymMatrixWeightedStiffIntegrator<ConfiguratorType, SemiNonlinearBendingEnergySubHessian_PartDiff1<MatOptConfType> > ( matOptConf._conf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _pf( matOptConf._confpf, pf ),
    _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ),
    _k(k), _l(l) {}

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, const int &QuadPoint, Matrix22 &Matrix ) const{

     const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
     const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );

     Matrix33 DxB_gB_Inv_DxBT =_xBStorage.getGradientGInv( El.getGlobalElementIdx(), QuadPoint ) * _xBStorage.getGradient( El.getGlobalElementIdx(), QuadPoint ).transpose();

     Matrix22 gAInv2ShapeTensorDiffMat = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * ( _xBStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) - _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) );

     TangentVecType gAInv2ShapeTensorDiffMatDdotHessianXb;
     for (int comp = 0; comp < 3; ++comp)
        gAInv2ShapeTensorDiffMatDdotHessianXb(comp) = pesopt::ddProd<RealType,Matrix22> ( gAInv2ShapeTensorDiffMat, _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[comp] );

     Tensor322Type gAInv2HessianB;
     for( int r=0; r<3; r++ ){
         Matrix22 tmpMat ( _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) );
         tmpMat *= (_xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint ))[r];
         for( int p=0; p<2; p++ )
           for( int q=0; q<2; q++ )
             gAInv2HessianB.set( r, p, q, tmpMat(p,q) );
      }

     Matrix33 gAInv2HessianxBddHessianxB;
      for( int r=0; r<3; ++r )
        for( int q=0; q<3; ++q )
         gAInv2HessianxBddHessianxB( r, q ) = pesopt::ddProd<RealType,Matrix22>( gAInv2HessianB[r], _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[q] );


     //Part 1.4 //this should be implemented without 4 loops
     Matrix22 Matrix_Part14; Matrix_Part14.setZero();
     for( int p=0; p<2; ++p )
        for( int q=0; q<2; ++q )
          for( int s=0; s<3; ++s )
            for( int v=0; v<3; ++v )
             Matrix_Part14(q,p) += gAInv2HessianxBddHessianxB(v,s) * _xBStorage.getGradientGInv( El.getGlobalElementIdx(), QuadPoint )(v,p) *  _xBStorage.getGradientGInv( El.getGlobalElementIdx(), QuadPoint )(s,q);
     Matrix_Part14 *= _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )(_k) * _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )(_l);

     //Part 2.3
     Matrix22 Matrix_Part23 (_xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ));
     Matrix_Part23 *= -1. * gAInv2ShapeTensorDiffMatDdotHessianXb(_l) *  _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )(_k);

     //Part 2.4.1
     Matrix22 Matrix_Part241; Matrix_Part241.setZero ();
     for (int p = 0; p < 2; ++p)
        for (int q = 0; q < 2; ++q)
          for (int s = 0; s < 3; ++s)
            Matrix_Part241(p,q) += gAInv2ShapeTensorDiffMatDdotHessianXb(s) * _xBStorage.getGradientGInv( El.getGlobalElementIdx(), QuadPoint )(s,q) *  _xBStorage.getGradientGInv( El.getGlobalElementIdx(), QuadPoint )(_l,p);
     Matrix_Part241 *= _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )(_k);

     //Part 2.4.2
     Matrix22 Matrix_Part242 ( _xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) );
     RealType tmpValue = 0.;
     for (int s = 0; s < 3; ++s) tmpValue += gAInv2ShapeTensorDiffMatDdotHessianXb(s) * DxB_gB_Inv_DxBT(s,_l);
     Matrix_Part242 *= _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )(_k) * tmpValue;

     //Part 2.5
     Matrix22 Matrix_Part25; Matrix_Part25.setZero ();
     for (int p = 0; p < 2; ++p)
        for (int q = 0; q < 2; ++q)
           for (int s = 0; s < 3; ++s)
             Matrix_Part25(p,q) +=  gAInv2ShapeTensorDiffMatDdotHessianXb(s) * _xBStorage.getGradientGInv( El.getGlobalElementIdx(), QuadPoint )(s,p) *  _xBStorage.getGradientGInv( El.getGlobalElementIdx(), QuadPoint )(_k,q);
     Matrix_Part25 *= _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )(_l);


      //! sum up
      Matrix.setZero();
      Matrix += Matrix_Part14;
      Matrix += Matrix_Part23;
      Matrix += Matrix_Part241;
      Matrix += Matrix_Part242;
      Matrix += Matrix_Part25;

      Matrix *= 1./12. * chi * _factorBendingEnergy * deltaSqr * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );
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

       Matrix = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * (_xBStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) - _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint ) );

       Vector = _xBStorage.getGradientGInv( El.getGlobalElementIdx(), QuadPoint ).row(_k );
       Vector *= _xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )(_l);
       Vector *=  -1. * 1./12. * chi * _factorBendingEnergy * deltaSqr * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );
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

        Tensor322Type gAInv2HessianB;
        for( int r=0; r<3; ++r ){
            Matrix22 tmpMat ( _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) );
            tmpMat *= (_xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint )[r]);
            for( int p=0; p<2; ++p )
              for( int q=0; q<2; ++q )
                gAInv2HessianB.set( r, p, q, tmpMat(p,q) );
        }
        for( int p=0; p <2 ; ++p )
          for( int q=0; q<2 ; ++q )
            for( int s=0; s <2 ; ++s ){
                RealType tmp = 0.0;
                for( int r=0; r <3 ; ++r ) tmp += gAInv2HessianB.get( r, p, q ) * (_xBStorage.getGradientGInv( El.getGlobalElementIdx(), QuadPoint )(r,s) );
                Tensor.set( p, q, s, tmp );
            }
        Tensor *= (_xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )(_l)) * (_xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )(_k) );
        Tensor *=  -1. * 1./12. * chi * _factorBendingEnergy * deltaSqr * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint ) ;
     }
 };




// Part 1.1
template<typename MatOptConfType>
class SemiNonlinearBendingEnergySubHessian_PartDiff2 :
public RefTriangleFELinMatrixDiff2Integrator<typename MatOptConfType::ConfiguratorType, SemiNonlinearBendingEnergySubHessian_PartDiff2<MatOptConfType> >
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
     RefTriangleFELinMatrixDiff2Integrator<ConfiguratorType, SemiNonlinearBendingEnergySubHessian_PartDiff2<MatOptConfType> > ( matOptConf._conf ),
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _pf( matOptConf._confpf, pf ),
    _HardMaterial ( matOptConf._materialInfo._HardMaterial ), _SoftMaterial ( matOptConf._materialInfo._SoftMaterial ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _thicknessHard( matOptConf._materialInfo._thicknessHard ), _thicknessSoft( matOptConf._materialInfo._thicknessSoft ),
    _k(k), _l(l) {}

   inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, const int QuadPoint,  Matrix22 &Matrix ) const{
      const RealType chi = _matOptConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _HardMaterial.getElastModulus(), _SoftMaterial.getElastModulus() );
      const RealType deltaSqr = _matOptConf.approxCharFct_thicknessSqr ( _pf.evaluateAtQuadPoint( El, QuadPoint ), _thicknessHard, _thicknessSoft );
      Matrix = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint ) * _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
      Matrix *= (_xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )(_k)) *  (_xBStorage.getNormal( El.getGlobalElementIdx(), QuadPoint )(_l));
      Matrix *= chi * 1./12. * _factorBendingEnergy * deltaSqr * _xAStorage.getArea( El.getGlobalElementIdx(), QuadPoint );
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
  const VectorType &_pf;
  const RealType _factorBendingEnergy;
  const int _k, _l;

 public:
  SemiNonlinearBendingEnergySubHessian( const MatOptConfType &matOptConf,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xAStorage,
                              const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> &xBStorage,
                              const VectorType &pf,
                              const RealType factorBendingEnergy,
                              const int k, const int l  ) :
     _matOptConf( matOptConf ),
     _xAStorage(xAStorage),
     _xBStorage ( xBStorage ),
     _pf( pf ),
    _factorBendingEnergy ( factorBendingEnergy ),
    _k(k), _l(l) {}

     void assembleTripletList ( std::vector<TripletType> &tripletList ) const{

       //Part 1.1
       std::vector<TripletType> tripletList_Part1;
       SemiNonlinearBendingEnergySubHessian_PartDiff2<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _pf, _factorBendingEnergy, _k, _l).assembleTripletList( tripletList_Part1, 1. );

       //Parts 1.4, 2.3, 2.4.1, 2.4.2, 2.5
       std::vector<TripletType> tripletList_Part2;
       SemiNonlinearBendingEnergySubHessian_PartDiff1<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _pf, _factorBendingEnergy, _k, _l).assembleTripletList( tripletList_Part2, 1. );

       //Part 1.2 = (1.3)^T
       std::vector<TripletType> tripletList_Part3;
       SemiNonlinearBendingEnergySubHessian_PartDiffMixed2<MatOptConfType> ( _matOptConf, _xAStorage, _xBStorage, _pf, _factorBendingEnergy, _k, _l).assembleTripletList( tripletList_Part3, 2. );

       //Part 2.1 = Part 2.2
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
      SemiNonlinearBendingEnergy<MatOptConfType> ( _matOptConf, _xAStorage, xBStorage, _pf, _factorBendingEnergy ).assembleAdd( Dest );
  }

  void evaluateStressOnElements( const VectorType &Displacement, VectorType &bendingStressVec ) const{
      bendingStressVec.setZero();
      VectorType xB ( Displacement.size() ); xB = _xA + Displacement;
      DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xBStorage ( _matOptConf._conf, xB, 3 );
      SemiNonlinearBendingEnergy<MatOptConfType> ( _matOptConf, _xAStorage, xBStorage, _pf, _factorBendingEnergy ).assembleOnElements( bendingStressVec );
  }

  void evaluateGradient(const VectorType& Displacement, VectorType& Dest) const {
      Dest.setZero();
      VectorType xB ( Displacement.size() ); xB = _xA + Displacement;
      DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xBStorage ( _matOptConf._conf, xB, 3 );

      MixedDiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xABStorage ( _matOptConf._conf, _xAStorage, xBStorage, 3 );

      SemiNonlinearBendingEnergyGradient<MatOptConfType> ( _matOptConf, _xAStorage, xBStorage, xABStorage, _pf, _factorBendingEnergy ).assembleAdd( Dest );
  }

  void evaluateMixedSecondDerivative(const VectorType& Displacement, const VectorType& adjointSol, VectorType& Dest) const {
      Dest.setZero();
      VectorType xB ( Displacement.size() ); xB = _xA + Displacement;
      DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xBStorage ( _matOptConf._conf, xB, 3 );
      SemiNonlinearBendingEnergyMixedSecondDerivative<MatOptConfType> ( _matOptConf, _xAStorage, xBStorage, _pf, adjointSol, _factorBendingEnergy ).assembleAdd( Dest );
  }

  void evaluateHessian( const VectorType& Displacement, SparseMatrixType& Hessian ) const {
    Hessian.setZero();
    std::vector<TripletType> tripletList;
//     tripletList.reserve( 9 * 9 * _matOptConf._conf.getInitializer().getNumElements() );  // TODO
    this->assembleTripletListHessian( Displacement, tripletList );
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
    Hessian.makeCompressed();
  }

  //
  void assembleTripletListHessian( const VectorType& Displacement, std::vector<TripletType>& tripletList, const RealType Fac = 1. ) const  {
      VectorType xB ( Displacement.size() ); xB = _xA + Displacement;
      DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xBStorage ( _matOptConf._conf, xB, 3 );
      const int numGlobalDofs = _matOptConf._conf.getNumGlobalDofs();
      for( int k=0; k<3; ++k )
          for( int l=0; l<3; ++l ){
              std::vector<TripletType> tripletList_kl;
              SemiNonlinearBendingEnergySubHessian<MatOptConfType> ( _matOptConf, _xAStorage, xBStorage, _pf, _factorBendingEnergy, k, l ).assembleTripletList( tripletList_kl );
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

  void evaluateMixedSecondDerivative(const VectorType& disp, const VectorType& adjointSol, VectorType& Deriv) const {
      NonlinearMembraneEnergyOp<MatOptConfType> membraneEnergyOp( this->_matOptConf, this->_xAStorage, this->_pf, _factorMembraneEnergy );
      VectorType membraneEnergyGrad ( Deriv.size() ); membraneEnergyOp.evaluateMixedSecondDerivative( disp, adjointSol, membraneEnergyGrad );
      SemiNonlinearBendingEnergyOp<MatOptConfType> bendingEnergyOp( this->_matOptConf, this->_xAStorage, this->_pf, _factorBendingEnergy );
      VectorType bendingEnergyGrad ( Deriv.size() ); bendingEnergyOp.evaluateMixedSecondDerivative( disp, adjointSol, bendingEnergyGrad );
      Deriv = membraneEnergyGrad + bendingEnergyGrad;
  }

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
