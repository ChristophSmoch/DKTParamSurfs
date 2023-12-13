#ifndef __DEFORMATIONENERGIESINTERFACES_H
#define __DEFORMATIONENERGIESINTERFACES_H

#include <pesopt_IO.h>
#include <pesopt_DKT.h>
#include <energyDefines.h>



enum ElastEnergyType {
  PlateLaplace = 1,
  KirchhoffLove = 3,
  NonLinMemLinPlateLaplace = 21,
  NonLinMemLinKiLo = 23,
  FullNonLinMem = 31,
};


string getElastEnergyString ( const int ElastEnergyType ){
    string ElastEnergyString;
        switch( ElastEnergyType ){
            case  1 : ElastEnergyString = "PlateLaplaceE"; break;
            case  3 : ElastEnergyString = "KiLoE"; break;
            case 21 : ElastEnergyString = "NonLinMemLinPlateLaplaceBendE"; break;
            case 23 : ElastEnergyString = "NonLinMemLinKiLoBendE"; break;
            case 31 : ElastEnergyString = "FullNonLinMemBendE"; break;
            default : throw std::invalid_argument ( pesopt::strprintf ( "wrong ElastEnergyType in File %s at line %d.", __FILE__, __LINE__ ).c_str() );  break;
        }
    return ElastEnergyString;
}


//*************************************************************************************
//*************************************************************************************
//************************************************************************************
//                                  Interfaces:
//*************************************************************************************
//*************************************************************************************
//*************************************************************************************



template <typename MatOptConfType, typename Imp>
class LinElastEnergy_EvaluationHelper{
public:
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
  typedef typename ConfiguratorType::InitType MeshType;

  const MatOptConfType &_matOptConf;
  const ConfiguratorType &_conf;
  const MeshType &_mesh;
  const Material<RealType> & _HardMaterial;  const Material<RealType> & _SoftMaterial;

public:
    LinElastEnergy_EvaluationHelper ( const MatOptConfType &matOpConf ) :
    _matOptConf( matOpConf ), _conf ( matOpConf._conf ), _mesh( _conf.getInitializer() ),
    _HardMaterial ( matOpConf._materialInfo._HardMaterial ), _SoftMaterial ( matOpConf._materialInfo._SoftMaterial ) { }


  virtual ~LinElastEnergy_EvaluationHelper( ) {}

  //-----------------------------------------------------------
  // functions to be provided in derived class
  //-----------------------------------------------------------
  RealType getFactorMembraneEnergy (  ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }

  RealType getFactorBendingEnergy (  ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }

  //! \brief This is for example used in the computation of D^2_{mat,disp} E(pf,u)(hat m, v) (the second partial derivative of the Energy wrt. both the displacment and the material)
  //  given     E(pf,u) = 1/2 int sqrt(det_g) W(u) chi(pf)
  //            -> D^2_{mat,disp} E(pf,u) (hat m, v) = int sqrt(det_g) DW(u)(v) Dchi(m) hat m
  //  here we provide DW(u)(v) Dchi(m)
  RealType evaluateMixedSecondDerivativeAtQuadPoint ( const ElementType &El, const int QuadPoint, const VectorType &u, const VectorType &v, const RealType pf ) const {
      throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }

  RealType evaluateMembraneStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }

  RealType evaluateMembraneStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }

  RealType evaluateBendingStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    return this->asImp().evaluateBendingStress ( El, RefCoord, material, displacement, xA );
  }

  RealType evaluateBendingStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    return this->asImp().evaluateBendingStressAtQuadPoint ( El, QuadPoint, material, displacement );
  }


  //-----------------------------------------------------------
  //
  //-----------------------------------------------------------
  RealType evaluateTotalStress( const ElementType &El, const DomVecType &RefCoord, const VectorType &material, const VectorType &displacement, const VectorType &xA ) const{
      return this->asImp().evaluateMembraneStress( El, RefCoord, material, displacement, xA )
           + this->asImp().evaluateBendingStress( El, RefCoord, material, displacement, xA );
  }

  RealType evaluateTotalStressAtQuadPoint( const ElementType &El, const int QuadPoint, const VectorType &material, const VectorType &displacement ) const{
      return this->asImp().evaluateMembraneStressAtQuadPoint( El, QuadPoint, material, displacement )
           + this->asImp().evaluateBendingStressAtQuadPoint( El, QuadPoint, material, displacement );
  }

  RealType evaluateMembraneStressOnElement( const ElementType &El, const VectorType &material, const VectorType &displacement ) const{
      RealType aux = 0.;
      const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );
      for( int q=0; q < numQuadPoints; ++q ){
         aux += this->asImp().evaluateMembraneStressAtQuadPoint( El, q, material, displacement ) * bfs.getWeight ( q );
      }
      return aux * El.getAreaOfRefTriangle();
  }

  RealType evaluateBendingStressOnElement( const ElementType &El, const VectorType &material, const VectorType &displacement ) const{
      RealType aux = 0.;
      const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );
      for( int q=0; q < numQuadPoints; ++q ){
         aux +=  this->asImp().evaluateBendingStressAtQuadPoint( El, q, material, displacement ) * bfs.getWeight ( q );
      }
      return aux * El.getAreaOfRefTriangle();
  }

  RealType evaluateTotalStressOnElement( const ElementType &El, const VectorType &material, const VectorType &displacement ) const{
      RealType aux = 0.;
      const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );
      for( int q=0; q < numQuadPoints; ++q ){
         aux +=  this->asImp().evaluateMembraneStressAtQuadPoint( El, q, material, displacement ) * bfs.getWeight ( q );
         aux +=  this->asImp().evaluateBendingStressAtQuadPoint( El, q, material, displacement ) * bfs.getWeight ( q );
      }
      return aux * El.getAreaOfRefTriangle();
  }

  void evaluateTotalStressOnElement( const ElementType &El, const VectorType &material, const VectorType &displacement,
                                    RealType &membraneStress, RealType &bendingStress, RealType &totalStress ) const{
      membraneStress = 0.; bendingStress = 0.;
      const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );
      for( int q=0; q < numQuadPoints; ++q ){
         membraneStress +=  this->asImp().evaluateMembraneStressAtQuadPoint( El, q, material, displacement ) * bfs.getWeight ( q );
         bendingStress +=  this->asImp().evaluateBendingStressAtQuadPoint( El, q, material, displacement ) * bfs.getWeight ( q );
      }
      membraneStress *= El.getAreaOfRefTriangle();
      bendingStress *= El.getAreaOfRefTriangle();
      totalStress = membraneStress + bendingStress;
  }

  //computes int_El stress
  void evaluateStressOnElements( const VectorType &material, const VectorType &displacement, const VectorType &xA,
                                 VectorType &membraneStressVec, VectorType &bendingStressVec, VectorType &totalStressVec ) const {
      for( int elementIdx=0; elementIdx<_mesh.getNumTriangs(); ++elementIdx ){
           const ElementType& El ( _mesh.getTriang( elementIdx ) );
           this->evaluateTotalStressOnElement( El, material, displacement, membraneStressVec[elementIdx], bendingStressVec[elementIdx], totalStressVec[elementIdx] );
           membraneStressVec[elementIdx];
           bendingStressVec[elementIdx];
           totalStressVec[elementIdx];
      }
   }

  //computes int_El stress / int_El 1
  void evaluateAveragedStressOnElements( const VectorType &material, const VectorType &displacement, const VectorType &xA,
                                         VectorType &membraneStressVec, VectorType &bendingStressVec, VectorType &totalStressVec ) const {
      DKTFEVectorFunctionEvaluator<ConfiguratorType> xADFD( _conf, xA, 3 );
      for( int elementIdx=0; elementIdx<_mesh.getNumTriangs(); ++elementIdx ){
           const ElementType& El ( _mesh.getTriang( elementIdx ) );
           const typename ConfiguratorType::BaseFuncSetType &bfs = _matOptConf._conf.getBaseFunctionSet ( El );
           const int numQuadPoints = bfs.numQuadPoints( );
           RealType areaElement = 0.;
           for ( int q = 0; q < numQuadPoints; ++q ){
               Matrix32 GradXA; xADFD.evaluateGradientAtQuadPoint( El, q, GradXA );
               Matrix22 gA; gA = GradXA.transpose() * GradXA;
               RealType weight = std::sqrt( gA.determinant() ) * bfs.getWeight ( q );
               areaElement += weight;
           }
           areaElement *= El.getAreaOfRefTriangle();
           this->evaluateTotalStressOnElement( El, material, displacement, membraneStressVec[elementIdx], bendingStressVec[elementIdx], totalStressVec[elementIdx] );
           membraneStressVec[elementIdx] /= areaElement;
           bendingStressVec[elementIdx] /= areaElement;
           totalStressVec[elementIdx] /= areaElement;
      }
   }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};




template <typename MatOptConfType, typename Imp>
class LinElastEnergyOpInterface {

 protected:
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::MaskType MaskType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  const MatOptConfType &_matOptConf;

 public:
    LinElastEnergyOpInterface ( const MatOptConfType &matOpConf ) : _matOptConf ( matOpConf ){}

  void evaluateElasticEnergies( const VectorType &disp, RealType &membraneEnergy, RealType &bendingEnergy, RealType &totalEnergy ) const {
      throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      this->asImp().evaluateElasticEnergies(disp, membraneEnergy, bendingEnergy, totalEnergy);
  }

  void evaluateEnergy( const VectorType &disp, RealType &totalEnergy ) const {
      RealType membraneEnergy, bendingEnergy;
      this->asImp().evaluateElasticEnergies( disp, membraneEnergy, bendingEnergy, totalEnergy );
  }

  void evaluateGradient( const VectorType &disp, VectorType &Deriv ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    this->asImp().evaluateGradient(disp,Deriv);
  }

  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
      throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      this->asImp().assembleTripletListDirichlet(tripletListMasked,boundaryMask,Factor);
  }

  void evaluateHessianDirichlet( SparseMatrixType &Hessian, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
     std::vector<TripletType> tripletListMasked;
     this->asImp().assembleTripletListDirichlet( tripletListMasked, boundaryMask, Factor );
     Hessian.setFromTriplets( tripletListMasked.cbegin(), tripletListMasked.cend() );
     Hessian.makeCompressed();
  }

  protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};















//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
//                                  Interfaces Nonlin:
//*************************************************************************************
//*************************************************************************************
//*************************************************************************************




template <typename MatOptConfType, typename Imp>
class NonLinElastEnergyOpInterface
: public pesopt::NonlinearEnergyOp<typename MatOptConfType::ConfiguratorType::DTContainer>
{

 protected:
  typedef typename MatOptConfType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::InitType MeshType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::MaskType MaskType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;

 protected:
  const MatOptConfType &_matOptConf;
  const ConfiguratorType &_conf;
  const MeshType &_mesh;
  const int _numDofsDisp;
  const int _numVertices;

 public:
    NonLinElastEnergyOpInterface ( const MatOptConfType &matOptConf ) :
    _matOptConf ( matOptConf ), _conf( _matOptConf._conf ), _mesh ( _conf.getInitializer() ),
    _numDofsDisp( 3 * _conf.getNumGlobalDofs() ),
    _numVertices( _mesh.getNumVertices() ) {}


  const int getNumDofs ( ) const override { return _numDofsDisp; }

  void evaluateElasticEnergies( const VectorType &disp, RealType &membraneEnergy, RealType &bendingEnergy, RealType &potEnergy, RealType &totalEnergy ) const {
      throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }

  void evaluateEnergy( const VectorType &disp, RealType &totalEnergy ) const override {
      RealType membraneEnergy, bendingEnergy, potEnergy;
      this->asImp().evaluateElasticEnergies( disp, membraneEnergy, bendingEnergy, potEnergy, totalEnergy );
  }

  void evaluateJacobian( const VectorType &disp, VectorType &Deriv ) const override {
      throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }

  void assembleTripletListHessian( const VectorType &Disp, std::vector<TripletType> & tripletListMasked, const RealType Factor = 1.0 ) const {
      throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }

  void evaluateHessian( const VectorType &Disp, SparseMatrixType &Hessian ) const override {
      std::vector<TripletType> tripletList;
      this->asImp().assembleTripletListHessian( Disp, tripletList, 1.0 );
      Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
      Hessian.makeCompressed();
  }

  void evaluateTripletListHessian( const VectorType &Disp, std::vector<TripletType> & tripletList ) const override {
    this->asImp().assembleTripletListHessian( Disp, tripletList );
  }

  void evaluateTripletListHessianSym( const VectorType &Disp, std::vector<TripletType> & tripletListSym ) const override {
    std::vector<TripletType> tripletListFull; this->asImp().assembleTripletListHessian( Disp, tripletListFull, 1. );
    tripletListSym.reserve( tripletListFull.size() );
    for( int i=0; i< tripletListFull.size(); ++i ){
      int row = tripletListFull[i].row();
      int col = tripletListFull[i].col();
      if( row >= col ) tripletListSym.push_back( TripletType( row, col, tripletListFull[i].value() ) );
    }
  }

  void evaluateStressOnElements( const VectorType &disp, VectorType &membraneStressVec, VectorType &bendingStressVec, VectorType &totalStressVec ) const{
      throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }

  protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};


#endif
