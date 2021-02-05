#ifndef __FEINTEGRATORS_H
#define __FEINTEGRATORS_H

#include <pesopt_IO.h>
 
  
//!===========================================================================================================================
//! FE OPERATOR (for any quadrature type specified in the template argument)
//!===========================================================================================================================


//!===========================================================================================================================
//! Scalar-Valued
//!===========================================================================================================================
  
//! Integrator to compute \f$ \int_\Omega f(...) dx \f$, where \f$ f \f$ is the argument of the operator.
template <typename ConfiguratorType, typename Imp>
class FEIntegrator{
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
protected:
  const ConfiguratorType &_conf;
public:

  FEIntegrator ( const ConfiguratorType & Config ) : _conf( Config ) {}
  
  virtual ~FEIntegrator( ) {}

  void assembleAdd ( RealType &Dest ) const {
    RealType res = 0.;
    for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
      const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
      const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );
      RealType a = 0.;
      for ( int q = 0; q < numQuadPoints; ++q )
        a += this->asImp().evaluateIntegrand ( El, q ) * bfs.getWeight ( q );
      res += El.getVolume() * a;
    }
    Dest += res;
  }
  

  //! base function, has to be provided in derived classes.
  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    return this->asImp().evaluateIntegrand ( El, QuadPoint );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};



//!===========================================================================================================================
//! Vector-Valued Intefaces
//!===========================================================================================================================


//! Integrator for \f$(\int_\Omega s(x) w_i(x) dx )_i\f$, of some scalar valued function \f$s\f$.
template <typename ConfiguratorType, typename Imp>
class FEIntegratorVec  {
    
  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;
  const ConfiguratorType &_conf;
  
public:
  FEIntegratorVec( const ConfiguratorType & Conf ) : _conf( Conf){ }

  virtual ~FEIntegratorVec( ) {}

  void assembleAdd ( VectorType &Dest, const RealType factor = 1.0 ) const {
    
    RealType *nl_cache = new RealType[ _conf.maxNumQuadPoints() ];

    for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
      const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
      const int numLocalDofs = _conf.getNumLocalDofs ( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q )
        nl_cache[q] = this->asImp().getNonlinearity ( El, q );

      RealType aux;
  
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux = 0.;   
        for ( int q = 0; q < numQuadPoints; ++q ) aux += nl_cache[q] * bfs.evaluate ( dof, q ) * bfs.getWeight ( q );
        Dest[ _conf.localToGlobal ( El, dof ) ] += factor * El.getVolume()  * aux;
      }
    }
    delete[] nl_cache;
  }
  
  void assemble ( VectorType &Dest, const RealType factor = 1.0 ) const {
      Dest.setZero();
      this->assembleAdd( Dest,  factor );
  }
  
  void assembleDirichlet ( VectorType &Dest, const MaskType &boundaryMask, const RealType factor = 1.0 ) const {
      this->assemble( Dest,  factor );
      for( int i = 0; i < Dest.size(); ++i ){
          if ( boundaryMask[i] ){
              Dest[i] = 0.0;
          }
      } 
  }

  //! base function, has to be provided in derived classes.
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};

//! Integrator for \f$ (\int_\Omega w_i(x) dx )_{i} \f$
template <typename ConfiguratorType>
class FEMassIntegratorVec {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
  const ConfiguratorType &_conf;
  
public:
  FEMassIntegratorVec( const ConfiguratorType & Conf ) : _conf(Conf){ }

  void assembleAdd ( VectorType &Dest ) const {
    for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
      const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
      const int numLocalDofs = _conf.getNumLocalDofs ( El );
      const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );
      RealType aux;
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux = 0.;   
        for ( int q = 0; q < numQuadPoints; ++q ) aux += bfs.evaluate ( dof, q ) * bfs.getWeight ( q );
        Dest[ _conf.localToGlobal ( El, dof ) ] += El.getVolume()  * aux;
      }
    }
  }

};



//! Integrator for \f$ (\int_\Omega A(x) \nabla w_i(x) dx )_{ij} \f$, of some vector valued function \f$ A \f$.
template <typename ConfiguratorType, typename Imp>
class FEDiffOpIntegratorVec  {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::RealVecChart RealVecChart;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
  const ConfiguratorType &_conf;
  
public:
  FEDiffOpIntegratorVec ( const ConfiguratorType & Conf ) : _conf( Conf ){ }

  virtual ~FEDiffOpIntegratorVec( ) {}

  void assembleAdd ( VectorType &Dest ) const {

    RealVecChart *nl_cache = new RealVecChart[ _conf.maxNumQuadPoints() ];

    for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
      const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
      const int numLocalDofs = _conf.getNumLocalDofs ( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q )
        this->asImp().getNonlinearity (  El, q, nl_cache[q] );

      RealType aux;
  
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux = 0.0;   
        
        for ( int q = 0; q < numQuadPoints; ++q )
          aux += bfs.getWeight ( q ) * ( nl_cache[q].dot(bfs.evaluateGradient( dof, q ) ) );

        Dest[ _conf.localToGlobal ( El, dof ) ] += El.getVolume()  * aux;
      }
    }
    delete[] nl_cache;
  }

  //! base function, has to be provided in derived classes.
  void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint,
                         RealVecChart &NL ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    this->asImp().getNonlinearity (  El, QuadPoint,  NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};




//! Integrator for \f$ (\int_\Omega A(x) : D^2 w_i(x) dx )_{ij} \f$, of some matrix valued function \f$ A \f$.
template <typename ConfiguratorType, typename Imp>
class FEDiff2OpIntegratorVec  {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::RealVecChart RealVecChart;
  typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
  const ConfiguratorType &_conf;
  
public:
  FEDiff2OpIntegratorVec ( const ConfiguratorType & Conf ) : _conf( Conf ){ }

  virtual ~FEDiff2OpIntegratorVec( ) {}

  void assembleAdd ( VectorType &Dest ) const {

    DerivativeVectorValuedType *nl_cache = new DerivativeVectorValuedType[ _conf.maxNumQuadPoints() ];

    for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
      const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
      const int numLocalDofs = _conf.getNumLocalDofs ( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q )
        this->asImp().getNonlinearity (  El, q, nl_cache[q] );

      RealType aux;
  
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux = 0.0;   
        
        for ( int q = 0; q < numQuadPoints; ++q )
          aux += bfs.getWeight ( q ) * pesopt::ddProd<RealType, DerivativeVectorValuedType> ( nl_cache[q], bfs.evaluateHessian( dof, q ) );

        Dest[ _conf.localToGlobal ( El, dof ) ] += El.getVolume()  * aux;
      }
    }
    delete[] nl_cache;
  }

  //! base function, has to be provided in derived classes.
  void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint,
                         DerivativeVectorValuedType &NL ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    this->asImp().getNonlinearity (  El, QuadPoint,  NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};

//!===========================================================================================================================
//! Multi-Vector-Valued Intefaces
//!===========================================================================================================================

//! Integrator for 
//! \f$ (\int_\Omega (v(x) \cdot w_i(x)) dx )_{i} \f$
//! of some vector valued function \f$ v \f$. 
//! \f$ v \f$  is supposed to be of dimension dimChartDomain? or dimOfWorld?
template <typename ConfiguratorType, typename Imp, int dimChartDomain = ConfiguratorType::dimChartDomain>
class FEVectorIntegratorMultiVec {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::PointType PointType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;
  const ConfiguratorType &_conf;
  
public:
  FEVectorIntegratorMultiVec ( const ConfiguratorType &conf ) : _conf(conf){ }

  void assembleAdd ( typename ConfiguratorType::VectorType &Dest ) const {
    
    PointType *nl_cache = new PointType[ _conf.maxNumQuadPoints() ];
    const int numGlobalDofs = _conf.getNumGlobalDofs();
    
    for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
      const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
      const int numLocalDofs = _conf.getNumLocalDofs ( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q )
        this->asImp().getNonlinearity ( El, q, nl_cache[q] );

      PointType aux;
  
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux.setZero();    

        for ( int q = 0; q < numQuadPoints; ++q )
          aux += bfs.evaluate ( dof, q ) * bfs.getWeight ( q ) * nl_cache[q];

        for ( int comp = 0; comp < dimChartDomain; ++comp )
          Dest[ _conf.localToGlobal ( El, dof ) + comp * numGlobalDofs ] += El.getVolume()  * aux[comp];
      }
    }
    delete[] nl_cache;
  }

  void assembleDirichlet ( typename ConfiguratorType::VectorType &Dest, const MaskType& boundaryMask ) const {
      Dest.setZero();
      this->assembleAdd( Dest );
      const int numGlobalDofs = _conf.getNumGlobalDofs();
      for( int i = 0; i < numGlobalDofs; ++i ){
          if ( boundaryMask[i] ){
            for( int comp=0; comp<dimChartDomain; ++comp ) Dest[i + comp * numGlobalDofs] = 0.0;
          }
      } 
  }
  
  //! base function, has to be provided in derived classes.
  void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint, PointType &NL ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    this->asImp().getNonlinearity (  El, QuadPoint, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};


//! Integrator for \f$ (\int_\Omega trace( f(x)^T w_i(x) ) dx )_{i} \f$, of some matrix valued function \f$ f \f$.
template <typename ConfiguratorType, typename Imp, int dimChartDomain = ConfiguratorType::dimChartDomain>
class FEVectorDiffOpIntegratorMultiVec {
    
  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::PointType PointType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
  typedef typename ConfiguratorType::MaskType MaskType;
  const ConfiguratorType &_conf;
  
public:
  FEVectorDiffOpIntegratorMultiVec ( const ConfiguratorType & Conf ) : _conf( Conf){ }

  virtual ~FEVectorDiffOpIntegratorMultiVec( ) {}

  
  void assembleAdd ( typename ConfiguratorType::VectorType &Dest ) const {
    
    DerivativeVectorValuedType *nl_cache = new DerivativeVectorValuedType[ _conf.maxNumQuadPoints() ];
    const int numGlobalDofs = _conf.getNumGlobalDofs();
    
    for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
      const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
      const int numLocalDofs = _conf.getNumLocalDofs ( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q ) this->asImp().getNonlinearity ( El, q, nl_cache[q] );

      PointType aux;
  
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux.setZero();    

        for ( int q = 0; q < numQuadPoints; ++q ) aux +=  bfs.getWeight ( q ) * nl_cache[q] * bfs.evaluateGradient ( dof, q );

        for ( int comp = 0; comp < dimChartDomain; ++comp )
          Dest[ _conf.localToGlobal ( El, dof ) + comp * numGlobalDofs ] += El.getVolume()  * aux[comp];
      }
    }
    delete[] nl_cache;
  }
  
  void assembleDirichlet ( typename ConfiguratorType::VectorType &Dest, const MaskType& boundaryMask ) const {
      Dest.setZero();
      this->assembleAdd( Dest );
      const int numGlobalDofs = _conf.getNumGlobalDofs();
      for( int i = 0; i < numGlobalDofs; ++i ){
          if ( boundaryMask[i] ){
            for( int comp=0; comp<dimChartDomain; ++comp ) Dest[i + comp * numGlobalDofs] = 0.0;
          }
      } 
  }

  //! base function, has to be provided in derived classes.
  void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint, DerivativeVectorValuedType &NL ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    this->asImp().getNonlinearity (  El, QuadPoint, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};







//!===========================================================================================================================
//! Matrix-Valued Intefaces
//!===========================================================================================================================


//! General Op for matrix valued integrators
template < typename ConfiguratorType, typename Imp >
class FEMatrixValuedIntegratorBase {
public:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;

  explicit FEMatrixValuedIntegratorBase ( const ConfiguratorType &conf ) : _conf ( conf ) {}

public:
  void assembleTripletList ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
    tripletList.reserve(_conf.getMesh().getNumElements() * pesopt::Sqr( _conf.getNumLocalDofs() ) );
    LocalMatrixType localMatrix;
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
      const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
      this->asImp().prepareLocalMatrix ( El, localMatrix );
      const int numLocalDofs = _conf.getNumLocalDofs ( El );

      for ( int i = 0; i < numLocalDofs; ++i ) globalDofs[ i ] = _conf.localToGlobal ( El, i );

      for ( int i = 0; i < numLocalDofs; ++i ) {
        int glob_i = globalDofs[ i ];
        for ( int j = 0; j < numLocalDofs; ++j ) {
          int glob_j = globalDofs[ j ];
          tripletList.push_back( TripletType( glob_i, glob_j, El.getVolume()  * Factor * localMatrix(i,j) ) );
        }
      }
    }
      
  }
  
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor ) const {
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );
      
    tripletListMasked.reserve(_conf.getMesh().getNumElements() * pesopt::Sqr( _conf.getNumLocalDofs() ) );

    for( unsigned iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row()]) || (boundaryMask[tripletList[iter].col()]) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _conf.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] )
         tripletListMasked.push_back( TripletType( i, i, 1.0 ) );
    }
    
  }

  template <typename SparseMatrixType>
  void assemble ( SparseMatrixType &Dest, const RealType Factor = 1.0 ) const {
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  template <typename SparseMatrixType>
  void assembleDirichlet ( SparseMatrixType &Dest, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
       
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );
    
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve(_conf.getMesh().getNumElements() * pesopt::Sqr( _conf.getNumLocalDofs() ) );

    for( unsigned iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row()]) || (boundaryMask[tripletList[iter].col()]) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _conf.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] )
         tripletListMasked.push_back( TripletType( i, i, 1.0 ) );
    }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
  const ConfiguratorType &_conf;
};



//! The corresponding matrix assembly yields 
//! \f$ \left( \int_\Omega  w(x) \phi_i \phi_j dx \right)_{ij} \f$ 
//! for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class FEWeightedMassIntegrator :
      public FEMatrixValuedIntegratorBase<  ConfiguratorType, FEWeightedMassIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_conf;
  
public:
  FEWeightedMassIntegrator ( const ConfiguratorType & Config ) : 
   FEMatrixValuedIntegratorBase<  ConfiguratorType, FEWeightedMassIntegrator<ConfiguratorType, Imp> > ( Config ),
  _conf ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                  LocalMatrixType &LocalMatrix ) const {
    
    const int numDofs = _conf.getNumLocalDofs ( El );	
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    RealType nonlinearity;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    //TODO cache b_i, b_j? use symmetry?
    for ( int q = 0; q < numQuadPoints; ++q ) {
      nonlinearity = getNonlinearity ( El, q );
      for ( int i = 0; i < numDofs; ++i ) {
          RealType b_i = bfs.evaluate( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          RealType b_j = bfs.evaluate( j, q );
          LocalMatrix(j,i) += nonlinearity * b_i * b_j * bfs.getWeight ( q );
        }
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};



//! The corresponding matrix assembly yields 
//! \f$ \left( \int_\Omega  \phi_i \phi_j dx \right)_{ij} \f$ 
//! for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType>
class FEMassMatrix 
: public FEWeightedMassIntegrator< ConfiguratorType, FEMassMatrix<ConfiguratorType> > {
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const ConfiguratorType &_conf;
    
  public:
  FEMassMatrix ( const ConfiguratorType &conf ) 
  : FEWeightedMassIntegrator< ConfiguratorType, FEMassMatrix<ConfiguratorType>  > ( conf ),
    _conf ( conf ){}
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      return 1.;
  }
};



//! The corresponding matrix assembly yields 
//! \f$ \left(\int_\Omega A(x) \phi_i' \phi_j' dx\right)_{ij} \f$ 
//! for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class FEWeightedStiffIntegrator :
      public FEMatrixValuedIntegratorBase<  ConfiguratorType, FEWeightedStiffIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::RealVecChart RealVecChart;
  const ConfiguratorType &_conf;
  
public:
  FEWeightedStiffIntegrator ( const ConfiguratorType & Config ) : 
  FEMatrixValuedIntegratorBase<  ConfiguratorType, FEWeightedStiffIntegrator<ConfiguratorType, Imp> > ( Config ), _conf ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, LocalMatrixType &LocalMatrix ) const {
    
    const int numDofs = _conf.getNumLocalDofs ( El );	
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    RealType nonlinearity;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    //TODO cache b_i, b_j? use symmetry?
    for ( int q = 0; q < numQuadPoints; ++q ) {
      nonlinearity = getNonlinearity ( El, q );
      for ( int i = 0; i < numDofs; ++i ) {
          RealVecChart grad_i = bfs.evaluateGradient( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          RealVecChart grad_j = bfs.evaluateGradient( j, q );
          LocalMatrix(j,i) += nonlinearity * grad_i.dot(grad_j) * bfs.getWeight ( q );
        }
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};


//! The corresponding matrix assembly yields 
//! \f$ \left(\int_\Omega  \nabla \phi_i \nabla \phi_j dx \right)_{ij} \f$ 
//! for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType>
class FEStiffMatrix 
: public FEWeightedStiffIntegrator< ConfiguratorType, FEStiffMatrix<ConfiguratorType> > {
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const ConfiguratorType &_conf;
    
  public:
  FEStiffMatrix ( const ConfiguratorType &conf ) 
  : FEWeightedStiffIntegrator< ConfiguratorType, FEStiffMatrix<ConfiguratorType>  > ( conf ),
    _conf ( conf ){}
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      return 1.;
  }
};


//! The corresponding matrix assembly yields 
//! \f$ \left(\int_\Omega s(x) D^2 \phi_i : D^2 \phi_j dx\right)_{ij} \f$ 
//! for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class FEWeightedBiStiffIntegrator :
      public FEMatrixValuedIntegratorBase<  ConfiguratorType, FEWeightedBiStiffIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::RealVecChart RealVecChart;
  typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
  const ConfiguratorType &_conf;
  
public:
  FEWeightedBiStiffIntegrator ( const ConfiguratorType & Config ) : 
  FEMatrixValuedIntegratorBase< ConfiguratorType, FEWeightedBiStiffIntegrator<ConfiguratorType, Imp> > ( Config ),
  _conf ( Config ) {}

  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, LocalMatrixType &LocalMatrix ) const {
    
    const int numDofs = _conf.getNumLocalDofs ( El );	
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    RealType nonlinearity;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    //TODO cache b_i, b_j? use symmetry?
    for ( int q = 0; q < numQuadPoints; ++q ) {
      nonlinearity = getNonlinearity ( El, q );
      for ( int i = 0; i < numDofs; ++i ) {
          DerivativeVectorValuedType hessian_i = bfs.evaluateHessian( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          DerivativeVectorValuedType hessian_j = bfs.evaluateHessian( j, q );
          LocalMatrix(j,i) += nonlinearity * bfs.getWeight ( q ) * pesopt::ddProd<RealType, DerivativeVectorValuedType>( hessian_i,  hessian_j );
        }
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};


//! The corresponding matrix assembly yields 
//! \f$ \left(\int_\Omega D^2 \phi_i : D^2 \nabla \phi_j dx \right)_{ij} \f$ 
//! for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType>
class FEBiStiffMatrix 
: public FEWeightedBiStiffIntegrator< ConfiguratorType, FEBiStiffMatrix<ConfiguratorType> > {
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const ConfiguratorType &_conf;
    
  public:
  FEBiStiffMatrix ( const ConfiguratorType &conf ) 
  : FEWeightedBiStiffIntegrator< ConfiguratorType, FEBiStiffMatrix<ConfiguratorType>  > ( conf ),
    _conf ( conf ){}
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      return 1.;
  }
};




//!===========================================================================================================================
//! BlockMatrix-Valued Intefaces
//!===========================================================================================================================

//! General Interface for matrix valued integrators
template < typename ConfiguratorType, typename Imp, int dimChartDomain = ConfiguratorType::dimChartDomain>
class FEBlockMatrixValuedIntegratorBase {
public:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::RealVecChart RealVecChart;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;
  typedef typename ConfiguratorType::TripletType TripletType;

  explicit FEBlockMatrixValuedIntegratorBase ( const ConfiguratorType &conf ): _conf ( conf ) {}

public:
    
    void assembleTripletList ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
        tripletList.reserve( dimChartDomain * dimChartDomain * pesopt::Sqr( _conf.getNumLocalDofs() ) *_conf.getMesh ().getNumElements ());
        LocalMatrixType localMatrix[dimChartDomain][dimChartDomain];
        const int numGlobalDofs = _conf.getNumGlobalDofs();
        int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
        for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
            const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
            // assemble the local matrix for the current element
            this->asImp().prepareLocalMatrix ( El, localMatrix );

            const int numLocalDofs = _conf.getNumLocalDofs ( El );

            for ( int i = 0; i < numLocalDofs; ++i )
                globalDofs[ i ] = _conf.localToGlobal ( El, i );
            
            for ( int argComp = 0; argComp < dimChartDomain; ++argComp )
                for ( int destComp = 0; destComp < dimChartDomain; ++destComp )
                 for ( int i = 0; i < numLocalDofs; ++i ) {
                    int glob_i = globalDofs[ i ];
                    for ( int j = 0; j < numLocalDofs; ++j ) {
                      int glob_j = globalDofs[ j ];
                      tripletList.push_back( TripletType( glob_i + destComp * numGlobalDofs, glob_j + argComp * numGlobalDofs, El.getVolume()  * Factor * localMatrix[argComp][destComp]( i, j ) ) );
                    }
                }
        }
    }
  


  template <typename BlockMatrixType>
  void assemble ( BlockMatrixType &Dest, const RealType Factor = 1.0 ) const {
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );
    Dest.setFromTriplets( tripletList.begin(), tripletList.end() ); 
  }
  
  template <typename BlockMatrixType>
  void assembleDirichlet ( BlockMatrixType &Dest, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );

    // Boundary Mask
    const int numGlobalDofs = _conf.getNumGlobalDofs();
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( dimChartDomain * dimChartDomain * pesopt::Sqr( _conf.getNumLocalDofs() ) *_conf.getMesh().getNumElements() );
    for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row() % numGlobalDofs] ) || (boundaryMask[tripletList[iter].col() % numGlobalDofs] ) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _conf.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] ){
        for ( int Comp = 0; Comp < dimChartDomain; ++Comp )
            tripletListMasked.push_back( TripletType( i + Comp * numGlobalDofs, i + Comp * numGlobalDofs, 1.0 ) );
      }
    }
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
  
    /*
     * Dest has additional size of dimChartDomain for constraints on \int u_i = 0
     * periodic Mask contains values 0 (no periodic node) and 1 (periodic node)
     * periodicIndicesMask contains for every periodic node the corresponding node on the boundary
     */
  template <typename BlockMatrixType>
  void assemblePeriodic ( BlockMatrixType &Dest, const MaskType& periodicMask, const std::vector<int> & periodicIndicesMask, const RealType Factor = 1.0 ) const {
    
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );

    // Boundary Mask
    const int numGlobalDofs = _conf.getNumGlobalDofs();
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( dimChartDomain * dimChartDomain * pesopt::Sqr( _conf.getNumLocalDofs() ) *_conf.getMesh().getNumElements ());

    for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
        const int colIndex = tripletList[iter].col(); 
        const int colNodeIdx = colIndex % numGlobalDofs;
        const int rowIndex = tripletList[iter].row(); 
        const int rowNodeIdx = rowIndex % numGlobalDofs;
      if( (periodicMask[rowNodeIdx] ) || (periodicMask[colNodeIdx] ) ){
       //Periodic Boundary node! 
        const int colZ = colIndex / numGlobalDofs;
        const int colIndexPeriodic = periodicIndicesMask[colNodeIdx] + colZ * numGlobalDofs;
        const int rowZ = rowIndex / numGlobalDofs;
        const int rowIndexPeriodic = periodicIndicesMask[rowNodeIdx] + rowZ * numGlobalDofs;
        if( (periodicMask[rowNodeIdx] ) && (periodicMask[colNodeIdx] ) ) tripletListMasked.push_back( TripletType( rowIndexPeriodic, colIndexPeriodic, tripletList[iter].value() ) );
        if( (!periodicMask[rowNodeIdx]) && (periodicMask[colNodeIdx]) )  tripletListMasked.push_back( TripletType( rowIndex, colIndexPeriodic, tripletList[iter].value() ) );
        if( (periodicMask[rowNodeIdx]) && (!periodicMask[colNodeIdx]) )  tripletListMasked.push_back( TripletType( rowIndexPeriodic, colIndex, tripletList[iter].value() ) );
      }else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    //diagonal
    for ( int i = 0; i < _conf.getNumGlobalDofs(); ++i ){
      if ( periodicMask[i] ){
        for ( int Comp = 0; Comp < dimChartDomain; ++Comp ) tripletListMasked.push_back( TripletType( i + Comp * numGlobalDofs, i + Comp * numGlobalDofs, 1.0 ) );
      }
    }
    
    //int u_i = 0
    typename ConfiguratorType::VectorType constraintVec ( numGlobalDofs ); constraintVec.setZero();
    FEMassIntegratorVec<ConfiguratorType> ( _conf ).assembleAdd( constraintVec );
    //colabse periodically
    for ( int nodeIdx=0; nodeIdx<numGlobalDofs; nodeIdx++ ) {
        if( periodicMask[nodeIdx] ){
            constraintVec[periodicIndicesMask[nodeIdx]] += constraintVec[nodeIdx];
            constraintVec[nodeIdx] = 0.0;
        }
    }
    //insert into matrix
    for( int nodeIdx=0; nodeIdx < numGlobalDofs; ++nodeIdx ){
        for( int comp=0; comp<dimChartDomain; ++comp){
            tripletListMasked.push_back( TripletType( nodeIdx + comp * numGlobalDofs,  dimChartDomain * numGlobalDofs + comp,    constraintVec[nodeIdx] ) );
            tripletListMasked.push_back( TripletType( dimChartDomain * numGlobalDofs + comp,  nodeIdx + comp * numGlobalDofs,    constraintVec[nodeIdx] ) );
        }
    }
    
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
  
  
  //TODO so far only for LagrangianOp
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );

    // Boundary Mask
    const int numGlobalDofs = _conf.getNumGlobalDofs();
    tripletListMasked.reserve( dimChartDomain * dimChartDomain * pesopt::Sqr( _conf.getNumLocalDofs() ) *_conf.getMesh ().getNumElements ());

    for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row() % numGlobalDofs] ) || (boundaryMask[tripletList[iter].col() % numGlobalDofs] ) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _conf.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] ){
        for ( int Comp = 0; Comp < dimChartDomain; ++Comp )
            tripletListMasked.push_back( TripletType( i + Comp * numGlobalDofs, i + Comp * numGlobalDofs, 1.0 ) );
      }
    }
  }
  

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
  const ConfiguratorType &_conf;
};







//! The corresponding matrix assembly yields 
//! \f$ \left( \int_\Omega  s(x) \sum_{k=1}^{dimDomain} D^2 phi^k_i : D^2 phi^k_j dx \right)_{ij} \f$ 
//! for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class FEWeightedVectorBiStiffIntegrator :
public FEBlockMatrixValuedIntegratorBase<ConfiguratorType, FEWeightedVectorBiStiffIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::RealVecChart RealVecChart;
  typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
  const ConfiguratorType &_conf;
  
public:
  FEWeightedVectorBiStiffIntegrator ( const ConfiguratorType & Config ) : 
  FEBlockMatrixValuedIntegratorBase< ConfiguratorType, FEWeightedVectorBiStiffIntegrator<ConfiguratorType, Imp> > ( Config ),
  _conf ( Config ) {}

  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }
  
  void prepareLocalMatrix( const typename ConfiguratorType::ElementType &El,
                           LocalMatrixType (&localMatrix)[ConfiguratorType::dimChartDomain][ConfiguratorType::dimChartDomain] ) const {
      
      for (int argComp = 0; argComp < ConfiguratorType::dimChartDomain; ++argComp)
          for (int destComp = 0; destComp < ConfiguratorType::dimChartDomain; ++destComp)
              localMatrix[argComp][destComp].setZero();
            
      const typename ConfiguratorType::BaseFuncSetType &bfs = _conf.getBaseFunctionSet(El);
      const int numDofs = _conf.getNumLocalDofs ( El );
      RealType nonlinearity ;    
      
      for (int quadPoint = 0; quadPoint < _conf.maxNumQuadPoints(); ++quadPoint) {
          nonlinearity = getNonlinearity ( El, quadPoint );
          for ( int i = 0; i < numDofs; ++i ) {
              DerivativeVectorValuedType hessian_i = bfs.evaluateHessian( i, quadPoint );
              for ( int j = 0; j < numDofs; ++j ) {
                  DerivativeVectorValuedType hessian_j = bfs.evaluateHessian( j, quadPoint );
                  for(int argComp=0; argComp < ConfiguratorType::dimChartDomain; ++argComp){
                      localMatrix[argComp][argComp](j,i) += nonlinearity *  bfs.getWeight ( quadPoint ) * pesopt::ddProd<RealType, DerivativeVectorValuedType>( hessian_i,  hessian_j );
//                       for(int destComp=0; destComp < ConfiguratorType::dimChartDomain; ++destComp ){
//                         localMatrix[argComp][destComp](j,i) += nonlinearity * bfs.getWeight ( quadPoint );
//                       }
                  }
              }
          }
      }
  }
  

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};





//! The corresponding matrix assembly yields 
//! \f$ \left(\int_\Omega \sum_{k=1}^{dimDomain} D^2 \phi^k_i : D^2 \nabla \phi^k_j dx \right)_{ij} \f$ 
//! for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType>
class FEVectorBiStiffMatrix 
: public FEWeightedVectorBiStiffIntegrator< ConfiguratorType, FEVectorBiStiffMatrix<ConfiguratorType> > {
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const ConfiguratorType &_conf;
    
  public:
  FEVectorBiStiffMatrix ( const ConfiguratorType &conf ) 
  : FEWeightedVectorBiStiffIntegrator< ConfiguratorType, FEVectorBiStiffMatrix<ConfiguratorType>  > ( conf ),
    _conf ( conf ){}
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      return 1.;
  }
};








//! ==========================================================================
//! base clases for integrals of functions u = u_FE + u_aff 
//! =========================================================================
    
    
    
//! 
template < typename ConfiguratorType, typename Imp>
class FEPlusAffineGradMatrixValuedIntegratorBase {
public:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::RealVecChart RealVecChart;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::GlobalAffineGradBaseFuncSet GlobalAffineGradBaseFuncSet;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::LocalMatrixTypeAffineGrad LocalMatrixTypeAffineGrad;
  typedef typename ConfiguratorType::LocalMatrixTypeMixedAffineGrad LocalMatrixTypeMixedAffineGrad;

  explicit FEPlusAffineGradMatrixValuedIntegratorBase ( const ConfiguratorType &conf ) : _conf ( conf ) {}

public:
    
    void assembleTripletList ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
        tripletList.reserve( pesopt::Sqr( _conf.getNumLocalDofs() ) * _conf.getMesh().getNumElements() );
        LocalMatrixType localMatrix;
        int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
        for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
            const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
            // assemble the local matrix for the current element
            this->asImp().prepareLocalMatrix ( El, localMatrix );
            const int numLocalDofs = _conf.getNumLocalDofs ( El );

            for ( int i = 0; i < numLocalDofs; ++i ) globalDofs[ i ] = _conf.localToGlobal ( El, i );
            
            for ( int i = 0; i < numLocalDofs; ++i ) {
                    int glob_i = globalDofs[ i ];
                    for ( int j = 0; j < numLocalDofs; ++j ) {
                      int glob_j = globalDofs[ j ];
                      tripletList.push_back( TripletType( glob_i, glob_j, El.getVolume() * Factor * localMatrix( i, j ) ) );
                    }
            }
        }
    }
    
  //assembles only dimAffine x dimVector Matrix (without offsets)
   void assembleTripletListMixedPart ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
       const int numGlobalDofs = _conf.getNumGlobalDofs();
        GlobalAffineGradBaseFuncSet globAffBfs; const int numAffineDofs = globAffBfs.numBaseFuncs;
        tripletList.reserve(  _conf.getNumLocalDofs() * numAffineDofs  *_conf.getMesh().getNumElements ());
        int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
        LocalMatrixTypeMixedAffineGrad localMatrix;
        for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
            const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
            // assemble the local matrix for the current element
            this->asImp().prepareLocalMatrixMixedPart ( El, localMatrix );
            const int numLocalDofs = _conf.getNumLocalDofs ( El );
            for ( int i = 0; i < numLocalDofs; ++i ) globalDofs[ i ] = _conf.localToGlobal ( El, i );
            for ( int i = 0; i < numLocalDofs; ++i ) {
                    int glob_i = globalDofs[ i ];
                    for ( int affIndex = 0; affIndex < globAffBfs.numBaseFuncs; ++affIndex ) {
                      tripletList.push_back( TripletType( affIndex, glob_i, El.getVolume() * Factor * localMatrix( affIndex, i ) ) );
                    }
            }
        }
    }
    
    //without offset
    void assembleTripletListAffinePart ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
        const int numDofsOtherDisp = _conf.getNumGlobalDofs() * _conf.dimChartDomain;
        GlobalAffineGradBaseFuncSet globAffBfs;
        LocalMatrixTypeAffineGrad localMatrix;
        for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
            const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
            this->asImp().prepareLocalMatrixAffinePart ( El, localMatrix );
            for( int affIndexArg = 0; affIndexArg < globAffBfs.numBaseFuncs; ++affIndexArg )
                for( int affIndexDest = 0; affIndexDest < globAffBfs.numBaseFuncs; ++affIndexDest ){
                   tripletList.push_back( TripletType( affIndexDest, affIndexArg, El.getVolume() * Factor * localMatrix( affIndexArg, affIndexDest ) ) );
                }
        }
    }
  
public:
  
    
    
  template <typename BlockMatrixType>
  void assemblePeriodic ( BlockMatrixType &Dest, const MaskType& periodicMask, const std::vector<int> & periodicIndicesMask, const RealType Factor = 1.0 ) const {
    const int numGlobalDofs = _conf.getNumGlobalDofs(); 
    const int numLocalDofs = _conf.getNumLocalDofs();
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( pesopt::Sqr( numLocalDofs * _conf.dimChartDomain ) * _conf.getMesh().getNumElements() + _conf.dimChartDomain * numGlobalDofs );
   
    //periodic part
        std::vector<TripletType> tripletList;
        assembleTripletList ( tripletList, Factor );

        // Boundary Mask
        for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
            const int colIndex = tripletList[iter].col(); 
            const int colNodeIdx = colIndex % numGlobalDofs;
            const int rowIndex = tripletList[iter].row(); 
            const int rowNodeIdx = rowIndex % numGlobalDofs;
            if( (periodicMask[rowNodeIdx] ) || (periodicMask[colNodeIdx] ) ){
            //Periodic Boundary node! 
                const int colZ = colIndex / numGlobalDofs;
                const int colIndexPeriodic = periodicIndicesMask[colNodeIdx] + colZ * numGlobalDofs;
                const int rowZ = rowIndex / numGlobalDofs;
                const int rowIndexPeriodic = periodicIndicesMask[rowNodeIdx] + rowZ * numGlobalDofs;
                if( (periodicMask[rowNodeIdx] ) && (periodicMask[colNodeIdx] ) ) tripletListMasked.push_back( TripletType( rowIndexPeriodic, colIndexPeriodic, tripletList[iter].value() ) );
                if( (!periodicMask[rowNodeIdx]) && (periodicMask[colNodeIdx]) )  tripletListMasked.push_back( TripletType( rowIndex, colIndexPeriodic, tripletList[iter].value() ) );
                if( (periodicMask[rowNodeIdx]) && (!periodicMask[colNodeIdx]) )  tripletListMasked.push_back( TripletType( rowIndexPeriodic, colIndex, tripletList[iter].value() ) );
            }else {
                tripletListMasked.push_back( tripletList[iter] );
            }
        }
        
        //diagonal
        for ( int i = 0; i < _conf.getNumGlobalDofs(); ++i ){
            if ( periodicMask[i] ){
               tripletListMasked.push_back( TripletType( i, i, 1.0 ) );
            }
        }
    
    
    //int u_i^per = 0
        typename ConfiguratorType::VectorType constraintVec ( numGlobalDofs ); constraintVec.setZero();
        FEMassIntegratorVec<ConfiguratorType> ( _conf ).assembleAdd( constraintVec );
        //colabse periodically
        for ( int nodeIdx=0; nodeIdx<numGlobalDofs; nodeIdx++ ) {
            if( periodicMask[nodeIdx] ){
                constraintVec[periodicIndicesMask[nodeIdx]] += constraintVec[nodeIdx];
                constraintVec[nodeIdx] = 0.0;
            }
        }
        //insert into matrix
        for( int nodeIdx=0; nodeIdx < numGlobalDofs; ++nodeIdx ){
                tripletListMasked.push_back( TripletType( nodeIdx,  numGlobalDofs,  constraintVec[nodeIdx] ) );
                tripletListMasked.push_back( TripletType( numGlobalDofs,  nodeIdx,  constraintVec[nodeIdx] ) );
        }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
    
  template <typename BlockMatrixType>
  void assembleMixedPeriodic ( BlockMatrixType &Dest, const MaskType& periodicMask, const std::vector<int> & periodicIndicesMask, const RealType Factor = 1.0 ) const {
    const int numGlobalDofs = _conf.getNumGlobalDofs();
    const int numLocalDofs = _conf.getNumLocalDofs();
    GlobalAffineGradBaseFuncSet globAffBfs; const int numAffineDofs = globAffBfs.numBaseFuncs;
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( ( pesopt::Sqr( numLocalDofs * _conf.dimChartDomain ) + 2 * numLocalDofs * _conf.dimChartDomain * numAffineDofs + pesopt::Sqr( numAffineDofs ) ) * _conf.getMesh().getNumElements() + _conf.dimChartDomain * numGlobalDofs );
    
    //Mixed Part
       std::vector<TripletType> tripletListMixed;
       assembleTripletListMixedPart( tripletListMixed, Factor );
       // Boundary Mask
       for( unsigned int iter=0; iter < tripletListMixed.size(); ++iter ){
            const int colIndex = tripletListMixed[iter].col(); 
            const int colNodeIdx = colIndex % numGlobalDofs;
            const int rowIndex = tripletListMixed[iter].row();
            if( periodicMask[colNodeIdx] ){
            //Periodic Boundary node! 
                const int colZ = colIndex / numGlobalDofs;
                const int colIndexPeriodic = periodicIndicesMask[colNodeIdx] + colZ * numGlobalDofs;
                tripletListMasked.push_back( TripletType( colIndexPeriodic, rowIndex, tripletListMixed[iter].value() ) );
            }else {
                tripletListMasked.push_back( TripletType( colIndex, rowIndex, tripletListMixed[iter].value() ) );
            }
       }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }

 template <typename BlockMatrixType>
  void assembleAffine ( BlockMatrixType &Dest, const RealType Factor = 1.0 ) const {
    //Affine Part
    std::vector<TripletType> tripletListAffine;
    assembleTripletListAffinePart( tripletListAffine, Factor );
    Dest.setFromTriplets( tripletListAffine.begin(), tripletListAffine.end() );
  }
  
protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
  const ConfiguratorType &_conf;
  
};







template < typename ConfiguratorType, typename Imp, int dimChartDomain = ConfiguratorType::dimChartDomain>
class FEPlusAffineSymGradBlockMatrixValuedIntegratorBase {
public:
  typedef typename ConfiguratorType::RealType                           RealType;
  typedef typename ConfiguratorType::RealVecChart                       RealVecChart;
  typedef typename ConfiguratorType::ElementType                        ElementType;
  typedef typename ConfiguratorType::MaskType                           MaskType;
  typedef typename ConfiguratorType::TripletType                        TripletType;
  typedef typename ConfiguratorType::GlobalAffineSymGradBaseFuncSet     GlobalAffineSymGradBaseFuncSet;
  typedef typename ConfiguratorType::LocalMatrixType                    LocalMatrixType;
  typedef typename ConfiguratorType::LocalMatrixTypeAffineSymGrad       LocalMatrixTypeAffineSymGrad;
  typedef typename ConfiguratorType::LocalMatrixTypeMixedAffineSymGrad  LocalMatrixTypeMixedAffineSymGrad;

  explicit FEPlusAffineSymGradBlockMatrixValuedIntegratorBase ( const ConfiguratorType &conf ): _conf ( conf ) {}

public:
    
    void assembleTripletList ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
        tripletList.reserve( pesopt::Sqr( dimChartDomain * _conf.getNumLocalDofs() ) * _conf.getMesh ().getNumElements ());
        LocalMatrixType localMatrix[dimChartDomain][dimChartDomain];
        const int numGlobalDofs = _conf.getNumGlobalDofs();
        int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
        for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
            const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
            // assemble the local matrix for the current element
            this->asImp().prepareLocalMatrix ( El, localMatrix );
            const int numLocalDofs = _conf.getNumLocalDofs ( El );

            for ( int i = 0; i < numLocalDofs; ++i ) globalDofs[ i ] = _conf.localToGlobal ( El, i );
            
            for ( int argComp = 0; argComp < dimChartDomain; ++argComp )
                for ( int destComp = 0; destComp < dimChartDomain; ++destComp )
                 for ( int i = 0; i < numLocalDofs; ++i ) {
                    int glob_i = globalDofs[ i ];
                    for ( int j = 0; j < numLocalDofs; ++j ) {
                      int glob_j = globalDofs[ j ];
                      tripletList.push_back( TripletType( glob_i + destComp * numGlobalDofs, glob_j + argComp * numGlobalDofs, El.getVolume() * Factor * localMatrix[argComp][destComp]( i, j ) ) );
                    }
                }
        }
    }
    
  //assembles only dimAffine x dimFEVector Matrix (without offsets)
   void assembleTripletListMixedPart ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
       const int numGlobalDofs = _conf.getNumGlobalDofs();
       const int numDofsOtherDisp = _conf.getNumGlobalDofs() * _conf.dimChartDomain;
        GlobalAffineSymGradBaseFuncSet globAffBfs; const int numAffineDofs = globAffBfs.numBaseFuncs;
        tripletList.reserve( dimChartDomain * _conf.getNumLocalDofs() * numAffineDofs  *_conf.getMesh ().getNumElements ());
        int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
        LocalMatrixTypeMixedAffineSymGrad localMatrix[dimChartDomain];
        for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
            const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
            // assemble the local matrix for the current element
            this->asImp().prepareLocalMatrixMixedPart ( El, localMatrix );
            const int numLocalDofs = _conf.getNumLocalDofs ( El );
            for ( int i = 0; i < numLocalDofs; ++i ) globalDofs[ i ] = _conf.localToGlobal ( El, i );
            for ( int comp = 0; comp < dimChartDomain; ++comp )
                 for ( int i = 0; i < numLocalDofs; ++i ) {
                    int glob_i = globalDofs[ i ];
                    for ( int affIndex = 0; affIndex < globAffBfs.numBaseFuncs; ++affIndex ) {
                      tripletList.push_back( TripletType( affIndex, glob_i + comp * numGlobalDofs, El.getVolume() * Factor * localMatrix[comp]( affIndex, i ) ) );
//                       cout << "triplet mixed at " <<  affIndex <<  ",  " <<  glob_i + comp * numGlobalDofs <<  " : " <<  El.getVolume() * Factor * localMatrix[comp]( affIndex, i ) <<  endl;
                    }
                }
        }
    }
    
    //without offset
    void assembleTripletListAffinePart ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
        const int numDofsOtherDisp = _conf.getNumGlobalDofs() * _conf.dimChartDomain;
        GlobalAffineSymGradBaseFuncSet globAffBfs;
        LocalMatrixTypeAffineSymGrad localMatrix;
        for ( int elementIdx = 0; elementIdx < _conf.getMesh().getNumElements(); ++elementIdx){
            const ElementType& El ( _conf.getMesh().getElement( elementIdx ) );
            this->asImp().prepareLocalMatrixAffinePart ( El, localMatrix );
            for( int affIndexArg = 0; affIndexArg < globAffBfs.numBaseFuncs; ++affIndexArg )
                for( int affIndexDest = 0; affIndexDest < globAffBfs.numBaseFuncs; ++affIndexDest ){
                   tripletList.push_back( TripletType( affIndexDest, affIndexArg, El.getVolume() * Factor * localMatrix( affIndexArg, affIndexDest ) ) );
                }
        }
    }
  
public:
    
    
  template <typename BlockMatrixType>
  void assemble ( BlockMatrixType &Dest, const RealType Factor = 1.0 ) const {
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );
    Dest.setFromTriplets( tripletList.begin(), tripletList.end() );
  }
    
  template <typename BlockMatrixType>
  void assembleDirichlet ( BlockMatrixType &Dest, 
                           const MaskType& boundaryMask, 
                           const RealType Factor = 1.0 ) const {
    
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );

    // Boundary Mask
    const int numGlobalDofs = _conf.getNumGlobalDofs();
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( pesopt::Sqr( dimChartDomain * _conf.getNumLocalDofs() ) *_conf.getMesh().getNumElements() );
    for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row() % numGlobalDofs] ) || (boundaryMask[tripletList[iter].col() % numGlobalDofs] ) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    //Diagonal
    for ( int i = 0; i < _conf.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] ){
        for ( int Comp = 0; Comp < dimChartDomain; ++Comp )
            tripletListMasked.push_back( TripletType( i + Comp * numGlobalDofs, i + Comp * numGlobalDofs, 1.0 ) );
      }
    }
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
    
  template <typename BlockMatrixType>
  void assembleNeumann ( BlockMatrixType &Dest, const RealType Factor = 1.0 ) const {
    const int numGlobalDofs = _conf.getNumGlobalDofs(); const int numNeumannDispDofs = dimChartDomain * numGlobalDofs;
    const int numLocalDofs = _conf.getNumLocalDofs();
    GlobalAffineSymGradBaseFuncSet globAffBfs; const int numAffineDofs = globAffBfs.numBaseFuncs;
    
    std::vector<TripletType> tripletList;
    tripletList.reserve( ( pesopt::Sqr( numLocalDofs * dimChartDomain ) + 2 * numLocalDofs * dimChartDomain * numAffineDofs + pesopt::Sqr( numAffineDofs ) ) * _conf.getMesh().getNumElements() + dimChartDomain * numNeumannDispDofs );
    
    //second derivative of energy part
    assembleTripletList ( tripletList, Factor );
    
    //constraint : int u_i^Neumann = 0
        typename ConfiguratorType::VectorType constraintVec ( numGlobalDofs ); constraintVec.setZero();
        FEMassIntegratorVec<ConfiguratorType> ( _conf ).assembleAdd( constraintVec );
        //insert into matrix
        for( int nodeIdx=0; nodeIdx < numGlobalDofs; ++nodeIdx ){
            for( int comp=0; comp<dimChartDomain; ++comp){
                tripletList.push_back( TripletType( nodeIdx + comp * numGlobalDofs,  numNeumannDispDofs + comp,   constraintVec[nodeIdx] ) );
                tripletList.push_back( TripletType( numNeumannDispDofs + comp,  nodeIdx + comp * numGlobalDofs,   constraintVec[nodeIdx] ) );
            }
        }
    Dest.setFromTriplets( tripletList.begin(), tripletList.end() );
  }
    
    
  template <typename BlockMatrixType>
  void assemblePeriodic ( BlockMatrixType &Dest, 
                          const MaskType& periodicMask, const std::vector<int> & periodicIndicesMask, 
                          const RealType Factor = 1.0 ) const {
                              
    const int numGlobalDofs = _conf.getNumGlobalDofs();
    const int numPeriodicDispDofs = dimChartDomain * numGlobalDofs;
    const int numLocalDofs = _conf.getNumLocalDofs();
    GlobalAffineSymGradBaseFuncSet globAffBfs;
    const int numAffineDofs = globAffBfs.numBaseFuncs;
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( ( pesopt::Sqr( numLocalDofs * dimChartDomain ) + 2 * numLocalDofs * dimChartDomain * numAffineDofs + pesopt::Sqr( numAffineDofs ) ) * _conf.getMesh().getNumElements() + dimChartDomain * numPeriodicDispDofs );
 
    //periodic part
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );

    // Boundary Mask
    for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
        const int colIndex = tripletList[iter].col(); 
        const int colNodeIdx = colIndex % numGlobalDofs;
        const int rowIndex = tripletList[iter].row(); 
        const int rowNodeIdx = rowIndex % numGlobalDofs;
        if( (periodicMask[rowNodeIdx] ) || (periodicMask[colNodeIdx] ) ){
            //Periodic Boundary node! 
            const int colZ = colIndex / numGlobalDofs;
            const int colIndexPeriodic = periodicIndicesMask[colNodeIdx] + colZ * numGlobalDofs;
            const int rowZ = rowIndex / numGlobalDofs;
            const int rowIndexPeriodic = periodicIndicesMask[rowNodeIdx] + rowZ * numGlobalDofs;
            if( (periodicMask[rowNodeIdx] ) && (periodicMask[colNodeIdx] ) ) tripletListMasked.push_back( TripletType( rowIndexPeriodic, colIndexPeriodic, tripletList[iter].value() ) );
            if( (!periodicMask[rowNodeIdx]) && (periodicMask[colNodeIdx]) )  tripletListMasked.push_back( TripletType( rowIndex, colIndexPeriodic, tripletList[iter].value() ) );
            if( (periodicMask[rowNodeIdx]) && (!periodicMask[colNodeIdx]) )  tripletListMasked.push_back( TripletType( rowIndexPeriodic, colIndex, tripletList[iter].value() ) );
        }else {
                tripletListMasked.push_back( tripletList[iter] );
        }
    }
        
    //diagonal
    for ( int i = 0; i < _conf.getNumGlobalDofs(); ++i ){
        if ( periodicMask[i] ){
            for ( int Comp = 0; Comp < dimChartDomain; ++Comp ) 
                tripletListMasked.push_back( TripletType( i + Comp * numGlobalDofs, i + Comp * numGlobalDofs, 1.0 ) );
        }
    }
    
    
    //int u_i^per = 0
    typename ConfiguratorType::VectorType constraintVec ( numGlobalDofs ); constraintVec.setZero();
    FEMassIntegratorVec<ConfiguratorType> ( _conf ).assembleAdd( constraintVec );
    //colabse periodically
    for ( int nodeIdx=0; nodeIdx<numGlobalDofs; nodeIdx++ ) {
            if( periodicMask[nodeIdx] ){
                constraintVec[periodicIndicesMask[nodeIdx]] += constraintVec[nodeIdx];
                constraintVec[nodeIdx] = 0.0;
            }
    }
    //insert into matrix
    for( int nodeIdx=0; nodeIdx < numGlobalDofs; ++nodeIdx ){
            for( int comp=0; comp<dimChartDomain; ++comp){
                tripletListMasked.push_back( TripletType( nodeIdx + comp * numGlobalDofs,  numPeriodicDispDofs + comp,   constraintVec[nodeIdx] ) );
                tripletListMasked.push_back( TripletType( numPeriodicDispDofs + comp,  nodeIdx + comp * numGlobalDofs,    constraintVec[nodeIdx] ) );
            }
    }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
  
  template <typename BlockMatrixType>
  void assembleMixedDirichlet ( BlockMatrixType &Dest, 
                                const MaskType& boundaryMask,
                                const RealType Factor = 1.0 ) const {
    const int numGlobalDofs = _conf.getNumGlobalDofs(); 
    const int numPeriodicDispDofs = dimChartDomain * numGlobalDofs;
    const int numLocalDofs = _conf.getNumLocalDofs();
    GlobalAffineSymGradBaseFuncSet globAffBfs; 
    const int numAffineDofs = globAffBfs.numBaseFuncs;
      
    //Mixed Part
    std::vector<TripletType> tripletListMixed;
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( ( pesopt::Sqr( numLocalDofs * dimChartDomain ) + 2 * numLocalDofs * dimChartDomain * numAffineDofs + pesopt::Sqr( numAffineDofs ) ) * _conf.getMesh().getNumElements() + dimChartDomain * numPeriodicDispDofs );
    assembleTripletListMixedPart( tripletListMixed, Factor );
      
    // Boundary Mask
    for( unsigned int iter=0; iter < tripletListMixed.size(); ++iter ){
        const int colIndex = tripletListMixed[iter].col(); 
        const int colNodeIdx = colIndex % numGlobalDofs;
        const int rowIndex = tripletListMixed[iter].row();
      if( (boundaryMask[colNodeIdx] ) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( TripletType( colIndex, rowIndex, tripletListMixed[iter].value() ) );
      }
    }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
  
  template <typename BlockMatrixType>
  void assembleMixedNeumann ( BlockMatrixType &Dest, const RealType Factor = 1.0 ) const {
    const int numGlobalDofs = _conf.getNumGlobalDofs(); 
    const int numNeumannDispDofs = dimChartDomain * numGlobalDofs;
    const int numLocalDofs = _conf.getNumLocalDofs();
    GlobalAffineSymGradBaseFuncSet globAffBfs; 
    const int numAffineDofs = globAffBfs.numBaseFuncs;
    
    //Mixed Part
    std::vector<TripletType> tripletListMixed;
    tripletListMixed.reserve( ( pesopt::Sqr( numLocalDofs * dimChartDomain ) + 2 * numLocalDofs * dimChartDomain * numAffineDofs + pesopt::Sqr( numAffineDofs ) ) * _conf.getMesh().getNumElements() + dimChartDomain * numNeumannDispDofs );
    assembleTripletListMixedPart( tripletListMixed, Factor );
    
    // transpose 
    std::vector<TripletType> tripletListTransposed;
    tripletListTransposed.reserve( tripletListMixed.size() );
    for( unsigned int iter=0; iter < tripletListMixed.size(); ++iter ){
        const int colIndex = tripletListMixed[iter].col(); 
        const int rowIndex = tripletListMixed[iter].row();
        tripletListTransposed.push_back( TripletType( colIndex, rowIndex, tripletListMixed[iter].value() ) );
    }
    
    Dest.setFromTriplets( tripletListTransposed.begin(), tripletListTransposed.end() );
  }
  
  
  template <typename BlockMatrixType>
  void assembleMixedPeriodic ( BlockMatrixType &Dest, 
                               const MaskType& periodicMask, const std::vector<int> & periodicIndicesMask, 
                               const RealType Factor = 1.0 ) const {
    const int numGlobalDofs = _conf.getNumGlobalDofs();
    const int numPeriodicDispDofs = dimChartDomain * numGlobalDofs;
    const int numLocalDofs = _conf.getNumLocalDofs();
    GlobalAffineSymGradBaseFuncSet globAffBfs; 
    const int numAffineDofs = globAffBfs.numBaseFuncs;
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( ( pesopt::Sqr( numLocalDofs * dimChartDomain ) + 2 * numLocalDofs * dimChartDomain * numAffineDofs + pesopt::Sqr( numAffineDofs ) ) * _conf.getMesh().getNumElements() + dimChartDomain * numPeriodicDispDofs );
    
    //Mixed Part
    std::vector<TripletType> tripletListMixed;
    assembleTripletListMixedPart( tripletListMixed, Factor );
    
    // Boundary Mask
    for( unsigned int iter=0; iter < tripletListMixed.size(); ++iter ){
            const int colIndex = tripletListMixed[iter].col(); 
            const int colNodeIdx = colIndex % numGlobalDofs;
            const int rowIndex = tripletListMixed[iter].row();
            if( periodicMask[colNodeIdx] ){
            //Periodic Boundary node! 
                const int colZ = colIndex / numGlobalDofs;
                const int colIndexPeriodic = periodicIndicesMask[colNodeIdx] + colZ * numGlobalDofs;
                tripletListMasked.push_back( TripletType( colIndexPeriodic, rowIndex, tripletListMixed[iter].value() ) );
            }else {
                tripletListMasked.push_back( TripletType( colIndex, rowIndex, tripletListMixed[iter].value() ) );
            }
    }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
  
  
 template <typename BlockMatrixType>
  void assembleAffine ( BlockMatrixType &Dest, const RealType Factor = 1.0 ) const {
    //Affine Part
    std::vector<TripletType> tripletListAffine;
    this->assembleTripletListAffinePart( tripletListAffine, Factor );
    Dest.setFromTriplets( tripletListAffine.begin(), tripletListAffine.end() );
  }
  
  
  
    // Dest has additional size of dimChartDomain for constraints on \int u_i = 0
  // periodic Mask contains values 0 (no periodic node) and 1 (periodic node)
  // periodicIndicesMask contains for every periodic node the corresponding node on the boundary
//   template <typename BlockMatrixType>
//   void assemblePeriodicPlusAffine ( BlockMatrixType &Dest, 
//                                     const MaskType& periodicMask, const std::vector<int> & periodicIndicesMask, 
//                                     const RealType Factor = 1.0 ) const {
//     const int numGlobalDofs = _conf.getNumGlobalDofs(); const int numPeriodicDispDofs = dimChartDomain * numGlobalDofs;
//     const int numLocalDofs = _conf.getNumLocalDofs();
//     GlobalAffineSymGradBaseFuncSet globAffBfs; const int numAffineDofs = globAffBfs.numBaseFuncs;
//     std::vector<TripletType> tripletListMasked;
//     tripletListMasked.reserve( ( pesopt::Sqr( numLocalDofs * dimChartDomain ) + 2 * numLocalDofs * dimChartDomain * numAffineDofs + pesopt::Sqr( numAffineDofs ) ) * _conf.getMesh().getNumElements() + dimChartDomain * numPeriodicDispDofs );
//  
//     //periodic part
//         std::vector<TripletType> tripletList;
//         assembleTripletList ( tripletList, Factor );
// 
//         // Boundary Mask
//         for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
//             const int colIndex = tripletList[iter].col(); 
//             const int colNodeIdx = colIndex % numGlobalDofs;
//             const int rowIndex = tripletList[iter].row(); 
//             const int rowNodeIdx = rowIndex % numGlobalDofs;
//             if( (periodicMask[rowNodeIdx] ) || (periodicMask[colNodeIdx] ) ){
//             //Periodic Boundary node! 
//                 const int colZ = colIndex / numGlobalDofs;
//                 const int colIndexPeriodic = periodicIndicesMask[colNodeIdx] + colZ * numGlobalDofs;
//                 const int rowZ = rowIndex / numGlobalDofs;
//                 const int rowIndexPeriodic = periodicIndicesMask[rowNodeIdx] + rowZ * numGlobalDofs;
//                 if( (periodicMask[rowNodeIdx] ) && (periodicMask[colNodeIdx] ) ) tripletListMasked.push_back( TripletType( rowIndexPeriodic, colIndexPeriodic, tripletList[iter].value() ) );
//                 if( (!periodicMask[rowNodeIdx]) && (periodicMask[colNodeIdx]) )  tripletListMasked.push_back( TripletType( rowIndex, colIndexPeriodic, tripletList[iter].value() ) );
//                 if( (periodicMask[rowNodeIdx]) && (!periodicMask[colNodeIdx]) )  tripletListMasked.push_back( TripletType( rowIndexPeriodic, colIndex, tripletList[iter].value() ) );
//             }else {
//                 tripletListMasked.push_back( tripletList[iter] );
//             }
//         }
//         
//         //diagonal
//         for ( int i = 0; i < _conf.getNumGlobalDofs(); ++i ){
//             if ( periodicMask[i] ){
//                 for ( int Comp = 0; Comp < dimChartDomain; ++Comp ) tripletListMasked.push_back( TripletType( i + Comp * numGlobalDofs, i + Comp * numGlobalDofs, 1.0 ) );
//             }
//         }
//     
//     
//     //Mixed Part
//        std::vector<TripletType> tripletListMixed;
//        assembleTripletListMixedPart( tripletListMixed, Factor );
//        // Boundary Mask
//        for( unsigned int iter=0; iter < tripletListMixed.size(); ++iter ){
//             const int colIndex = tripletListMixed[iter].col(); 
//             const int colNodeIdx = colIndex % numGlobalDofs;
//             const int rowIndex = tripletListMixed[iter].row();
//             if( periodicMask[colNodeIdx] ){
//             //Periodic Boundary node! 
//                 const int colZ = colIndex / numGlobalDofs;
//                 const int colIndexPeriodic = periodicIndicesMask[colNodeIdx] + colZ * numGlobalDofs;
//                 tripletListMasked.push_back( TripletType( rowIndex + numPeriodicDispDofs, colIndexPeriodic, tripletListMixed[iter].value() ) );
//                 tripletListMasked.push_back( TripletType( colIndexPeriodic, rowIndex + numPeriodicDispDofs, tripletListMixed[iter].value() ) );
//             }else {
//                 tripletListMasked.push_back( TripletType( rowIndex + numPeriodicDispDofs, colIndex, tripletListMixed[iter].value() ) );
//                 tripletListMasked.push_back( TripletType( colIndex, rowIndex + numPeriodicDispDofs, tripletListMixed[iter].value() ) );
//             }
//        }
//     
//     //Affine Part
//        std::vector<TripletType> tripletListAffine;
//        this->assembleTripletListAffinePart( tripletListAffine, Factor );
//        for( unsigned int iter=0; iter < tripletListAffine.size(); ++iter ){
//             const int colIndex = tripletListAffine[iter].col(); 
//             const int rowIndex = tripletListAffine[iter].row();
//             tripletListMasked.push_back( TripletType( rowIndex + numPeriodicDispDofs, colIndex + numPeriodicDispDofs, tripletListAffine[iter].value() ) );
//        }
//     
//     
//     //int u_i^per = 0
//         typename ConfiguratorType::VectorType constraintVec ( numGlobalDofs ); constraintVec.setZero();
//         FEMassIntegratorVec<ConfiguratorType> ( _conf ).assembleAdd( constraintVec );
//         //colabse periodically
//         for ( int nodeIdx=0; nodeIdx<numGlobalDofs; nodeIdx++ ) {
//             if( periodicMask[nodeIdx] ){
//                 constraintVec[periodicIndicesMask[nodeIdx]] += constraintVec[nodeIdx];
//                 constraintVec[nodeIdx] = 0.0;
//             }
//         }
//         //insert into matrix
//         for( int nodeIdx=0; nodeIdx < numGlobalDofs; ++nodeIdx ){
//             for( int comp=0; comp<dimChartDomain; ++comp){
//                 tripletListMasked.push_back( TripletType( nodeIdx + comp * numGlobalDofs,  numPeriodicDispDofs + numAffineDofs  + comp,   constraintVec[nodeIdx] ) );
//                 tripletListMasked.push_back( TripletType( numPeriodicDispDofs + numAffineDofs + comp,  nodeIdx + comp * numGlobalDofs,    constraintVec[nodeIdx] ) );
//             }
//         }
//     
//     Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
//   }
  
protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
  const ConfiguratorType &_conf;
  
};





#endif
