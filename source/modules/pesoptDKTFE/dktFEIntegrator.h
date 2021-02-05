#ifndef __DKTFEINTEGRATOR_H
#define __DKTFEINTEGRATOR_H

#include <dktFEFunctionEvaluator.h>
  
  
//!===========================================================================================================================
//! FE OPERATOR INTERFACES (for any quadrature type specified in the template argument)
//!===========================================================================================================================


//!===========================================================================================================================
//! Scalar-Valued Intefaces
//!===========================================================================================================================
  
//! Integrator to compute \f$\int_\Omega f(...) dx\f$, where \f$\f$ is the argument of the operator.
template <typename ConfiguratorType, typename Imp>
class RefTriangleIntegrator{
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
protected:
  const ConfiguratorType &_config;
public:

  RefTriangleIntegrator ( const ConfiguratorType & Config ) : _config( Config ) {}
  
  virtual ~RefTriangleIntegrator( ) {}

  void assembleAdd ( RealType &Dest ) const {
    RealType res = 0.;
    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );
      RealType a = 0.;
      for ( int q = 0; q < numQuadPoints; ++q ) a += this->asImp().evaluateIntegrand ( El, q ) * bfs.getWeight ( q );
      res += El.getAreaOfRefTriangle() * a;
    }
    Dest += res;
  }
  
  //assembles vector with \int_El f(...)
  template<typename VectorType>
  void assembleOnElements ( VectorType &Dest ) const {
    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );
      RealType a = 0.;
      for ( int q = 0; q < numQuadPoints; ++q ) a += this->asImp().evaluateIntegrand ( El, q ) * bfs.getWeight ( q );
      Dest[elementIdx] = El.getAreaOfRefTriangle() * a;
    }
  }

  //! interface function, has to be provided in derived classes.
  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
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


//! General Interface for vector valued integrators
// template < typename ConfiguratorType, typename Imp >
// class VectorValuedIntegratorBase {
// public:
//   typedef typename ConfiguratorType::RealType   RealType;
//   typedef typename ConfiguratorType::ElementType ElementType;
// 
//   explicit VectorValuedIntegratorBase ( const ConfiguratorType &conf ) : _conf ( conf ) {}
// 
//   template <typename VectorType>
//   void assembleAdd ( VectorType &Dest, const RealType Factor = 1.0 ) const {
//     VectorType localVec (ConfiguratorType::maxNumLocalDofs);
//     for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
//       const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
//       this->asImp().prepareLocalVec ( El, localVec );
//       for ( int i = 0; i < _conf.getNumLocalDofs ( El ); ++i )
//           Dest[_config.localToGlobal ( El, dof )] += localVec[dof] * Factor;
//     }
//   }
// 
// protected:
//   // barton-nackman
//   inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
//   inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
//   
//   const ConfiguratorType &_conf;
// };


//! Integrator for \f$ (\int_\Omega s(x)  w_i(x) \da )_{i} \f$, of some scalar valued function \f$ s\f$.
template <typename ConfiguratorType, typename Imp>
class RefTriangleFENonlinOpIntegrator  {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
  const ConfiguratorType & _config;
  
public:
  RefTriangleFENonlinOpIntegrator( const ConfiguratorType & Conf ) : _config( Conf){ }

  virtual ~RefTriangleFENonlinOpIntegrator( ) {}

  void assembleAdd ( VectorType &Dest ) const {
    RealType *nl_cache = new RealType[ _config.maxNumQuadPoints() ];
    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
      const int numLocalDofs = _config.getNumLocalDofs ( El );
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );
      for ( int q = 0; q < numQuadPoints; ++q ) nl_cache[q] = this->asImp().getNonlinearity ( El, q );
      RealType aux;
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux = 0.;   
        for ( int q = 0; q < numQuadPoints; ++q ) aux += nl_cache[q] * bfs.evaluateOnRefTriang( dof, q ) * bfs.getWeight ( q );
        Dest[ _config.localToGlobal ( El, dof ) ] += El.getAreaOfRefTriangle() * aux;
      }
    }
    delete[] nl_cache;
  }
  
  void assembleDirichlet ( VectorType &Dest, const typename ConfiguratorType::MaskType &boundaryMask ) const {
     Dest.setZero();
     this->assembleAdd(Dest);
     for( int i=0; i<Dest.size(); ++i ){ if( boundaryMask[i] ) Dest[i] = 0.0; }
  }

  //! interface function, has to be provided in derived classes.
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};



//! Integrator for \f$ (\int_\Omega A(x) \nabla w_i(x) \da )_{ij} \f$, of some matrix valued function \f$ A\f$.
template <typename ConfiguratorType, typename Imp>
class RefTriangleFENonlinDiffOpIntegrator  {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  const ConfiguratorType & _config;
  
public:
  RefTriangleFENonlinDiffOpIntegrator ( const ConfiguratorType & Conf ) : _config( Conf ){ }

  virtual ~RefTriangleFENonlinDiffOpIntegrator( ) {}

  void assembleAdd ( VectorType &Dest ) const {

    DomVecType *nl_cache = new DomVecType[ _config.maxNumQuadPoints() ];

    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
      const int numLocalDofs = _config.getNumLocalDofs ( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q )
        this->asImp().getNonlinearity (  El, q, nl_cache[q] );

      RealType aux;
  
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux = 0.0;   
        
        for ( int q = 0; q < numQuadPoints; ++q )
          aux += bfs.getWeight ( q ) * ( nl_cache[q].dot(bfs.evaluateGradientOnRefTriang( dof, q ) ) );

        Dest[ _config.localToGlobal ( El, dof ) ] += El.getAreaOfRefTriangle() * aux;
      }
    }
    delete[] nl_cache;
  }
  
  void assembleDirichlet ( VectorType &Dest, const typename ConfiguratorType::MaskType &boundaryMask ) const {
     Dest.setZero();
     this->assembleAdd(Dest);
     for( int i=0; i<Dest.size(); ++i ){
       if( boundaryMask[i] ) Dest[i] = 0.0;
     }
  }

  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint, DomVecType &NL ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
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

//! General Interface for vector valued integrators
// template < typename ConfiguratorType, typename Imp >
// class MultiVectorValuedIntegratorBase {
// public:
//   typedef typename ConfiguratorType::RealType   RealType;
//   typedef typename ConfiguratorType::ElementType ElementType;
// 
//   explicit MultiVectorValuedIntegratorBase ( const ConfiguratorType &conf ) : _conf ( conf ) {}
// 
//   template <typename MultiVectorType>
//   void assembleAdd ( MultiVectorType &Dest, const RealType Factor = 1.0 ) const {
//     MultiVectorType localMultiVec ( ConfiguratorType::maxNumLocalDofs, 3 );
//     for ( ElementIteratorType it = _conf.begin(); it != _conf.end(); ++it ) {
//       this->asImp().prepareLocalMultiVec ( El, localMultiVec );
//       for( int comp=0; comp<3; ++comp )
//         for ( int dof = 0; dof < _conf.getNumLocalDofs ( El ); ++dof )
//           Dest[comp][_conf.localToGlobal ( El, dof )] += localMultiVec[dof][comp] * Factor;
//     }
//   }
// 
// protected:
//   // barton-nackman
//   inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
//   inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
//   
//   const ConfiguratorType &_conf;
// };

//! Integrator for \f$ (\int_\Omega (v(x) \cdot w_i(x)) \da )_{i} \f$, of some vector valued function \f$v\f$. \f$v\f$  is supposed to be of dimension 3.
template <typename ConfiguratorType, typename Imp>
class RefTriangleFENonlinVectorOpIntegrator  {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;
  const ConfiguratorType & _config;
  
public:
  RefTriangleFENonlinVectorOpIntegrator ( const ConfiguratorType & Conf ) : _config( Conf){ }

  virtual ~RefTriangleFENonlinVectorOpIntegrator( ) {}

  void assembleAdd ( typename ConfiguratorType::VectorType &Dest ) const {
    
    TangentVecType *nl_cache = new TangentVecType[ _config.maxNumQuadPoints() ];
    const int numGlobalDofs = _config.getNumGlobalDofs();
    
    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
      const int numLocalDofs = _config.getNumLocalDofs ( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q )
        this->asImp().getNonlinearity ( El, q, nl_cache[q] );

      TangentVecType aux;
  
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux.setZero();    

        for ( int q = 0; q < numQuadPoints; ++q )
          aux += bfs.evaluateOnRefTriang( dof, q ) * bfs.getWeight ( q ) * nl_cache[q];

        for ( int comp = 0; comp < 3; ++comp )
          Dest[ _config.localToGlobal ( El, dof ) + comp * numGlobalDofs ] += El.getAreaOfRefTriangle() * aux[comp];
      }
    }
    delete[] nl_cache;
  }

  void assembleDirichlet ( typename ConfiguratorType::VectorType &Dest, const MaskType& boundaryMask ) const {
      Dest.setZero();
      this->assembleAdd( Dest );
      const int numGlobalDofs = _config.getNumGlobalDofs();
      for( int i = 0; i < numGlobalDofs; ++i ){
          if ( boundaryMask[i] ){
            for( int comp=0; comp<3; ++comp ) Dest[i + comp * numGlobalDofs] = 0.0;
          }
      } 
  }
  
  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint, TangentVecType &NL ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    this->asImp().getNonlinearity (  El, QuadPoint, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};



//! Integrator for \f$ (\int_\Omega (A(x) : D w_i(x)) \da )_{ij} \f$, of some matrix valued function \f$A\f$ supposed to be a 3x2 - matrix.
template <typename ConfiguratorType, typename Imp>
class RefTriangleMVDiffOpIntegrator  {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  const ConfiguratorType & _config;
  
public:
  RefTriangleMVDiffOpIntegrator ( const ConfiguratorType & Conf ) :  _config( Conf ) { }

  virtual ~RefTriangleMVDiffOpIntegrator( ) {}

  void assembleAdd ( VectorType &Dest ) const {       
    
    const int numGlobalDofs = _config.getNumGlobalDofs();
    Matrix32 *nl_cache = new Matrix32[ _config.maxNumQuadPoints() ];

    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
      const int numLocalDofs = _config.getNumLocalDofs( El );
      
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      // pre-assemble for all quadrature points (for efficiency)
      for ( int q = 0; q < numQuadPoints; ++q )
        this->asImp().getNonlinearity ( El, q, nl_cache[q] );
      
      TangentVecType aux;
      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        aux.setZero();
        for ( int q = 0; q < numQuadPoints; ++q ) {
          DomVecType gradbf = bfs.evaluateGradientOnRefTriang( dof, q );
          Matrix32 nl = bfs.getWeight ( q ) * nl_cache[q];
          aux += nl * gradbf;
        }

        for ( int d = 0; d < 3; ++d )
          Dest[ _config.localToGlobal ( El, dof ) + d * numGlobalDofs ] += El.getAreaOfRefTriangle() * aux[d];
      }
    }
    delete[] nl_cache;

  }

  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint, Matrix32 &NL ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    this->asImp().getNonlinearity (  El, QuadPoint,  NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};



//! Integrator for \f$ (\int_\Omega (A(x) : D^2 w_i(x))_j \da )_{ij} \f$, \$ j = 0,1,2 \$ where \f$ A \f is supposed to be a 3x2x2 tensor.
template <typename ConfiguratorType, typename Imp>
class RefTriangleMVDiff2OpIntegrator  {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::VectorType VectorType;
  const ConfiguratorType & _config;

public:
  RefTriangleMVDiff2OpIntegrator ( const ConfiguratorType & Conf ) :  _config( Conf ) { }

  virtual ~RefTriangleMVDiff2OpIntegrator( ) {}

  void assembleAdd ( VectorType &Dest ) const {   
    const int numGlobalDofs = _config.getNumGlobalDofs();
    Tensor322Type *nl_cache = new Tensor322Type[ _config.maxNumQuadPoints() ];

    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
      const int numLocalDofs = _config.getNumLocalDofs( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const typename ConfiguratorType::ApproxGradientBaseFuncSetType &approxBfs = _config.getApproxGradientBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      // pre-assemble for all quadrature points (for efficiency)
      for ( int q = 0; q < numQuadPoints; ++q )
        this->asImp().getNonlinearity ( El, q, nl_cache[q] );
      
      TangentVecType aux;
      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        aux.setZero();
        for ( int q = 0; q < numQuadPoints; ++q ) {
          Matrix22 hessian; hessian = approxBfs.evaluateApproxHessianSymOnRefTriang( dof, q ); 
          Tensor322Type nl ( nl_cache[q] );
          
          TangentVecType tmp;
          for( int i=0; i<3; i++ ) tmp[i] = pesopt::ddProd<RealType,Matrix22>(nl[i], hessian );
           
          tmp *= bfs.getWeight ( q );
          aux += tmp;
        }
        for ( int d = 0; d < 3; ++d )
          Dest[ _config.localToGlobal ( El, dof ) + d * numGlobalDofs ] += El.getAreaOfRefTriangle() * aux[d];

      }
    }
    delete[] nl_cache;

  }

  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint,  Tensor322Type &NL ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    this->asImp().getNonlinearity ( El, QuadPoint,  NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};





//! provides an interface for Finite Element operators of the form \f$ A(x) : (\phi_1 \phi_3 | \phi_3 \phi_2 ) \f$, where \f$A\f$ is a SYMMETRIC coefficient matrix.
// template <typename ConfiguratorType, typename Imp >
// class RefTriangleFELinSymMatrixWeightedVecIntegrator {
// protected:
//   typedef typename ConfiguratorType::RealType RealType;
//   typedef typename ConfiguratorType::Matrix22 Matrix22;
//   typedef typename ConfiguratorType::TangentVecType TangentVecType;
//   typedef typename ConfiguratorType::ElementType ElementType;
//   const ConfiguratorType &_config;
//   
// public:
//   RefTriangleFELinSymMatrixWeightedVecIntegrator ( const ConfiguratorType & Config ) :
//   _config ( Config ) {}
// 
//    virtual ~RefTriangleFELinSymMatrixWeightedVecIntegrator( ) {}
// 
//   void assembleAdd ( pesopt::MultiVector<RealType, typename ConfiguratorType::VectorType> &Dest ) const {       
//     
//     Matrix22 *nl_cache = new Matrix22[ _config.maxNumQuadPoints() ];
// 
//     for ( ElementIteratorType it = _config.begin(); it != _config.end(); ++it ) {
//       const int numLocalDofs = _config.getNumLocalDofs( El );
//       
//       const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
//       const int numQuadPoints = bfs.numQuadPoints( );
// 
//       // pre-assemble for all quadrature points (for efficiency)
//       for ( int q = 0; q < numQuadPoints; ++q )
//         this->asImp().getNonlinearity ( El, q, nl_cache[q] );
//       
//       TangentVecType aux;
//       for ( int dof = 0; dof < numLocalDofs; dof++ ) {
//         aux.setZero();
//         for ( int q = 0; q < numQuadPoints; ++q ) {
//           RealType valuebf = bfs.evaluateOnRefTriang( dof, q ); 
//           Matrix22 nl;
//           nl = nl_cache[q];
//           
//           TangentVecType tmp;
//           tmp[0] = valuebf * nl[0][0];
//           tmp[1] = valuebf * nl[1][1];
//           tmp[2] = 2.0 * valuebf * nl[1][0];
//            
//           tmp *= bfs.getWeight ( q );
//           aux += tmp;
//         }
// 
//         for ( int d = 0; d < 3; ++d )
//           Dest[d][ _config.localToGlobal ( El, dof ) ] += 0.5 * aux[d];
//       }
//     }
//     delete[] nl_cache;
// 
//   }
// 
//   //! interface function, has to be provided in derived classes.
//   void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint,
//                          Matrix22 &NL ) const {
//     throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
//     this->asImp().getNonlinearity ( El, QuadPoint,  NL );
//   }
// 
// protected:
//   inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
//   inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }
// 
// };



//!===========================================================================================================================
//! Matrix-Valued Intefaces
//!===========================================================================================================================


//! General Interface for matrix valued integrators
template < typename ConfiguratorType, typename Imp >
class MatrixValuedIntegratorBase {
public:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;

  explicit MatrixValuedIntegratorBase ( const ConfiguratorType &conf ) : _config ( conf ) {}

public:
  void assembleTripletList ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
    tripletList.reserve(_config.getInitializer().getNumTriangs() * pesopt::Sqr( _config.getNumLocalDofs() ) );
    LocalMatrixType localMatrix;
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
      this->asImp().prepareLocalMatrix ( El, localMatrix );
      const int numLocalDofs = _config.getNumLocalDofs ( El );

      for ( int i = 0; i < numLocalDofs; ++i ) globalDofs[ i ] = _config.localToGlobal ( El, i );

      for ( int i = 0; i < numLocalDofs; ++i ) {
        int glob_i = globalDofs[ i ];
        for ( int j = 0; j < numLocalDofs; ++j ) {
          int glob_j = globalDofs[ j ];
          tripletList.push_back( TripletType( glob_i, glob_j, El.getAreaOfRefTriangle() * Factor * localMatrix(i,j) ) );
        }
      }
    }
  }
  
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor ) const {
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );
      
    tripletListMasked.reserve(_config.getInitializer().getNumTriangs() * pesopt::Sqr( _config.getNumLocalDofs() ) );

    for( unsigned iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row()]) || (boundaryMask[tripletList[iter].col()]) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _config.getNumGlobalDofs(); ++i ){
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
    tripletListMasked.reserve(_config.getInitializer().getNumTriangs() * pesopt::Sqr( _config.getNumLocalDofs() ) );

    for( unsigned iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row()]) || (boundaryMask[tripletList[iter].col()]) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _config.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] )
         tripletListMasked.push_back( TripletType( i, i, 1.0 ) );
    }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
  const ConfiguratorType &_config;
};




//! \brief provides an interface to Finite Element operators of the form
//! The corresponding matrix assembly yields \f$ \left(\int_\Omega  A(x)\nabla\phi_j\cdot \nabla\phi_i dx\right)_{ij} \f$ 
//! for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinAsymMatrixWeightedStiffIntegrator :
      public MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinAsymMatrixWeightedStiffIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinAsymMatrixWeightedStiffIntegrator ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinAsymMatrixWeightedStiffIntegrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, Matrix22 &Matrix ) const {
    this->asImp().getCoeffMatrix ( El, QuadPoint, Matrix );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                  LocalMatrixType &LocalMatrix ) const {
    const int numLocalDofs = _config.getNumLocalDofs ( El );
    
    for ( int i = 0; i < numLocalDofs; ++i )
      for ( int j = 0; j < numLocalDofs; ++j )
        LocalMatrix(i,j) = 0.;

    Matrix22 mat;
    DomVecType MatGrad1;
    DomVecType gradbf [numLocalDofs];

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffMatrix ( El, q, mat );
      
      for( int locIdx = 0; locIdx < numLocalDofs; ++locIdx ) 
          gradbf[locIdx] = bfs.evaluateGradientOnRefTriang( locIdx, q );
      
      for ( int i = 0; i < numLocalDofs; ++i ) {
          MatGrad1 = mat * gradbf[i];
        for ( int j = 0; j < numLocalDofs; ++j ) {
          LocalMatrix(j,i) += MatGrad1.dot( gradbf[j] ) * bfs.getWeight ( q );
        }
      }
    }
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};




//! \brief provides an interface to Finite Element operators of the form
//! The corresponding matrix assembly yields \f$ \left(\int_\Omega  A(x)\nabla_h \phi_j\cdot \nabla_h \phi_i dx\right)_{ij} \f$ 
//! for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinAsymMatrixWeightedApproxStiffIntegrator :
      public MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinAsymMatrixWeightedApproxStiffIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinAsymMatrixWeightedApproxStiffIntegrator ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinAsymMatrixWeightedApproxStiffIntegrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                               Matrix22 &Matrix ) const {
    this->asImp().getCoeffMatrix ( El, QuadPoint, Matrix );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                  LocalMatrixType &LocalMatrix ) const {
    const int numLocalDofs = _config.getNumLocalDofs ( El );
    
    for ( int i = 0; i < numLocalDofs; ++i )
      for ( int j = 0; j < numLocalDofs; ++j )
        LocalMatrix(i,j) = 0.;

    Matrix22 mat;
    DomVecType MatGrad1;
    DomVecType gradbf [numLocalDofs];

    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &approxBfs = _config.getApproxGradientBaseFunctionSet( El );
    const int numQuadPoints = approxBfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffMatrix ( El, q, mat );
      
      for( int locIdx = 0; locIdx < numLocalDofs; ++locIdx ) gradbf[locIdx] = approxBfs.evaluateApproxGradientOnRefTriang( locIdx, q );
      
      for ( int i = 0; i < numLocalDofs; ++i ) {
          MatGrad1 = mat * gradbf[i];
        for ( int j = 0; j < numLocalDofs; ++j ) {
          LocalMatrix(j,i) += MatGrad1.dot( gradbf[j] ) * approxBfs.getWeight ( q );
        }
      }
    }
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};



//! The corresponding matrix assembly yields \f$ \left(\int_\Omega  w(x) \phi_i \phi_j dx\right)_{ij} \f$ for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinWeightedMassIntegrator :
      public MatrixValuedIntegratorBase<  ConfiguratorType, RefTriangleFELinWeightedMassIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinWeightedMassIntegrator ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase<  ConfiguratorType, RefTriangleFELinWeightedMassIntegrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                  LocalMatrixType &LocalMatrix ) const {
    const int numDofs = _config.getNumLocalDofs ( El );	
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    RealType nonlinearity;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    //TODO cache b_i, b_j? use symmetry?
    for ( int q = 0; q < numQuadPoints; ++q ) {
      nonlinearity = getNonlinearity ( El, q );
      for ( int i = 0; i < numDofs; ++i ) {
          RealType b_i = bfs.evaluateOnRefTriang( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          RealType b_j = bfs.evaluateOnRefTriang( j, q );
          LocalMatrix(j,i) += nonlinearity * b_i * b_j * bfs.getWeight ( q );
        }
      }
    }
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};


//! The corresponding matrix assembly yields \f$ \left(\int_\Omega  w(x) nabla\phi_i nabla\phi_j dx\right)_{ij} \f$ for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinScalarWeightedStiffIntegrator :
      public MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinScalarWeightedStiffIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinScalarWeightedStiffIntegrator ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase<  ConfiguratorType, RefTriangleFELinScalarWeightedStiffIntegrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                  LocalMatrixType &LocalMatrix ) const {
    const int numDofs = _config.getNumLocalDofs ( El );	
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    RealType nonlinearity;;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    //TODO cache b_i, b_j? use symmetry?
    for ( int q = 0; q < numQuadPoints; ++q ) {
      nonlinearity = getNonlinearity ( El, q );
      for ( int i = 0; i < numDofs; ++i ) {
          DomVecType grad_i = bfs.evaluateGradientOnRefTriang( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          DomVecType grad_j = bfs.evaluateGradientOnRefTriang( j, q );
          LocalMatrix(j,i) += nonlinearity * grad_i.dot(grad_j) * bfs.getWeight ( q );
        }
      }
    }
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};









// The corresponding matrix assembly yields \f$ \left(\int_\Omega v(x) \nabla\phi_i \cdot \phi_j dx\right)_{ij} \f$ for FE basis functions \f$ \phi_i,\phi_j \f$, where v is a vector in R^2
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinScalarWeightedSemiDiffIntegrator :
      public MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinScalarWeightedSemiDiffIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinScalarWeightedSemiDiffIntegrator ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinScalarWeightedSemiDiffIntegrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffVector ( const typename ConfiguratorType::ElementType &El, int QuadPoint, DomVecType &Vector ) const {
    this->asImp().getCoeffVector ( El, QuadPoint, Vector );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                   LocalMatrixType &LocalMatrix ) const {
    
    const int numDofs = _config.getNumLocalDofs ( El ); 
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    DomVecType vec;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffVector ( El, q, vec );
      for ( int i = 0; i < numDofs; ++i ) {
        DomVecType grad_i = bfs.evaluateGradientOnRefTriang( i, q );
        RealType vec_gradi = vec.dot(grad_i);
        for ( int j = 0; j < numDofs; ++j )
          LocalMatrix(j,i) += vec_gradi * bfs.evaluateOnRefTriang( j, q ) * bfs.getWeight ( q );
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};



//! The corresponding matrix assembly yields \f$ \left(\int_\Omega v1(x) \cdot \nabla\phi_j  v2(x) \cdot \nabla\phi_i dx\right)_{ij} \f$ for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinSeperatedVectorWeightedStiffIntegrator :
      public MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinSeperatedVectorWeightedStiffIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinSeperatedVectorWeightedStiffIntegrator ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase<  ConfiguratorType, RefTriangleFELinSeperatedVectorWeightedStiffIntegrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffVectors ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                                DomVecType &vec1,  DomVecType &vec2 ) const {
    this->asImp().getCoeffVectors ( El, QuadPoint, vec1, vec2 );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                   LocalMatrixType &LocalMatrix ) const {
    
    const int numDofs = _config.getNumLocalDofs ( El ); 
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    DomVecType vec1, vec2;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffVectors ( El, q, vec1, vec2 );
      for ( int i = 0; i < numDofs; ++i ) {
        DomVecType grad_i = bfs.evaluateGradientOnRefTriang( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          DomVecType grad_j = bfs.evaluateGradientOnRefTriang( j, q );
          LocalMatrix(j,i) += ( vec2.dot(grad_j) ) * ( vec1.dot(grad_i) ) * bfs.getWeight ( q );
        }
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};





//! provides an interface to Finite Element operators of the form \f$ (S(x) : D^2 u_j) v(x) \cdot D u_i \f$, where \f$S\f$ is 2x2 Matrix, v is a vector in R^2.
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator :
      public MatrixValuedIntegratorBase<  ConfiguratorType, RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffMatrixAndVector ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                                        DomVecType &Vector,  Matrix22 &Matrix ) const {
    this->asImp().getCoeffMatrixAndVector ( El, QuadPoint, Vector, Matrix );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                   LocalMatrixType &LocalMatrix ) const {
    const int numDofs = _config.getNumLocalDofs ( El );    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;
    DomVecType vec;
    Matrix22 mat;
    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &approxBfs = _config.getApproxGradientBaseFunctionSet( El );
    const int numQuadPoints = bfs.numQuadPoints( );
    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffMatrixAndVector ( El, q, vec, mat );
      for ( int i = 0; i < numDofs; ++i ) {
        DomVecType grad_i = bfs.evaluateGradientOnRefTriang( i, q );
        RealType vecgradi = grad_i.dot( vec );
        for ( int j = 0; j < numDofs; ++j ) {
          Matrix22 hessian_j = approxBfs.evaluateApproxHessianSymOnRefTriang( j, q );
          RealType mathessian_j = pesopt::ddProd<RealType, Matrix22> ( mat, hessian_j );
          LocalMatrix(j,i) += vecgradi * mathessian_j * bfs.getWeight ( q );
        }
      }
    }
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};


//! provides an interface to Finite Element operators of the form \f$ (T(x) : D^2 u_i \otimes D u_j \f$, where \f$T\f$ is 2x2x2 Tensor
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator_Tensor :
      public MatrixValuedIntegratorBase<  ConfiguratorType, RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator_Tensor<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Tensor222Type Tensor222Type;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator_Tensor ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase<  ConfiguratorType, RefTriangleFELinMatrixMixedFirstSecondDiffIntegrator_Tensor<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}

  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffTensor ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                               Tensor222Type &Tensor ) const {
    this->asImp().getCoeffTensor( El, QuadPoint, Tensor );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                   LocalMatrixType &LocalMatrix ) const {
    
    const int numDofs = _config.getNumLocalDofs ( El );      
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    Tensor222Type tensor;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &approxBfs = _config.getApproxGradientBaseFunctionSet( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffTensor ( El, q, tensor );
      for ( int i = 0; i < numDofs; ++i ) {
        DomVecType grad_i = bfs.evaluateGradientOnRefTriang( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          Matrix22 hessian_j = approxBfs.evaluateApproxHessianSymOnRefTriang( j, q );
          RealType aux = 0.0;
          for( int a=0; a<2; a++ )
              for( int b=0; b<2; b++ )
                  for( int c=0; c<2; c++ )
                      aux += tensor.get(a,b,c) * hessian_j(a,b) * grad_i( c );    
          LocalMatrix(j,i) += aux * bfs.getWeight ( q );
        }
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};




//! provides an interface to Finite Element operators of the form \f$ (s(x) D^2 u_j : D^2 u_i) \f$, where \f$s\f$ is scalar valued function
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinWeightedDiff2Integrator :
      public MatrixValuedIntegratorBase<  ConfiguratorType, RefTriangleFELinWeightedDiff2Integrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinWeightedDiff2Integrator ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinWeightedDiff2Integrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return this->asImp().getCoeff ( El, QuadPoint );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                   LocalMatrixType &LocalMatrix ) const {
    
    
    const int numDofs = _config.getNumLocalDofs ( El );    
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    RealType scalar;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &approxBfs = _config.getApproxGradientBaseFunctionSet( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      scalar = getCoeff( El, q );
      for ( int i = 0; i < numDofs; ++i ) {
        Matrix22 hessian_i = approxBfs.evaluateApproxHessianSymOnRefTriang( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          Matrix22 hessian_j = approxBfs.evaluateApproxHessianSymOnRefTriang( j, q );
          RealType aux = pesopt::ddProd<RealType, Matrix22> ( hessian_i, hessian_j );
         
          LocalMatrix(j,i) += scalar * aux * bfs.getWeight ( q );
        }
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};


//! provides an interface to Finite Element operators of the form \f$ (A(x) D^2 u_j : D^2 u_i) \f$, where \f$A\f$ is 2x2 Matrix
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinMatrixDiff2Integrator :
      public MatrixValuedIntegratorBase<  ConfiguratorType, RefTriangleFELinMatrixDiff2Integrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinMatrixDiff2Integrator ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase<  ConfiguratorType, RefTriangleFELinMatrixDiff2Integrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, Matrix22 &Matrix ) const {
    this->asImp().getCoeffMatrix ( El, QuadPoint, Matrix );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El,  LocalMatrixType &LocalMatrix ) const {
    
    const int numDofs = _config.getNumLocalDofs ( El );        
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    Matrix22 mat;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &approxBfs = _config.getApproxGradientBaseFunctionSet( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffMatrix( El, q, mat );
      for ( int i = 0; i < numDofs; ++i ) {
        Matrix22 hessian_i = approxBfs.evaluateApproxHessianSymOnRefTriang( i, q );
        Matrix22 mathessiani ( mat );
        mathessiani *= hessian_i;
        for ( int j = 0; j < numDofs; ++j ) {
          Matrix22 hessian_j = approxBfs.evaluateApproxHessianSymOnRefTriang( j, q );
          RealType aux = pesopt::ddProd<RealType, Matrix22> ( mathessiani, hessian_j );
         
          LocalMatrix(j,i) += aux * bfs.getWeight ( q );
        }
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};



//! provides an interface to Finite Element operators of the form \f$ (A(x) : D^2 u_j  B(x) D^2 u_i) \f$, where \f$A,B\f$ are 2x2 Matrix
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinMatrixSeperatedDiff2Integrator :
      public MatrixValuedIntegratorBase<  ConfiguratorType, RefTriangleFELinMatrixSeperatedDiff2Integrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinMatrixSeperatedDiff2Integrator ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinMatrixSeperatedDiff2Integrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}

  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffMatrices ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                                 Matrix22 &MatrixA, Matrix22 &MatrixB ) const {
    this->asImp().getCoeffMatrices ( El, QuadPoint, MatrixA, MatrixB );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                   LocalMatrixType &LocalMatrix ) const {
   
    const int numDofs = _config.getNumLocalDofs ( El );        
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    Matrix22 matA, matB;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &approxBfs = _config.getApproxGradientBaseFunctionSet( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffMatrices( El, q, matA, matB );
      for ( int i = 0; i < numDofs; ++i ) {
        Matrix22 hessian_i = approxBfs.evaluateApproxHessianSymOnRefTriang( i, q );
        RealType auxA = pesopt::ddProd<RealType, Matrix22> ( matA, hessian_i );
        for ( int j = 0; j < numDofs; ++j ) {
          Matrix22 hessian_j = approxBfs.evaluateApproxHessianSymOnRefTriang( j, q );
          RealType auxB = pesopt::ddProd<RealType, Matrix22>( matB, hessian_j );
         
          LocalMatrix(j,i) += auxA * auxB * bfs.getWeight ( q );
        }
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};



//! provides an interface to Finite Element operators of the form \f$ \int w(x) (<v_1(x), D^2 u_j>  <v_2(x), D^2 u_i) \f$, where D^2 u = (d_11, d_22, d_12 + d_21)
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinMatrixSeperatedDiff2AsVecIntegrator :
      public MatrixValuedIntegratorBase<  ConfiguratorType, RefTriangleFELinMatrixSeperatedDiff2AsVecIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinMatrixSeperatedDiff2AsVecIntegrator ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinMatrixSeperatedDiff2AsVecIntegrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}

  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffVecs ( const typename ConfiguratorType::ElementType &El, int QuadPoint, RealType &w, TangentVecType &vecA, TangentVecType &vecB ) const {
    this->asImp().getCoeffVecs ( El, QuadPoint, w, vecA, vecB );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, LocalMatrixType &LocalMatrix ) const {
   
    const int numDofs = _config.getNumLocalDofs ( El );        
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    TangentVecType vecA, vecB;
    RealType w;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &approxBfs = _config.getApproxGradientBaseFunctionSet( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffVecs( El, q, w, vecA, vecB );
      for ( int i = 0; i < numDofs; ++i ) {
        TangentVecType hessian_i = approxBfs.evaluateApproxHessianAsVecOnRefTriang( i, q );
        RealType auxA = vecA.dot(hessian_i);
        for ( int j = 0; j < numDofs; ++j ) {
          TangentVecType hessian_j = approxBfs.evaluateApproxHessianAsVecOnRefTriang( j, q );
          RealType auxB = vecB.dot(hessian_j);
          LocalMatrix(j,i) += w * auxA * auxB * bfs.getWeight ( q );
        }
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};


//! provides an interface to Finite Element operators of the form \f$ \int A(x) D^2 u_j \cdot D^2 u_i \f$, where D^2 u = (d_11, d_22, d_12 + d_21), A is a 3x3 matrix
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinMatrixDiff2AsVecIntegrator :
      public MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinMatrixDiff2AsVecIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::Point3DType Point3DType;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinMatrixDiff2AsVecIntegrator ( const ConfiguratorType & Config ) : 
   MatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinMatrixDiff2AsVecIntegrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}

  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, Matrix33 &A ) const {
    this->asImp().getCoeffMatrix ( El, QuadPoint, A );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, LocalMatrixType &LocalMatrix ) const {
   
    const int numDofs = _config.getNumLocalDofs ( El );        
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    Matrix33 A;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const typename ConfiguratorType::ApproxGradientBaseFuncSetType &approxBfs = _config.getApproxGradientBaseFunctionSet( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffMatrix( El, q, A );
      for ( int i = 0; i < numDofs; ++i ) {
        Point3DType hessian_i = approxBfs.evaluateApproxHessianAsVecOnRefTriang( i, q );
        Point3DType A_hessian_i = A * hessian_i;
        for ( int j = 0; j < numDofs; ++j ) {
          Point3DType hessian_j = approxBfs.evaluateApproxHessianAsVecOnRefTriang( j, q );
          LocalMatrix(j,i) += approxBfs.getWeight ( q ) * A_hessian_i.dot(hessian_j);
        }
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};






//!===========================================================================================================================
//! BlockMatrix-Valued Intefaces
//!===========================================================================================================================

//! General Interface for matrix valued integrators
template < typename ConfiguratorType, typename Imp, int NumVecCompsArg = 3, int NumVecCompsDest = 3>
class BlockMatrixValuedIntegratorBase {
public:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;
  typedef typename ConfiguratorType::TripletType TripletType;

  explicit BlockMatrixValuedIntegratorBase ( const ConfiguratorType &conf ): _config ( conf ) {}

protected:
    void assembleTripletList ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
        tripletList.reserve( NumVecCompsArg * NumVecCompsDest * pesopt::Sqr( _config.getNumLocalDofs() ) *_config.getInitializer ().getNumTriangs ());
        LocalMatrixType localMatrix[NumVecCompsArg][NumVecCompsDest];
        const int numGlobalDofs = _config.getNumGlobalDofs();
        int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
        for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumTriangs(); ++elementIdx){
            const ElementType& El ( _config.getInitializer().getTriang( elementIdx ) );
            // assemble the local matrix for the current element
            this->asImp().prepareLocalMatrix ( El, localMatrix );

            const int numLocalDofs = _config.getNumLocalDofs ( El );

            for ( int i = 0; i < numLocalDofs; ++i )
                globalDofs[ i ] = _config.localToGlobal ( El, i );
            
            for ( int argComp = 0; argComp < NumVecCompsArg; ++argComp )
                for ( int destComp = 0; destComp < NumVecCompsDest; ++destComp )
                 for ( int i = 0; i < numLocalDofs; ++i ) {
                    int glob_i = globalDofs[ i ];
                    for ( int j = 0; j < numLocalDofs; ++j ) {
                      int glob_j = globalDofs[ j ];
                      tripletList.push_back( TripletType( glob_i + destComp * numGlobalDofs, glob_j + argComp * numGlobalDofs, El.getAreaOfRefTriangle() * Factor * localMatrix[argComp][destComp]( i, j ) ) );
                    }
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
  void assembleDirichlet ( BlockMatrixType &Dest, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );

    // Boundary Mask
    const int numGlobalDofs = _config.getNumGlobalDofs();
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( NumVecCompsArg * NumVecCompsDest * pesopt::Sqr( _config.getNumLocalDofs() ) *_config.getInitializer ().getNumTriangs ());

    for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row() % numGlobalDofs] ) || (boundaryMask[tripletList[iter].col() % numGlobalDofs] ) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _config.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] ){
        for ( int Comp = 0; Comp < 3; ++Comp )
            tripletListMasked.push_back( TripletType( i + Comp * numGlobalDofs, i + Comp * numGlobalDofs, 1.0 ) );
      }
    }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
    
  }
  
  //TODO so far only for LagrangianOp
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );

    // Boundary Mask
    const int numGlobalDofs = _config.getNumGlobalDofs();
    tripletListMasked.reserve( NumVecCompsArg * NumVecCompsDest * pesopt::Sqr( _config.getNumLocalDofs() ) *_config.getInitializer ().getNumTriangs ());

    for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row() % numGlobalDofs] ) || (boundaryMask[tripletList[iter].col() % numGlobalDofs] ) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _config.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] ){
        for ( int Comp = 0; Comp < 3; ++Comp )
            tripletListMasked.push_back( TripletType( i + Comp * numGlobalDofs, i + Comp * numGlobalDofs, 1.0 ) );
      }
    }
    
  }
  

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
  const ConfiguratorType &_config;
};




//! The corresponding matrix assembly yields \f$ \left(\int_\Omega  w(x) \phi_i \phi_j dx\right)_{ij} \f$ for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class RefTriangleFELinWeightedBlockMassIntegrator :
      public BlockMatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinWeightedBlockMassIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_config;
  
public:
  RefTriangleFELinWeightedBlockMassIntegrator ( const ConfiguratorType & Config ) : 
   BlockMatrixValuedIntegratorBase< ConfiguratorType, RefTriangleFELinWeightedBlockMassIntegrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }

  void prepareLocalMatrix( const typename ConfiguratorType::ElementType &El, LocalMatrixType (&localMatrix)[3][3]) const {
      
      for (int argComp = 0; argComp < 3; ++argComp)
          for (int destComp = 0; destComp < 3; ++destComp)
              localMatrix[argComp][destComp].setZero();
            
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet(El);
      const int numDofs = _config.getNumLocalDofs ( El );
      RealType nonlinearity = 0.;    
      
      for (int quadPoint = 0; quadPoint < _config.maxNumQuadPoints(); ++quadPoint) {
          nonlinearity = getNonlinearity ( El, quadPoint );
          for ( int i = 0; i < numDofs; ++i ) {
              RealType b_i = bfs.evaluateOnRefTriang( i, quadPoint );
              for ( int j = 0; j < numDofs; ++j ) {
                  RealType b_j = bfs.evaluateOnRefTriang( j, quadPoint );
                  for(int argComp=0; argComp < 3; ++argComp)
                      localMatrix[argComp][argComp](j,i) += nonlinearity * b_i * b_j * bfs.getWeight ( quadPoint );
              }
          }
      }
        
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};


#endif
