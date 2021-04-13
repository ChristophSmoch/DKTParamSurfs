#ifndef __OPTIMALDEFORMISOLAGRANGIAN_H
#define __OPTIMALDEFORMISOLAGRANGIAN_H

#include <pesopt_IO.h>
#include <pesopt_DKT.h>

#include <energyDefines.h>
#include <solverInfo.h>





//-----------------------------------------------------------------
//-----------------------------------------------------------------
//              Pointwise Isometry Constraint Op and Derivatives
//
//              I(w, lambda ) = sum_nodes (g_B - g_A) : lambda
//              where lambda = ( l_1 l_3 | 
//                               l_3 l_2 )
//-----------------------------------------------------------------
//-----------------------------------------------------------------
template< typename ConfiguratorType >
class ShellIsometryHandler_Pointwise{
  
    typedef typename ConfiguratorType::DTContainer DataTypeContainer;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    typedef typename ConfiguratorType::MaskType MaskType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const VectorType& _xA;
    const MaskType & _mask;
    const int _numVertices, _numGlobalDofs, _numBoundaryNodes, _numInteriorNodes;
    mutable std::vector<Eigen::Ref<const VectorType> > _refsD1X, _refsD2X;
    
protected :
    void getRefsToDisp( const VectorType & x,
                        std::vector<Eigen::Ref<const VectorType> > & refsD1, 
                        std::vector<Eigen::Ref<const VectorType> > & refsD2 ) const {
         for( int comp=0; comp<3; ++comp ){
           refsD1.push_back ( x.segment(     _numVertices + comp * _numGlobalDofs, _numVertices) );
           refsD2.push_back ( x.segment( 2 * _numVertices + comp * _numGlobalDofs, _numVertices) );
         }
    }
    
    void getRefsToMult( const VectorType & x, std::vector<Eigen::Ref<const VectorType> > & refsMult ) const {
         for( int comp=0; comp<3; ++comp ) refsMult.push_back ( x.segment( comp * _numInteriorNodes, _numInteriorNodes) );
    }
  
public:
  
    ShellIsometryHandler_Pointwise ( const ConfiguratorType &conf,
                                     const VectorType &xA,
                                     const MaskType &mask,
                                     const int numBoundaryNodes   ) : 
    _xA( xA ), _mask( mask ),
    _numVertices ( conf.getInitializer().getNumVertices() ), 
    _numGlobalDofs ( conf.getNumGlobalDofs() ),
    _numBoundaryNodes ( numBoundaryNodes ),
    _numInteriorNodes ( _numVertices - _numBoundaryNodes ){ 
        this->getRefsToDisp (  _xA, _refsD1X, _refsD2X );
    }
  

  void evaluateConstraint ( const VectorType & Disp, VectorType &ConstraintVec ) const{

      ConstraintVec.setZero();
      std::vector<Eigen::Ref<const VectorType> > refsD1U, refsD2U; 
      this->getRefsToDisp ( Disp, refsD1U, refsD2U );

      //Part 1
      unsigned int bdrIt = 0;
      for( int i=0; i < _numVertices; ++i ){
          if( _mask[i] ) bdrIt++;        
          else{
            for( int comp = 0; comp <3; ++comp )  ConstraintVec[i-bdrIt] += ( 2. * _refsD1X[comp][i] + refsD1U[comp][i] ) * refsD1U[comp][i];
          }
      }
        
      //Part 2
      bdrIt = 0;
      for( int i=0; i < _numVertices; ++i ){
          if( _mask[i] ) bdrIt++;        
          else{
            for( int comp = 0; comp <3; ++comp ) ConstraintVec[i-bdrIt+ 1*(_numVertices - _numBoundaryNodes)] += ( 2. * _refsD2X[comp][i] + refsD2U[comp][i] ) * refsD2U[comp][i];
          }
      }
       
      //Part 12
      bdrIt = 0;
      for( int i=0; i < _numVertices; ++i ){
          if( _mask[i] ) bdrIt++;        
          else{ 
            for( int comp = 0; comp <3; ++comp )
                ConstraintVec[i-bdrIt + 2*(_numVertices - _numBoundaryNodes )] += 2. * (_refsD1X[comp][i] * refsD2U[comp][i] + _refsD2X[comp][i] * refsD1U[comp][i] + refsD1U[comp][i] * refsD2U[comp][i] );
          }
      }
  }
    
  // Mult^T D Iso
  void applyJacobian( const VectorType & Disp, const VectorType & Multiplier, VectorType & jacConstraintMult ) const {
      
      jacConstraintMult.setZero();
      std::vector<Eigen::Ref<const VectorType> > refsD1U, refsD2U, refsMult; 
      this->getRefsToDisp (  Disp, refsD1U, refsD2U );
      this->getRefsToMult ( Multiplier, refsMult );
       
        unsigned bdrIt = 0;
        for( int i=0; i<_numVertices; ++i ){
            if( _mask[i] ) bdrIt++;
            else{
                for( int comp = 0; comp <3; ++comp ){
                  jacConstraintMult[i +     _numVertices + comp * 3 * _numVertices] += 2. * refsMult[0][i-bdrIt]  * ( _refsD1X[comp][i] + refsD1U[comp][i] );
                  jacConstraintMult[i + 2 * _numVertices + comp * 3 * _numVertices] += 2. * refsMult[1][i-bdrIt]  * ( _refsD2X[comp][i] + refsD2U[comp][i] );
                  jacConstraintMult[i +     _numVertices + comp * 3 * _numVertices] += 2. * refsMult[2][i-bdrIt]  * ( _refsD2X[comp][i] + refsD2U[comp][i] );
                  jacConstraintMult[i + 2 * _numVertices + comp * 3 * _numVertices] += 2. * refsMult[2][i-bdrIt]  * ( _refsD1X[comp][i] + refsD1U[comp][i] );
                }
            }
        }
  }
  
  void assembleAddTripletListHessian( const VectorType & Disp, std::vector<TripletType> & tripletList, const VectorType & Multiplier, const RealType fac = 1. ) const {
      
      std::vector<Eigen::Ref<const VectorType> > refsD1U, refsD2U, refsMult; 
      this->getRefsToDisp ( Disp, refsD1U, refsD2U );
      this->getRefsToMult ( Multiplier, refsMult );
      
      
      // Mult^T D^2 Iso
       unsigned bdrIt = 0;
       for( int i=0; i<_numVertices; ++i){
         if( _mask[i] ) bdrIt++;        
         else{  
           for( int comp=0; comp<3; ++comp ){
               tripletList.push_back( TripletType ( i +   _numVertices + comp*3*_numVertices, i +   _numVertices + comp*3*_numVertices, fac*2.*refsMult[0][i - bdrIt] ) );
               tripletList.push_back( TripletType ( i + 2*_numVertices + comp*3*_numVertices, i + 2*_numVertices + comp*3*_numVertices, fac*2.*refsMult[1][i - bdrIt] ) );
               tripletList.push_back( TripletType ( i +   _numVertices + comp*3*_numVertices, i + 2*_numVertices + comp*3*_numVertices, fac*2.*refsMult[2][i - bdrIt] ) );
               tripletList.push_back( TripletType ( i + 2*_numVertices + comp*3*_numVertices, i +   _numVertices + comp*3*_numVertices, fac*2.*refsMult[2][i - bdrIt] ) );
           }
         }
       }
       
     // D Iso
     bdrIt = 0;
     const int offSet = 3 * _numGlobalDofs;
     for( int i=0; i<_numVertices; ++i){
       if( _mask[i] )  bdrIt++;        
       else{  
        for( int comp = 0; comp <3; ++comp ){
          tripletList.push_back( TripletType ( i - bdrIt +                       offSet, i +   _numVertices + comp*3*_numVertices, fac *2.*( _refsD1X[comp][i] + refsD1U[comp][i] ) ) );
          tripletList.push_back( TripletType ( i - bdrIt +   _numInteriorNodes + offSet, i + 2*_numVertices + comp*3*_numVertices, fac *2.*( _refsD2X[comp][i] + refsD2U[comp][i] ) ) );
          tripletList.push_back( TripletType ( i - bdrIt + 2*_numInteriorNodes + offSet, i +   _numVertices + comp*3*_numVertices, fac *2.*( _refsD2X[comp][i] + refsD2U[comp][i] ) ) );
          tripletList.push_back( TripletType ( i - bdrIt + 2*_numInteriorNodes + offSet, i + 2*_numVertices + comp*3*_numVertices, fac *2.*( _refsD1X[comp][i] + refsD1U[comp][i] ) ) );
          
          //symmetric
          tripletList.push_back( TripletType ( i +   _numVertices + comp*3*_numVertices, i - bdrIt +                       offSet, fac*2.*( _refsD1X[comp][i] + refsD1U[comp][i] ) ) );
          tripletList.push_back( TripletType ( i + 2*_numVertices + comp*3*_numVertices, i - bdrIt +   _numInteriorNodes + offSet, fac*2.*( _refsD2X[comp][i] + refsD2U[comp][i] ) ) );
          tripletList.push_back( TripletType ( i +   _numVertices + comp*3*_numVertices, i - bdrIt + 2*_numInteriorNodes + offSet, fac*2.*( _refsD2X[comp][i] + refsD2U[comp][i] ) ) );
          tripletList.push_back( TripletType ( i + 2*_numVertices + comp*3*_numVertices, i - bdrIt + 2*_numInteriorNodes + offSet, fac*2.*( _refsD1X[comp][i] + refsD1U[comp][i] ) ) );
          
        }
       }
     }
       
  }
  
    int getNumVertices( ) const { return _numVertices; }
    int getNumGlobalDofs( ) const { return _numGlobalDofs; }
    int getNumDofsMult( ) const { return 3 * _numInteriorNodes; }
    const MaskType & getDirichletMask() const { return _mask; }

};





template< typename ConfiguratorType >
class LagrangianIsometryPointwise 
: public pesopt::NonlinearEnergyOp<typename ConfiguratorType::DTContainer> {    
protected :

    typedef typename ConfiguratorType::DTContainer DataTypeContainer;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    typedef typename ConfiguratorType::MaskType MaskType;
    typedef typename ConfiguratorType::TripletType TripletType;
    
    const pesopt::QuadraticEnergyOp<DataTypeContainer> & _ElasticEnergy;
    const SparseMatrixType & _HessianElastic;  // const std::vector<TripletType> & _ElasticHessianTriplet;
    const std::vector<TripletType> &_ElasticHessianTriplet;
    const ShellIsometryHandler_Pointwise<ConfiguratorType> & _isometryOp;
    const MaskType & _mask;
    
    const int _numVertices, _numDofsDisp, _numDofsMultiplier, _numDofsTotal;
    const RealType _facIsoConstr;
    
    mutable RealType _energyElast;
    mutable VectorType _constraintVecEnergy;
public : 
    
    LagrangianIsometryPointwise ( const pesopt::QuadraticEnergyOp<DataTypeContainer> & ElasticEnergy, 
                                  const std::vector<TripletType> &ElasticHessianTriplet,
                                  const ShellIsometryHandler_Pointwise<ConfiguratorType> &isoOp,
                                  const RealType facIsoConstr = 1.0  ) : 
//     _ElasticHessianTriplet ( ElasticEnergy.getTripletList() ), //! \todo getTripletList
    _ElasticEnergy ( ElasticEnergy ),
    _HessianElastic ( ElasticEnergy.getHessian() ),
    _ElasticHessianTriplet ( ElasticHessianTriplet ),
    _isometryOp ( isoOp ), _mask ( isoOp.getDirichletMask() ),
    _numVertices( isoOp.getNumVertices() ), _numDofsDisp( 3 * isoOp.getNumGlobalDofs() ), _numDofsMultiplier( isoOp.getNumDofsMult() ), _numDofsTotal( _numDofsDisp + _numDofsMultiplier ),
    _facIsoConstr ( facIsoConstr ),
    _constraintVecEnergy( _numDofsMultiplier ){}

    
    const int getNumDofs( ) const override { return _numDofsTotal; }
    
    void evaluateEnergy ( const VectorType &DispAndMultiplier, RealType & L ) const override{
        
        const Eigen::Ref<const VectorType> Disp = DispAndMultiplier.segment( 0, _numDofsDisp );
        const Eigen::Ref<const VectorType> Multiplier = DispAndMultiplier.segment( _numDofsDisp, _numDofsMultiplier );
        
        _ElasticEnergy.evaluateEnergy( Disp, _energyElast );
        _isometryOp.evaluateConstraint( Disp,  _constraintVecEnergy );
        
        L = _energyElast + _facIsoConstr * Multiplier.dot( _constraintVecEnergy );
    }
    
    RealType getLastElasticEnergy( ) const { return _energyElast; }
    void getLastConstraintVecEnergy( VectorType & constraintVec ) const { constraintVec = _constraintVecEnergy; }
    
    void evaluateJacobian ( const VectorType &DispAndMultiplier, VectorType & jacL ) const override{    
        const Eigen::Ref<const VectorType> Disp = DispAndMultiplier.segment( 0, _numDofsDisp );
        const Eigen::Ref<const VectorType> Multiplier = DispAndMultiplier.segment( _numDofsDisp, _numDofsMultiplier );
        
        Eigen::Ref<VectorType> jacL_Disp = jacL.segment( 0, _numDofsDisp );
        Eigen::Ref<VectorType> jacL_Mult = jacL.segment( _numDofsDisp, _numDofsMultiplier );
        
        VectorType jacElast ( _numDofsDisp );
        _ElasticEnergy.evaluateJacobian( Disp, jacElast  );
        
        VectorType jacConstraint( _numDofsDisp );
        _isometryOp.applyJacobian( Disp, Multiplier, jacConstraint );
        jacL_Disp = jacElast + _facIsoConstr * jacConstraint;
        
        //bc 
        for( int i = 0; i < 3*_numVertices; ++i ){
            if ( _mask[i] ){
              for( int comp=0; comp<3; ++comp )
                 jacL_Disp[i + comp * 3 * _numVertices] = 0.0;
            }
         } 
        
        VectorType constraintVec ( _numDofsMultiplier );
        _isometryOp.evaluateConstraint( Disp,  constraintVec );
        jacL_Mult = _facIsoConstr * constraintVec;
        
    }
    
    void evaluateHessian ( const VectorType &DispAndMultiplier, SparseMatrixType & hessL ) const override{
    
        const int _numGlobalDofs = _numDofsDisp / 3;
        
        const Eigen::Ref<const VectorType> Disp = DispAndMultiplier.segment( 0, _numDofsDisp );
        const Eigen::Ref<const VectorType> Multiplier = DispAndMultiplier.segment( _numDofsDisp, _numDofsMultiplier );
        
        std::vector<TripletType> tripletListLagrangianHessian ( _ElasticHessianTriplet ); //COPY of EnergyHessian
        tripletListLagrangianHessian.reserve( _ElasticHessianTriplet.size() + (4+3*8)* _numGlobalDofs );
        _isometryOp.assembleAddTripletListHessian( Disp, tripletListLagrangianHessian, Multiplier, _facIsoConstr );

        hessL.setFromTriplets( tripletListLagrangianHessian.begin(), tripletListLagrangianHessian.end() ); 
    }
    
    void evaluateTripletListHessian( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const override{
        throw std::invalid_argument( pesopt::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
    void evaluateTripletListHessianSym( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const override{
        throw std::invalid_argument( pesopt::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
    void print( const VectorType &DispAndMultiplier ) const {
        RealType energyL;
        this->evaluateEnergy( DispAndMultiplier, energyL );
        cout << "energyL = " << energyL << endl;
        cout << "elastic Energy = " <<  this->getLastElasticEnergy( ) << endl;
        VectorType constraintVec ( _numDofsMultiplier );
        this->getLastConstraintVecEnergy( constraintVec );
        cout << "constraintVec.normSqr = " << constraintVec.squaredNorm() << endl;
    }
      
};







//-----------------------------------------------------------------
// used for IPOPT
//-----------------------------------------------------------------
template< typename ConfiguratorType >
class ShellIsometryConstraint : public pesopt::SparseJacobianNonlinearConstraintOps<typename ConfiguratorType::DTContainer>{
  
    typedef typename ConfiguratorType::DTContainer DataTypeContainer;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    typedef typename ConfiguratorType::MaskType MaskType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const VectorType &_xA;
    const RealType _facIsoConstr;
    const MaskType & _mask;
    const int _numVertices, _numGlobalDofs, _numBoundaryNodes, _numInteriorNodes;
    mutable std::vector<Eigen::Ref<const VectorType> > _refsD1X, _refsD2X;
    
protected :
    void getRefsToDisp( const VectorType & x, std::vector<Eigen::Ref<const VectorType> > & refsD1, std::vector<Eigen::Ref<const VectorType> > & refsD2 ) const {
         for( int comp=0; comp<3; ++comp ){
           refsD1.push_back ( x.segment(     _numVertices + comp * _numGlobalDofs, _numVertices) );
           refsD2.push_back ( x.segment( 2 * _numVertices + comp * _numGlobalDofs, _numVertices) );
         }
    }
  
public:
    
    ShellIsometryConstraint ( const ConfiguratorType &conf,
                              const VectorType &xA,
                              const MaskType &mask,
                              const int numBoundaryNodes, 
                              const RealType facIsoConstr = 1.0 ) : 
    pesopt::SparseJacobianNonlinearConstraintOps< typename ConfiguratorType::DTContainer  > ( 3 * ( conf.getInitializer().getNumVertices() - numBoundaryNodes ) ),
    _xA( xA ),
    _facIsoConstr ( facIsoConstr ),
    _mask ( mask ),
    _numVertices ( conf.getInitializer().getNumVertices() ), _numGlobalDofs ( conf.getNumGlobalDofs() ),
    _numBoundaryNodes ( numBoundaryNodes ),
    _numInteriorNodes ( _numVertices - _numBoundaryNodes )
     {
         this->getRefsToDisp (  _xA, _refsD1X, _refsD2X );
     }
  
  const RealType getLowerBound( const int numConstraint ) const override { return 0.; };
  const RealType getUpperBound( const int numConstraint ) const override { return 0.; };
  
  void evaluate ( const VectorType & Disp, VectorType & ConstraintVec ) const override{

      ConstraintVec.setZero();
      std::vector<Eigen::Ref<const VectorType> > refsD1U, refsD2U; 
      this->getRefsToDisp ( Disp, refsD1U, refsD2U );

      //Part 1
      unsigned int bdrIt = 0;
      for( int nodeIndex=0; nodeIndex < _numVertices; ++nodeIndex ){
          if( _mask[nodeIndex] ) bdrIt++;        
          else{
            for( int comp = 0; comp <3; ++comp )  
                ConstraintVec[nodeIndex - bdrIt] += ( 2. * _refsD1X[comp][nodeIndex] + refsD1U[comp][nodeIndex] ) * refsD1U[comp][nodeIndex];
          }
      }
        
      //Part 2
      bdrIt = 0;
      for( int nodeIndex=0; nodeIndex < _numVertices; ++nodeIndex ){
          if( _mask[nodeIndex] ) bdrIt++;        
          else{
            for( int comp = 0; comp <3; ++comp ) 
                ConstraintVec[nodeIndex - bdrIt + 1 * _numInteriorNodes] += ( 2. * _refsD2X[comp][nodeIndex] + refsD2U[comp][nodeIndex] ) * refsD2U[comp][nodeIndex];
          }
      }
       
      //Part 12
      bdrIt = 0;
      for( int nodeIndex=0; nodeIndex < _numVertices; ++nodeIndex ){
          if( _mask[nodeIndex] ) bdrIt++;        
          else{ 
            for( int comp = 0; comp <3; ++comp )
                ConstraintVec[nodeIndex - bdrIt + 2 * _numInteriorNodes] += 2. * (_refsD1X[comp][nodeIndex] * refsD2U[comp][nodeIndex] + _refsD2X[comp][nodeIndex] * refsD1U[comp][nodeIndex] + refsD1U[comp][nodeIndex] * refsD2U[comp][nodeIndex] );
          }
      }
      
      ConstraintVec *= _facIsoConstr;
  }
    

  const int sizeJacobian ( ) const override { return 4 * this->getNumConstraints(); }
    
  void evaluateJacobian( const VectorType & Disp, std::vector<TripletType> & tripletListJac ) const override{

      tripletListJac.reserve( this->sizeJacobian() );
      
      std::vector<Eigen::Ref<const VectorType> > refsD1U, refsD2U, refsMult; 
      this->getRefsToDisp (  Disp, refsD1U, refsD2U );
       
      unsigned bdrIt = 0;
      for( int nodeIndex=0; nodeIndex<_numVertices; ++nodeIndex ){
            if( _mask[nodeIndex] ) bdrIt++;
            else{
                for( int comp = 0; comp <3; ++comp ){
                  //Part1
                  tripletListJac.push_back( 
                    TripletType(   nodeIndex-bdrIt,
                                   nodeIndex + _numVertices + comp * 3 * _numVertices,
                                   _facIsoConstr * 2. * ( _refsD1X[comp][nodeIndex] + refsD1U[comp][nodeIndex] ) ) 
                  );
                  //Part2
                  tripletListJac.push_back( 
                    TripletType(   nodeIndex-bdrIt + 1 * _numInteriorNodes,
                                   nodeIndex + 2 * _numVertices + comp * 3 * _numVertices,
                                   _facIsoConstr * 2. * ( _refsD2X[comp][nodeIndex] + refsD2U[comp][nodeIndex] ) ) 
                  );
                  //Part3
                  tripletListJac.push_back( 
                    TripletType(   nodeIndex-bdrIt + 2 * _numInteriorNodes,
                                   nodeIndex +     _numVertices + comp * 3 * _numVertices,
                                   _facIsoConstr * 2. * ( _refsD2X[comp][nodeIndex] + refsD2U[comp][nodeIndex] ) ) 
                  );
                  tripletListJac.push_back( 
                    TripletType(   nodeIndex-bdrIt + 2 * _numInteriorNodes,
                                   nodeIndex + 2 * _numVertices + comp * 3 * _numVertices,
                                   _facIsoConstr * 2. * ( _refsD1X[comp][nodeIndex] + refsD1U[comp][nodeIndex] )) 
                  );
                }
                
            }
     }
  }
  
  const int sizeHessian ( ) const override { return 3 * this->getNumConstraints(); }
  void evaluateHessian( const VectorType &/*Arg*/, std::vector<std::vector<TripletType>> &hessian ) const override{
       unsigned bdrIt = 0;
       for( int i=0; i<_numVertices; ++i){
         if( _mask[i] ) bdrIt++;        
         else{  
           for( int comp=0; comp<3; ++comp ){
               hessian[i - bdrIt].push_back( TripletType ( i +   _numVertices + comp*3*_numVertices, i +   _numVertices + comp*3*_numVertices, 2. * _facIsoConstr ) );
               hessian[i - bdrIt + 1 * _numInteriorNodes].push_back( TripletType ( i + 2* _numVertices + comp*3*_numVertices, i + 2*_numVertices + comp*3*_numVertices, 2. * _facIsoConstr ) );
               hessian[i - bdrIt + 2 * _numInteriorNodes].push_back( TripletType ( i + 2* _numVertices + comp*3*_numVertices, i +   _numVertices + comp*3*_numVertices, 2. * _facIsoConstr ) );
           }
         }
       }
       
  }
  
    int getNumVertices( ) const { return _numVertices; }
    int getNumGlobalDofs( ) const { return _numGlobalDofs; }
    const MaskType & getDirichletMask() const { return _mask; }

};







//-----------------------------------------------------------------
//-----------------------------------------------------------------
//              Integrated Isometry Constraint
//               \int_chart sqrt det(g_A) |g_A - g_B|^2
//-----------------------------------------------------------------
//-----------------------------------------------------------------
template <typename MatOptConfiguratorType>
class IsometryConstraintIntegrated :
public RefTriangleIntegrator < typename MatOptConfiguratorType::ConfiguratorType, IsometryConstraintIntegrated<MatOptConfiguratorType> >{
  
protected :
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::ElementType ElementType;
  
  const ConfiguratorType & _conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
  mutable DKTFEVectorFunctionEvaluator<ConfiguratorType> *_dispPtr;
  
  public:
    IsometryConstraintIntegrated ( const MatOptConfiguratorType & matOptConf, const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage  ) : 
     RefTriangleIntegrator<ConfiguratorType, IsometryConstraintIntegrated<MatOptConfiguratorType> > ( matOptConf._conf ),
     _conf ( matOptConf._conf ), _xAStorage ( xAStorage ), _dispPtr ( NULL ) {}
     
   IsometryConstraintIntegrated ( const ConfiguratorType & conf, const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage  ) : 
     RefTriangleIntegrator<ConfiguratorType, IsometryConstraintIntegrated<MatOptConfiguratorType> > ( conf ),
     _conf ( conf ), _xAStorage ( xAStorage ), _dispPtr ( NULL ) {}

    ~IsometryConstraintIntegrated() {delete _dispPtr; };
     
  void apply ( const VectorType &Disp, RealType &Dest ) const {
    Dest = 0.0;
    delete _dispPtr;
    _dispPtr = new DKTFEVectorFunctionEvaluator<ConfiguratorType> ( _conf, Disp, 3);
    RefTriangleIntegrator< ConfiguratorType, IsometryConstraintIntegrated<MatOptConfiguratorType> >::assembleAdd( Dest );
  }
  
  void applyOnElements ( const VectorType &Disp, VectorType &Dest ) const {
    delete _dispPtr;
    _dispPtr = new DKTFEVectorFunctionEvaluator<ConfiguratorType> ( _conf, Disp, 3);
    RefTriangleIntegrator< ConfiguratorType, IsometryConstraintIntegrated<MatOptConfiguratorType> >::assembleOnElements( Dest );
  } 

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
    Matrix32 DxA, Du;
//     DxA = _xAStorage.getGradientInTangentSpace(El.getGlobalElementIdx(),QuadPoint);
    DxA = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint);
//     _dispPtr->evaluateGradientInTangentSpaceAtQuadPoint( El, QuadPoint, Du );
    _dispPtr->evaluateGradientAtQuadPoint( El, QuadPoint, Du );
    Matrix22 gAMinusgB = Du.transpose() * DxA + DxA.transpose() * Du + Du.transpose() * Du;
     return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * pesopt::ddProd<RealType,Matrix22>(gAMinusgB,gAMinusgB);
  }
};




//-----------------------------------------------------------------
//-----------------------------------------------------------------
//              Integrated Isometry Constraint
//               \int_chart sqrt det(g_A) |g_A - g_B|
//-----------------------------------------------------------------
//-----------------------------------------------------------------
template <typename MatOptConfiguratorType>
class IsometryConstraintL1 :
public RefTriangleIntegrator < typename MatOptConfiguratorType::ConfiguratorType, IsometryConstraintL1<MatOptConfiguratorType> >{
  
protected :
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::ElementType ElementType;
  
  const ConfiguratorType & _conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
  mutable DKTFEVectorFunctionEvaluator<ConfiguratorType> *_dispPtr;
  
  public:
    IsometryConstraintL1 ( const MatOptConfiguratorType & matOptConf, const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage  ) : 
     RefTriangleIntegrator<ConfiguratorType, IsometryConstraintL1<MatOptConfiguratorType> > ( matOptConf._conf ),
     _conf ( matOptConf._conf ), _xAStorage ( xAStorage ), _dispPtr ( NULL ) {}
     
   IsometryConstraintL1 ( const ConfiguratorType & conf, const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage  ) : 
     RefTriangleIntegrator<ConfiguratorType, IsometryConstraintL1<MatOptConfiguratorType> > ( conf ),
     _conf ( conf ), _xAStorage ( xAStorage ), _dispPtr ( NULL ) {}

    ~IsometryConstraintL1() {delete _dispPtr; };
     
  void apply ( const VectorType &Disp, RealType &Dest ) const {
    Dest = 0.0;
    delete _dispPtr;
    _dispPtr = new DKTFEVectorFunctionEvaluator<ConfiguratorType> ( _conf, Disp, 3);
    RefTriangleIntegrator< ConfiguratorType, IsometryConstraintL1<MatOptConfiguratorType> >::assembleAdd( Dest );
  }
  
  void applyOnElements ( const VectorType &Disp, VectorType &Dest ) const {
    delete _dispPtr;
    _dispPtr = new DKTFEVectorFunctionEvaluator<ConfiguratorType> ( _conf, Disp, 3);
    RefTriangleIntegrator< ConfiguratorType, IsometryConstraintL1<MatOptConfiguratorType> >::assembleOnElements( Dest );
  } 

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
    Matrix32 DxA, Du;
    //     DxA = _xAStorage.getGradientInTangentSpace(El.getGlobalElementIdx(),QuadPoint);
    DxA = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint);
    //     _dispPtr->evaluateGradientInTangentSpaceAtQuadPoint( El, QuadPoint, Du );
    _dispPtr->evaluateGradientAtQuadPoint( El, QuadPoint, Du );
    Matrix22 gAMinusgB = Du.transpose() * DxA + DxA.transpose() * Du + Du.transpose() * Du;
     return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * sqrt( pesopt::ddProd<RealType,Matrix22>(gAMinusgB,gAMinusgB) );
  }
};




//-----------------------------------------------------------------
//-----------------------------------------------------------------
//              Integrated Isometry Constraint
//               \int_chart sqrt det(g_A) (tr( g_A^-1 g_B ) - 2.)
//-----------------------------------------------------------------
//-----------------------------------------------------------------
template <typename MatOptConfiguratorType>
class IsometryConstraint_TraceIntegrated :
public RefTriangleIntegrator < typename MatOptConfiguratorType::ConfiguratorType, IsometryConstraint_TraceIntegrated<MatOptConfiguratorType> >{
  
protected :
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::ElementType ElementType;
  
  const ConfiguratorType & _conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
  mutable DKTFEVectorFunctionEvaluator<ConfiguratorType> *_dispPtr;
  
  public:
    IsometryConstraint_TraceIntegrated ( const MatOptConfiguratorType & matOptConf, const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage  ) : 
     RefTriangleIntegrator<ConfiguratorType, IsometryConstraint_TraceIntegrated<MatOptConfiguratorType> > ( matOptConf._conf ),
     _conf ( matOptConf._conf ), _xAStorage ( xAStorage ), _dispPtr ( NULL ) {}
     
   IsometryConstraint_TraceIntegrated ( const ConfiguratorType & conf, const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage  ) : 
     RefTriangleIntegrator<ConfiguratorType, IsometryConstraint_TraceIntegrated<MatOptConfiguratorType> > ( conf ),
     _conf ( conf ), _xAStorage ( xAStorage ), _dispPtr ( NULL ) {}

    ~IsometryConstraint_TraceIntegrated() {delete _dispPtr; };
     
  void apply ( const VectorType &Disp, RealType &Dest ) const {
    Dest = 0.0;
    delete _dispPtr;
    _dispPtr = new DKTFEVectorFunctionEvaluator<ConfiguratorType> ( _conf, Disp, 3);
    RefTriangleIntegrator< ConfiguratorType, IsometryConstraint_TraceIntegrated<MatOptConfiguratorType> >::assembleAdd( Dest );
  }
  
  void applyOnElements ( const VectorType &Disp, VectorType &Dest ) const {
    delete _dispPtr;
    _dispPtr = new DKTFEVectorFunctionEvaluator<ConfiguratorType> ( _conf, Disp, 3);
    RefTriangleIntegrator< ConfiguratorType, IsometryConstraint_TraceIntegrated<MatOptConfiguratorType> >::assembleOnElements( Dest );
  } 

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
    Matrix32 DxA, Du;
//     DxA = _xAStorage.getGradientInTangentSpace(El.getGlobalElementIdx(),QuadPoint);
    DxA = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint);
//     _dispPtr->evaluateGradientInTangentSpaceAtQuadPoint( El, QuadPoint, Du );
    _dispPtr->evaluateGradientAtQuadPoint( El, QuadPoint, Du );
    Matrix22 gA = DxA.transpose() * DxA;
    Matrix22 gAinv; gAinv = gA.inverse();
    Matrix22 gB = Du.transpose() * DxA + DxA.transpose() * Du + Du.transpose() * Du + gA;
    Matrix22 gAinvGB = gAinv * gB;
     return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * ( gAinvGB.trace() - 2. );
  }
};





//-----------------------------------------------------------------
//-----------------------------------------------------------------
//              Integrated Isometry Constraint
//               \int_chart sqrt det(g_A) | g_A^-1 g_B - id |^2
//-----------------------------------------------------------------
//-----------------------------------------------------------------
template <typename MatOptConfiguratorType>
class IsometryConstraint_L2 :
public RefTriangleIntegrator < typename MatOptConfiguratorType::ConfiguratorType, IsometryConstraint_L2<MatOptConfiguratorType> >{
  
protected :
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::Matrix22 Matrix22;
  typedef typename ConfiguratorType::Matrix32 Matrix32;
  typedef typename ConfiguratorType::Matrix33 Matrix33;
  typedef typename ConfiguratorType::ElementType ElementType;
  
  const ConfiguratorType & _conf;
  const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &_xAStorage;
  mutable DKTFEVectorFunctionEvaluator<ConfiguratorType> *_dispPtr;
  
  public:
    IsometryConstraint_L2 ( const MatOptConfiguratorType & matOptConf, const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage  ) : 
     RefTriangleIntegrator<ConfiguratorType, IsometryConstraint_L2<MatOptConfiguratorType> > ( matOptConf._conf ),
     _conf ( matOptConf._conf ), _xAStorage ( xAStorage ), _dispPtr ( NULL ) {}
     
   IsometryConstraint_L2 ( const ConfiguratorType & conf, const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage  ) : 
     RefTriangleIntegrator<ConfiguratorType, IsometryConstraint_L2<MatOptConfiguratorType> > ( conf ),
     _conf ( conf ), _xAStorage ( xAStorage ), _dispPtr ( NULL ) {}

    ~IsometryConstraint_L2() {delete _dispPtr; };
     
  void apply ( const VectorType &Disp, RealType &Dest ) const {
    Dest = 0.0;
    delete _dispPtr;
    _dispPtr = new DKTFEVectorFunctionEvaluator<ConfiguratorType> ( _conf, Disp, 3);
    RefTriangleIntegrator< ConfiguratorType, IsometryConstraint_L2<MatOptConfiguratorType> >::assembleAdd( Dest );
  }
  
  void applyOnElements ( const VectorType &Disp, VectorType &Dest ) const {
    delete _dispPtr;
    _dispPtr = new DKTFEVectorFunctionEvaluator<ConfiguratorType> ( _conf, Disp, 3);
    RefTriangleIntegrator< ConfiguratorType, IsometryConstraint_L2<MatOptConfiguratorType> >::assembleOnElements( Dest );
  } 

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
    Matrix32 DxA, Du;
//     DxA = _xAStorage.getGradientInTangentSpace(El.getGlobalElementIdx(),QuadPoint);
    DxA = _xAStorage.getGradient(El.getGlobalElementIdx(),QuadPoint);
//     _dispPtr->evaluateGradientInTangentSpaceAtQuadPoint( El, QuadPoint, Du );
    _dispPtr->evaluateGradientAtQuadPoint( El, QuadPoint, Du );
    Matrix22 gA = DxA.transpose() * DxA;
    Matrix22 gAinv; gAinv = gA.inverse();
    Matrix22 gBminGA = Du.transpose() * DxA + DxA.transpose() * Du + Du.transpose() * Du;
    Matrix22 gAinvGBminGA = gAinv * gBminGA;
     return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * pesopt::ddProd<RealType,Matrix22>( gAinvGBminGA, gAinvGBminGA );
  }
};

#endif
