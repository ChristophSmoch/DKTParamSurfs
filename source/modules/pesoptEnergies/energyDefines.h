#ifndef __ENERGYDEFINES_H
#define __ENERGYDEFINES_H

namespace pesopt {


template <typename DataTypeContainer >
class NonlinearEnergyOp {

protected:
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  typedef pesopt::BoostParser ParameterParserType;
  
public:
  
  NonlinearEnergyOp( ) { }

  // Destroy polymorphic Ops correctly, important!
  virtual ~NonlinearEnergyOp () {}
  
  virtual void evaluateEnergy( const VectorType &Arg, RealType & energy ) const = 0;
  RealType evaluateEnergy ( VectorType &Arg ) { RealType energy = 0; evaluateEnergy( Arg, energy ); return energy; }
  virtual void evaluateJacobian( const VectorType &Arg, VectorType & jacobian ) const = 0;
  virtual void evaluateHessian( const VectorType &Arg, SparseMatrixType & hessian ) const = 0;
  virtual void evaluateTripletListHessian( const VectorType &Arg, std::vector<typename DataTypeContainer::TripletType> & tripletListHessian ) const = 0;
  virtual void evaluateTripletListHessianSym( const VectorType &Arg, std::vector<typename DataTypeContainer::TripletType> & tripletListHessian ) const = 0;
  
  virtual const int getNumDofs ( ) const = 0;
  
  virtual void getEnergyInfo ( const VectorType &Arg, ParameterParserType &energyInfo) const {};
  
};


// Energy is quadratic, i.e of type 
// E(u) = 0.5 Mu \cdot u - F \cdot u,
// for a constant matrix M and a constant vector F
// Thus, DE(u) = M u - F and D^2E(u) = M
template <typename DataTypeContainer >
class QuadraticEnergyOp {
protected :
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  
  const SparseMatrixType &_Hessian;
  const VectorType &_rhs;
  
public:
  QuadraticEnergyOp( const SparseMatrixType & Hessian, const VectorType & rhs ) : _Hessian( Hessian ), _rhs( rhs ){ }

  RealType evaluateEnergy ( const VectorType &Arg ) const { return 0.5 * (_Hessian * Arg ).dot( Arg ) - _rhs.dot( Arg );}
  void evaluateEnergy( const VectorType &Arg, RealType & energy ) const { energy = this->evaluateEnergy(Arg);}
  void evaluateJacobian( const VectorType &Arg, VectorType & jacobian ) const { jacobian = _Hessian * Arg - _rhs;}
  void evaluateHessian( const VectorType &/*Arg*/, SparseMatrixType & hessian ) const { hessian = _Hessian;}
  const SparseMatrixType & getHessian( ) const {return _Hessian;}
};



template <typename DataTypeContainer >
class NonlinearConstraintOps {

protected:
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  typedef pesopt::BoostParser ParameterParserType;
  
  mutable int _numConstraints;
  
public:
  
  NonlinearConstraintOps( const int numConstraints) : _numConstraints ( numConstraints ) { }

  // Destroy polymorphic Ops correctly, important!
  virtual ~NonlinearConstraintOps () {}
  
  virtual void evaluateEnergy( const int numConstraint, const VectorType &Arg, RealType & energy ) const = 0;
  RealType evaluateEnergy ( const int numConstraint, VectorType &Arg ) { RealType energy = 0; evaluateEnergy( numConstraint, Arg, energy ); return energy; }
  virtual void evaluateJacobian( const int numConstraint, const VectorType &Arg, VectorType & jacobian ) const = 0;
  virtual void evaluateHessian( const int numConstraint, const VectorType &Arg, SparseMatrixType & hessian ) const = 0;
  
  virtual const RealType getLowerBound( const int numConstraint ) const = 0;
  virtual const RealType getUpperBound( const int numConstraint ) const = 0;
  
  const int getNumConstraints( ) const { return _numConstraints;}
  
  virtual void getConstraintInfo ( const VectorType &Arg, ParameterParserType & energyInfo) const {};
  
};




//! version with pointer
template <typename DataTypeContainer >
class CombinedNonlinearConstraintOps :
public NonlinearConstraintOps<DataTypeContainer>{

protected:
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  typedef pesopt::BoostParser ParameterParserType;

  //! version I - pointer
  mutable std::vector<const pesopt::NonlinearConstraintOps<DataTypeContainer>*> _nonlinearConstraintOps;
//   //! version II - shared_pointer
//   mutable std::vector< std::shared_ptr<pesopt::NonlinearConstraintOps<DataTypeContainer> > > _nonlinearConstraintOps;
  
  // map: numconstraint -> (numOp,numConstraintInOp)
  mutable std::vector<std::pair<int,int>> _numOpAndNumConstraintVec;
  
public:
  
  CombinedNonlinearConstraintOps( ) : 
  NonlinearConstraintOps<DataTypeContainer>( 0 ) { }

 
  
  
  // Destroy polymorphic Ops correctly, important!
  virtual ~CombinedNonlinearConstraintOps () {}
  
  void push_back( const NonlinearConstraintOps<DataTypeContainer>& constraintOps ) const {
      
      const int oldNumOps = _nonlinearConstraintOps.size();
      const int newNumConstraints = constraintOps.getNumConstraints();
      this->_numConstraints += newNumConstraints;
     
      //! version I - pointer
     _nonlinearConstraintOps.push_back( &constraintOps );
     //       //! version II - shared pointer
     //       _nonlinearConstraintOps.push_back( std::make_shared<NonlinearConstraintOps<DataTypeContainer>> (&constraintOps) );
      
      
      for( int j=0; j< newNumConstraints; ++j )
          _numOpAndNumConstraintVec.push_back(std::make_pair(oldNumOps,j));
  }
 
  void getNumOpAndNumConstraint( const int numConstraint, int& numOp, int& numConstraintOp ) const {
      numOp = _numOpAndNumConstraintVec[numConstraint].first;
      numConstraintOp = _numOpAndNumConstraintVec[numConstraint].second;
  }
  
  void evaluateEnergy( const int numConstraint, const VectorType &Arg, RealType & energy ) const override{
      int numOp, numConstraintOp;
      this->getNumOpAndNumConstraint( numConstraint, numOp, numConstraintOp );
      _nonlinearConstraintOps[numOp]->evaluateEnergy(numConstraintOp, Arg, energy );
  }

  void evaluateJacobian( const int numConstraint, const VectorType &Arg, VectorType & jacobian ) const override{
      int numOp, numConstraintOp;
      this->getNumOpAndNumConstraint( numConstraint, numOp, numConstraintOp );
      _nonlinearConstraintOps[numOp]->evaluateJacobian(numConstraintOp, Arg, jacobian );
  }
  
  void evaluateHessian( const int numConstraint, const VectorType &Arg, SparseMatrixType & hessian ) const override{
      int numOp, numConstraintOp;
      this->getNumOpAndNumConstraint( numConstraint, numOp, numConstraintOp );
      _nonlinearConstraintOps[numOp]->evaluateHessian(numConstraintOp, Arg, hessian );
  }
  
  const RealType getLowerBound( const int numConstraint ) const {
      int numOp, numConstraintOp;
      this->getNumOpAndNumConstraint( numConstraint, numOp, numConstraintOp );
      return _nonlinearConstraintOps[numOp]->getLowerBound(numConstraintOp );
  }
  const RealType getUpperBound( const int numConstraint ) const {
      int numOp, numConstraintOp;
      this->getNumOpAndNumConstraint( numConstraint, numOp, numConstraintOp );
      return _nonlinearConstraintOps[numOp]->getUpperBound(numConstraintOp );
  }
  
  void getConstraintInfo ( const VectorType &Arg, ParameterParserType & energyInfo) const override {
      for( int i=0; i<_nonlinearConstraintOps.size(); ++i )
          _nonlinearConstraintOps[i]->getConstraintInfo( Arg, energyInfo );
  };
  
};


template <typename DataTypeContainer >
class SparseJacobianNonlinearConstraintOps {

protected:
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  typedef typename DataTypeContainer::TripletType TripletType;
  typedef pesopt::BoostParser ParameterParserType;
  
  const int _numConstraints;
  
public:
  
  SparseJacobianNonlinearConstraintOps( const int numConstraints) : _numConstraints ( numConstraints ) { }

  // Destroy polymorphic Ops correctly, important!
  virtual ~SparseJacobianNonlinearConstraintOps () {}
  
  virtual void evaluate( const VectorType &Arg, VectorType & constraintVec ) const = 0;
  virtual void evaluateJacobian( const VectorType &Arg, std::vector<TripletType> & jacobian ) const = 0;
  virtual void evaluateHessian( const VectorType &Arg, std::vector<std::vector<TripletType>> &hessian ) const = 0;
  
  virtual const RealType getLowerBound( const int numConstraint ) const = 0;
  virtual const RealType getUpperBound( const int numConstraint ) const = 0;
  
  virtual void getConstraintInfo ( const VectorType &Arg, ParameterParserType & energyInfo) const {};
  
  const int getNumConstraints( ) const { return _numConstraints;}
  virtual const int sizeJacobian( ) const = 0;
  virtual const int sizeHessian( ) const = 0;
};


// template <typename DataTypeContainer >
// class MultipleNonlinearConstraintOps 
// : public NonlinearConstraintOps<DataTypeContainer> {
// 
// protected:
//   typedef DataTypeContainer DTContainer;
//   typedef typename DataTypeContainer::RealType RealType;
//   typedef typename DataTypeContainer::VectorType VectorType;
//   typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
//   
//   const std::vector< std::reference_wrapper<pesopt::NonlinearConstraintOps<DataTypeContainer>> > &_constraintOps;
//   
// public:
//   
//   MultipleNonlinearConstraintOpsr( const std::vector< std::reference_wrapper<pesopt::NonlinearConstraintOps<DataTypeContainer>> > &constraintOps) : 
//   NonlinearConstraintOps<DataTypeContainer> ( constraintOps.size() ),
//   _constraintOps ( constraintOps ) { }
// 
//   // Destroy polymorphic Ops correctly, important!
//   ~NonlinearConstraintOps_FromVector () {}
//   
//   void evaluateEnergy( const int numConstraint, const VectorType &Arg, RealType & energy ) const {
//       _constraintOps[numConstraint].get().evaluateEnergy( Arg, energy );
//   };
//   void evaluateJacobian( const int numConstraint, const VectorType &Arg, VectorType & jacobian ) const {
//       _constraintOps[numConstraint].get().evaluateJacobian( Arg, jacobian );
//   };
//   void evaluateHessian( const int numConstraint, const VectorType &Arg, SparseMatrixType & hessian ) const {
//       _constraintOps[numConstraint].get().evaluateHessian( Arg, hessian );
//   };
//   
// };



// template <typename DataTypeContainer >
// class NonlinearConstraintOps_FromVector : public NonlinearConstraintOps<DataTypeContainer> {
// 
// protected:
//   typedef DataTypeContainer DTContainer;
//   typedef typename DataTypeContainer::RealType RealType;
//   typedef typename DataTypeContainer::VectorType VectorType;
//   typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
//   
//   const std::vector< std::reference_wrapper<pesopt::NonlinearEnergyOp<DataTypeContainer>> > &_constraintOps;
//   
// public:
//   
//   NonlinearConstraintOps_FromVector( const std::vector< std::reference_wrapper<pesopt::NonlinearEnergyOp<DataTypeContainer>> > &constraintOps) : 
//   NonlinearConstraintOps<DataTypeContainer> ( constraintOps.size() ),
//   _constraintOps ( constraintOps ) { }
// 
//   // Destroy polymorphic Ops correctly, important!
//   ~NonlinearConstraintOps_FromVector () {}
//   
//   void evaluateEnergy( const int numConstraint, const VectorType &Arg, RealType & energy ) const {
//       _constraintOps[numConstraint].get().evaluateEnergy( Arg, energy );
//   };
//   void evaluateJacobian( const int numConstraint, const VectorType &Arg, VectorType & jacobian ) const {
//       _constraintOps[numConstraint].get().evaluateJacobian( Arg, jacobian );
//   };
//   void evaluateHessian( const int numConstraint, const VectorType &Arg, SparseMatrixType & hessian ) const {
//       _constraintOps[numConstraint].get().evaluateHessian( Arg, hessian );
//   };
//   
// };








}//end namespace


#endif
