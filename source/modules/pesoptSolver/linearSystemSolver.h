#ifndef __LINEARSYSTEMSOLVER_H
#define __LINEARSYSTEMSOLVER_H

#include <pesopt_IO.h>
#include <pesopt_Eigen.h>
// #include <unsupported/Eigen/IterativeSolvers>
// #include <EigenBiCGSTAB_New.h>
// #include <EigenGMRES_New.h>


#define _DEBUGLINEARSOLVER

#ifdef _DEBUGLINEARSOLVER
 int debugCounterIterativeSolver = 0;
#endif

 
template<typename DataTypeContainer>
class LinearSystemSolver{
    
protected:
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  
public:
  
  LinearSystemSolver( ) { }

  // Destroy polymorphic Ops correctly, important!
  virtual ~LinearSystemSolver () {}
  
  virtual void solve( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs ) const = 0;
    
    
};
 
 
enum IterativeLinearSolverMethod {
    EigenConjugateGradient = 11,
    EigenConjugateGradient_NEW = 12,
//     PCGDiagonalPre = 13, 
//     EigenBiCGSTAB = 21,
//     EigenBiCGSTAB_WithGuess = 22,
//     EigenBiCGSTAB_ILU = 23,
//     EigenBiCGSTAB_NEW = 24,
//     EigenGMRES = 31,
//     EigenGMRES_WithGuess = 32,
//     EigenGMRES_ILU = 33,
//     EigenGMRES_NEW = 34
// //     EigenDGMRES = 15,
// //     EigenDGMRES_ILU = 16
// //     EigenSparseLU = 21
};
 
 
// template<typename DataTypeContainer, typename IterativeLinearSolverType>
// class IterativeLinearSystemSolver :
// public LinearSystemSolver<DataTypeContainer> {
//     
// protected:
//     typedef typename DataTypeContainer::RealType RealType;
//     typedef typename DataTypeContainer::VectorType  VectorType;
//     typedef typename DataTypeContainer::SparseMatrixType  SparseMatrixType;
//     
//     mutable IterativeLinearSolverType _linearSolver;
//     
//     mutable RealType _tolerance;
//     mutable RealType _maxItersFac;
//     const string _errorMessage;
//     const string _saveDirectory;
//     const bool _debug;
//     
// public:
//      
//   IterativeLinearSystemSolver( const RealType tolerance = Eigen::NumTraits<RealType>::epsilon(),
//                                const RealType maxItersFac = 1.,
//                                const string &errorMessage = "", 
//                                const string &saveDirectory = "", 
//                                const int outputLevel = 0, 
//                                const bool debug = true  ) :
//   _tolerance ( tolerance ), _maxItersFac ( maxItersFac ),
//   _errorMessage( errorMessage ), _saveDirectory ( saveDirectory ), _debug ( debug ) { }   
//   
//   void setTolerance ( const RealType tolerance ) const { _tolerance = tolerance; }
//   void setMaxItersFac ( const RealType maxItersFac ) const { _maxItersFac = maxItersFac; }
//   
//   //
//   bool prepareSolver( const SparseMatrixType& systemMatrix ) const{
//       _linearSolver.compute( systemMatrix );
//       if(_linearSolver.info()!=Eigen::Success) {
//           cout << "prepare linear solver failed" << endl;
//           return false;
//       }
//       return true;
//   }
//   
//   //
//   bool solveForPreparedMatrix( VectorType& solution, const VectorType& rhs  ) const{
//     solution = _linearSolver.solve( rhs );
//     if(_linearSolver.info()!=Eigen::Success ) {
//         cerr <<  endl
//              << "=================================================" << endl 
//              << "error comes from solve in " << _errorMessage << endl
//              << "#iterations:     " << _linearSolver.iterations() << endl
//              << "estimated error: " << std::setprecision( 32 ) << _linearSolver.error() << endl
//              << "but tolerance is: " <<  _tolerance << endl
//              << "=================================================" << endl;
//         return false;  
//     }
//             
//     return true;
//   }
//   
//   
//   bool solveForPreparedMatrixWithGuess( VectorType& solution, const VectorType& rhs, const VectorType& guess ) const{
//       
//     solution = _linearSolver.solveWithGuess( rhs, guess );
//     if(_linearSolver.info()!=Eigen::Success ) {
//         cerr <<  endl
//              << "=================================================" << endl 
//              << "error comes from solve in " << _errorMessage << endl
//              << "#iterations:     " << _linearSolver.iterations() << endl
//              << "estimated error: " << std::setprecision( 32 ) << _linearSolver.error() << endl
//              << "but tolerance is: " <<  _tolerance << endl
//              << "=================================================" << endl;
//         return false;  
//     }
//     return true;
//   }
//   
//   
//   // solve Ax = b for x
//   void solve( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs ) const override{
//       bool succesPrepareSolver = prepareSolver( systemMatrix );
//       const RealType newTolerance = sqrt( _tolerance * _tolerance * systemMatrix.cols() );
//       _linearSolver.setTolerance( newTolerance );
//       const int maxIters = static_cast<int> ( _maxItersFac * static_cast<RealType> ( systemMatrix.cols() ) );
//       _linearSolver.setMaxIterations(maxIters);
//       bool success = solveForPreparedMatrix( solution, rhs );
// #ifdef _DEBUGLINEARSOLVER
//       if( !success && _debug ) {
//         cout << endl 
//              << "======================================================" << endl
//              << "linear solver failed" << endl;
//         cout << "norm of rhs = " << rhs.norm() << endl;
//         if(!(_saveDirectory.empty())){
//             pesopt::printVectorToFile<VectorType> ( rhs, pesopt::strprintf( "%s/rhs_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver ).c_str(), 32 ); 
//             pesopt::printSparseMatrixToFile<SparseMatrixType> ( systemMatrix, pesopt::strprintf( "%s/systemmatrix_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver ).c_str(), 32 ); 
//             std::ofstream outSolver ( pesopt::strprintf ( "%s/iterativeSolverOutput_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver  ).c_str ()  );
//             _linearSolver.setOutputStream( outSolver );
//         }
//         solveForPreparedMatrix( solution, rhs );
//         debugCounterIterativeSolver++;
//         throw invalid_argument ( pesopt::strprintf ( " solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
//       }
// #endif
//   }
//   
//   void solveWithGuess( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs,
//                        const VectorType &guess  ) const{
//       bool succesPrepareSolver = prepareSolver( systemMatrix );
//       const RealType newTolerance = sqrt( _tolerance * _tolerance * systemMatrix.cols() );
//       _linearSolver.setTolerance( newTolerance );
//       const int maxIters = static_cast<int> ( _maxItersFac * static_cast<RealType> ( systemMatrix.cols() ) );
//       _linearSolver.setMaxIterations(maxIters);
//       bool success = solveForPreparedMatrixWithGuess( solution, rhs, guess );
// #ifdef _DEBUGLINEARSOLVER
//       if( !success && _debug ) {
//         cout << endl 
//              << "======================================================" << endl
//              << "linear solver failed" << endl;
//         cout << "norm of rhs = " << rhs.norm() << endl;
//         if(!(_saveDirectory.empty())){
//             pesopt::printVectorToFile<VectorType> ( rhs, pesopt::strprintf( "%s/rhs_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver ).c_str(), 32 ); 
//             pesopt::printSparseMatrixToFile<SparseMatrixType> ( systemMatrix, pesopt::strprintf( "%s/systemmatrix_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver ).c_str(), 32 ); 
//             std::ofstream outSolver ( pesopt::strprintf ( "%s/iterativeSolverOutput_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver  ).c_str()  );
//             _linearSolver.setOutputStream( outSolver );
//         }
//         solveForPreparedMatrixWithGuess( solution, rhs, guess );
//         debugCounterIterativeSolver++;
//         throw invalid_argument ( pesopt::strprintf ( " solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
//       }
// #endif
//   }
//    
// };


template<typename DataTypeContainer, IterativeLinearSolverMethod iterativeLSMethod>
class IterativeLinearSystemSolver {};


template<typename DataTypeContainer>
class IterativeLinearSystemSolver<DataTypeContainer, EigenConjugateGradient> :
public LinearSystemSolver<DataTypeContainer> {
    
protected:
    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType  VectorType;
    typedef typename DataTypeContainer::SparseMatrixType  SparseMatrixType;
    
    mutable Eigen::ConjugateGradient<SparseMatrixType,  Eigen::Lower|Eigen::Upper> _linearSolver;
    
    mutable RealType _tolerance;
    mutable RealType _maxItersFac;
    const string _errorMessage;
    const string _saveDirectory;
    const bool _debug;
    
public:
     
  IterativeLinearSystemSolver( const RealType tolerance = Eigen::NumTraits<RealType>::epsilon(),
                               const RealType maxItersFac = 1.,
                               const string &errorMessage = "", 
                               const string &saveDirectory = "", 
                               const int outputLevel = 0, 
                               const bool debug = true  ) :
  _tolerance ( tolerance ), _maxItersFac ( maxItersFac ),
  _errorMessage( errorMessage ), _saveDirectory ( saveDirectory ), _debug ( debug ) { }   
  
  void setTolerance ( const RealType tolerance ) const { _tolerance = tolerance; }
  void setMaxItersFac ( const RealType maxItersFac ) const { _maxItersFac = maxItersFac; }
  
  //
  bool prepareSolver( const SparseMatrixType& systemMatrix ) const{
      _linearSolver.compute( systemMatrix );
      if(_linearSolver.info()!=Eigen::Success) {
          cout << "prepare linear solver failed" << endl;
          return false;
      }
      return true;
  }
  
  //
  bool solveForPreparedMatrix( VectorType& solution, const VectorType& rhs  ) const{
    solution = _linearSolver.solve( rhs );
    if(_linearSolver.info()!=Eigen::Success ) {
        cerr <<  endl
             << "=================================================" << endl 
             << "error comes from solve in " << _errorMessage << endl
             << "#iterations:     " << _linearSolver.iterations() << endl
             << "estimated error: " << std::setprecision( 32 ) << _linearSolver.error() << endl
             << "but tolerance is: " <<  _tolerance << endl
             << "=================================================" << endl;
        return false;  
    }
            
    return true;
  }
  
  
  bool solveForPreparedMatrixWithGuess( VectorType& solution, const VectorType& rhs, const VectorType& guess ) const{
      
    solution = _linearSolver.solveWithGuess( rhs, guess );
    if(_linearSolver.info()!=Eigen::Success ) {
        cerr <<  endl
             << "=================================================" << endl 
             << "error comes from solve in " << _errorMessage << endl
             << "#iterations:     " << _linearSolver.iterations() << endl
             << "estimated error: " << std::setprecision( 32 ) << _linearSolver.error() << endl
             << "but tolerance is: " <<  _tolerance << endl
             << "=================================================" << endl;
        return false;  
    }
    return true;
  }
  
  
  // solve Ax = b for x
  void solve( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs ) const override{
      bool succesPrepareSolver = prepareSolver( systemMatrix );
      const RealType newTolerance = sqrt( _tolerance * _tolerance * systemMatrix.cols() );
      _linearSolver.setTolerance( newTolerance );
      const int maxIters = static_cast<int> ( _maxItersFac * static_cast<RealType> ( systemMatrix.cols() ) );
      _linearSolver.setMaxIterations(maxIters);
      bool success = solveForPreparedMatrix( solution, rhs );
#ifdef _DEBUGLINEARSOLVER
      if( !success && _debug ) {
        cout << endl 
             << "======================================================" << endl
             << "linear solver failed" << endl;
        cout << "norm of rhs = " << rhs.norm() << endl;
        if(!(_saveDirectory.empty())){
            pesopt::printVectorToFile<VectorType> ( rhs, pesopt::strprintf( "%s/rhs_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver ).c_str(), 32 ); 
            pesopt::printSparseMatrixToFile<SparseMatrixType> ( systemMatrix, pesopt::strprintf( "%s/systemmatrix_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver ).c_str(), 32 ); 
            std::ofstream outSolver ( pesopt::strprintf ( "%s/iterativeSolverOutput_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver  ).c_str ()  );
//             _linearSolver.setOutputStream( outSolver );
        }
        solveForPreparedMatrix( solution, rhs );
        debugCounterIterativeSolver++;
        throw invalid_argument ( pesopt::strprintf ( " solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      }
#endif
  }
  
  void solveWithGuess( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs,
                       const VectorType &guess  ) const{
      bool succesPrepareSolver = prepareSolver( systemMatrix );
      const RealType newTolerance = sqrt( _tolerance * _tolerance * systemMatrix.cols() );
      _linearSolver.setTolerance( newTolerance );
      const int maxIters = static_cast<int> ( _maxItersFac * static_cast<RealType> ( systemMatrix.cols() ) );
      _linearSolver.setMaxIterations(maxIters);
      bool success = solveForPreparedMatrixWithGuess( solution, rhs, guess );
#ifdef _DEBUGLINEARSOLVER
      if( !success && _debug ) {
        cout << endl 
             << "======================================================" << endl
             << "linear solver failed" << endl;
        cout << "norm of rhs = " << rhs.norm() << endl;
        if(!(_saveDirectory.empty())){
            pesopt::printVectorToFile<VectorType> ( rhs, pesopt::strprintf( "%s/rhs_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver ).c_str(), 32 ); 
            pesopt::printSparseMatrixToFile<SparseMatrixType> ( systemMatrix, pesopt::strprintf( "%s/systemmatrix_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver ).c_str(), 32 ); 
            std::ofstream outSolver ( pesopt::strprintf ( "%s/iterativeSolverOutput_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver  ).c_str()  );
            _linearSolver.setOutputStream( outSolver );
        }
        solveForPreparedMatrixWithGuess( solution, rhs, guess );
        debugCounterIterativeSolver++;
        throw invalid_argument ( pesopt::strprintf ( " solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      }
#endif
  }
   
};



template<typename DataTypeContainer>
class IterativeLinearSystemSolver<DataTypeContainer, EigenConjugateGradient_NEW> :
public LinearSystemSolver<DataTypeContainer> {
    
protected:
    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType  VectorType;
    typedef typename DataTypeContainer::SparseMatrixType  SparseMatrixType;
    
    mutable Eigen::ConjugateGradient_NEW<SparseMatrixType,  Eigen::Lower|Eigen::Upper> _linearSolver;
    
    mutable RealType _tolerance;
    mutable RealType _maxItersFac;
    const string _errorMessage;
    const string _saveDirectory;
    const bool _debug;
    
public:
     
  IterativeLinearSystemSolver( const RealType tolerance = Eigen::NumTraits<RealType>::epsilon(),
                               const RealType maxItersFac = 1.,
                               const string &errorMessage = "", 
                               const string &saveDirectory = "", 
                               const int outputLevel = 0, 
                               const bool debug = true  ) :
  _tolerance ( tolerance ), _maxItersFac ( maxItersFac ),
  _errorMessage( errorMessage ), _saveDirectory ( saveDirectory ), _debug ( debug ) { }   
  
  void setTolerance ( const RealType tolerance ) const { _tolerance = tolerance; }
  void setMaxItersFac ( const RealType maxItersFac ) const { _maxItersFac = maxItersFac; }
  
  //
  bool prepareSolver( const SparseMatrixType& systemMatrix ) const{
      _linearSolver.compute( systemMatrix );
      if(_linearSolver.info()!=Eigen::Success) {
          cout << "prepare linear solver failed" << endl;
          return false;
      }
      return true;
  }
  
  //
  bool solveForPreparedMatrix( VectorType& solution, const VectorType& rhs  ) const{
    solution = _linearSolver.solve( rhs );
    if(_linearSolver.info()!=Eigen::Success ) {
        cerr <<  endl
             << "=================================================" << endl 
             << "error comes from solve in " << _errorMessage << endl
             << "#iterations:     " << _linearSolver.iterations() << endl
             << "estimated error: " << std::setprecision( 32 ) << _linearSolver.error() << endl
             << "but tolerance is: " <<  _tolerance << endl
             << "=================================================" << endl;
        return false;  
    }
            
    return true;
  }
  
  
  bool solveForPreparedMatrixWithGuess( VectorType& solution, const VectorType& rhs, const VectorType& guess ) const{
      
    solution = _linearSolver.solveWithGuess( rhs, guess );
    if(_linearSolver.info()!=Eigen::Success ) {
        cerr <<  endl
             << "=================================================" << endl 
             << "error comes from solve in " << _errorMessage << endl
             << "#iterations:     " << _linearSolver.iterations() << endl
             << "estimated error: " << std::setprecision( 32 ) << _linearSolver.error() << endl
             << "but tolerance is: " <<  _tolerance << endl
             << "=================================================" << endl;
        return false;  
    }
    return true;
  }
  
  
  // solve Ax = b for x
  void solve( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs ) const override{
      bool succesPrepareSolver = prepareSolver( systemMatrix );
      const RealType newTolerance = sqrt( _tolerance * _tolerance * systemMatrix.cols() );
      _linearSolver.setTolerance( newTolerance );
      const int maxIters = static_cast<int> ( _maxItersFac * static_cast<RealType> ( systemMatrix.cols() ) );
      _linearSolver.setMaxIterations(maxIters);
      bool success = solveForPreparedMatrix( solution, rhs );
#ifdef _DEBUGLINEARSOLVER
      if( !success && _debug ) {
        cout << endl 
             << "======================================================" << endl
             << "linear solver failed" << endl;
        cout << "norm of rhs = " << rhs.norm() << endl;
        if(!(_saveDirectory.empty())){
            pesopt::printVectorToFile<VectorType> ( rhs, pesopt::strprintf( "%s/rhs_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver ).c_str(), 32 ); 
            pesopt::printSparseMatrixToFile<SparseMatrixType> ( systemMatrix, pesopt::strprintf( "%s/systemmatrix_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver ).c_str(), 32 ); 
            std::ofstream outSolver ( pesopt::strprintf ( "%s/iterativeSolverOutput_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver  ).c_str ()  );
            _linearSolver.setOutputStream( outSolver );
        }
        solveForPreparedMatrix( solution, rhs );
        debugCounterIterativeSolver++;
        throw invalid_argument ( pesopt::strprintf ( " solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      }
#endif
  }
  
  void solveWithGuess( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs,
                       const VectorType &guess  ) const{
      bool succesPrepareSolver = prepareSolver( systemMatrix );
      const RealType newTolerance = sqrt( _tolerance * _tolerance * systemMatrix.cols() );
      _linearSolver.setTolerance( newTolerance );
      const int maxIters = static_cast<int> ( _maxItersFac * static_cast<RealType> ( systemMatrix.cols() ) );
      _linearSolver.setMaxIterations(maxIters);
      bool success = solveForPreparedMatrixWithGuess( solution, rhs, guess );
#ifdef _DEBUGLINEARSOLVER
      if( !success && _debug ) {
        cout << endl 
             << "======================================================" << endl
             << "linear solver failed" << endl;
        cout << "norm of rhs = " << rhs.norm() << endl;
        if(!(_saveDirectory.empty())){
            pesopt::printVectorToFile<VectorType> ( rhs, pesopt::strprintf( "%s/rhs_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver ).c_str(), 32 ); 
            pesopt::printSparseMatrixToFile<SparseMatrixType> ( systemMatrix, pesopt::strprintf( "%s/systemmatrix_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver ).c_str(), 32 ); 
            std::ofstream outSolver ( pesopt::strprintf ( "%s/iterativeSolverOutput_%d.txt", _saveDirectory.c_str(), debugCounterIterativeSolver  ).c_str()  );
            _linearSolver.setOutputStream( outSolver );
        }
        solveForPreparedMatrixWithGuess( solution, rhs, guess );
        debugCounterIterativeSolver++;
        throw invalid_argument ( pesopt::strprintf ( " solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      }
#endif
  }
   
};




// template< typename RealType,  typename VectorType >
// class DiagonalPreconditioner {
// protected:
//   pesopt::DiagonalMatrix< RealType > _diagInv;
// 
// public:
//   template< typename MatrixType >
//   explicit DiagonalPreconditioner ( const MatrixType &Matrix ) : 
//   _diagInv ( Matrix.getNumRows() ) {
//     _diagInv.setToAbsInverseDiagonalOf ( Matrix );
//   }
// 
//   void apply ( const VectorType &Arg, VectorType &Dest ) const {
//     _diagInv.apply ( Arg, Dest );
//   }
// 
//   template< typename MatrixType >
//   void setForMatrix ( const MatrixType &Matrix ) {
//     _diagInv.reallocate ( Matrix.getNumRows() );
//     _diagInv.setToAbsInverseDiagonalOf ( Matrix );
//   }
// };
// 
// 
// template<typename DataTypeContainer>
// class IterativeLinearSystemSolver<DataTypeContainer, PCGDiagonalPre> :
// public LinearSystemSolver<DataTypeContainer> {
//     
// protected:
//     typedef typename DataTypeContainer::RealType RealType;
//     typedef typename DataTypeContainer::VectorType  VectorType;
//     typedef typename DataTypeContainer::SparseMatrixType  SparseMatrixType;
//     
//     mutable RealType _tolerance;
//     mutable RealType _maxItersFac;
//     const string _errorMessage;
//     const string _saveDirectory;
//     const bool _debug;
//     
// public:
//      
//   IterativeLinearSystemSolver( const RealType tolerance = Eigen::NumTraits<RealType>::epsilon(),
//                                const RealType maxItersFac = 1.,
//                                const string &errorMessage = "", 
//                                const string &saveDirectory = "", 
//                                const int outputLevel = 0, 
//                                const bool debug = true  ) :
//   _tolerance ( tolerance ), _maxItersFac ( maxItersFac ),
//   _errorMessage( errorMessage ), _saveDirectory ( saveDirectory ), _debug ( debug ) { }   
// 
//   
//   void setTolerance ( const RealType tolerance ) const { _tolerance = tolerance; }
//   void setMaxItersFac ( const RealType maxItersFac ) const { _maxItersFac = maxItersFac; }
//   
//   // solve Ax = b for x
//   void solve( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs ) const override{
// 
//     RealType alpha_numer, alpha_denom, beta_numer, beta_denom, spn;
// 
//     VectorType g ( rhs.size() );
//     VectorType d ( rhs.size() );
//     VectorType h ( rhs.size() );
// 
//     h = systemMatrix * solution;
// 
//     g = h;
//     g -= rhs;
// 
//     h.setZero();
//     _DiagonalPreconditioner.apply ( g, h );
// 
//     spn = g.dot(g);
// 
//     d -= h;
// 
//     this->_infoPtr->startIterations ( rhs.normSqr(), spn, "p-cg", "l_2 norm ^2" );
// 
//     while ( ! ( this->_infoPtr->stoppingCriterionIsFulfilled() ) 
//          && ! ( this->_infoPtr->maxIterIsReached() ) 
//          && ! ( this->_infoPtr->currentResidualIsNaN() ) ) {
//       this->_infoPtr->startStep();
// 
//       beta_denom = alpha_numer = g.dot( h );
// 
//       h = systemMatrix * d;
//       alpha_denom = d * h;
// 
//       solution +=  (alpha_numer / alpha_denom ) * d;
// 
//       g +=  (alpha_numer / alpha_denom) * h;
// 
//       h.setZero();
//       _DiagonalPreconditioner.apply ( g, h );
//       beta_numer = g.dot( h );
// 
//       d *= ( beta_numer / beta_denom );
//       d -= h;
// 
//       spn = g.dot( g );
// 
//       this->_infoPtr->finishStep ( spn );
//     }
// 
//     h = systemMatrix * solution;
//     g = h;
//     g -= rhs;
//     spn = g.dot(g);
// 
//     this->_infoPtr->finishIterations ( spn );
//   }
// 
// };




// template<typename DataTypeContainer, class DirectLinearSolverType>
// class DirectLinearSystemSolver :
// public LinearSystemSolver<DataTypeContainer> {
//     
// protected:
//     typedef typename DataTypeContainer::RealType RealType;
//     typedef typename DataTypeContainer::VectorType  VectorType;
//     typedef typename DataTypeContainer::SparseMatrixType  SparseMatrixType;
//     
//     mutable DirectLinearSolverType _linearSolver;
//     
//     const string _errorMessage;
//     const string _saveDirectory;
//     
// public:
//      
//   DirectLinearSystemSolver( const string &errorMessage = "", const string &saveDirectory = "" )
//   : _errorMessage( errorMessage ), _saveDirectory ( saveDirectory ) { }   
//   
//   
//   void analyzePattern( const SparseMatrixType& systemMatrix ) const{_linearSolver.analyzePattern(systemMatrix);}
//   
//   //
//   bool prepareSolver( const SparseMatrixType& systemMatrix ) const{
//       _linearSolver.compute(systemMatrix);
//       if(_linearSolver.info()!=Eigen::Success) {
//           cout << "prepare linear solver failed" << endl;
//           return false;
//       }
//       return true;
//   }
//   
//   //
//   bool solve( VectorType& solution, const VectorType& rhs ) const{
//     bool succes = true;
//     solution = _linearSolver.solve( rhs );
//     return succes;
//   }
//   
//   // solve Ax = b for x
//   void solve( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs ) const override {
//       prepareSolver( systemMatrix );
//       bool success = solve( solution, rhs );
// // #ifdef _DEBUGLINEARSOLVER
// //       if( !success ) {
// //         cout << "failed" << endl;
// //         cout << "norm of rhs = " << rhs.norm() << endl;
// //         pesopt::printVectorToFile<VectorType> ( rhs, pesopt::strprintf( "%s/rhs_%d.txt", _saveDirectory.c_str(), debugCounter ).c_str(), 10 ); 
// //         pesopt::printSparseMatrixToFile<SparseMatrixType> ( systemMatrix, pesopt::strprintf( "%s/systemmatrix_%d.txt", _saveDirectory.c_str(), debugCounter ).c_str(), 10 ); 
// //         debugCounter++;
// //         throw std::invalid_argument ( pesopt::strprintf ( " solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
// //       }
// // #endif
//   }
//    
// };  



enum DirectLinearSolverMethod {
#ifdef PESOPT_WITH_SUITESPARSE
    EigenUmfPackLU = 1,
    EigenCholmodSupernodalLLT = 2,
    PESOPTUmfPackLU = 3,
//  PESOPTCholmod = 4, 
#endif   //PESOPT_WITH_SUITESPARSE
   EigenSparseLU = 5, 
};



template<typename DataTypeContainer, DirectLinearSolverMethod directLSType>
class DirectLinearSystemSolver {};

#ifdef PESOPT_WITH_SUITESPARSE

template<typename DataTypeContainer>
class DirectLinearSystemSolver<DataTypeContainer, EigenUmfPackLU> :
public LinearSystemSolver<DataTypeContainer> {
    
protected:
    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType  VectorType;
    typedef typename DataTypeContainer::SparseMatrixType  SparseMatrixType;
    
    mutable Eigen::UmfPackLU<SparseMatrixType>  _linearSolver;
    
    const string _errorMessage;
    const string _saveDirectory;
    
public:
     
  DirectLinearSystemSolver( const string &errorMessage = "", const string &saveDirectory = "" )
  : _errorMessage( errorMessage ), _saveDirectory ( saveDirectory ) { }   
  
  
  void analyzePattern( const SparseMatrixType& systemMatrix ) const{_linearSolver.analyzePattern(systemMatrix);}
  
  //
  bool prepareSolver( const SparseMatrixType& systemMatrix ) const{
      _linearSolver.compute(systemMatrix);
      if(_linearSolver.info()!=Eigen::Success) {
          cout << "prepare linear solver failed" << endl;
          return false;
      }
      return true;
  }
  
  //
  bool solve( VectorType& solution, const VectorType& rhs ) const{
    bool succes = true;
    solution = _linearSolver.solve( rhs );
    return succes;
  }
  
  // solve Ax = b for x
  void solve( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs ) const override {
      prepareSolver( systemMatrix );
      bool success = solve( solution, rhs );
  }
   
};  



template<typename DataTypeContainer>
class DirectLinearSystemSolver<DataTypeContainer, EigenCholmodSupernodalLLT> :
public LinearSystemSolver<DataTypeContainer> {
    
protected:
    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType  VectorType;
    typedef typename DataTypeContainer::SparseMatrixType  SparseMatrixType;
    
    mutable Eigen::CholmodSupernodalLLT<SparseMatrixType> _linearSolver;

    const string _errorMessage;
    const string _saveDirectory;
    
public:
     
  DirectLinearSystemSolver( const string &errorMessage = "", const string &saveDirectory = "" )
  : _errorMessage( errorMessage ), _saveDirectory ( saveDirectory ) { }   
  
  
  void analyzePattern( const SparseMatrixType& systemMatrix ) const{_linearSolver.analyzePattern(systemMatrix);}
  
  //
  bool prepareSolver( const SparseMatrixType& systemMatrix ) const{
      _linearSolver.compute(systemMatrix);
      if(_linearSolver.info()!=Eigen::Success) {
          cout << "prepare linear solver failed" << endl;
          return false;
      }
      return true;
  }
  
  //
  bool solve( VectorType& solution, const VectorType& rhs ) const{
    bool succes = true;
    solution = _linearSolver.solve( rhs );
    return succes;
  }
  
  // solve Ax = b for x
  void solve( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs ) const override {
      prepareSolver( systemMatrix );
      bool success = solve( solution, rhs );
  }
   
};  


typedef long int SuiteSparseIntType;

template<typename DataTypeContainer>
class DirectLinearSystemSolver<DataTypeContainer, DirectLinearSolverMethod::PESOPTUmfPackLU> :
public LinearSystemSolver<DataTypeContainer> {
    
    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType VectorType;
    typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;

private:
    
    const string _errorMessage;
    const string _saveDirectory;
    
  // the system matrix in umfpack format (number of rows and columns, indices and values of compressed column format)
  mutable SuiteSparseIntType _n;
  mutable std::vector<SuiteSparseIntType> _Ap, _Ai;
  mutable std::vector<RealType> _Ax; 
  mutable std::vector<SuiteSparseIntType> _rowInd, _colInd;
  mutable std::vector<RealType> _val;
  // umfpack control objects (the LU factors of the matrix, control and info variables)
  mutable void *_umfpackNumeric;
  mutable RealType _umfpackControl[UMFPACK_CONTROL];
  mutable RealType _umfpackInfo[UMFPACK_INFO];

public:
  DirectLinearSystemSolver( const string &errorMessage = "", const string &saveDirectory = "" ) :
    _errorMessage( errorMessage ), _saveDirectory ( saveDirectory ), 
    _umfpackNumeric( NULL ) {
    umfpack_dl_defaults( _umfpackControl );
    umfpack_dl_report_control ( _umfpackControl ) ;
  }
  
  ~DirectLinearSystemSolver() {  if ( _umfpackNumeric != NULL ) umfpack_dl_free_numeric( &_umfpackNumeric );}

  
  void analyzePattern( const SparseMatrixType& systemMatrix ) const {
     cout <<  "WARNING: analyzePattern is so for not implemented for DirectLinearSystemSolver<DataTypeContainer, PESOPTUmfPackLU>" <<  endl;   
  }
  
  void prepareSolver ( const SparseMatrixType &systemMatrix ) const {
    // read in all values and put them into the triplet vectors
    _n = static_cast<SuiteSparseIntType>( systemMatrix.rows() );
    
    SuiteSparseIntType numEntries = 0;
    for (int k=0; k<systemMatrix.outerSize(); ++k)
        for ( typename SparseMatrixType::InnerIterator it(systemMatrix,k); it; ++it){
                numEntries++;
        }
        
    _colInd.resize( numEntries );
    _rowInd.resize( numEntries );
    _val.resize( numEntries  );  
    
    SuiteSparseIntType idx = 0;
    for (int k=0; k<systemMatrix.outerSize(); ++k)
        for ( typename SparseMatrixType::InnerIterator it(systemMatrix,k); it; ++it){
                _rowInd[idx] = it.row(); 
                _colInd[idx] = it.col();
                _val[idx] = it.value();
                idx++;
        }

    // convert triplet format into required format  
    _Ap.resize( _n + 1 );
    _Ai.resize( static_cast<SuiteSparseIntType>(_val.size()) );
    _Ax.resize( static_cast<SuiteSparseIntType>(_val.size()) );
    umfpack_dl_triplet_to_col( _n, _n, _val.size(), _rowInd.data(), _colInd.data(), _val.data(), _Ap.data(), _Ai.data(), _Ax.data(), NULL );

    // perform symbolic and numeric factorization
    if ( _umfpackNumeric != NULL ) umfpack_dl_free_numeric( &_umfpackNumeric );
    void *umfpackSymbolic;

    // check factorization status
    int status;

    // symbolic factorization
    status = umfpack_dl_symbolic( _n, _n, _Ap.data(), _Ai.data(), _Ax.data(), &umfpackSymbolic, _umfpackControl, _umfpackInfo );
    if ( status != UMFPACK_OK ){
      umfpack_dl_report_status ( _umfpackControl, status) ;
      cerr << "umfpack_dl_symbolic failed!" << endl;
    }
    // numeric factorization
    status = umfpack_dl_numeric( _Ap.data(), _Ai.data(), _Ax.data(), umfpackSymbolic, &_umfpackNumeric, _umfpackControl, _umfpackInfo );
    if ( status != UMFPACK_OK ){
      umfpack_dl_report_status ( _umfpackControl, status) ;
      cerr << "umfpack_dl_numeric failed!" << endl;
    }
    // free symbolic factorization
    umfpack_dl_free_symbolic( &umfpackSymbolic );
  }
  
  
  
  VectorType solve ( const VectorType &rhs ) const {
    VectorType sol( rhs.size() );
    // solve and copy solution into result vector
    umfpack_dl_solve( UMFPACK_A, _Ap.data(), _Ai.data(), _Ax.data(), sol.data(), rhs.data(), _umfpackNumeric, _umfpackControl, _umfpackInfo );
    // if umfpackControl [UMFPACK_PRL] is set, then this we give you al the information you need
    umfpack_dl_report_info (_umfpackControl, _umfpackInfo);
    return sol;
  }
  
  //
  bool solve( VectorType& solution, const VectorType& rhs ) const{
    bool succes = true;
    solution = solve( rhs );
    return succes;
  }
  
  // solve Ax = b for x
  void solve( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs ) const override {
      prepareSolver( systemMatrix );
      bool success = solve( solution, rhs );
  }
  
};



# endif   //PESOPT_WITH_SUITESPARSE






template<typename DataTypeContainer>
class DirectLinearSystemSolver<DataTypeContainer, EigenSparseLU> :
public LinearSystemSolver<DataTypeContainer> {
    
protected: 
    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType  VectorType;
    typedef typename DataTypeContainer::SparseMatrixType  SparseMatrixType;
    
    mutable Eigen::SparseLU<SparseMatrixType>  _linearSolver;
    
    const string _errorMessage;
    const string _saveDirectory;
    
public:
     
  DirectLinearSystemSolver( const string &errorMessage = "", const string &saveDirectory = "" )
  : _errorMessage( errorMessage ), _saveDirectory ( saveDirectory ) { }   
  
  
  void analyzePattern( const SparseMatrixType& systemMatrix ) const{_linearSolver.analyzePattern(systemMatrix);}
  
  //
  bool prepareSolver( const SparseMatrixType& systemMatrix ) const{
      _linearSolver.compute(systemMatrix);
      if(_linearSolver.info()!=Eigen::Success) {
          cout << "prepare linear solver failed" << endl;
          return false;
      }
      return true;
  }
  
  //
  bool solve( VectorType& solution, const VectorType& rhs ) const{
    bool succes = true;
    solution = _linearSolver.solve( rhs );
    return succes;
  }
  
  // solve Ax = b for x
  void solve( const SparseMatrixType& systemMatrix, VectorType& solution, const VectorType& rhs ) const override {
      prepareSolver( systemMatrix );
      bool success = solve( solution, rhs );
  }
   
};  




#endif //__LINEARSYSTEMSOLVER_H
