#ifndef __IPOPTNONLINEARCONSTRAINTSOLVER_H
#define __IPOPTNONLINEARCONSTRAINTSOLVER_H


#ifdef PESOPT_WITH_IPOPT

#include <includesIpopt.h>
#include <solverInfo.h>

/* use IPOPT to solve a constraint minimization problem of type 
 *  minimize f(x) 
 * over all x s.t. 
 * x_l \leq x \leq x_u
 * g_l(x) \leq g_i(x) \leq g_u(x)
 */


template<typename DataTypeContainer >
class IpoptConstraintFirstOrderSolverBase : public Ipopt::TNLP {

protected:
  
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
    
  const pesopt::NonlinearEnergyOp<DataTypeContainer> &_energyOp;  
  const pesopt::NonlinearConstraintOps<DataTypeContainer> &_constraintOps; 

  const VectorType &_startingPoint;
  VectorType &_nextIter;
 
  unsigned int _numTotalDofs;
  unsigned int _numConstraints;
  
  const VectorType &_x_l, &_x_u;
  RealType _nlpError;
  int _iter;

public:
    
  IpoptConstraintFirstOrderSolverBase ( 
    const pesopt::NonlinearEnergyOp<DataTypeContainer> &energyOp,
    const pesopt::NonlinearConstraintOps<DataTypeContainer> &constraintOps,
    const VectorType &startingPoint, VectorType &nextIter,
    const VectorType &x_l, const VectorType &x_u )
  : _energyOp ( energyOp ),  _constraintOps ( constraintOps ),
    _startingPoint ( startingPoint ), 
    _nextIter ( nextIter ),
    _numTotalDofs ( startingPoint.size () ), 
    _numConstraints ( constraintOps.getNumConstraints() ),
    _x_l ( x_l ), _x_u ( x_u )
    { }
    

  //returns info about the nlp
  virtual bool get_nlp_info ( Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style ) {
    // The number of DOFs.
    n = _numTotalDofs;
    // The number of constraints.
    m = _numConstraints;
    // The number of nonzeros in the jacobian: Possibly every entry.
    nnz_jac_g = _numConstraints * _numTotalDofs;
    //Nonzeros in hessian: 0 because not implemented
    nnz_h_lag = 0;
    //Use C index style
    index_style = C_STYLE;

    return true;
  }

  virtual bool get_bounds_info ( Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, 
                                 Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u ) {
      for ( Ipopt::Index i = 0; i < n; ++i ) {
          x_l[i] = _x_l[i]; x_u[i] = _x_u[i];
      }
      for ( Ipopt::Index i = 0; i < _numConstraints; ++i ) {
          g_l[i] = _constraintOps.getLowerBound(i); g_u[i] = _constraintOps.getUpperBound(i);
      }
      return true;
  }

  virtual bool get_starting_point ( Ipopt::Index n, bool init_x,  Ipopt::Number* x, 
                                    bool init_z, Ipopt::Number* /*z_L*/, Ipopt::Number* /*z_U*/, 
                                    Ipopt::Index /*m*/, bool init_lambda, Ipopt::Number* /*lambda*/ ) {
    init_x = true; init_z = false; init_lambda = false;
    for ( int j = 0; j < _startingPoint.size (); ++j ) x[j] = _startingPoint[j];
    return true;
  }

  virtual bool eval_f ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number& obj_value ) {
    // Convert
    VectorType v ( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    RealType s;
    _energyOp.evaluateEnergy ( v, s );
    obj_value = s;
    
    return true;
  }

  virtual bool eval_grad_f ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f ) {
    // Convert
    VectorType v ( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

   VectorType res ( v.size() );
   _energyOp.evaluateJacobian( v, res );
   for ( int i = 0; i < n; ++i ) grad_f[i] = res[i];
   
   return true;
  }

  virtual bool eval_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number *g ) {
    // Convert
    VectorType v ( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];
    
    for( int j=0; j<_numConstraints; ++j ){
        RealType s;
        _constraintOps.evaluateEnergy( j, v, s );
        g[j] = s;
    }

    return true;
  }

  virtual bool eval_jac_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                            Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values ) {
      
    if (values == NULL) {
        // return the structure of the Jacobian
        // this particular Jacobian is dense
        Ipopt::Index index = 0;
        for( Ipopt::Index j=0; j<m; ++j )
            for( Ipopt::Index i = 0; i < n; ++i ){
                iRow[index] = j; jCol[index] = i;
                index++;
            }
    }
    else {
        // return the values of the Jacobian of the constraints
        // Convert
        VectorType v ( n );
        for ( int i = 0; i < n; ++i ) v[i] = x[i];

        Ipopt::Index index = 0;
        for( Ipopt::Index j=0; j<m; ++j ){
            VectorType res ( v.size() );
            _constraintOps.evaluateJacobian( j, v, res );
            for( Ipopt::Index i = 0; i < n; ++i ){
                values[index] = res[i];
                index++;
            }
        }
    }
      
    return true;
  }
  
  virtual bool intermediate_callback ( Ipopt::AlgorithmMode /*mode*/, Ipopt::Index iter, 
                                       Ipopt::Number /*obj_value*/,  Ipopt::Number /*inf_pr*/,  Ipopt::Number /*inf_du*/, Ipopt::Number /*mu*/,  Ipopt::Number /*d_norm*/, Ipopt::Number /*regularization_size*/, Ipopt::Number /*d_du*/, Ipopt::Number /*d_pr*/,  Ipopt::Index /*ls_trials*/, 
                                       const Ipopt::IpoptData* ip_data,
                                       Ipopt::IpoptCalculatedQuantities* ip_cq ) {

    // Retrieve primal variables (code similar to code from the IPOPT documentation):
    Ipopt::TNLPAdapter* tnlpAdapter = NULL;

    if ( ip_cq != 0 ) {
      Ipopt::OrigIpoptNLP* origNLP;
      origNLP = dynamic_cast< Ipopt::OrigIpoptNLP*> ( Ipopt::GetRawPtr ( ip_cq->GetIpoptNLP () ) );
      // If in restoration mode, origNLP will be NULL. Quit method in this case
      if ( origNLP == NULL )  return true;
      tnlpAdapter = dynamic_cast < Ipopt::TNLPAdapter* > ( GetRawPtr ( origNLP->nlp () ) );
    }
    else throw std::invalid_argument ( pesopt::strprintf ( "ip_cq == NULL in FILE %s at line %d.", __FILE__, __LINE__ ).c_str() );
    double *x = new double[_startingPoint.size ()];   // n == _numDofs
    tnlpAdapter->ResortX ( *ip_data->curr ()->x (), x );
    delete[] x;

    // Plot solution.
     _iter = iter;
    
    return true;
  }

  virtual void finalize_solution ( Ipopt::SolverReturn /*status*/, Ipopt::Index /*n*/,  const Ipopt::Number* x,
                                   const Ipopt::Number* /*z_L*/,  const Ipopt::Number* /*z_U*/, 
                                   Ipopt::Index /*m*/,  const Ipopt::Number* /*g*/, const Ipopt::Number* /*lambda*/, 
                                   Ipopt::Number /*obj_value*/,  const Ipopt::IpoptData* /*ip_data*/,  Ipopt::IpoptCalculatedQuantities* ip_cq ) {
    for ( int i = 0; i < _nextIter.size (); ++i ) _nextIter[i] = x[i];
    _nlpError = ip_cq->curr_nlp_error();
  }
  
  public : 
    RealType getNLPError( ) const { return _nlpError; }
    int getNumIterations() const { return _iter; }
  
};




template <typename DataTypeContainer >
class IpoptConstraintFirstOrderSolver {
protected:
 
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
    
  const pesopt::NonlinearEnergyOp<DataTypeContainer> &_energyOp;  
  const pesopt::NonlinearConstraintOps<DataTypeContainer> &_constraintOp;  
  const RealType _ipoptTol;
  const int _MaxIterations;
  mutable VectorType _x_l, _x_u;
  const int _linearSolverTypeIpopt;
  const RealType _boundRelaxFactorIpopt;
  const int _ipoptPrintLevel;
  const string _ipoptOutputFile;

public:
  IpoptConstraintFirstOrderSolver ( const pesopt::NonlinearEnergyOp<DataTypeContainer> &E,
                                    const pesopt::NonlinearConstraintOps<DataTypeContainer> &constraintOp,
                                    const int MaxIterations, const RealType Tolerance,
                                    const RealType x_l, const RealType x_u,
                                    const int linearSolverTypeIpopt = 0,
                                    const int ipoptPrintLevel = 5,
                                    const string ipoptOutputFile = ""
                            )
  : _energyOp ( E ), _constraintOp ( constraintOp ),
    _ipoptTol ( Tolerance ),
    _MaxIterations ( MaxIterations ),
    _x_l( E.getNumDofs() ), _x_u( E.getNumDofs() ),
    _linearSolverTypeIpopt ( linearSolverTypeIpopt ),
    _ipoptPrintLevel ( ipoptPrintLevel ), _ipoptOutputFile ( ipoptOutputFile ) {
        for( int i=0; i<E.getNumDofs(); ++i ){
         _x_l[i] = x_l; _x_u[i] = x_u;   
        }
    }
    
  IpoptConstraintFirstOrderSolver ( const pesopt::NonlinearEnergyOp<DataTypeContainer> &E,
                              const pesopt::NonlinearConstraintOps<DataTypeContainer> &constraintOp,
                              const int MaxIterations, const RealType Tolerance,
                              const VectorType & x_l, const VectorType & x_u,
                              const int linearSolverTypeIpopt = 0,
                              const RealType boundRelaxFactorIpopt = 1.e-8, 
                              const int ipoptPrintLevel = 5,
                              const string ipoptOutputFile = ""
                            )
  : _energyOp ( E ), _constraintOp ( constraintOp ),
    _ipoptTol ( Tolerance ),
    _MaxIterations ( MaxIterations ),
    _x_l( x_l ), _x_u( x_u ),
    _linearSolverTypeIpopt ( linearSolverTypeIpopt ),
    _boundRelaxFactorIpopt ( boundRelaxFactorIpopt ), 
     _ipoptPrintLevel ( ipoptPrintLevel ), _ipoptOutputFile ( ipoptOutputFile ) { }

  virtual ~IpoptConstraintFirstOrderSolver() {}

  void solve ( const VectorType &StartPosition, VectorType &Solution, SolverInfo<DataTypeContainer> &solverInfo ) const {

    // Set up masking operators:
    Ipopt::SmartPtr<IpoptConstraintFirstOrderSolverBase<DataTypeContainer> > tsOpt 
    = new IpoptConstraintFirstOrderSolverBase<DataTypeContainer> ( _energyOp, _constraintOp, StartPosition, Solution, _x_l, _x_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp = new Ipopt::IpoptApplication ();
    
    setToleranceOptionsIpopt<RealType>( _ipoptTol, _MaxIterations, ipoptApp );
    switchLinearSolverTypeIpopt ( _linearSolverTypeIpopt, ipoptApp );
    setPrintOptionsIpopt( _ipoptPrintLevel, _ipoptOutputFile, ipoptApp );
    
    // Enable quasi-Newton hessian approximation
    ipoptApp->Options()->SetStringValue ( "hessian_approximation", "limited-memory" );
    
    ipoptApp->Options()->SetNumericValue ( "bound_relax_factor", _boundRelaxFactorIpopt );
    
    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();

    Ipopt::SmartPtr<Ipopt::TNLP> nlp;
    
    //! Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
    outputIpoptStatus ( ipoptStatus, true );
    
    solverInfo.setSolverStatus( getIpoptStatus( ipoptStatus ).c_str() );
    solverInfo.setError( tsOpt->getNLPError() );
    solverInfo.setNumIterations( tsOpt->getNumIterations() );
    
  }
  
  void testDerivative ( const VectorType &StartPosition ) const {

    // Set up masking operators:
    VectorType Solution ( StartPosition.size() );
    Ipopt::SmartPtr<IpoptConstraintFirstOrderSolverBase<DataTypeContainer> > tsOpt 
    = new IpoptConstraintFirstOrderSolverBase<DataTypeContainer> ( _energyOp, _constraintOp, StartPosition, Solution, _x_l, _x_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp = new Ipopt::IpoptApplication ();
    
    setPrintOptionsIpopt( _ipoptPrintLevel, _ipoptOutputFile, ipoptApp );
    ipoptApp->Options ()->SetStringValue  ( "derivative_test", "first-order" );
    setDerivateTestOptionsIpopt ( ipoptApp );

    ipoptApp->Options()->SetStringValue ( "hessian_approximation", "limited-memory" );

    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();
    
    Ipopt::SmartPtr<Ipopt::TNLP> nlp;
    
    //! Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
  }
  
  
};












template<typename DataTypeContainer >
class IpoptConstraintSecondOrderSolverBase : public Ipopt::TNLP {

protected:
  
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  typedef typename DataTypeContainer::TripletType TripletType;
    
  const pesopt::NonlinearEnergyOp<DataTypeContainer> &_energyOp;
  const pesopt::SparseJacobianNonlinearConstraintOps<DataTypeContainer> &_constraintOps;  

  const VectorType &_startingPoint;
  VectorType &_nextIter;
  const VectorType &_startingLambda;
  VectorType &_nextLambda;
 
  unsigned int _numTotalDofs;
  unsigned int _numConstraints;
  
  const VectorType& _x_l_vec; const VectorType& _x_u_vec;
  RealType _nlpError;
  int _iter;

public:
  IpoptConstraintSecondOrderSolverBase ( const pesopt::NonlinearEnergyOp<DataTypeContainer> &energyOp,
                                              const pesopt::SparseJacobianNonlinearConstraintOps<DataTypeContainer> &constraintOps,
                                              const VectorType &startingPoint, VectorType &nextIter,
                                              const VectorType &startingLambda, VectorType &nextLambda,
                                              const VectorType& x_l_vec, const VectorType& x_u_vec
                                            )
  : _energyOp ( energyOp ), _constraintOps ( constraintOps ),
    _startingPoint ( startingPoint ), _nextIter ( nextIter ),
    _startingLambda ( startingLambda ), _nextLambda ( nextLambda ),
    _numTotalDofs ( startingPoint.size () ), 
    _numConstraints ( constraintOps.getNumConstraints() ),
    _x_l_vec( x_l_vec ), _x_u_vec( x_u_vec )
    { }
    

  //returns info about the nlp
  virtual bool get_nlp_info ( Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style ) {
    // The number of DOFs.
    n = _numTotalDofs;
    // The number of constraints.
    m = _numConstraints;
    // The number of nonzeros in the jacobian: Possibly every entry.
    nnz_jac_g = _constraintOps.sizeJacobian();
    //Nonzeros in hessian
    std::vector<TripletType> tripletListHessian;
    VectorType tmpVec ( n ); tmpVec.setZero();
    _energyOp.evaluateTripletListHessianSym( tmpVec, tripletListHessian );
    nnz_h_lag = tripletListHessian.size() + _constraintOps.sizeHessian();
    //Use C index style
    index_style = C_STYLE;

    return true;
  }

  virtual bool get_bounds_info ( Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, 
                                 Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u ) {
    for ( Ipopt::Index i = 0; i < n; ++i ) {
          x_l[i] = _x_l_vec[i]; x_u[i] = _x_u_vec[i];
    }
    for ( Ipopt::Index i = 0; i < _numConstraints; ++i ) {
          g_l[i] = _constraintOps.getLowerBound(i); g_u[i] = _constraintOps.getUpperBound(i);
    }
    return true;
  }

  virtual bool get_starting_point ( Ipopt::Index n, 
                                    bool init_x,  Ipopt::Number* x, 
                                    bool init_z, Ipopt::Number* /*z_L*/, Ipopt::Number* /*z_U*/, 
                                    Ipopt::Index /*m*/, 
                                    bool init_lambda, Ipopt::Number* /*lambda*/ ) {
    init_x = true; init_z = false;
    for ( int j = 0; j < _startingPoint.size (); ++j ) x[j] = _startingPoint[j];
    init_lambda = false;
  
    return true;
  }

  virtual bool eval_f ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number& obj_value ) {
    // Convert
    VectorType v ( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    RealType s;
    _energyOp.evaluateEnergy ( v, s );
    obj_value = s;
    
    return true;
  }

  virtual bool eval_grad_f ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f ) {
    // Convert
    VectorType v ( n ); for ( int i = 0; i < n; ++i ) v[i] = x[i];
    VectorType res ( v.size() ); _energyOp.evaluateJacobian( v, res );  for ( int i = 0; i < n; ++i ) grad_f[i] = res[i];
    return true;
  }

  virtual bool eval_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number *g ) {
    // Convert
    VectorType v ( n ); for ( int i = 0; i < n; ++i ) v[i] = x[i];
    
    VectorType constraintVec ( _numConstraints );
    _constraintOps.evaluate( v, constraintVec );

    for( int j=0; j<_numConstraints; ++j ) g[j] = constraintVec[j];
    
    return true;
  }

  virtual bool eval_jac_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                            Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values ) {
      
    if (values == NULL) {
        
        VectorType tmpVec ( n ); tmpVec.setZero();
        std::vector<TripletType> tripletListJac;
        _constraintOps.evaluateJacobian( tmpVec, tripletListJac );
        for( int index = 0; index < tripletListJac.size(); ++index ){
            iRow[index] = tripletListJac[index].row();
            jCol[index] = tripletListJac[index].col();
        }
    }
    else {
        
        // return the values of the Jacobian of the constraints
        VectorType v ( n ); for ( int i = 0; i < n; ++i ) v[i] = x[i];

        std::vector<TripletType> tripletListJac;
        _constraintOps.evaluateJacobian( v, tripletListJac );
        for( int index = 0; index < tripletListJac.size(); ++index ){
            values[index] = tripletListJac[index].value();
        }
    }
    
    return true;
  }
  
  
    // Hessian of Lagrangian: L = sigma * D^2 E + lambda D^2 Constraint
    // since Hessian is symmetric, only fill lower left part (otherwise IPOPT does not understand the matrix correctly)
  virtual bool eval_h( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor,
                       Ipopt::Index m, const Ipopt::Number* lambda, bool new_lambda, 
                       Ipopt::Index nele_hess, Ipopt::Index* iRow,Ipopt::Index* jCol, Ipopt::Number* values){
    if (values == NULL) {
        // return the structure.
        Ipopt::Index idx=0;
        VectorType tmpVec ( n ); tmpVec.setZero();
        
        std::vector<TripletType> tripletListHessian;
        _energyOp.evaluateTripletListHessianSym( tmpVec, tripletListHessian );
        for( int index = 0; index < tripletListHessian.size(); ++index ){
              iRow[idx] = tripletListHessian[index].row();
              jCol[idx] = tripletListHessian[index].col();
              idx++;
        }
    
        std::vector<std::vector<TripletType>> tripletListConstraintHessian ( _numConstraints );
        _constraintOps.evaluateHessian( tmpVec, tripletListConstraintHessian );
        for( int constraint = 0; constraint < _numConstraints; ++ constraint ){
            for( int index = 0; index < tripletListConstraintHessian[constraint].size(); ++index ){
              iRow[idx] = tripletListConstraintHessian[constraint][index].row();
              jCol[idx] = tripletListConstraintHessian[constraint][index].col();
              idx++;
            }
        } 
    }
    else {
        // Convert
        VectorType v ( n ); for ( int i = 0; i < n; ++i ) v[i] = x[i];

        Ipopt::Index idx=0;
        
        std::vector<TripletType> tripletListHessian;     
        _energyOp.evaluateTripletListHessianSym( v, tripletListHessian );
        for( int index = 0; index < tripletListHessian.size(); ++index ){
              values[idx] = obj_factor * tripletListHessian[index].value();
              idx++;
        }
        
        std::vector<std::vector<TripletType>> tripletListConstraintHessian ( _numConstraints );
        _constraintOps.evaluateHessian( v, tripletListConstraintHessian );
        for( int constraint = 0; constraint < _numConstraints; ++ constraint ){
            for( int index = 0; index < tripletListConstraintHessian[constraint].size(); ++index ){
              values[idx] = lambda[constraint] * tripletListConstraintHessian[constraint][index].value();
              idx++;
            }
        }
    }

    return true;
 }
  
  
  virtual bool intermediate_callback ( Ipopt::AlgorithmMode /*mode*/, Ipopt::Index iter, 
                                       Ipopt::Number /*obj_value*/,  Ipopt::Number /*inf_pr*/,  Ipopt::Number /*inf_du*/, Ipopt::Number /*mu*/,  Ipopt::Number /*d_norm*/, Ipopt::Number /*regularization_size*/, Ipopt::Number /*d_du*/, Ipopt::Number /*d_pr*/,  Ipopt::Index /*ls_trials*/, 
                                       const Ipopt::IpoptData* ip_data,
                                       Ipopt::IpoptCalculatedQuantities* ip_cq ) {

    // Retrieve primal variables (code similar to code from the IPOPT documentation):
    Ipopt::TNLPAdapter* tnlpAdapter = NULL;

    if ( ip_cq != 0 ) {
      Ipopt::OrigIpoptNLP* origNLP;
      origNLP = dynamic_cast< Ipopt::OrigIpoptNLP*> ( Ipopt::GetRawPtr ( ip_cq->GetIpoptNLP () ) );
      // If in restoration mode, origNLP will be NULL. Quit method in this case
      if ( origNLP == NULL )  return true;
      tnlpAdapter = dynamic_cast < Ipopt::TNLPAdapter* > ( GetRawPtr ( origNLP->nlp () ) );
    }
    else throw std::invalid_argument ( pesopt::strprintf ( "ip_cq == NULL in FILE %s at line %d.", __FILE__, __LINE__ ).c_str() );
    double *x = new double[_startingPoint.size ()];   // n == _numDofs
    tnlpAdapter->ResortX ( *ip_data->curr ()->x (), x );
    delete[] x;

    // Plot solution.
     _iter = iter;
    
    return true;
  }

  virtual void finalize_solution ( Ipopt::SolverReturn /*status*/, Ipopt::Index /*n*/,  const Ipopt::Number* x,
                                   const Ipopt::Number* /*z_L*/,  const Ipopt::Number* /*z_U*/, 
                                   Ipopt::Index /*m*/,  const Ipopt::Number* /*g*/, const Ipopt::Number* lambda, 
                                   Ipopt::Number /*obj_value*/,  const Ipopt::IpoptData* /*ip_data*/,  Ipopt::IpoptCalculatedQuantities* ip_cq ) {
    for ( int i = 0; i < _nextIter.size (); ++i ) _nextIter[i] = x[i];
    for ( int i = 0; i < _nextLambda.size (); ++i ) _nextLambda[i] = lambda[i];
    _nlpError = ip_cq->curr_nlp_error();
  }
  
  public : 
    RealType getNLPError( ) const { return _nlpError; }
    int getNumIterations() const { return _iter; }
  
};




template <typename DataTypeContainer >
class IpoptConstraintSecondOrderSolver {
protected:
 
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  typedef typename DataTypeContainer::TripletType TripletType;  
  
  const pesopt::NonlinearEnergyOp<DataTypeContainer> &_energyOp;  
  const pesopt::SparseJacobianNonlinearConstraintOps<DataTypeContainer> &_constraintOps;  
  const RealType _ipoptTol;
  const int _MaxIterations;
  mutable VectorType _x_l, _x_u;
  const int _linearSolverTypeIpopt;
  const RealType _boundRelaxFactorIpopt;
  const int _ipoptPrintLevel;
  const string _ipoptOutputFile;

public:

   IpoptConstraintSecondOrderSolver ( const pesopt::NonlinearEnergyOp<DataTypeContainer> &E,
                                      const pesopt::SparseJacobianNonlinearConstraintOps<DataTypeContainer> &constraintOps,
                                      const int MaxIterations, const RealType Tolerance,
                                      const RealType x_l, const RealType x_u,
                                      const int linearSolverTypeIpopt = 0,
                                      const RealType boundRelaxFactorIpopt = 1.e-8, 
                                      const int ipoptPrintLevel = 5,
                                      const string ipoptOutputFile = ""
                            )
  : _energyOp ( E ), _constraintOps ( constraintOps ),
    _ipoptTol ( Tolerance ),
    _MaxIterations ( MaxIterations ),
    _x_l( E.getNumDofs() ), _x_u( E.getNumDofs() ),
    _linearSolverTypeIpopt ( linearSolverTypeIpopt ),
    _boundRelaxFactorIpopt ( boundRelaxFactorIpopt ), 
    _ipoptPrintLevel ( ipoptPrintLevel ), _ipoptOutputFile ( ipoptOutputFile ) {
        for( int i=0; i<E.getNumDofs(); ++i ){
         _x_l[i] = x_l; _x_u[i] = x_u;   
        }
    }
    
    IpoptConstraintSecondOrderSolver (  const pesopt::NonlinearEnergyOp<DataTypeContainer> &energyOp,
                                        const pesopt::SparseJacobianNonlinearConstraintOps<DataTypeContainer> &constraintOps,
                                        const int MaxIterations, const RealType Tolerance,
                                        const VectorType& x_l, const VectorType& x_u,
                                        const int linearSolverTypeIpopt = 0,
                                        const RealType boundRelaxFactorIpopt = 1.e-8, 
                                        const int ipoptPrintLevel = 5, 
                                        const string ipoptOutputFile = "" )
  : _energyOp ( energyOp ), _constraintOps ( constraintOps ),
    _ipoptTol ( Tolerance ),
    _MaxIterations ( MaxIterations ),
    _x_l( x_l ), _x_u( x_u ),
    _linearSolverTypeIpopt ( linearSolverTypeIpopt ),
    _boundRelaxFactorIpopt ( boundRelaxFactorIpopt ), 
    _ipoptPrintLevel ( ipoptPrintLevel ), _ipoptOutputFile ( ipoptOutputFile ){   }

  virtual ~IpoptConstraintSecondOrderSolver() {}

  void solve ( const VectorType &StartPosition, VectorType &Solution, 
               const VectorType &StartLambda, VectorType &SolLambda, 
               SolverInfo<DataTypeContainer> &solverInfo ) const {

    // Set up masking operators:
    Ipopt::SmartPtr<IpoptConstraintSecondOrderSolverBase<DataTypeContainer> > tsOpt 
    = new IpoptConstraintSecondOrderSolverBase<DataTypeContainer> ( _energyOp, _constraintOps, StartPosition, Solution, StartLambda, SolLambda, _x_l, _x_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp = new Ipopt::IpoptApplication ();
    
    setToleranceOptionsIpopt<RealType>( _ipoptTol, _MaxIterations, ipoptApp );
    switchLinearSolverTypeIpopt ( _linearSolverTypeIpopt, ipoptApp );
    setPrintOptionsIpopt ( _ipoptPrintLevel, _ipoptOutputFile, ipoptApp );
    
    ipoptApp->Options()->SetStringValue ( "hessian_approximation", "exact" );
    
    ipoptApp->Options()->SetNumericValue ( "bound_relax_factor", _boundRelaxFactorIpopt );

    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();

    //! Optimization step
    Ipopt::SmartPtr<Ipopt::TNLP> nlp;
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
    outputIpoptStatus ( ipoptStatus, true );
    
    solverInfo.setSolverStatus( getIpoptStatus( ipoptStatus ).c_str() );
    solverInfo.setError( tsOpt->getNLPError() );
    solverInfo.setNumIterations( tsOpt->getNumIterations() );
    
  }
  
  
  void testDerivative ( const VectorType &StartPosition, const VectorType &StartLambda, const int order ) const {

    // Set up masking operators:
    VectorType Solution ( StartPosition.size() ); VectorType SolLambda ( StartLambda.size() );
    Ipopt::SmartPtr<IpoptConstraintSecondOrderSolverBase<DataTypeContainer> > tsOpt 
    = new IpoptConstraintSecondOrderSolverBase<DataTypeContainer> ( _energyOp, _constraintOps, StartPosition, Solution, StartLambda, SolLambda, _x_l, _x_u  );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp = new Ipopt::IpoptApplication ();
    
    setPrintOptionsIpopt ( _ipoptPrintLevel, _ipoptOutputFile, ipoptApp );
    
    //Derivative Test
    if( order == 1 ) ipoptApp->Options()->SetStringValue  ( "derivative_test", "first-order" );
    if( order == 2 ) ipoptApp->Options()->SetStringValue  ( "derivative_test", "second-order" );
    setDerivateTestOptionsIpopt ( ipoptApp );
    ipoptApp->Options()->SetStringValue ( "hessian_approximation", "exact" );

    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();
    
    Ipopt::SmartPtr<Ipopt::TNLP> nlp;
    
    //! Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
    
  }
  
};


















template<typename DataTypeContainer >
class IpoptConstraintSecondOrderSolverBase_QuadraticEnergy : public Ipopt::TNLP {

protected:
  
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  typedef typename DataTypeContainer::TripletType TripletType;
    
  const pesopt::QuadraticEnergyOp<DataTypeContainer> &_energyOp;  
  const std::vector<TripletType> &_HessianEnergyTripletList;
  const pesopt::SparseJacobianNonlinearConstraintOps<DataTypeContainer> &_constraintOps;  

  const VectorType &_startingPoint;
  VectorType &_nextIter;
  const VectorType &_startingLambda;
  VectorType &_nextLambda;
 
  unsigned int _numTotalDofs;
  unsigned int _numConstraints;
  
  const RealType _x_l, _x_u;
  RealType _nlpError;
  int _iter;

public:
  IpoptConstraintSecondOrderSolverBase_QuadraticEnergy ( const pesopt::QuadraticEnergyOp<DataTypeContainer> &energyOp,
                                             const std::vector<TripletType> &HessianEnergyTripletList,
                                             const pesopt::SparseJacobianNonlinearConstraintOps<DataTypeContainer> &constraintOps,
                                             const VectorType &startingPoint, VectorType &nextIter,
                                             const VectorType &startingLambda, VectorType &nextLambda,
                                             const RealType x_l = -2.e+19, const RealType x_u = 2.e+19
)
  : _energyOp ( energyOp ), _HessianEnergyTripletList ( HessianEnergyTripletList ),  _constraintOps ( constraintOps ),
    _startingPoint ( startingPoint ), _nextIter ( nextIter ),
    _startingLambda ( startingLambda ), _nextLambda ( nextLambda ),
    _numTotalDofs ( startingPoint.size () ), 
    _numConstraints ( constraintOps.getNumConstraints() ),
    _x_l ( x_l ), _x_u ( x_u )
    { }
    

  //returns info about the nlp
  virtual bool get_nlp_info ( Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style ) {
    // The number of DOFs.
    n = _numTotalDofs;
    // The number of constraints.
    m = _numConstraints;
    // The number of nonzeros in the jacobian: Possibly every entry.
    nnz_jac_g = _constraintOps.sizeJacobian();
    //Nonzeros in hessian
    nnz_h_lag = _HessianEnergyTripletList.size() + _constraintOps.sizeHessian();
    //Use C index style
    index_style = C_STYLE;

    return true;
  }

  virtual bool get_bounds_info ( Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, 
                                 Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u ) {
    for ( Ipopt::Index i = 0; i < n; ++i ) {
      x_l[i] = _x_l; x_u[i] = _x_u;
    }
    for ( Ipopt::Index i = 0; i < _numConstraints; ++i ) {
          g_l[i] = _constraintOps.getLowerBound(i); g_u[i] = _constraintOps.getUpperBound(i);
    }
    return true;
  }

  virtual bool get_starting_point ( Ipopt::Index n, 
                                    bool init_x,  Ipopt::Number* x, 
                                    bool init_z, Ipopt::Number* /*z_L*/, Ipopt::Number* /*z_U*/, 
                                    Ipopt::Index /*m*/, 
                                    bool init_lambda, Ipopt::Number* /*lambda*/ ) {
    init_x = true; init_z = false;
    for ( int j = 0; j < _startingPoint.size (); ++j ) x[j] = _startingPoint[j];
    init_lambda = false;

    return true;
  }

  virtual bool eval_f ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number& obj_value ) {
    // Convert
    VectorType v ( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

    RealType s;
    _energyOp.evaluateEnergy ( v, s );
    obj_value = s;
    
    return true;
  }

  virtual bool eval_grad_f ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f ) {
    // Convert
    VectorType v ( n ); for ( int i = 0; i < n; ++i ) v[i] = x[i];
    VectorType res ( v.size() ); _energyOp.evaluateJacobian( v, res );  for ( int i = 0; i < n; ++i ) grad_f[i] = res[i];
    return true;
  }

  virtual bool eval_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number *g ) {
    // Convert
    VectorType v ( n ); for ( int i = 0; i < n; ++i ) v[i] = x[i];
    
    VectorType constraintVec ( _numConstraints );
    _constraintOps.evaluate( v, constraintVec );

    for( int j=0; j<_numConstraints; ++j ) g[j] = constraintVec[j];
    
    return true;
  }

  virtual bool eval_jac_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                            Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values ) {
      
    if (values == NULL) {
        
        VectorType tmpVec ( n ); tmpVec.setZero();
        std::vector<TripletType> tripletListJac;
        _constraintOps.evaluateJacobian( tmpVec, tripletListJac );
        for( int index = 0; index < tripletListJac.size(); ++index ){
            iRow[index] = tripletListJac[index].row();
            jCol[index] = tripletListJac[index].col();
        }
    }
    else {
        
        // return the values of the Jacobian of the constraints
        VectorType v ( n ); for ( int i = 0; i < n; ++i ) v[i] = x[i];

        std::vector<TripletType> tripletListJac;
        _constraintOps.evaluateJacobian( v, tripletListJac );
        for( int index = 0; index < tripletListJac.size(); ++index ){
            values[index] = tripletListJac[index].value();
        }
    }
    
    return true;
  }
  
  
    // Hessian of Lagrangian: L = sigma * D^2 E + lambda D^2 Constraint
    // since Hessian is symmetric, only fill lower left part (otherwise IPOPT does not understand the matrix correctly)
  virtual bool eval_h( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor,
                       Ipopt::Index m, const Ipopt::Number* lambda, bool new_lambda, 
                       Ipopt::Index nele_hess, Ipopt::Index* iRow,Ipopt::Index* jCol, Ipopt::Number* values){
    if (values == NULL) {
        
        Ipopt::Index idx=0;
        
        for( int index = 0; index < _HessianEnergyTripletList.size(); ++index ){
            iRow[idx] = _HessianEnergyTripletList[index].row();
            jCol[idx] = _HessianEnergyTripletList[index].col();
            idx++;
        }
        
        std::vector<std::vector<TripletType>> tripletListHessian ( _numConstraints );
        VectorType tmpVec ( n ); tmpVec.setZero();
        _constraintOps.evaluateHessian( tmpVec, tripletListHessian );
        for( int constraint = 0; constraint < _numConstraints; ++ constraint ){
            for( int index = 0; index < tripletListHessian[constraint].size(); ++index ){
              iRow[idx] = tripletListHessian[constraint][index].row();
              jCol[idx] = tripletListHessian[constraint][index].col();
              idx++;
            }
        }
        
        // return the structure.
//         Ipopt::Index idx=0;
//         for (int k=0; k<_StructureOfHessian.outerSize(); ++k)
//             for ( typename SparseMatrixType::InnerIterator it(_StructureOfHessian,k); it; ++it){
//                 iRow[idx] = it.row(); 
//                 jCol[idx] = it.col();
//                 idx++;
//             }
//         assert(idx == nele_hess);
    }
    else {
        // Convert
        VectorType v ( n ); for ( int i = 0; i < n; ++i ) v[i] = x[i];
        
        Ipopt::Index idx=0;
        
        for( int index = 0; index < _HessianEnergyTripletList.size(); ++index ){
            values[idx] = obj_factor * _HessianEnergyTripletList[index].value();
            idx++;
        }
        
        std::vector<std::vector<TripletType>> tripletListHessian ( _numConstraints );
        _constraintOps.evaluateHessian( v, tripletListHessian );
        for( int constraint = 0; constraint < _numConstraints; ++ constraint ){
            for( int index = 0; index < tripletListHessian[constraint].size(); ++index ){
              values[idx] = lambda[constraint] * tripletListHessian[constraint][index].value();
              idx++;
            }
        }

//         Ipopt::Index idx=0;
//         for (int k=0; k<Hessian.outerSize(); ++k)
//             for ( typename SparseMatrixType::InnerIterator it(Hessian,k); it; ++it){
//                 values[idx] = it.value();
//                 idx++;
//             }
//         assert(idx == nele_hess);
    }

    return true;
 }
  
  
  virtual bool intermediate_callback ( Ipopt::AlgorithmMode /*mode*/, Ipopt::Index iter, 
                                       Ipopt::Number /*obj_value*/,  Ipopt::Number /*inf_pr*/,  Ipopt::Number /*inf_du*/, Ipopt::Number /*mu*/,  Ipopt::Number /*d_norm*/, Ipopt::Number /*regularization_size*/, Ipopt::Number /*d_du*/, Ipopt::Number /*d_pr*/,  Ipopt::Index /*ls_trials*/, 
                                       const Ipopt::IpoptData* ip_data,
                                       Ipopt::IpoptCalculatedQuantities* ip_cq ) {

    // Retrieve primal variables (code similar to code from the IPOPT documentation):
    Ipopt::TNLPAdapter* tnlpAdapter = NULL;

    if ( ip_cq != 0 ) {
      Ipopt::OrigIpoptNLP* origNLP;
      origNLP = dynamic_cast< Ipopt::OrigIpoptNLP*> ( Ipopt::GetRawPtr ( ip_cq->GetIpoptNLP () ) );
      // If in restoration mode, origNLP will be NULL. Quit method in this case
      if ( origNLP == NULL )  return true;
      tnlpAdapter = dynamic_cast < Ipopt::TNLPAdapter* > ( GetRawPtr ( origNLP->nlp () ) );
    }
    else throw std::invalid_argument ( pesopt::strprintf ( "ip_cq == NULL in FILE %s at line %d.", __FILE__, __LINE__ ).c_str() );
    double *x = new double[_startingPoint.size ()];   // n == _numDofs
    tnlpAdapter->ResortX ( *ip_data->curr ()->x (), x );
    delete[] x;

    // Plot solution.
     _iter = iter;
    
    return true;
  }

  virtual void finalize_solution ( Ipopt::SolverReturn /*status*/, Ipopt::Index /*n*/,  const Ipopt::Number* x,
                                   const Ipopt::Number* /*z_L*/,  const Ipopt::Number* /*z_U*/, 
                                   Ipopt::Index /*m*/,  const Ipopt::Number* /*g*/, const Ipopt::Number* lambda, 
                                   Ipopt::Number /*obj_value*/,  const Ipopt::IpoptData* /*ip_data*/,  Ipopt::IpoptCalculatedQuantities* ip_cq ) {
    for ( int i = 0; i < _nextIter.size (); ++i ) _nextIter[i] = x[i];
    for ( int i = 0; i < _nextLambda.size (); ++i ) _nextLambda[i] = lambda[i];
    _nlpError = ip_cq->curr_nlp_error();
  }
  
  public : 
    RealType getNLPError( ) const { return _nlpError; }
    int getNumIterations() const { return _iter; }
  
};




template <typename DataTypeContainer >
class IpoptConstraintSecondOrderSolver_QuadraticEnergy {
protected:
 
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
  typedef typename DataTypeContainer::TripletType TripletType;  
  
  const pesopt::QuadraticEnergyOp<DataTypeContainer> &_energyOp;  
  const std::vector<TripletType> &_HessianEnergyTripletList;
  const pesopt::SparseJacobianNonlinearConstraintOps<DataTypeContainer> &_constraintOp;  
  const RealType _ipoptTol;
  const int _MaxIterations;
  const RealType _x_l, _x_u;
  const int _linearSolverTypeIpopt;
  const RealType _boundRelaxFactorIpopt;
  const int _ipoptPrintLevel;
  const string _ipoptOutputFile;

public:
  IpoptConstraintSecondOrderSolver_QuadraticEnergy ( 
    const pesopt::QuadraticEnergyOp<DataTypeContainer> &energyOp,
    const std::vector<TripletType> &HessianEnergyTripletList,
    const pesopt::SparseJacobianNonlinearConstraintOps<DataTypeContainer> &constraintOps,
    const int MaxIterations, const RealType Tolerance,
    const RealType x_l = - 2.e+19, const RealType x_u = 2.e+19,
    const int linearSolverTypeIpopt = 0,
    const RealType boundRelaxFactorIpopt = 1.e-8, 
    const int ipoptPrintLevel = 5,
    const string ipoptOutputFile = ""
  )
  : _energyOp ( energyOp ), _HessianEnergyTripletList( HessianEnergyTripletList ), _constraintOp ( constraintOps ),
    _ipoptTol ( Tolerance ),
    _MaxIterations ( MaxIterations ),
    _x_l( x_l), _x_u( x_u ),
    _linearSolverTypeIpopt ( linearSolverTypeIpopt ),
    _boundRelaxFactorIpopt ( boundRelaxFactorIpopt ), 
    _ipoptPrintLevel ( ipoptPrintLevel ), _ipoptOutputFile ( ipoptOutputFile ) { }

  virtual ~IpoptConstraintSecondOrderSolver_QuadraticEnergy() {}

  void solve ( const VectorType &StartPosition, VectorType &Solution, 
               const VectorType &StartLambda, VectorType &SolLambda, 
               SolverInfo<DataTypeContainer> &solverInfo ) const {

    // Set up masking operators:
    Ipopt::SmartPtr<IpoptConstraintSecondOrderSolverBase_QuadraticEnergy<DataTypeContainer> > tsOpt 
    = new IpoptConstraintSecondOrderSolverBase_QuadraticEnergy<DataTypeContainer> ( _energyOp, _HessianEnergyTripletList, _constraintOp, StartPosition, Solution, StartLambda, SolLambda, _x_l, _x_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp = new Ipopt::IpoptApplication ();
    
    setToleranceOptionsIpopt<RealType>( _ipoptTol, _MaxIterations, ipoptApp );
    switchLinearSolverTypeIpopt ( _linearSolverTypeIpopt, ipoptApp );
    setPrintOptionsIpopt ( _ipoptPrintLevel, _ipoptOutputFile, ipoptApp );
    
    ipoptApp->Options()->SetStringValue ( "hessian_approximation", "exact" );
    
    ipoptApp->Options()->SetNumericValue ( "bound_relax_factor", _boundRelaxFactorIpopt );

    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();

    //! Optimization step
    Ipopt::SmartPtr<Ipopt::TNLP> nlp;
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
    outputIpoptStatus ( ipoptStatus, true );
    
    solverInfo.setSolverStatus( getIpoptStatus( ipoptStatus ).c_str() );
    solverInfo.setError( tsOpt->getNLPError() );
    solverInfo.setNumIterations( tsOpt->getNumIterations() );
    
  }
  
};


#endif //PESOPT_WITH_IPOPT


#endif //__IPOPTNONLINEARCONSTRAINTSOLVER_H
