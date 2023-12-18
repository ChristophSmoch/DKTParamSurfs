#ifndef __IPOPTBOXCONSTRAINTSOLVER_H
#define __IPOPTBOXCONSTRAINTSOLVER_H

#include <includesIpopt.h>
#include <solverInfo.h>


#ifdef PESOPT_WITH_IPOPT




/* use IPOPT to solve a constraint minimization problem of type 
 *  minimize f(x) 
 * over all x s.t. 
 * x_l \leq x \leq x_u
 */
template<typename DataTypeContainer>
class IpoptBoxConstraintFirstOrderSolverBase : public Ipopt::TNLP {

protected:
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  
  const pesopt::NonlinearEnergyOp<DataTypeContainer> &_energyOp;

  const VectorType &_startingPoint;
  VectorType &_nextIter;
 
  unsigned int _numTotalDofs;
  unsigned int _numConstraints;
  
  RealType _nlpError;
  int _iter;

  const VectorType &_x_l, &_x_u;

public:
   
    IpoptBoxConstraintFirstOrderSolverBase ( const pesopt::NonlinearEnergyOp<DataTypeContainer> &energyOp,
                                                  const VectorType &startingPoint, VectorType &nextIter,
                                                  const VectorType & x_l, const VectorType & x_u  )
  : _energyOp ( energyOp ),
    _startingPoint ( startingPoint ), _nextIter ( nextIter ),
    _numTotalDofs ( startingPoint.size () ), _numConstraints ( 0 ),
    _x_l ( x_l ), _x_u ( x_u ) {   }
    
    
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
        x_l[i] = _x_l[i];
        x_u[i] = _x_u[i];
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
    _energyOp.evaluateEnergy( v, s );
    obj_value = s;
    return true;

  }

  virtual bool eval_grad_f ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f ) {

    // Convert
    VectorType v ( n );
    for ( int i = 0; i < n; ++i ) v[i] = x[i];

   VectorType res ( v.size() );
   _energyOp.evaluateJacobian ( v, res );
   for ( int i = 0; i < n; ++i )
     grad_f[i] = res[i];
   return true;

  }

  virtual bool eval_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m, Ipopt::Number *g ) {
    return true;
  }

  virtual bool eval_jac_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                            Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values ) {
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
                                   Ipopt::Number /*obj_value*/,  const Ipopt::IpoptData* /*ip_data*/,
                                   Ipopt::IpoptCalculatedQuantities* ip_cq ) {
    for ( int i = 0; i < _nextIter.size (); ++i ) _nextIter[i] = x[i];
    _nlpError = ip_cq->curr_nlp_error();
  }
  
public : 
    RealType getNLPError( ) const { return _nlpError; }
    int getNumIterations() const { return _iter; }
  
};




template <typename DataTypeContainer>
class IpoptBoxConstraintFirstOrderSolver {
protected:
 
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
    
  const pesopt::NonlinearEnergyOp<DataTypeContainer> &_energyOp;
  RealType _ipoptTol;
  RealType _MaxIterations;
  const int _linearSolverTypeIpopt;
  const RealType _boundRelaxFactorIpopt;
  const int _ipoptPrintLevel;
  const string _ipoptOutputFile;

  mutable VectorType _x_l, _x_u;
  
public:
  IpoptBoxConstraintFirstOrderSolver ( const pesopt::NonlinearEnergyOp<DataTypeContainer> &E,
                                       const int MaxIterations, const RealType Tolerance,
                                       const RealType x_l = - 2.e+19, const RealType x_u = 2.e+19,
                                       const int linearSolverTypeIpopt = 0,
                                       const RealType boundRelaxFactorIpopt = 1.e-8, 
                                       const int ipoptPrintLevel = 5,
                                       const string ipoptOutputFile = ""   )
  : _energyOp ( E ),
    _ipoptTol ( Tolerance ), _MaxIterations ( MaxIterations ),
    _linearSolverTypeIpopt ( linearSolverTypeIpopt ),
    _boundRelaxFactorIpopt ( boundRelaxFactorIpopt ), 
    _ipoptPrintLevel ( ipoptPrintLevel ), _ipoptOutputFile ( ipoptOutputFile ),
    _x_l ( E.getNumDofs() ), _x_u ( E.getNumDofs() ) { 
        for( int i=0; i<E.getNumDofs(); ++i ){
          _x_l[i] = x_l; _x_u[i] = x_u;   
        }
     }
    
  IpoptBoxConstraintFirstOrderSolver ( const pesopt::NonlinearEnergyOp<DataTypeContainer> &E,
                                       const int MaxIterations, const RealType Tolerance,
                                       const VectorType & x_l, const VectorType & x_u,
                                       const int linearSolverTypeIpopt,
                                       const RealType boundRelaxFactorIpopt = 1.e-8, 
                                       const int ipoptPrintLevel = 5,
                                       const string ipoptOutputFile = "" )
  : _energyOp ( E ),
    _ipoptTol ( Tolerance ), _MaxIterations ( MaxIterations ),
    _linearSolverTypeIpopt ( linearSolverTypeIpopt ),
    _boundRelaxFactorIpopt ( boundRelaxFactorIpopt ), 
    _ipoptPrintLevel ( ipoptPrintLevel ),
    _x_l( x_l ), _x_u ( x_u ) {}
  
  void testDerivative ( const VectorType &StartPosition, const int order ) const {

    // Set up masking operators:
    VectorType Solution ( StartPosition.size() );
    Ipopt::SmartPtr<IpoptBoxConstraintFirstOrderSolverBase<DataTypeContainer> > tsOpt 
    = new IpoptBoxConstraintFirstOrderSolverBase<DataTypeContainer> ( _energyOp, StartPosition, Solution, _x_l, _x_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp = new Ipopt::IpoptApplication ();
    
    //Derivative Test
    if( order == 1 ) ipoptApp->Options()->SetStringValue  ( "derivative_test", "first-order" );
    if( order == 2 ) throw std::invalid_argument( pesopt::strprintf ( "derivative of order 2 not implmented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    
    setDerivateTestOptionsIpopt ( ipoptApp ); 
    ipoptApp->Options()->SetStringValue ( "hessian_approximation", "exact" );    
    setPrintOptionsIpopt( _ipoptPrintLevel, _ipoptOutputFile, ipoptApp );

    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();
    
    Ipopt::SmartPtr<Ipopt::TNLP> nlp;
    
    //! Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
  }
  
  
  void solve ( const VectorType &StartPosition, VectorType &Solution, SolverInfo<DataTypeContainer> &solverInfo ) const {

    // Set up masking operators:
    Ipopt::SmartPtr<IpoptBoxConstraintFirstOrderSolverBase<DataTypeContainer> > tsOpt 
    = new IpoptBoxConstraintFirstOrderSolverBase<DataTypeContainer> ( _energyOp, StartPosition, Solution, _x_l, _x_u );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp = new Ipopt::IpoptApplication ();
    
    setToleranceOptionsIpopt<RealType>( _ipoptTol, _MaxIterations, ipoptApp );
    switchLinearSolverTypeIpopt ( _linearSolverTypeIpopt, ipoptApp );
    setPrintOptionsIpopt( _ipoptPrintLevel, _ipoptOutputFile, ipoptApp );
    
    // Enable quasi-Newton hessian approximation
    ipoptApp->Options ()->SetStringValue ( "hessian_approximation", "limited-memory" );

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
  
  void solve ( const VectorType &StartPosition, VectorType &Solution ) const {
    SolverInfo<DataTypeContainer> solverInfo;
    this->solve( StartPosition, Solution, solverInfo );   
  }

};







template<typename DataTypeContainer >
class IpoptBoxConstraintSecondOrderSolverBase : public Ipopt::TNLP {

protected:
    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType VectorType;
    typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
    typedef typename DataTypeContainer::MaskType MaskType;
    typedef typename DataTypeContainer::TripletType TripletType;
    
  const pesopt::NonlinearEnergyOp<DataTypeContainer> & _energyOp;

  const VectorType &_startingPoint;
  VectorType &_nextIter;
  unsigned int _numTotalDofs;
  const VectorType& _x_l_vec; const VectorType& _x_u_vec;
  RealType _nlpError;
  int _iter;


public:
    
  IpoptBoxConstraintSecondOrderSolverBase ( const pesopt::NonlinearEnergyOp<DataTypeContainer> &Energy,
                                                 const VectorType &startingPoint, VectorType &nextIter,
                                                 const VectorType &x_l_vec, const VectorType &x_u_vec )
  : _energyOp ( Energy ),
    _startingPoint ( startingPoint ), 
    _nextIter ( nextIter ),
    _numTotalDofs ( startingPoint.size () ),
    _x_l_vec ( x_l_vec ), _x_u_vec ( x_u_vec ) { }  
  

  //returns info about the nlp
  virtual bool get_nlp_info ( Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style ) {
    // The number of DOFs.
    n = _numTotalDofs;
    // The number of constraints.
    m = 0;
    // The number of nonzeros in the jacobian: Possibly every entry. 
    nnz_jac_g = 0;
    
    std::vector<TripletType> tripletListHessian;
    VectorType tmpVec ( n ); tmpVec.setZero();
    _energyOp.evaluateTripletListHessianSym( tmpVec, tripletListHessian );
    nnz_h_lag = tripletListHessian.size();
    
    //Use C index style
    index_style = C_STYLE;

    return true;
  }

  virtual bool get_bounds_info ( Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, 
                                 Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u ) {
   for ( Ipopt::Index i = 0; i < n; ++i ) {
      x_l[i] = _x_l_vec[i];
      x_u[i] = _x_u_vec[i];
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
    RealType energy = 0.;
    _energyOp.evaluateEnergy( v, energy );
    obj_value = energy;
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
    return true;
  }

  virtual bool eval_jac_g ( Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Index m,
                            Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values ) {
    return true;
  }
  
  // Hessian of Lagrangian: L = sigma * D^2 E + lambda D^2 Constraint
  virtual bool eval_h( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor,
                       Ipopt::Index m, const Ipopt::Number* lambda, bool new_lambda, 
                       Ipopt::Index nele_hess, Ipopt::Index* iRow,Ipopt::Index* jCol, Ipopt::Number* values){
      
    if (values == NULL) {
        // return the structure.
        Ipopt::Index idx=0;
        std::vector<TripletType> tripletListHessian;
        VectorType tmpVec ( n ); tmpVec.setZero();
        _energyOp.evaluateTripletListHessianSym( tmpVec, tripletListHessian );
        for( int index = 0; index < tripletListHessian.size(); ++index ){
              iRow[idx] = tripletListHessian[index].row();
              jCol[idx] = tripletListHessian[index].col();
              idx++;
        }
    }else {
        // return the values
        VectorType v ( n ); for ( int i = 0; i < n; ++i ) v[i] = x[i];
        Ipopt::Index idx=0;
        std::vector<TripletType> tripletListHessian;
        _energyOp.evaluateTripletListHessianSym( v, tripletListHessian );
        for( int index = 0; index < tripletListHessian.size(); ++index ){
              values[idx] = obj_factor * tripletListHessian[index].value();
              idx++;
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
class IpoptBoxConstraintSecondOrderSolver {
protected:
 
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
    
  const pesopt::NonlinearEnergyOp<DataTypeContainer> &_energyOp;
  const RealType _ipoptTol;
  const RealType _MaxIterations;
  mutable VectorType _x_l_vec, _x_u_vec;
  const int _linearSolverTypeIpopt;
  const RealType _boundRelaxFactorIpopt;
  const int _ipoptPrintLevel;
  const string _ipoptOutputFile;
  
public:
  IpoptBoxConstraintSecondOrderSolver ( const pesopt::NonlinearEnergyOp<DataTypeContainer> &energyOp,
                                        const int MaxIterations, const RealType Tolerance,
                                        const int numDofs, const RealType x_l = - 2.e+19, const RealType x_u = 2.e+19,
                                        const int linearSolverTypeIpopt = 0,
                                        const int ipoptPrintLevel = 5,
                                        const RealType boundRelaxFactorIpopt = 1.e-8, 
                                        const string ipoptOutputFile = "" )
  : _energyOp( energyOp ),
    _ipoptTol ( Tolerance ),
    _MaxIterations ( MaxIterations ),
    _x_l_vec( numDofs ), _x_u_vec( numDofs ),
    _linearSolverTypeIpopt ( linearSolverTypeIpopt ),
    _boundRelaxFactorIpopt ( boundRelaxFactorIpopt ), 
    _ipoptPrintLevel ( ipoptPrintLevel ), _ipoptOutputFile ( ipoptOutputFile ) {
        for( int i=0; i<numDofs; ++i ){
          _x_l_vec[i] = x_l; _x_u_vec[i] = x_u;
        }
    }

    
  IpoptBoxConstraintSecondOrderSolver ( const pesopt::NonlinearEnergyOp<DataTypeContainer> &energyOp,
                                        const int MaxIterations, const RealType Tolerance,
                                        const VectorType &x_l_vec, const VectorType &x_u_vec,
                                        const int linearSolverTypeIpopt,
                                        const RealType boundRelaxFactorIpopt = 1.e-8, 
                                        const int ipoptPrintLevel = 5,
                                        const string ipoptOutputFile = ""  )
  : _energyOp( energyOp ),
    _ipoptTol ( Tolerance ),
    _MaxIterations ( MaxIterations ),
    _x_l_vec( x_l_vec ), _x_u_vec( x_u_vec ),
    _linearSolverTypeIpopt ( linearSolverTypeIpopt ),
    _boundRelaxFactorIpopt ( boundRelaxFactorIpopt ), 
    _ipoptPrintLevel ( ipoptPrintLevel ), _ipoptOutputFile ( ipoptOutputFile ) { }
  
  void testDerivative ( const VectorType &StartPosition, const int order ) const {

    // Set up masking operators:
    VectorType Solution ( StartPosition.size() );
    Ipopt::SmartPtr<IpoptBoxConstraintSecondOrderSolverBase<DataTypeContainer> > tsOpt 
    = new IpoptBoxConstraintSecondOrderSolverBase<DataTypeContainer> ( _energyOp, StartPosition, Solution, _x_l_vec, _x_u_vec );

    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp = new Ipopt::IpoptApplication ();
    
    //Derivative Test
    if( order == 1 ) ipoptApp->Options()->SetStringValue  ( "derivative_test", "first-order" );
    if( order == 2 ) ipoptApp->Options()->SetStringValue  ( "derivative_test", "second-order" );
    setDerivateTestOptionsIpopt ( ipoptApp );
    ipoptApp->Options()->SetStringValue ( "hessian_approximation", "exact" );
    
    setPrintOptionsIpopt( _ipoptPrintLevel, _ipoptOutputFile, ipoptApp );

    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();
    
    Ipopt::SmartPtr<Ipopt::TNLP> nlp;
    
    //! Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
    
  }
  
  void solve ( const VectorType &StartPosition, VectorType &Solution, SolverInfo<DataTypeContainer> &solverInfo ) const {

    // Set up masking operators:
    Ipopt::SmartPtr<IpoptBoxConstraintSecondOrderSolverBase<DataTypeContainer> > tsOpt 
    = new IpoptBoxConstraintSecondOrderSolverBase<DataTypeContainer> ( _energyOp, StartPosition, Solution, _x_l_vec, _x_u_vec );
    
    // Create an instance of the IpoptApplication
    Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp = new Ipopt::IpoptApplication();
    
    setToleranceOptionsIpopt<RealType>( _ipoptTol, _MaxIterations, ipoptApp );
    switchLinearSolverTypeIpopt ( _linearSolverTypeIpopt, ipoptApp );
    setPrintOptionsIpopt( _ipoptPrintLevel, _ipoptOutputFile, ipoptApp );
    
    ipoptApp->Options()->SetStringValue ( "hessian_approximation", "exact" );

    ipoptApp->Options()->SetNumericValue ( "bound_relax_factor", _boundRelaxFactorIpopt );
    
    
    // Initialize the IpoptApplication
    Ipopt::ApplicationReturnStatus ipoptStatus;
    ipoptStatus = ipoptApp->Initialize();

    Ipopt::SmartPtr<Ipopt::TNLP> nlp;
    
    //! Optimization step
    ipoptStatus = ipoptApp->OptimizeTNLP( tsOpt );
    // true means, output only if solver failed
    outputIpoptStatus ( ipoptStatus, true );
    
    solverInfo.setSolverStatus( getIpoptStatus( ipoptStatus ).c_str() );
    solverInfo.setError( tsOpt->getNLPError() );
    solverInfo.setNumIterations( tsOpt->getNumIterations() );
  }
  
  void solve ( const VectorType &StartPosition, VectorType &Solution ) const {
    SolverInfo<DataTypeContainer> solverInfo;
    this->solve( StartPosition, Solution, solverInfo );   
  }

};


#endif //PESOPT_WITH_IPOPT

#endif //__IPOPTBOXCONSTRAINTSOLVER_H
