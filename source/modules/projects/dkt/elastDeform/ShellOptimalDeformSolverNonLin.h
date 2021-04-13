#ifndef __OPTIMALDEFORMSOLVERNONLINEAR_H
#define __OPTIMALDEFORMSOLVERNONLINEAR_H


#include <pesopt_Ipopt.h>

#include <linearSystemSolver.h>
#include <lineSearchMethods.h>
#include <derivativeTest.h>

#include "ShellForces.h"
#include "ShellDeformationEnergiesNonLin.h"
#include "ShellOptimalDeformSolverInterface.h"


#define __DEBUGOPTDEFORMSOLVERNL


template <typename MatOptConfigurator, typename NonLinElastEnergyType>
class optimalDeformSolver<MatOptConfigurator, false, NonLinElastEnergyType, false> 
: public optimalDeformSolverInterface<MatOptConfigurator> {

  typedef typename MatOptConfigurator::ConfiguratorType         ConfiguratorType;
  typedef typename MatOptConfigurator::ConfiguratorTypePf       ConfiguratorTypePf;
  typedef typename ConfiguratorType::DTContainer                DataTypeContainer;
  typedef typename ConfiguratorType::RealType                   RealType;
  typedef typename ConfiguratorType::MaskType                   MaskType;
  typedef typename ConfiguratorType::InitType                   MeshType;
  typedef typename ConfiguratorType::VectorType                 VectorType;
  typedef typename ConfiguratorType::SparseMatrixType           SparseMatrixType;
  typedef pesopt::BoostParser ParameterParserType;
  typedef typename DataTypeContainer::TripletType               TripletType;
  
protected :

  const ShellHandler<ConfiguratorType,NonLinElastEnergyType::_DiscreteFunctionCacheType> &_shellHandler;
  mutable VectorType _initDisp;
  mutable SparseMatrixType _SystemMatAdjoint;
  
public:
      
      optimalDeformSolver( const ParameterParserType &parser,
                           const MatOptConfigurator & matOpConf,
                           const VectorType & pf,
                           const ShellHandler<ConfiguratorType,NonLinElastEnergyType::_DiscreteFunctionCacheType> & shellHandler,
                           const VectorType & rhs_force,
                           const VectorType & initDisp  ) : 
         optimalDeformSolverInterface<MatOptConfigurator> ( parser, matOpConf, pf, shellHandler.getDirichletMask(), rhs_force ),
         _shellHandler ( shellHandler ),
         _initDisp ( initDisp ), 
          _SystemMatAdjoint( 3*matOpConf._conf.getNumGlobalDofs(), 3*matOpConf._conf.getNumGlobalDofs() )
        {
          solve();
        }

    const ShellHandler<ConfiguratorType,NonLinElastEnergyType::_DiscreteFunctionCacheType> & getShellHandler ( ) const { return _shellHandler; }

    
    RealType getLastResidual( ) const{
      VectorType grad ( this->_solDisp.size() );
      NonLinElastEnergyType ( this->_matOpConf, _shellHandler.getDirichletMask(), _shellHandler.getChartToUndeformedShell_Cache(), this->_pf, this->_rhs ).evaluateJacobian( this->_solDisp, grad );
      return grad.norm();
    }
    
    void getLastStressOnElements( VectorType &membraneStressVec, VectorType &bendingStressVec, VectorType &totalStressVec ) const{
      NonLinElastEnergyType ( this->_matOpConf, _shellHandler.getDirichletMask(), _shellHandler.getChartToUndeformedShell_Cache(), this->_pf, this->_rhs ).evaluateStressOnElements( this->_solDisp, membraneStressVec, bendingStressVec, totalStressVec );
    }
    
   void getEnergyInfo(  DeformationOptimizationShellEnergyInfo<RealType> & energyInfo ) const {
        this->template getPartialEnergyInfo<false> ( energyInfo );
        energyInfo.setResidualConstraint( this->getLastResidual() );
   }
    
    void updatePhasefield( const VectorType &pf ) const {
      const RealType diff = ( this->_pf - pf ).squaredNorm();
      if ( diff > 1.e-15 ){ this->_pf = pf; solve(); }
    }
    
    const SparseMatrixType & getSystemMatForAdjointProblem () const { return _SystemMatAdjoint; }

protected:
 
  void solve() const { 
      
    const string nonlinearSolver = this->_parser.template get<string> ( "ConstraintProblem.nonlinearSolver" );
    if( nonlinearSolver == "NEWTON" ) this->solveWithNewton();
    if( nonlinearSolver == "IPOPT" ) this->solveWithIPOPT();
    if( nonlinearSolver == "FirstTryNEWTONelseIPOPT" ){
        const bool NewtonSuccess = this->solveWithNewton();
        if( !NewtonSuccess ) const bool IPOPTSuccess = this->solveWithIPOPT();
    }

    RealType energyTotal = 0.;
    NonLinElastEnergyType ( this->_matOpConf, _shellHandler.getDirichletMask(), _shellHandler.getChartToUndeformedShell_Cache(), this->_pf, this->_rhs ).evaluateElasticEnergies( this->_solDisp, this->_energy_Membrane, this->_energy_Bending, this->_energy_Pot, energyTotal ); 
    
    //set initial values to solution (for next computation)
    if( this->_parser.template get<bool>("ConstraintProblem.updateInitialization") ) _initDisp = this->_solDisp;
  }
 
  bool solveWithNewton() const {    
#ifdef __DEBUGOPTDEFORMSOLVERNL
   cout << endl 
        << "===========================================" << endl
        << "start newton solver for constraint" << endl; 
#endif
      
     NonLinElastEnergyType NonLinearEnergyOp( this->_matOpConf, _shellHandler.getDirichletMask(), _shellHandler.getChartToUndeformedShell_Cache(), this->_pf, this->_rhs );
    
//BEGIN DerivativeTest
    if( this->_parser.template get<int>( "DerivativeTestConstraint.order" ) > 0 ){
        pesopt::consoleOutput( "test Derivative for opt deform nonlin" );
        VectorType testPoint = VectorType::Random( this->_solDisp.size() );
        for( int i = 0; i < this->_numGlobalDofsDisp; ++i ){
            if ( this->_mask[i] ){
                for( int comp=0; comp<3; ++comp ) testPoint[i + comp * this->_numGlobalDofsDisp] = 0.0;
            }
        } 
        const int outputLevelDerivativeTest = 0;
        DerivativeTester<DataTypeContainer> derivTester ( NonLinearEnergyOp, this->_parser.template get<RealType>( "DerivativeTestConstraint.stepSize" ), testPoint, outputLevelDerivativeTest, this->_parser.template get<RealType>( "DerivativeTestConstraint.tolerance" ) );
        if( this->_parser.template get<int>( "DerivativeTestConstraint.order" ) == 1 ) derivTester.testFirstDerivative_AllDirections( );
        if( this->_parser.template get<int>( "DerivativeTestConstraint.order" ) == 2 ) derivTester.testSecondDerivative_AllDirections( );
        if( this->_parser.template get<int>( "DerivativeTestConstraint.order" ) == 12 ) {
            derivTester.testFirstDerivative_AllDirections( );
            derivTester.testSecondDerivative_AllDirections( );
            derivTester.testSecondDerivative_AllDirections( this->_parser.template get<string>("saving.saveDirectory"), "testHessian" );
        }
        return true;
    }
//END DerivativeTest

    VectorType startPoint ( 3 * this->_numGlobalDofsDisp ); startPoint = _initDisp;
    
    const int maxIterations       = this->_parser.template get<int> ( "ConstraintProblem.maxIterations" );
    const RealType breakCondition = this->_parser.template get<RealType> ( "ConstraintProblem.breakCondition" );
     
    //const int linearSolverType    = this->_parser.template get<int> ( "ConstraintProblem.linearSolver" );
# ifdef PESOPT_WITH_SUITESPARSE
    DirectLinearSystemSolver<DataTypeContainer,PESOPTUmfPackLU> directLinearSystemSolver( "opt disp", this->_parser.template get<string> ("saving.saveDirectory" ).c_str() );
#   else 
    DirectLinearSystemSolver<DataTypeContainer,EigenSparseLU> directLinearSystemSolver( "opt disp", this->_parser.template get<string> ("saving.saveDirectory" ).c_str() );
#   endif 
    
    
    NewtonMethod< DataTypeContainer, NewtonOptimalStepsizeControl > NewtonSolver ( NonLinearEnergyOp, directLinearSystemSolver, 3 * this->_numGlobalDofsDisp, maxIterations, breakCondition, this->_outputLevel );
    //TODO setAcceptTolerance???
    NewtonSolver.setAcceptTolerance( breakCondition * 1000. * startPoint.size() );
    const bool solvedToAcceptTol = NewtonSolver.solve( startPoint, this->_solDisp );
    NonLinearEnergyOp.evaluateHessian( this->_solDisp, _SystemMatAdjoint );
     
     
#ifdef __DEBUGOPTDEFORMSOLVERNL 
   cout << "finished newton solver for constraint" << endl
        << "===========================================" << endl;
#endif
     
     return solvedToAcceptTol;
  }


  bool solveWithIPOPT() const {    
      
#ifdef __DEBUGOPTDEFORMSOLVERNL
   cout << endl << "===========================================" << endl  << "start ipopt solver for constraint" << endl; 
   auto startTime_Solve = std::chrono::high_resolution_clock::now();
#endif
      
    //For IPOPT do not use boundary conditions in matrix, instead set them as constraint
    MaskType zeroMask ( this->_mask.size() ); for( int i=0; i<zeroMask.size(); ++i ) zeroMask[i] = false;
    NonLinElastEnergyType NonLinearEnergyOp( this->_matOpConf, zeroMask, _shellHandler.getChartToUndeformedShell_Cache(), this->_pf, this->_rhs );
    
    VectorType startPoint ( 3 * this->_numGlobalDofsDisp ); startPoint = _initDisp;

    const int linearSolverTypeIpopt    = this->_parser.template get<int> ( "ConstraintProblem.linearSolverTypeIpopt" );
    const int maxIterationsIpopt       = this->_parser.template get<int> ( "ConstraintProblem.maxIterationsIpopt" );
    const RealType breakConditionIpopt = this->_parser.template get<RealType> ( "ConstraintProblem.breakConditionIpopt" );
    const int outputLevelIpopt         = this->_parser.template get<int> ( "ConstraintProblem.outputLevelIpopt" );

    VectorType x_l_vec ( startPoint.size() ), x_u_vec ( startPoint.size() );
    for( int globalDof=0; globalDof<this->_numGlobalDofsDisp; ++globalDof ){
        if( _shellHandler.getDirichletMask()[globalDof] ){
            for(int comp=0; comp<3; ++comp ){
                x_l_vec[globalDof + comp * this->_numGlobalDofsDisp]  = 0.;
                x_u_vec[globalDof + comp * this->_numGlobalDofsDisp]  = 0.;
            }
        }else{
            for(int comp=0; comp<3; ++comp ){
                x_l_vec[globalDof + comp * this->_numGlobalDofsDisp]  = - 2.e+19;
                x_u_vec[globalDof + comp * this->_numGlobalDofsDisp]  =   2.e+19;
            }  
        }
    }
    
    IpoptBoxConstraintSecondOrderSolver<DataTypeContainer> ipoptSolver ( NonLinearEnergyOp, maxIterationsIpopt, breakConditionIpopt, x_l_vec , x_u_vec, linearSolverTypeIpopt, outputLevelIpopt );
    
//BEGIN DerivativeTest
    if( this->_parser.template get<int>( "DerivativeTestConstraint.order" ) > 0 ){
        pesopt::consoleOutput( "test Derivative for opt deform nonlin" );
        VectorType testPoint = VectorType::Random( this->_solDisp.size() );
        for( int i = 0; i < this->_numGlobalDofsDisp; ++i ){
            if ( this->_mask[i] ){
                for( int comp=0; comp<3; ++comp ) testPoint[i + comp * this->_numGlobalDofsDisp] = 0.0;
            }
        } 
        if( this->_parser.template get<int>( "DerivativeTestConstraint.order" ) ==  1 ) ipoptSolver.testDerivative( testPoint, 1 );
        if( this->_parser.template get<int>( "DerivativeTestConstraint.order" ) ==  2 ) ipoptSolver.testDerivative( testPoint, 2 );
        if( this->_parser.template get<int>( "DerivativeTestConstraint.order" ) == 12 ) ipoptSolver.testDerivative( testPoint, 2 );

        return true;
    }
//END DerivativeTest
    

    SolverInfo<DataTypeContainer> solverInfo;
    ipoptSolver.solve( startPoint, this->_solDisp, solverInfo );
    
    //compute system matrix for adjoint problem
    NonLinElastEnergyType ( this->_matOpConf, _shellHandler.getDirichletMask(), _shellHandler.getChartToUndeformedShell_Cache(), this->_pf, this->_rhs ).evaluateHessian( this->_solDisp, _SystemMatAdjoint );
    
    bool solvedToAcceptTol = true;
    if( solverInfo.getError() > breakConditionIpopt ) solvedToAcceptTol = false;
    if( (this->_outputLevel > 0) && !solvedToAcceptTol ){
     cout << "=======================================================================" << endl;
     cout << "WARNING: IPOPT failed to solve with tolerance = " << breakConditionIpopt << "." << endl;
     cout << "                                        error = " << solverInfo.getError() << "." << endl;
     cout << "=======================================================================" << endl;
    }
    
    
#ifdef __DEBUGOPTDEFORMSOLVERNL 
   //! print elapsed time
    std::chrono::duration<RealType> diff_Solve = std::chrono::high_resolution_clock::now() - startTime_Solve;
   cout << "finished ipopt solver for constraint in " << diff_Solve.count() << " sec"  << endl << "===========================================" << endl;
#endif
    
    return solvedToAcceptTol;
  }
  
  
};


#endif
