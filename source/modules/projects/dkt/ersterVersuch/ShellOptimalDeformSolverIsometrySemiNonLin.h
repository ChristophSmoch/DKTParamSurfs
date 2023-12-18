#ifndef __OPTIMALDEFORMSOLVERISOMETRYSEMINONLIN_H
#define __OPTIMALDEFORMSOLVERISOMETRYSEMINONLIN_H

#include <ipoptNonlinearConstraintSolver.h>

#include <linearSystemSolver.h>
#include <lineSearchMethods.h>
#include <derivativeTest.h>

#include "../elastDeform/ShellOptimalDeformSolverInterface.h"
#include "../elastDeform/ShellForces.h"
#include "ShellDeformationEnergiesSemiNonLin.h"
#include "../elastDeform/ShellIsometryConstraint.h"


template <typename MatOptConfigurator, typename NonLinElastEnergyType>
class optimalDeformSolver<MatOptConfigurator, false, NonLinElastEnergyType, true> 
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
  
  const ShellHandler<ConfiguratorType,NonLinElastEnergyType::_DiscreteFunctionCacheType> &_shellHandler;
  const ShellIsometryHandler_Pointwise<ConfiguratorType> _isoOp;
  
protected :

  const int _numBoundaryNodes, _numInteriorNodes;
  const RealType _factorIsoConstraint;
  mutable VectorType _initDisp, _initMultiplierIso;
  mutable VectorType _solMultiplierIso;
  mutable RealType _energy_Iso;
  
public:
      
      optimalDeformSolver( const ParameterParserType &parser,
                           const MatOptConfigurator & matOpConf,
                           const VectorType & pf,
                           const ShellHandler<ConfiguratorType,NonLinElastEnergyType::_DiscreteFunctionCacheType> & shellHandler,
                           const VectorType & rhs_force,
                           const VectorType & initDisp  ) : 
         optimalDeformSolverInterface<MatOptConfigurator> ( parser, matOpConf, pf, shellHandler.getDirichletMask(), rhs_force ),
         _shellHandler( shellHandler ),
         _isoOp ( matOpConf._conf, shellHandler.getChartToUndeformedShell(), shellHandler.getDirichletMask(), shellHandler.getNumBoundaryNodes() ),
         _numBoundaryNodes ( shellHandler.getNumBoundaryNodes() ), _numInteriorNodes ( this->_numVertices - _numBoundaryNodes ),
         _factorIsoConstraint ( parser.template get<RealType> ( "ConstraintProblem.factorIsoConstraint" ) / static_cast<RealType> ( this->_numVertices ) ),
         _initDisp ( initDisp ), _initMultiplierIso( 3 * (this->_numVertices - _numBoundaryNodes ) ),
         _solMultiplierIso( 3 * (this->_numVertices - _numBoundaryNodes) )
        {
          for( int i=0; i<_initMultiplierIso.size(); ++i ) _initMultiplierIso[i] = 1.0;  
          solve();
        }
        
    const int getNumInteriorNodes( ) const { return _numInteriorNodes; };
    const ShellIsometryHandler_Pointwise<ConfiguratorType> & getIsometryOp ( ) const { return _isoOp; };
    
    void updatePhasefield( const VectorType &pf ) const {
      const RealType diff = ( this->_pf - pf ).squaredNorm();
      if ( diff > 1.e-15 ){ this->_pf = pf; solve(); }
    }
    

    const SparseMatrixType getSystemMatForAdjointProblem () const { 
        std::vector<TripletType> tripletList ( this->getHessianLinElastTriplet() );
        tripletList.reserve( tripletList.size() + (4+3*8) * this->_numGlobalDofsDisp );
//         this->getIsometryOp().assembleAddTripletListHessian( this->_solDisp, tripletList, _solMultiplierIso, this->_factorIsoConstraint );
        SparseMatrixType systemMat( 3 * this->_numGlobalDofsDisp + 3 * _numInteriorNodes, 3 * this->_numGlobalDofsDisp + 3 * _numInteriorNodes );
        systemMat.setFromTriplets( tripletList.begin(), tripletList.end() ); 
        return systemMat; 
    }
    
    
    const ShellHandler<ConfiguratorType,NonLinElastEnergyType::_DiscreteFunctionCacheType> & getShellHandler ( ) const { return _shellHandler; }

    const VectorType& getSolutionMultiplierIso() const { return _solMultiplierIso;}
    

    RealType getLastResidual( ) const{
//       QuadraticEnergyOp<DataTypeContainer> ElasticEnergy ( _HessianLinElast, this->_rhs );
//       LagrangianIsometryPointwise<ConfiguratorType> Lagrangian ( ElasticEnergy, _HessianLinElastTriplet, _isoOp, _factorIsoConstraint );
//      
//       VectorType sol ( 3 * _numGlobalDofsDisp + 3 * _numInteriorNodes );
//       Eigen::Ref<VectorType> sol_Disp = sol.segment( 0, 3 * _numGlobalDofsDisp ); sol_Disp = this->_solDisp;
//       Eigen::Ref<VectorType> sol_Mult = sol.segment( 3 * _numGlobalDofsDisp, 3 * _numInteriorNodes ); sol_Mult = _solMultiplierIso;
//      
//       VectorType grad ( sol.size() );
//       Lagrangian.evaluateJacobian( sol, grad );
//       return grad.norm();
        return 0.;
    }
    
    
   void getLastStressOnElements( VectorType &membraneStressVec, VectorType &bendingStressVec, VectorType &totalStressVec ) const{
      NonLinElastEnergyType ( this->_matOpConf, _shellHandler.getDirichletMask(), _shellHandler.getChartToUndeformedShell_Cache(), this->_pf, this->_rhs ).evaluateStressOnElements( this->_solDisp, membraneStressVec, bendingStressVec, totalStressVec );
    }
    
    
   void getEnergyInfo(  DeformationOptimizationShellEnergyInfo<RealType> & energyInfo ) const {
        this->template getPartialEnergyInfo<true> ( energyInfo );
        energyInfo.setResidualConstraint( this->getLastResidual() );
   }

protected:
 
  void solve() const { 
    const string nonlinearSolver = this->_parser.template get<string> ( "ConstraintProblem.nonlinearSolver" );
//     if( nonlinearSolver == "NEWTON" ) this->solveWithNewton();
    if( nonlinearSolver == "IPOPT" ) this->solveWithIPOPT();
//     if( nonlinearSolver == "FirstTryNEWTONelseIPOPT" ){
//         const bool NewtonSuccess = this->solveWithNewton();
//         if( !NewtonSuccess ) const bool IPOPTSuccess = this->solveWithIPOPT();
//     }
    
    
    RealType energyTotal = 0.;
    NonLinElastEnergyType ( this->_matOpConf, _shellHandler.getDirichletMask(), _shellHandler.getChartToUndeformedShell_Cache(), this->_pf, this->_rhs ).evaluateElasticEnergies( this->_solDisp, this->_energy_Membrane, this->_energy_Bending, this->_energy_Pot, energyTotal ); 
    
    //set initial values to solution (for next computation)
    if( this->_parser.template get<bool>("ConstraintProblem.updateInitialization") ){ 
      _initDisp = this->_solDisp; _initMultiplierIso = _solMultiplierIso;
    }
  }
  
//   bool solveWithNewton() const {     
//     
//     QuadraticEnergyOp<DataTypeContainer> ElasticEnergy ( _HessianLinElast, this->_rhs );
//     LagrangianIsometryPointwise<ConfiguratorType> Lagrangian ( ElasticEnergy, _HessianLinElastTriplet, _isoOp, _factorIsoConstraint );
//     
// // //BEGIN DerivativeTest
// //     if( this->_parser.template get<int>( "DerivativeTest.order" ) > 0 ){
// //         //     VectorType testPoint = VectorType::Random( _solDisp.size() + 3 * _numInteriorNodes );
// //         VectorType testPoint ( 3 * _numGlobalDofsDisp + 3 * _numInteriorNodes  );
// //         Eigen::Ref<VectorType> testPoint_Disp = testPoint.segment( 0, 3 * _numGlobalDofsDisp );
// //         Eigen::Ref<VectorType> testPoint_Mult = testPoint.segment( 3 * _numGlobalDofsDisp, 3 * _numInteriorNodes );
// //         testPoint_Disp = VectorType::Random( 3 * _numGlobalDofsDisp );
// //         for( int i = 0; i < _numGlobalDofsDisp; ++i ){
// //             if ( this->_mask[i] ){
// //                 for( int comp=0; comp<3; ++comp ) testPoint_Disp[i + comp * _numGlobalDofsDisp] = 0.0;
// //             }
// //         } 
// //         testPoint_Mult = VectorType::Random( 3 * _numInteriorNodes );
// //         
// //         DerivativeTester<DataTypeContainer> derivTester ( Lagrangian, this->_parser.template get<RealType>( "DerivativeTest.stepSize" ), testPoint, 0, this->_parser.template get<RealType>( "DerivativeTest.tolerance" ) );
// //         if( this->_parser.template get<int>( "DerivativeTest.order" ) == 1 ) derivTester.testFirstDerivative_AllDirections( );
// //         if( this->_parser.template get<int>( "DerivativeTest.order" ) == 2 ) derivTester.testSecondDerivative_AllDirections( );
// //         return true;
// //     }
// // //END DerivativeTest
// 
//     VectorType startPoint ( 3 * _numGlobalDofsDisp + 3 * _numInteriorNodes ), solPoint ( 3 * _numGlobalDofsDisp + 3 * _numInteriorNodes  );
//     Eigen::Ref<VectorType> startPoint_Disp = startPoint.segment( 0, 3 * _numGlobalDofsDisp ); startPoint_Disp = _initDisp;
//     Eigen::Ref<VectorType> startPoint_Mult = startPoint.segment( 3 * _numGlobalDofsDisp, 3 * _numInteriorNodes ); startPoint_Mult = _initMultiplierIso;
//      
//      // Newton: use LU-Solver since Matrix is not spd
        //DirectLinearSystemSolver<DataTypeContainer,PESOPTUmfPackLU> directLinearSystemSolver( "opt disp", this->_parser.template get<string> ("saving.saveDirectory" ).c_str() );
//      NewtonMethod< DataTypeContainer, NewtonOptimalStepsizeControl > NewtonSolver ( Lagrangian, directLinearSystemSolver, 3 * _numGlobalDofsDisp + 3 * _numInteriorNodes, this->_parser.template get<int>("ConstraintProblem.maxIterations"), this->_parser.template get<RealType>("ConstraintProblem.breakCondition"), this->_outputLevel );
//      const bool solvedToAcceptTol = NewtonSolver.solve( startPoint, solPoint );
//      this->_solDisp = solPoint.segment( 0, 3 * _numGlobalDofsDisp );
//      _solMultiplierIso = solPoint.segment( 3 * _numGlobalDofsDisp, 3 * _numInteriorNodes );  
//      
//      return solvedToAcceptTol;
//   }
  
  bool solveWithIPOPT() const {    
      
    //For IPOPT do not use boundary conditions in matrix, instead set them as constraint
    MaskType zeroMask ( this->_mask.size() ); for( int i=0; i<zeroMask.size(); ++i ) zeroMask[i] = false;
    NonLinElastEnergyType NonLinearEnergyOp( this->_matOpConf, zeroMask, _shellHandler.getChartToUndeformedShell_Cache(), this->_pf, this->_rhs );
    ShellIsometryConstraint<ConfiguratorType> isoOp ( this->_conf, this->_shellHandler.getChartToUndeformedShell(), this->_shellHandler.getDirichletMask(), this->_shellHandler.getNumBoundaryNodes(), _factorIsoConstraint );
    
    VectorType startPoint ( 3 * this->_numGlobalDofsDisp + 3 * _numInteriorNodes ), solPoint ( 3 * this->_numGlobalDofsDisp + 3 * _numInteriorNodes  );
    Eigen::Ref<VectorType> startPoint_Disp = startPoint.segment( 0, 3 * this->_numGlobalDofsDisp ); startPoint_Disp = _initDisp;
    Eigen::Ref<VectorType> startPoint_Mult = startPoint.segment( 3 * this->_numGlobalDofsDisp, 3 * _numInteriorNodes ); startPoint_Mult = _initMultiplierIso;
    
    const RealType boundRelaxFactorIpopt = this->_parser.template get<RealType> ("ConstraintProblem.boundRelaxFactorIpopt" );
    const int linearSolverTypeIpopt    = this->_parser.template get<int> ( "ConstraintProblem.linearSolverTypeIpopt" );
    const int maxIterationsIpopt       = this->_parser.template get<int> ( "ConstraintProblem.maxIterationsIpopt" );
    const RealType breakConditionIpopt = this->_parser.template get<RealType> ( "ConstraintProblem.breakConditionIpopt" );
    const int outputLevelIpopt         = this->_parser.template get<int> ( "ConstraintProblem.outputLevelIpopt" );
    
    VectorType x_l_vec ( startPoint_Disp.size() ), x_u_vec ( startPoint_Disp.size() );
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
    
    IpoptConstraintSecondOrderSolver<DataTypeContainer> 
    ipoptSolver( NonLinearEnergyOp, isoOp, maxIterationsIpopt, breakConditionIpopt, x_l_vec , x_u_vec, linearSolverTypeIpopt, boundRelaxFactorIpopt, outputLevelIpopt );
    
//     //BEGIN DerivativeTest
//     if( this->_parser.template get<int>( "DerivativeTestConstraint.order" ) > 0 ){
//         pesopt::consoleOutput( "test Derivative for opt deform nonlin" );
//         VectorType testPoint = VectorType::Random( this->_solDisp.size() );
//         for( int i = 0; i < this->_numGlobalDofsDisp; ++i ){
//             if ( this->_mask[i] ){
//                 for( int comp=0; comp<3; ++comp ) testPoint[i + comp * this->_numGlobalDofsDisp] = 0.0;
//             }
//         } 
//         if( this->_parser.template get<int>( "DerivativeTestConstraint.order" ) ==  1 ) ipoptSolver.testDerivative( testPoint, 1 );
//         if( this->_parser.template get<int>( "DerivativeTestConstraint.order" ) ==  2 ) ipoptSolver.testDerivative( testPoint, 2 );
//         if( this->_parser.template get<int>( "DerivativeTestConstraint.order" ) == 12 ) ipoptSolver.testDerivative( testPoint, 2 );
// 
//         return true;
//     }
//     //END DerivativeTest


//BEGIN DerivativeTest
    if( this->_parser.template get<int>( "DerivativeTest.order" ) > 0 ){
        VectorType testPoint_Disp = VectorType::Random( 3 * this->_numGlobalDofsDisp );
        for( int i = 0; i < this->_numGlobalDofsDisp; ++i ){
            if ( this->_mask[i] ){
                for( int comp=0; comp<3; ++comp ) testPoint_Disp[i + comp * this->_numGlobalDofsDisp] = 0.0;
            }
        } 
        DerivativeTester<DataTypeContainer> derivTester ( NonLinearEnergyOp, this->_parser.template get<RealType>( "DerivativeTest.stepSize" ), testPoint_Disp, 1, this->_parser.template get<RealType>( "DerivativeTest.tolerance" ) );
        if( this->_parser.template get<int>( "DerivativeTest.order" ) == 1 ) derivTester.testFirstDerivative_AllDirections( );
        if( this->_parser.template get<int>( "DerivativeTest.order" ) == 2 ) derivTester.testSecondDerivative_AllDirections( );
        return true;
    }
//END DerivativeTest
    
    
    SolverInfo<DataTypeContainer> solverInfo;
    ipoptSolver.solve( startPoint_Disp, this->_solDisp, startPoint_Mult, this->_solMultiplierIso, solverInfo );
    
    
    bool solvedToAcceptTol = true;
    if( solverInfo.getError() > breakConditionIpopt ) solvedToAcceptTol = false;
    if( (this->_outputLevel > 0) && !solvedToAcceptTol ){
     std::cout << "=======================================================================" << std::endl;
     cout << "WARNING: IPOPT failed to solve with tolerance = " << breakConditionIpopt << "." << endl;
     cout << "                                        error = " << solverInfo.getError() << "." << endl;
     std::cout << "=======================================================================" << std::endl;
    }
    return solvedToAcceptTol;
  }
  
};


#endif
