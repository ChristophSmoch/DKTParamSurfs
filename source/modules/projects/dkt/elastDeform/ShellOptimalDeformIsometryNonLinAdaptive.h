#ifndef __SHELLBUCKLINGOPTIMALDEFORMISOMETRYNONLINADADAPTIVE_H
#define __SHELLBUCKLINGOPTIMALDEFORMISOMETRYNONLINADADAPTIVE_H

#include "ShellOptimalDeformSolverInterfaceAdaptive.h"
#include "ShellOptimalDeformSolverIsometryNonLin.h"
#include "ShellCurvature.h"


// TEST 
# include "ShellDeformationEnergiesNonLin.h"
# include "../paramSurfaces/ShellDeformationEnergiesSemiNonLin.h"

template<typename MatOptConfigurator, typename NonLinElastEnergyType>
class optimalDeformSolverIsometryNonLinAdaptive 
: public optimalDeformSolverInterfaceAdaptive<MatOptConfigurator>{
  
  typedef typename MatOptConfigurator::ConfiguratorType         ConfiguratorType;
  typedef typename MatOptConfigurator::ConfiguratorTypePf       ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType                   RealType;
  typedef typename ConfiguratorType::MaskType                   MaskType;
  typedef typename ConfiguratorType::InitType                   MeshType;
  typedef typename ConfiguratorTypePf::DTContainer              DataTypeContainer;
  typedef typename ConfiguratorType::VectorType                 VectorType;
  typedef pesopt::BoostParser ParameterParserType;
  
public:
  optimalDeformSolverIsometryNonLinAdaptive ( const ParameterParserType &Parser ) :
  optimalDeformSolverInterfaceAdaptive<MatOptConfigurator> ( Parser ) { } 

  void computeOptimalDeformation_OnCurrentMesh( const MeshType & mesh, 
                                           const VectorType& material, 
                                           const VectorType &initDisp, VectorType &solDisp,
                                           MaskType &DirichletMask,
                                           VectorType &membraneStressVec, VectorType &bendingStressVec, VectorType &totalStressVec,
                                           DeformationOptimizationShellEnergyInfo<RealType> &energyInfo,
                                           const string saveDirectory,
                                           const int refinementStep = 0,
                                           const RealType factorForce = 1.     ) const {
       
    ParameterParserType parser; parser = this->_parser;
    parser.set ( "saving.saveDirectory", saveDirectory );
    parser.saveToFile( pesopt::strprintf ( "ParameterParser.ini" ) );
    mesh.printBasicInfo();
    auto startTime = std::chrono::high_resolution_clock::now(); 
    
    //initialize configurators on current level
    ConfiguratorType conf ( mesh ); ConfiguratorTypePf confpf ( mesh );
    MatOptConfigurator matOptConf ( parser, conf, confpf );
      
    ShellHandler<ConfiguratorType,NonLinElastEnergyType::_DiscreteFunctionCacheType> shellHandler ( parser, conf );
    if( parser.template get<bool> ( "ConstraintProblemAdaptiveProlongationTypes.constructProlongatedBoundaryByShellHandler" ) )  DirichletMask = shellHandler.getDirichletMask();
    else shellHandler.setDirichletMask( DirichletMask ); 
    
    //! construct rhs 
    VectorType rhs_force ( 3 * conf.getNumGlobalDofs() );
    this->constructForceRHS( parser, matOptConf, shellHandler.getChartToUndeformedShell_Cache(), shellHandler.getDirichletMask(), rhs_force, factorForce );
    
    //! Compute Optimal Deformation
    optimalDeformSolver<MatOptConfigurator, false, NonLinElastEnergyType, true> OptDeformFinder ( parser, matOptConf, material, shellHandler, rhs_force, initDisp );
    solDisp.resize( OptDeformFinder.getSolutionDisplacement().size() ); solDisp = OptDeformFinder.getSolutionDisplacement();  
    OptDeformFinder.getEnergyInfo( energyInfo );
    OptDeformFinder.getLastStressOnElements( membraneStressVec, bendingStressVec, totalStressVec );
    
    
    VectorType xB ( solDisp.size() ); xB = shellHandler.getChartToUndeformedShell() + solDisp;
    DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xBStorage ( matOptConf._conf, xB, 3 );    
    // seminonlin
    SemiNonlinearBendingEnergy<MatOptConfigurator> testSemiNonlinearEnergyOp( matOptConf,
                             shellHandler.getChartToUndeformedShell_Cache(),
                             xBStorage,
                             material,
                             matOptConf._materialInfo._factorBendingEnergy  );
    RealType testSemiNonlinearEnergy;
    testSemiNonlinearEnergyOp.assembleAdd ( testSemiNonlinearEnergy );
    // full nonlin
    NonlinearBendingEnergy<MatOptConfigurator> testNonlinearEnergyOp( matOptConf,
                             shellHandler.getChartToUndeformedShell_Cache(),
                             xBStorage,
                             material,
                             matOptConf._materialInfo._factorBendingEnergy  );
    RealType testNonlinearEnergy;
    testNonlinearEnergyOp.assembleAdd ( testNonlinearEnergy );
    cout << "==========================" <<  endl
         << "seminonlin energy after optimization = " <<  testSemiNonlinearEnergy  <<  endl
         << "nonlin energy after optimization     = " <<  testNonlinearEnergy  <<  endl
         << "==========================" <<  endl;
    
    
    //! print elapsed time
    std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << endl << "duration of refinementStep = " << diff.count() << " sec" << endl;
    //solverInfo.setRunningTime( diff.count() );
    pesopt::consoleOutput( pesopt::strprintf( "finished computation optimal design for refinementStep %d", refinementStep ).c_str() );
    
    
    //! save results
    VectorType areaOnElementsVec ( mesh.getNumElements() );
    evaluateAreaOfElements<ConfiguratorType,NonLinElastEnergyType::_DiscreteFunctionCacheType>( conf, shellHandler, areaOnElementsVec ); 
    VectorType averagedMembraneStressVec ( membraneStressVec ), averagedBendingStressVec( bendingStressVec ), averagedTotalStressVec( totalStressVec );
    for( int elIdx=0; elIdx < mesh.getNumElements(); ++elIdx ){
        averagedMembraneStressVec[elIdx] = membraneStressVec[elIdx] / areaOnElementsVec[elIdx];
        averagedBendingStressVec[elIdx] = bendingStressVec[elIdx] / areaOnElementsVec[elIdx];
        averagedTotalStressVec[elIdx] = totalStressVec[elIdx] / areaOnElementsVec[elIdx];  
    }
    const std::vector<VectorType> StressOnElementVec = {membraneStressVec, bendingStressVec, totalStressVec,
                                                        averagedMembraneStressVec, averagedBendingStressVec, averagedTotalStressVec,
                                                        areaOnElementsVec };
    const std::vector<string> nameStressOnElementVec = {"membraneStressOnElements", "bendingStressOnElements", "totalStressOnElements",
                                                        "averagedMembraneStressOnElements", "averagedBendingStressOnElements", "averagedTotalStressOnElements", 
                                                        "areaOnElements"};
    this->template saveAllResults<NonLinElastEnergyType::_DiscreteFunctionCacheType,true> ( saveDirectory, matOptConf, shellHandler, material, solDisp, StressOnElementVec, nameStressOnElementVec, totalStressVec, areaOnElementsVec, energyInfo );
    
    //!plot results
    if( parser.template get<int> ( "saving.plotResults" ) == 1 ){
        cout << endl << "plot results" << endl;
        this->plotResults( saveDirectory, nameStressOnElementVec );
        runPdflatexWithBashFile( "AllResults.tex", saveDirectory );
    }
    
    
  }//end computeOptimalDeformation_OnCurrentMesh


  
  

  
  
  
  void computeOptimalDeformation_Adaptive ( MeshType &mesh,
                                            VectorType& material,
                                            VectorType &initDisp,
                                            const int numAdaptiveRefinementSteps ) const {
                                                
    pesopt::consoleOutput( "compute optimal design with adaptive refinement" );
    
    //! informations to store for all refinementsteps
    std::vector<VectorType> SolutionDisplamentVec;
    std::vector<RealType> gridSizeVec;
    std::vector<VectorType> membraneStressVector, bendingStressVector, totalStressVector;
    std::vector<DeformationOptimizationShellEnergyInfo<RealType>> energyInfoVec;
    std::vector<string> saveDirectoryVec;
    
    
    ConfiguratorType conf ( mesh );
    ShellHandlerInterface<ConfiguratorType> shellHandler( conf, this->_parser.template get<string>( "InputMesh.chartXAType" ), static_cast<ShellBoundaryType>( this->_parser.template get<int>( "InputMesh.ShellType" ) ), this->_parser.template get<bool> ( "InputMesh.ClampedBoundaryCondition" ) );
    MaskType DirichletMask ( shellHandler.getDirichletMask() );
    if( this->_parser.template get<int> ( "InputMesh.tangentSpaceType" ) == 2 ){
        mesh.generateApproximativeTangentSpaceAtNodes( DirichletMask ); //! \todo should be done directly in triangleMesh?
//         mesh.updateAllProjectionCoefficients(); 
    }
    
    
    //! use adaptive refinement
    for( int refinementStep=0; refinementStep <= numAdaptiveRefinementSteps; ++refinementStep ) {
        
      const string saveDirectoryRefinementStep = this->_parser.createSubDirectory( "refinementStep" + std::to_string(refinementStep) );
      saveDirectoryVec.push_back( saveDirectoryRefinementStep );
        
      VectorType InitDisplacement, SolutionDisplacement;
      MeshType oldMesh( mesh ); VectorType markedElements;  
      
      //! MARK and REFINE
      if( refinementStep > 0 ){
           cout << endl << "start mark and refine" << endl;
           auto startTime = std::chrono::high_resolution_clock::now();         
           InitDisplacement.resize( SolutionDisplamentVec[refinementStep - 1].size() ); InitDisplacement = SolutionDisplamentVec[refinementStep - 1];
           SolutionDisplacement.resize( InitDisplacement.size() ); SolutionDisplacement = InitDisplacement;
          if( this->_parser.template get<bool> ( "ConstraintProblem.adaptiveMarkingUseTotalStressVec") ){
              this->template markAndRefineForGivenElementErrorVector<NonLinElastEnergyType::_DiscreteFunctionCacheType> ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask, totalStressVector[refinementStep-1], markedElements  ); 
          }else{
              this->template markAndRefine<NonLinElastEnergyType::_DiscreteFunctionCacheType> ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask, markedElements );
          }
          std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
          std::cout << endl << "duration for markAndRefine = " << diff.count() << " sec" << endl;
      }else{
          InitDisplacement.resize( initDisp.size() ); InitDisplacement = initDisp;
          SolutionDisplacement.resize( InitDisplacement.size() ); SolutionDisplacement = InitDisplacement;
      }

      //! set initialmaterial on refined mesh
      VectorType initMaterial ( material.size() );
      if( refinementStep == 0 ) this->setInitialMaterial( this->_parser.template get<int> ( "ConstraintProblem.firstStepMaterialInitializationType" ), mesh, initMaterial, material  );
      else                      this->setInitialMaterial( this->_parser.template get<int> ( "ConstraintProblem.adaptiveMaterialInitializationType" ), mesh, initMaterial, material  );
      material = initMaterial;
      
      //! SOLVE
      VectorType membraneStressVec ( mesh.getNumElements() ), bendingStressVec ( mesh.getNumElements() ), totalStressVec( mesh.getNumElements() );
      DeformationOptimizationShellEnergyInfo<RealType> energyInfo;
      computeOptimalDeformation_OnCurrentMesh ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask, 
                                                membraneStressVec, bendingStressVec, totalStressVec,
                                                energyInfo, saveDirectoryRefinementStep, refinementStep ); 
      
      //! SAVE
      SolutionDisplamentVec.push_back( SolutionDisplacement );
      gridSizeVec.push_back( DKTTriangMeshInfo<MeshType> ( mesh ).getMinAreaSqrt() ); 
      membraneStressVector.push_back( membraneStressVec ); bendingStressVector.push_back( bendingStressVec ); totalStressVector.push_back( totalStressVec );
      energyInfoVec.push_back( energyInfo );
      
    }// end for refinement step
    
    
    //
    initDisp.resize( SolutionDisplamentVec[numAdaptiveRefinementSteps].size() ); initDisp = SolutionDisplamentVec[numAdaptiveRefinementSteps];
    
    this->template plotConvergenceOfApdaptiveRefinement<true>( energyInfoVec );
    
    if( this->_parser.template get<bool> ("saving.removeOldRefinementSteps") ){
        for( int refinementStep = 0; refinementStep < saveDirectoryVec.size() - 1; ++refinementStep ) pesopt::deleteDirectory( saveDirectoryVec[refinementStep] );
        const string saveDirectoryFinal = this->_parser.template get<string> ("saving.saveDirectory") + "/FinalResult";
        pesopt::moveDirectory( saveDirectoryVec[saveDirectoryVec.size()-1], saveDirectoryFinal );
    }else{
        const string saveDirectoryFinal = this->_parser.template get<string> ("saving.saveDirectory") + "/FinalResult";
        pesopt::copyDirectory( saveDirectoryVec[saveDirectoryVec.size()-1], saveDirectoryFinal );
    }
    
  }//end computeOptimalDeformation_Adaptive
    
};
    


#endif
