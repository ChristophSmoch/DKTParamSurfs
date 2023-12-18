#ifndef __DKTSHELLBUCKLINGOPTIMALDEFORMNLADADAPTIVE_H
#define __DKTSHELLBUCKLINGOPTIMALDEFORMNLADADAPTIVE_H

#include "ShellOptimalDeformSolverInterfaceAdaptive.h"
#include "ShellOptimalDeformSolverNonLin.h"


template<typename MatOptConfigurator, typename NonLinElastEnergyType>
class optimalDeformSolverNonLinAdaptive : public optimalDeformSolverInterfaceAdaptive<MatOptConfigurator>{
  
  typedef typename MatOptConfigurator::ConfiguratorType         ConfiguratorType;
  typedef typename MatOptConfigurator::ConfiguratorTypePf       ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType                   RealType;
  typedef typename ConfiguratorType::MaskType                   MaskType;
  typedef typename ConfiguratorType::InitType                   MeshType;
  typedef typename ConfiguratorTypePf::DTContainer              DataTypeContainer;
  typedef typename ConfiguratorType::VectorType                 VectorType;
  typedef pesopt::BoostParser ParameterParserType;
  
public:
  optimalDeformSolverNonLinAdaptive ( const ParameterParserType &Parser ) :
  optimalDeformSolverInterfaceAdaptive<MatOptConfigurator> ( Parser ) { } 

  void computeOptimalDeformation_OnCurrentMesh( const MeshType & mesh, 
                                           const VectorType& material, 
                                           const VectorType &initDisp, VectorType &solDisp,
                                           MaskType &DirichletMask,
                                           VectorType &membraneStressVec, VectorType &bendingStressVec, VectorType &totalStressVec,
                                           DeformationOptimizationShellEnergyInfo<RealType> &energyInfo,
                                           const string saveDirectory,
                                           const int refinementStep = 0,
                                           const RealType factorForce = 1.
                                         ) const {
                                             
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
    optimalDeformSolver<MatOptConfigurator, false, NonLinElastEnergyType, false> OptDeformFinder ( parser, matOptConf, material, shellHandler, rhs_force, initDisp );
    solDisp.resize( OptDeformFinder.getSolutionDisplacement().size() ); solDisp = OptDeformFinder.getSolutionDisplacement();  
    OptDeformFinder.getEnergyInfo( energyInfo );
    OptDeformFinder.getLastStressOnElements( membraneStressVec, bendingStressVec, totalStressVec );


    
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
    const std::vector<string> nameStressOnElementVec = {"membraneStressOnElements", 
                                                        "bendingStressOnElements", 
                                                        "totalStressOnElements",
                                                        "averagedMembraneStressOnElements", "averagedBendingStressOnElements", 
                                                        "averagedTotalStressOnElements",
                                                        "areaOnElements"};
    this->template saveAllResults<NonLinElastEnergyType::_DiscreteFunctionCacheType,false> ( saveDirectory,  matOptConf, shellHandler, material, solDisp, StressOnElementVec, nameStressOnElementVec, totalStressVec, areaOnElementsVec, energyInfo );
        
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
        mesh.generateApproximativeTangentSpaceAtNodes( DirichletMask ); 
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
           InitDisplacement.resize( SolutionDisplamentVec[refinementStep - 1].size() ); InitDisplacement = SolutionDisplamentVec[refinementStep - 1];
           SolutionDisplacement.resize( InitDisplacement.size() ); SolutionDisplacement = InitDisplacement;
          
          if( this->_parser.template get<bool> ( "ConstraintProblem.adaptiveMarkingUseTotalStressVec") ){
              this->template markAndRefineForGivenElementErrorVector<NonLinElastEnergyType::_DiscreteFunctionCacheType> ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask, totalStressVector[refinementStep-1], markedElements  ); 
          }else{
              this->template markAndRefine<NonLinElastEnergyType::_DiscreteFunctionCacheType> ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask, markedElements );
          }
      }else{
          InitDisplacement.resize( initDisp.size() ); InitDisplacement = initDisp;
          SolutionDisplacement.resize( InitDisplacement.size() ); SolutionDisplacement = InitDisplacement;
      }

      
      //! SOLVE
      VectorType membraneStressVec ( mesh.getNumElements() ), bendingStressVec ( mesh.getNumElements() ), totalStressVec( mesh.getNumElements() );
      DeformationOptimizationShellEnergyInfo<RealType> energyInfo;
      computeOptimalDeformation_OnCurrentMesh ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask, 
                                                membraneStressVec, bendingStressVec, totalStressVec,
                                                energyInfo, saveDirectoryRefinementStep, refinementStep  ); 
      
      //! SAVE
      SolutionDisplamentVec.push_back( SolutionDisplacement );
      gridSizeVec.push_back( DKTTriangMeshInfo<MeshType> ( mesh ).getMinAreaSqrt() ); 
      membraneStressVector.push_back( membraneStressVec ); bendingStressVector.push_back( bendingStressVec ); totalStressVector.push_back( totalStressVec );
      energyInfoVec.push_back( energyInfo );
      
      //!save and plot marked elements 
//       if( refinementStep > 0 ){
//         cout << endl << "save marked elements" << endl;
//         ConfiguratorType oldConf ( oldMesh );
//         ShellHandler<ConfiguratorType,NonLinElastEnergyType::_DiscreteFunctionCacheType> oldShellHandler ( this->_parser, oldConf );
//         oldShellHandler.clearData();
//         oldShellHandler.updateMesh( "undeform", oldShellHandler.getChartToUndeformedShell() );
//         oldShellHandler.addScalarData ( markedElements, "markedElements", FACE_DATA );
//         oldShellHandler.plotShell( "markedElements_Ref" + to_string(refinementStep) );
//         cout << "finished save marked elements" << endl;
//         if( this->_parser.template get<int> ( "saving.plotResults" ) == 1 ){
//             cout << endl << "plot marked elements" << endl;
//             ParameterParserType parserVTKUndeformed ( this->_parser.template get<string> ( "InputMesh.parserFileVTKPlotUndeformed" ) );
//             ShellVTKPlotter<DataTypeContainer> shellVtkPlotter( this->_parser.template get<string> ("saving.saveDirectory"), this->_VTKFileType );
//             shellVtkPlotter.plotMarkedElementsOnUndeformedShell( "markedElements_Ref" + to_string(refinementStep), parserVTKUndeformed, "markedElements");
//             cout << endl << "finished plot marked elements" << endl;
//         }
//       }
      
    }// end for refinement step
    
    
    //
    initDisp.resize( SolutionDisplamentVec[numAdaptiveRefinementSteps].size() ); initDisp = SolutionDisplamentVec[numAdaptiveRefinementSteps];
    
    
    this->template plotConvergenceOfApdaptiveRefinement<false>( energyInfoVec );
    
    if( this->_parser.template get<bool> ("saving.removeOldRefinementSteps") ){
        for( int refinementStep = 0; refinementStep < saveDirectoryVec.size() - 1; ++refinementStep ) pesopt::deleteDirectory( saveDirectoryVec[refinementStep] );
        const string saveDirectoryFinal = this->_parser.template get<string> ("saving.saveDirectory") + "/FinalResult";
        pesopt::moveDirectory( saveDirectoryVec[saveDirectoryVec.size()-1], saveDirectoryFinal );
    }else{
        const string saveDirectoryFinal = this->_parser.template get<string> ("saving.saveDirectory") + "/FinalResult";
        pesopt::copyDirectory( saveDirectoryVec[saveDirectoryVec.size()-1], saveDirectoryFinal );
    }
    
  }//end computeOptimalDeformation_Adaptive
  
  
  void compareOptimalDeformationForScaledForces ( MeshType &mesh,  VectorType& material,  VectorType &initDisp ) const {                                   
        std::vector<double> checkForcesVec;
        this->_parser.template getFreeSizeVector<double, std::vector<double>> ("Compare.checkScaledForceForFinalSolutionVec", checkForcesVec );
        ConfiguratorType conf ( mesh );
        ShellHandlerInterface<ConfiguratorType> shellHandler( conf, this->_parser.template get<string>( "InputMesh.chartXAType" ), 
                                                                   static_cast<ShellBoundaryType>( this->_parser.template get<int>( "InputMesh.ShellType" ) ), 
                                                                   this->_parser.template get<bool> ( "InputMesh.ClampedBoundaryCondition" ) );
        MaskType DirichletMask ( shellHandler.getDirichletMask() ); 
        VectorType initDispCheck ( initDisp ), solDispCheck( initDisp.size() );
       
        const int numAdaptiveRefinementSteps = 0;
        for( int i=0; i<checkForcesVec.size(); ++i ){
            
            VectorType membraneStressVec ( mesh.getNumElements() ), bendingStressVec ( mesh.getNumElements() ), totalStressVec( mesh.getNumElements() );
            DeformationOptimizationShellEnergyInfo<RealType> energyInfo;
            
            RealType factorScaleForce = checkForcesVec[i];
            if( i > 0 ) initDispCheck *= checkForcesVec[i] / checkForcesVec[i-1];
            else        initDispCheck *= checkForcesVec[i];
            computeOptimalDeformation_OnCurrentMesh ( mesh, material, initDispCheck, solDispCheck, DirichletMask, 
                                                      membraneStressVec, bendingStressVec, totalStressVec, energyInfo,
                                                      numAdaptiveRefinementSteps,
                                                      "-CheckForceScale" + std::to_string( factorScaleForce ), 
                                                      factorScaleForce ); 
            initDispCheck = solDispCheck;
        }                             
  }
    
};
    


#endif
