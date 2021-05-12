#ifndef __SHELLBUCKLINGOPTIMALDEFORMISOMETRYSEMINONLINADADAPTIVE_H
#define __SHELLBUCKLINGOPTIMALDEFORMISOMETRYSEMINONLINADADAPTIVE_H

#include <pesopt_VTK.h>
#include <energyDefines.h>

#include "../elastDeform/ShellOptimalDeformSolverInterfaceAdaptive.h"
#include "ShellOptimalDeformSolverIsometrySemiNonLin.h"
#include "../elastDeform/ShellCurvature.h"


template<typename MatOptConfigurator, typename NonLinElastEnergyType>
class optimalDeformSolverIsometrySemiNonLinAdaptive
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
  optimalDeformSolverIsometrySemiNonLinAdaptive ( const ParameterParserType &Parser ) :
  optimalDeformSolverInterfaceAdaptive<MatOptConfigurator> ( Parser ) { }

  void computeOptimalDeformation_OnCurrentMesh( const MeshType & mesh,
                                           const VectorType& material,
                                           const VectorType &initDisp, VectorType &solDisp,
                                           MaskType &DirichletMask,
                                           VectorType &membraneStressVec, VectorType &bendingStressVec, VectorType &totalStressVec,
                                           VectorType &GaussCurvVec,
                                           DeformationOptimizationShellEnergyInfo<RealType> &energyInfo,
                                           IsometryInfo<RealType> &isometryInfo,
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

    const int numVertices = mesh.getNumVertices();
    const int numGlobalDofsDeform = conf.getNumGlobalDofs();
    VectorType xB ( solDisp.size() ); xB = solDisp + shellHandler.getChartToUndeformedShell();
    DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xBStorage ( matOptConf._conf, xB, 3 );

    //Gauss Curvature
    VectorType GaussCurvVector ( totalStressVec.size() );
    GaussCurvatureL1<ConfiguratorType> ( matOptConf._conf, xBStorage ).assembleOnElements( GaussCurvVector );
    RealType integralGauss = 0.;
    RealType integralGauss1 = 0.;
    GaussCurvatureL1<ConfiguratorType> ( matOptConf._conf, xBStorage ).assembleAdd( integralGauss );
    isometryInfo.setGaussCurvatureL1( integralGauss );
    GaussCurvatureInt<ConfiguratorType> ( matOptConf._conf, xBStorage ).assembleAdd( integralGauss1 );
    isometryInfo.setGaussCurvatureInt( integralGauss1 );
    GaussCurvVec = GaussCurvVector;


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
    std::vector<VectorType> GaussCurvVector;
    std::vector<DeformationOptimizationShellEnergyInfo<RealType>> energyInfoVec;
    std::vector<IsometryInfo<RealType>> isometryInfoVec;
    std::vector<string> saveDirectoryVec;


    ConfiguratorType conf ( mesh );
    ShellHandlerInterface<ConfiguratorType> shellHandler( conf, this->_parser.template get<string>( "InputMesh.chartXAType" ), static_cast<ShellBoundaryType>( this->_parser.template get<int>( "InputMesh.ShellType" ) ), this->_parser.template get<bool> ( "InputMesh.ClampedBoundaryCondition" ) );
    MaskType DirichletMask ( shellHandler.getDirichletMask() );
    if( this->_parser.template get<int> ( "InputMesh.tangentSpaceType" ) == 2 ){
        mesh.generateApproximativeTangentSpaceAtNodes( DirichletMask ); //! \todo should be done directly in triangleMesh?
        mesh.updateAllProjectionCoefficients();
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
          if( this->_parser.template get<bool> ( "ConstraintProblem.adaptiveMarkingUseGaussCurvatureVec") ){
              this->template markAndRefineForGivenElementErrorVector<NonLinElastEnergyType::_DiscreteFunctionCacheType> ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask, GaussCurvVector[refinementStep-1], markedElements  );
          }else{
              this->template markAndRefine<NonLinElastEnergyType::_DiscreteFunctionCacheType> ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask, markedElements );
          }
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
      VectorType GaussCurvVec ( mesh.getNumElements() );
      DeformationOptimizationShellEnergyInfo<RealType> energyInfo;
      IsometryInfo<RealType> isometryInfo;
        computeOptimalDeformation_OnCurrentMesh ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask,
                                                membraneStressVec, bendingStressVec, totalStressVec, GaussCurvVec,
                                                energyInfo, isometryInfo, saveDirectoryRefinementStep, refinementStep );

      //! SAVE
      // SolutionDisplamentVec.push_back( SolutionDisplacement );
      // gridSizeVec.push_back( DKTTriangMeshInfo<MeshType> ( mesh ).getMinAreaSqrt() );
      // membraneStressVector.push_back( membraneStressVec ); bendingStressVector.push_back( bendingStressVec ); totalStressVector.push_back( totalStressVec );
      // energyInfoVec.push_back( energyInfo );
      // isometryInfoVec.push_back( isometryInfo );

      SolutionDisplamentVec.push_back( SolutionDisplacement );
      gridSizeVec.push_back( DKTTriangMeshInfo<MeshType> ( mesh ).getMinAreaSqrt() );
      membraneStressVector.push_back( membraneStressVec ); bendingStressVector.push_back( bendingStressVec ); totalStressVector.push_back( totalStressVec );
      GaussCurvVector.push_back( GaussCurvVec );
      energyInfoVec.push_back( energyInfo );
      isometryInfo.setGridSize( gridSizeVec[refinementStep] );
      isometryInfoVec.push_back( isometryInfo );
      isometryInfo.template saveToFile<ParameterParserType>( "isometryInfo", saveDirectoryRefinementStep );

    }// end for refinement step


    //
    initDisp.resize( SolutionDisplamentVec[numAdaptiveRefinementSteps].size() ); initDisp = SolutionDisplamentVec[numAdaptiveRefinementSteps];

    // this->template plotConvergenceOfApdaptiveRefinement<true>( energyInfoVec );
    this->template plotConvergenceOfApdaptiveRefinement<true>( energyInfoVec );

    ShellHandler<ConfiguratorType,NonLinElastEnergyType::_DiscreteFunctionCacheType> shellHandlerFine ( this->_parser, conf );


    VectorType fineSolution = SolutionDisplamentVec[numAdaptiveRefinementSteps];
    std::vector<VectorType> SolutionDisplamentVecProlongated;
    std::vector<RealType> ErrorD2L2_FineSolutionVec(numAdaptiveRefinementSteps+1), EOCD2L2_FineSolutionVec(numAdaptiveRefinementSteps), ErrorGaussCurvL1_FineSolutionVec(numAdaptiveRefinementSteps+1);
    for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
            VectorType currentSolutionProlongated = SolutionDisplamentVec[refinementStep];
            mesh.prolongateVectorFct_DKTLinearly( currentSolutionProlongated );
            SolutionDisplamentVecProlongated.push_back( currentSolutionProlongated );
    }
    for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
        DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> fineSolutionDFD ( conf, shellHandlerFine.getChartToUndeformedShell (  ) + fineSolution, 3);
        DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> solDFD( conf, shellHandlerFine.getChartToUndeformedShell (  ) + SolutionDisplamentVecProlongated[refinementStep], 3 );
        RealType tmp = 0.;
        SecondDerivativeEnergy<ConfiguratorType> ( conf, fineSolutionDFD, solDFD ).assembleAdd( tmp );
        ErrorD2L2_FineSolutionVec[refinementStep] = sqrt( tmp );
        isometryInfoVec[refinementStep].setErrorApproxD2uToFineSolutionL2( sqrt(tmp) );

        RealType tmp1 = 0.;
        RealType tmp2 = 0.;
        GaussCurvatureL1Diff<ConfiguratorType> ( conf, fineSolutionDFD, solDFD ).assembleAdd( tmp1 );
        GaussCurvatureL1<ConfiguratorType> ( conf,  solDFD ).assembleAdd( tmp2 );
        ErrorGaussCurvL1_FineSolutionVec[refinementStep] =  tmp1;
        isometryInfoVec[refinementStep].setGaussCurvatureL1Diff( tmp1 );
        // isometryInfoVec[refinementStep].setGaussCurvatureL1( tmp2 );

        // ShellPlotter<ConfiguratorType> shellPlotter( conf, shellHandlerFine.getChartToUndeformedShell(), shellHandlerFine.getDirichletMask(), saveDirectoryVec[refinementStep],   this->_parser.template get<string>("saving.VTKFileType")  );
        // shellPlotter.saveShellToFile ( "disp", SolutionDisplamentVecProlongated[refinementStep], pesopt::strprintf( "prolongatedDeformedShell" ).c_str() );


        // VectorType GaussCurvVector ( totalStressVec.size() );
        // GaussCurvatureL1<ConfiguratorType> ( matOptConf._conf, xBStorage ).assembleOnElements( GaussCurvVector );

    }

    this->plotConvergenceIsometryOfApdaptiveRefinement( isometryInfoVec );


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
