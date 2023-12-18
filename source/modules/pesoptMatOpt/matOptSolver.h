#ifndef __MATOPTSOLVER_H
#define __MATOPTSOLVER_H

#include <pesopt_IO.h>
#include "matOptDefines.h"
#include "matOptVTKPlotter.h"


template<typename ConfiguratorType>
class MaterialOptimizationMacroSolverBase {
  
  typedef typename ConfiguratorType::RealType                   RealType;
  typedef typename ConfiguratorType::MaskType                   MaskType;
  typedef typename ConfiguratorType::InitType                   MeshType;
  typedef typename ConfiguratorType::DTContainer                DataTypeContainer;
  typedef typename ConfiguratorType::VectorType                 VectorType;  
  typedef pesopt::BoostParser ParameterParserType;

  const pesopt::BoostParser _parser;

public:
  MaterialOptimizationMacroSolverBase ( const ParameterParserType & parser ) :
  _parser( parser ) { } 
  
  virtual void computeOptimalDesign_RefinementStep( const ParameterParserType & parserRefinementStep,
                                                    const MeshType & mesh, VectorType &initMaterial, VectorType &solutionMaterial,
                                                    const int refinementStep = 0, const bool onlyComputeInitialEnergy = false ) const {
    throw std::logic_error( pesopt::strprintf( "MaterialOptimizationMultipleLoadSolverBase::computeOptimalDesign():  virtual function! Should be provided in derived class!" ).c_str() );
  }
  
  void computeOptimalDesign_MultiLevel( const MeshType &meshStartLevel, 
                                        const int numRefinementSteps, 
                                        const int StartLevel = 0,
                                        const bool onlyComputeInitialEnergy = false ) const {
                                           
    pesopt::consoleOutput( "compute optimal design with multilevel scheme" );
    ConfiguratorType confStartLevel ( meshStartLevel );
    FEMaterialGenerator<ConfiguratorType> FEMaterialGeneratorStartLevel ( _parser, confStartLevel );
    //! initialize material
    VectorType initMaterialStartLevel ( confStartLevel.getNumGlobalDofs() );
    FEMaterialGeneratorStartLevel.switchMaterialType ( initMaterialStartLevel );
    
    std::vector<string> saveDirectoryVec;
    
    FEAdaptiveMesh<MeshType> adaptiveMesh( meshStartLevel );
    VectorType initMaterialCurrentLevel ( initMaterialStartLevel );
    VectorType solMaterialCurrentLevel;
    
    for( int level = StartLevel; level <= numRefinementSteps; ++level ){
        
       auto startTime = std::chrono::high_resolution_clock::now(); 
       if( onlyComputeInitialEnergy ) 
           pesopt::consoleOutput( pesopt::strprintf( "compute initial energy for refinementStep %d", level ).c_str() );
       else                          
           pesopt::consoleOutput( pesopt::strprintf( "compute optimal design for refinementStep %d", level ).c_str() );
        
        pesopt::BoostParser parserRefinementStep; parserRefinementStep = _parser;
        const string saveDirectoryRefinementStep = parserRefinementStep.createSubDirectory( "refinementStep" + std::to_string(level) );
        parserRefinementStep.set ( "saving.saveDirectory", saveDirectoryRefinementStep );
        saveDirectoryVec.push_back( saveDirectoryRefinementStep );
       
        if( level > StartLevel ){  
            //refinement option
            const string adaptiveRefinementType = parserRefinementStep.template get<string> ( "MaterialOptimization.adaptiveRefinementType" );
            if( adaptiveRefinementType == "refineAll" )
                adaptiveMesh.refineAll( );
            else if( adaptiveRefinementType == "refineByFactor" )
                adaptiveMesh.refineByFactor( _parser.template get<RealType> ( "MaterialOptimization.factorMultiLevel" ) );
            else
                throw std::invalid_argument( pesopt::strprintf ( "Wrong adaptiveRefinementType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
            
            adaptiveMesh.template prolongateScalarFunctionAtNodes<ConfiguratorType> ( solMaterialCurrentLevel, initMaterialCurrentLevel );
        }
        
        //! save mesh and parser to file 
        VTKMeshSaver<MeshType> meshSaver( adaptiveMesh.getCurrentMesh() );
        const string fileNameMesh = pesopt::strprintf( "%s/mesh.vtk", saveDirectoryRefinementStep.c_str() );
        meshSaver.save( fileNameMesh, MeshType::_VTKDATATYPEUNDEFORMED );
        adaptiveMesh.getCurrentMesh().saveToParser( parserRefinementStep, fileNameMesh );
        parserRefinementStep.saveToFile( pesopt::strprintf ( "ParameterParser.ini" ) );
        
        //! solve for refinement step
        solMaterialCurrentLevel.resize( initMaterialCurrentLevel.size() );
        this->computeOptimalDesign_RefinementStep ( parserRefinementStep, adaptiveMesh.getCurrentMesh(), initMaterialCurrentLevel, solMaterialCurrentLevel, level, onlyComputeInitialEnergy ); 
        
        //! print elapsed time
        std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
        std::cout << endl << "duration of refinementStep = " << diff.count() << " sec" << endl;
        pesopt::consoleOutput( pesopt::strprintf( "finished computation optimal design for refinementStep %d", level ).c_str() );

    }

    const string saveDirectoryFinal = _parser.createSubDirectory( "FinalResult" ); 
    pesopt::moveDirectory( saveDirectoryVec[saveDirectoryVec.size()-1], saveDirectoryFinal );
    if( _parser.template get<bool> ("saving.removeOldRefinementSteps") ){
        for( int refinementStep = 0; refinementStep < saveDirectoryVec.size() - 1; ++refinementStep ) 
            pesopt::deleteDirectory( saveDirectoryVec[refinementStep] );
    }
    
  }
  
private :
    
  void saveInfo( const pesopt::NonlinearEnergyOp<DataTypeContainer> &totalEnergyOp,
                 const pesopt::NonlinearConstraintOps<DataTypeContainer> &constraintOps,
                 const VectorType &material, 
                 ParameterParserType & info
               ) const{

    totalEnergyOp.getEnergyInfo( material, info );
    constraintOps.getConstraintInfo( material, info );
    info.saveToFile( pesopt::strprintf ( "energyInfo.ini" ) );
  }
    
    
  template <typename OptimalDeformSolver>
  void saveResultsToFile( const ParameterParserType & parserRefinementStep,
                          const OptimalDeformSolver &OptDeformFinder,
                          const VectorType &material
                        ) const{
      
    const string saveDirectoryRefinementStep = parserRefinementStep.template get<string> ( "saving.saveDirectory" );

    //! save material vector
    pesopt::printVectorToFile<VectorType> ( material, pesopt::strprintf( "%s/materialCollabsed.txt", saveDirectoryRefinementStep.c_str() ).c_str(), 10 ); 
    
    //! save extended material vector
    VectorType materialExtended ( material );
    OptDeformFinder.getBoundaryHandlerMaterial().extendVector( materialExtended );
    pesopt::printVectorToFile<VectorType> ( materialExtended, pesopt::strprintf( "%s/materialExtended.txt", saveDirectoryRefinementStep.c_str() ).c_str(), 10 );
    
    //! save deformation vectors 
    const string saveDirectoryDeformation = parserRefinementStep.createSubDirectory( "Deformation" );
    for( int loadIdx=0; loadIdx< OptDeformFinder.getNumLoads(); ++loadIdx ){
         pesopt::printVectorToFile<VectorType> ( OptDeformFinder.getSolutionDisplacement(loadIdx),                                                 
                                              pesopt::strprintf( "%s/Displacement_Dir%d.txt", saveDirectoryDeformation.c_str(), loadIdx ).c_str(), 10 );
    }
  }
    
  void saveResultsVTK( const ParameterParserType & parserRefinementStep ) const{
      
//     MaterialOptimizationPeriodicCellVTKSaver<ConfiguratorType> vtkSaver( parserRefinementStep );
    MaterialOptimizationMacroVTKSaver<ConfiguratorType> vtkSaver( parserRefinementStep );
    vtkSaver.saveMaterial(  );
//     vtkSaver.saveSingleCell( );

//     vtkSaver.saveBlockOfCells(  parserRefinementStep.template get<int> ("saving.numBlocks") );
    if( parserRefinementStep.template get<bool> ( "saving.saveDeformations" ) )
        vtkSaver.saveDeformations( );
//     if( parserRefinementStep.template get<bool> ("saving.saveInterface") ) 
//         vtkSaver.saveInterface( );
//     if( parserRefinementStep.template get<bool> ("saving.saveVonMisesOnInterface") )
//         vtkSaver.saveVonMisesOnInterface();
//     if( parserRefinementStep.template get<bool>("saving.saveVonMises") ) 
//         vtkSaver.saveVonMisesStresses( );
  }
 
  //NOTE can only plot to png if it was saved to vtk before, i.e. in function saveResultsVTK
  void plotResultsVTKtoPNG( const ParameterParserType & parserRefinementStep ) const{    
    
    //! plot vtk files to png
    MaterialOptimizationVTKtoPNGPlotter<DataTypeContainer> VTKtoPNGPlotter( parserRefinementStep );
    const int numLoads = 1;
    
    if( parserRefinementStep.template get<bool> ("saving.plotPNGUndeformed") )
        VTKtoPNGPlotter.plotUndeformedToPNG( );  
//     if( parserRefinementStep.template get<bool> ("saving.plotPNGBlockOfCells") )
//         VTKtoPNGPlotter.plotBlockOfCellsToPNG( );  
    if( parserRefinementStep.template get<bool> ("saving.plotPNGDeformed") )
        VTKtoPNGPlotter.plotDeformationToPNG( numLoads ); 
//     if( parserRefinementStep.template get<bool> ("saving.plotPNGInterface") ){
//             VTKtoPNGPlotter.plotInterfaceToPNG( );
//             if( parserRefinementStep.template get<bool> ("saving.plotPNGVonMisesOnInterface") )
//                VTKtoPNGPlotter.plotVonMisesStressesOnInterfaceToPNG( numLoads );
//     }
    
  }

protected :
  
  template <typename OptimalDeformSolver>
  void saveAndPlotAllResults_RefinementStep( const string name, 
                                             const ParameterParserType & parserRefinementStep,
                                             const pesopt::NonlinearEnergyOp<DataTypeContainer> &totalEnergyOp,
                                             const pesopt::NonlinearConstraintOps<DataTypeContainer> &constraintOps,
                                             const OptimalDeformSolver &OptDeformFinder,
                                             const VectorType &material,
                                             const bool saveVTK,
                                             const bool plotVTKtoPNG
                                           ) const {
 
    pesopt::BoostParser parserMaterial ( parserRefinementStep );
    const string subDirectoryMaterial = parserMaterial.createSubDirectory( name.c_str() );
    parserMaterial.set ( "saving.saveDirectory", subDirectoryMaterial );
    
    pesopt::BoostParser energyInfo; energyInfo.set ( "saving.saveDirectory", subDirectoryMaterial );
    this->saveInfo( totalEnergyOp, constraintOps, material, energyInfo );
   
    this->template saveResultsToFile<OptimalDeformSolver>( parserMaterial, OptDeformFinder, material );
    
    //! save results as vtk
    if( saveVTK ) this->saveResultsVTK( parserMaterial );
    
    //!plot results to png
    if( plotVTKtoPNG ) this->plotResultsVTKtoPNG( parserMaterial );
                                                    
   }
  
};



template<typename ConfiguratorType>
class MaterialOptimizationFineScaleMultipleLoadSolverBase {
  
  typedef typename ConfiguratorType::RealType                   RealType;
  typedef typename ConfiguratorType::MaskType                   MaskType;
  typedef typename ConfiguratorType::InitType                   MeshType;
  typedef typename ConfiguratorType::DTContainer                DataTypeContainer;
  typedef typename DataTypeContainer::PointType                 PointType;
  typedef typename ConfiguratorType::VectorType                 VectorType;
  typedef typename ConfiguratorType::VoigtTensorVec             VoigtTensorVec;
  typedef pesopt::BoostParser ParameterParserType;

  const pesopt::BoostParser _parser;

public:
  MaterialOptimizationFineScaleMultipleLoadSolverBase ( const ParameterParserType & parser ) :
  _parser( parser ) { } 
  
  virtual void computeOptimalDesign_RefinementStep( const ParameterParserType & parserRefinementStep,
                                                    const MeshType & mesh, VectorType &initMaterial, VectorType &solutionMaterial,
                                                    const int refinementStep = 0, const bool onlyComputeInitialEnergy = false ) const {
    throw std::logic_error( "MaterialOptimizationMultipleLoadSolverBase::computeOptimalDesign(): Unimplemented function! Should be provided in derived class!");
  }
 
  void computeOptimalDesign_MultiLevel( const MeshType &meshStartLevel, 
                                        const int numRefinementSteps, 
                                        const int StartLevel = 0,
                                        const bool onlyComputeInitialEnergy = false ) const {
                                           
    pesopt::consoleOutput( "compute optimal design with multilevel scheme" );
    ConfiguratorType confStartLevel ( meshStartLevel );
    FEMaterialGenerator<ConfiguratorType> FEMaterialGeneratorStartLevel ( _parser, confStartLevel );
    //! initialize material
    VectorType initMaterialStartLevel ( confStartLevel.getNumGlobalDofs() );
    FEMaterialGeneratorStartLevel.switchMaterialType ( initMaterialStartLevel );
    
    std::vector<string> saveDirectoryVec;
    
    FEAdaptiveMesh<MeshType> adaptiveMesh( meshStartLevel );
    VectorType initMaterialCurrentLevel ( initMaterialStartLevel );
    VectorType solMaterialCurrentLevel;
    
    for( int level = StartLevel; level <= numRefinementSteps; ++level ){
        
       auto startTime = std::chrono::high_resolution_clock::now(); 
       if( onlyComputeInitialEnergy ) 
           pesopt::consoleOutput( pesopt::strprintf( "compute initial energy for refinementStep %d", level ).c_str() );
       else                          
           pesopt::consoleOutput( pesopt::strprintf( "compute optimal design for refinementStep %d", level ).c_str() );
        
        pesopt::BoostParser parserRefinementStep; parserRefinementStep = _parser;
        const string saveDirectoryRefinementStep = parserRefinementStep.createSubDirectory( "refinementStep" + std::to_string(level) );
        parserRefinementStep.set ( "saving.saveDirectory", saveDirectoryRefinementStep );
        saveDirectoryVec.push_back( saveDirectoryRefinementStep );
        
        if( level > StartLevel ){  
            
            //refinement option
            const string adaptiveRefinementType = parserRefinementStep.template get<string> ( "MaterialOptimization.adaptiveRefinementType" );
            if( adaptiveRefinementType == "refineAll" )
                adaptiveMesh.refineAll( );
            else if( adaptiveRefinementType == "refineByFactor" )
                adaptiveMesh.refineByFactor( _parser.template get<RealType> ( "MaterialOptimization.factorMultiLevel" ) );

            else
                throw std::invalid_argument( pesopt::strprintf ( "Wrong adaptiveRefinementType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
            
            adaptiveMesh.template prolongateScalarFunctionAtNodes<ConfiguratorType> ( solMaterialCurrentLevel, initMaterialCurrentLevel );
        }
        
        //! save mesh and parser to file 
        VTKMeshSaver<MeshType> meshSaver( adaptiveMesh.getCurrentMesh() );
        const string fileNameMesh = pesopt::strprintf( "%s/mesh.vtk", saveDirectoryRefinementStep.c_str() );
        meshSaver.save( fileNameMesh, MeshType::_VTKDATATYPEUNDEFORMED );
        adaptiveMesh.getCurrentMesh().saveToParser( parserRefinementStep, fileNameMesh );
        parserRefinementStep.saveToFile( pesopt::strprintf ( "ParameterParser.ini" ) );
        
        //! solve for refinement step
        solMaterialCurrentLevel.resize( initMaterialCurrentLevel.size() );
        this->computeOptimalDesign_RefinementStep ( parserRefinementStep, adaptiveMesh.getCurrentMesh(), initMaterialCurrentLevel, solMaterialCurrentLevel, level, onlyComputeInitialEnergy ); 
       
        //! print elapsed time
        std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
        std::cout << endl << "duration of refinementStep = " << diff.count() << " sec" << endl;
        pesopt::consoleOutput( pesopt::strprintf( "finished computation optimal design for refinementStep %d", level ).c_str() );

    }
    
    const string saveDirectoryFinal = _parser.createSubDirectory( "FinalResult" ); 
    pesopt::moveDirectory( saveDirectoryVec[saveDirectoryVec.size()-1], saveDirectoryFinal );
    if( _parser.template get<bool> ("saving.removeOldRefinementSteps") ){
        for( int refinementStep = 0; refinementStep < saveDirectoryVec.size() - 1; ++refinementStep ) 
            pesopt::deleteDirectory( saveDirectoryVec[refinementStep] );
    }
    
  }
  
private :
    
  void saveInfo( const pesopt::NonlinearEnergyOp<DataTypeContainer> &totalEnergyOp,
                 const pesopt::NonlinearConstraintOps<DataTypeContainer> &constraintOps,
                 const VectorType &material, 
                 ParameterParserType & info
               ) const{

    totalEnergyOp.getEnergyInfo( material, info );
    constraintOps.getConstraintInfo( material, info );
    info.saveToFile( pesopt::strprintf ( "energyInfo.ini" ) );
  }
    
    
  template <typename OptimalDeformSolver>
  void saveResultsToFile( const ParameterParserType & parserRefinementStep,
                          const OptimalDeformSolver &OptDeformFinder,
                          const VectorType &material
                        ) const{
      
    const string saveDirectoryRefinementStep = parserRefinementStep.template get<string> ( "saving.saveDirectory" );
    
    //! save material vector
    pesopt::printVectorToFile<VectorType> ( material, pesopt::strprintf( "%s/materialCollabsed.txt", saveDirectoryRefinementStep.c_str() ).c_str(), 10 ); 
    
    //! save extended material vector
    VectorType materialExtended ( material );
    OptDeformFinder.getBoundaryHandlerMaterial().extendVector( materialExtended );
    pesopt::printVectorToFile<VectorType> ( materialExtended, pesopt::strprintf( "%s/materialExtended.txt", saveDirectoryRefinementStep.c_str() ).c_str(), 10 );
    
    
    //!save deformation vectors 
    const string saveDirectoryDeformation = parserRefinementStep.createSubDirectory( "Deformation" );
    for( int loadIdx=0; loadIdx< OptDeformFinder.getNumLoads(); ++loadIdx ){
        pesopt::printVectorToFile<VectorType> ( OptDeformFinder.getSolutionDisplacementCollabsed(loadIdx), 
                                                pesopt::strprintf( "%s/Displacement_Dir%d.txt", saveDirectoryDeformation.c_str(), loadIdx ).c_str(), 10 ); 
        pesopt::printVectorToFile<VectorType> ( OptDeformFinder.getAffineDisplacement(loadIdx), 
                                                pesopt::strprintf( "%s/DisplacementAffine_Dir%d.txt", saveDirectoryDeformation.c_str(), loadIdx ).c_str(), 10 ); 
    }
    
    //!save boundary conditions
    pesopt::printVectorToFile<MaskType> ( OptDeformFinder.getBoundaryHandlerMaterial( ).getDirichletMask(), pesopt::strprintf( "%s/DirichletBoundaryMaterial.txt", saveDirectoryRefinementStep.c_str() ).c_str(), 10 ); 
    
    pesopt::printVectorToFile<MaskType> ( OptDeformFinder.getBoundaryHandlerDisplacement( ).getDirichletMask(), pesopt::strprintf( "%s/DirichletBoundaryDisplacement.txt", saveDirectoryRefinementStep.c_str() ).c_str(), 10 ); 
    
  }
    
  void saveResultsVTK( const ParameterParserType & parserRefinementStep ) const{
      
    MaterialOptimizationFineScaleVTKSaver<ConfiguratorType> vtkSaver( parserRefinementStep );
    vtkSaver.saveMaterial(  );
    vtkSaver.saveBoundaryMask( );
    if( parserRefinementStep.template get<bool> ( "saving.saveDeformations" ) )
        vtkSaver.saveDeformations( );
  }
 
  //NOTE can only plot to png if it was saved to vtk before, i.e. in function saveResultsVTK
  void plotResultsVTKtoPNG( const ParameterParserType & parserRefinementStep ) const{    
    
    //! plot vtk files to png
    MaterialOptimizationVTKtoPNGPlotter<DataTypeContainer> VTKtoPNGPlotter( parserRefinementStep );

    const int numLoads = 1;
    
    if( parserRefinementStep.template get<bool> ("saving.plotPNGUndeformed") )
        VTKtoPNGPlotter.plotUndeformedToPNG( );  
    if( parserRefinementStep.template get<bool> ("saving.plotPNGDeformed") )
        VTKtoPNGPlotter.plotDeformationToPNG( numLoads ); 
    
  }

protected :
  
  template <typename OptimalDeformSolver>
  void saveAndPlotAllResults_RefinementStep( const string name, 
                                             const ParameterParserType & parserRefinementStep,
                                             const pesopt::NonlinearEnergyOp<DataTypeContainer> &totalEnergyOp,
                                             const pesopt::NonlinearConstraintOps<DataTypeContainer> &constraintOps,
                                             const OptimalDeformSolver &OptDeformFinder,
                                             const VectorType &material,
                                             const bool saveVTK,
                                             const bool plotVTKtoPNG
                                           ) const {
 
    pesopt::BoostParser parserMaterial ( parserRefinementStep );
    const string subDirectoryMaterial = parserMaterial.createSubDirectory( name.c_str() );
    parserMaterial.set ( "saving.saveDirectory", subDirectoryMaterial );
    
    pesopt::BoostParser energyInfo; energyInfo.set ( "saving.saveDirectory", subDirectoryMaterial );
    this->saveInfo( totalEnergyOp, constraintOps, material, energyInfo );
   
    this->template saveResultsToFile<OptimalDeformSolver>( parserMaterial, OptDeformFinder, material );
    
    //! save results as vtk
    if( saveVTK ) this->saveResultsVTK( parserMaterial );
    
    //!plot results to png
    if( plotVTKtoPNG ) this->plotResultsVTKtoPNG( parserMaterial );
                                                    
   }
  
      
};







template<typename ConfiguratorType>
class MaterialOptimizationCellMultipleLoadSolverBase {
  
  typedef typename ConfiguratorType::RealType                   RealType;
  typedef typename ConfiguratorType::MaskType                   MaskType;
  typedef typename ConfiguratorType::InitType                   MeshType;
  typedef typename ConfiguratorType::DTContainer                DataTypeContainer;
  typedef typename DataTypeContainer::PointType                 PointType;
  typedef typename ConfiguratorType::VectorType                 VectorType;
  typedef typename ConfiguratorType::VoigtTensorVec             VoigtTensorVec; 
  typedef pesopt::BoostParser ParameterParserType;

  const ParameterParserType _parser;

public:
  MaterialOptimizationCellMultipleLoadSolverBase ( const ParameterParserType & parser ) :
  _parser( parser ) { } 
  
  virtual void computeOptimalDesign_RefinementStep( const ParameterParserType & parserRefinementStep,
                                                    const MeshType & mesh, VectorType &initMaterial, VectorType &solutionMaterial,
                                                    const int refinementStep = 0, const bool onlyComputeInitialEnergy = false ) const {
    throw std::logic_error( "MaterialOptimizationMultipleLoadSolverBase::computeOptimalDesign(): Unimplemented function! Should be provided in derived class!");
  }
 
  void computeOptimalDesign_MultiLevel( const MeshType &meshStartLevel, 
                                        const int numRefinementSteps, 
                                        const int StartLevel = 0,
                                        const bool onlyComputeInitialEnergy = false ) const {
                                           
    pesopt::consoleOutput( "compute optimal design with multilevel scheme" );
    ConfiguratorType confStartLevel ( meshStartLevel );
    FEMaterialGenerator<ConfiguratorType> FEMaterialGeneratorStartLevel ( _parser, confStartLevel );
    //! initialize material
    VectorType initMaterialStartLevel ( confStartLevel.getNumGlobalDofs() );
    FEMaterialGeneratorStartLevel.switchMaterialType ( initMaterialStartLevel );
    
    std::vector<string> saveDirectoryVec;
    
    FEAdaptiveMesh<MeshType> adaptiveMesh( meshStartLevel );
    VectorType initMaterialCurrentLevel ( initMaterialStartLevel );
    VectorType solMaterialCurrentLevel;
    
    for( int level = StartLevel; level <= numRefinementSteps; ++level ){
        
       auto startTime = std::chrono::high_resolution_clock::now(); 
       if( onlyComputeInitialEnergy ) 
           pesopt::consoleOutput( pesopt::strprintf( "compute initial energy for refinementStep %d", level ).c_str() );
       else                          
           pesopt::consoleOutput( pesopt::strprintf( "compute optimal design for refinementStep %d", level ).c_str() );
        
        pesopt::BoostParser parserRefinementStep; parserRefinementStep = _parser;
        const string saveDirectoryRefinementStep = parserRefinementStep.createSubDirectory( "refinementStep" + std::to_string(level) );
        parserRefinementStep.set ( "saving.saveDirectory", saveDirectoryRefinementStep );
        saveDirectoryVec.push_back( saveDirectoryRefinementStep );
        
        if( level > StartLevel ){  
            
            const ConfiguratorType confOldLevel ( adaptiveMesh.getCurrentMesh() );
//             cout <<  "meshOldLevel numVertices = " <<  confOldLevel.getMesh().getNumVertices() <<  endl;
//             cout <<  "boundary handler in computeOptimalDesign_MultiLevel: start" <<  endl;
            QuocDirichletAndPeriodicBoundaryConditionHandler<ConfiguratorType> bdryHandlerMaterialOldLevel( parserRefinementStep, confOldLevel );
            bdryHandlerMaterialOldLevel.extendVector( solMaterialCurrentLevel );
//             cout <<  "boundary handler in computeOptimalDesign_MultiLevel: end" <<  endl;
            
            //refinement option
            const string adaptiveRefinementType = parserRefinementStep.template get<string> ( "MaterialOptimization.adaptiveRefinementType" );
            if( adaptiveRefinementType == "refineAll" )
                adaptiveMesh.refineAll( );
            else if( adaptiveRefinementType == "refineByFactor" )
                adaptiveMesh.refineByFactor( _parser.template get<RealType> ( "MaterialOptimization.factorMultiLevel" ) );

            else
                throw std::invalid_argument( pesopt::strprintf ( "Wrong adaptiveRefinementType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
            
            adaptiveMesh.template prolongateScalarFunctionAtNodes<ConfiguratorType> ( solMaterialCurrentLevel, initMaterialCurrentLevel );
        }
        
        //! save mesh and parser to file 
        VTKMeshSaver<MeshType> meshSaver( adaptiveMesh.getCurrentMesh() );
        const string fileNameMesh = pesopt::strprintf( "%s/mesh.vtk", saveDirectoryRefinementStep.c_str() );
        meshSaver.save( fileNameMesh, MeshType::_VTKDATATYPEUNDEFORMED );
        adaptiveMesh.getCurrentMesh().saveToParser( parserRefinementStep, fileNameMesh );
        parserRefinementStep.saveToFile( pesopt::strprintf ( "ParameterParser.ini" ) );
        
        //! solve for refinement step
        solMaterialCurrentLevel.resize( initMaterialCurrentLevel.size() );
        this->computeOptimalDesign_RefinementStep ( parserRefinementStep, adaptiveMesh.getCurrentMesh(), initMaterialCurrentLevel, solMaterialCurrentLevel, level, onlyComputeInitialEnergy ); 
       
        //! print elapsed time
        std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
        std::cout << endl << "duration of refinementStep = " << diff.count() << " sec" << endl;
        pesopt::consoleOutput( pesopt::strprintf( "finished computation optimal design for refinementStep %d", level ).c_str() );

    }
    
    const string saveDirectoryFinal = _parser.createSubDirectory( "FinalResult" ); 
    pesopt::moveDirectory( saveDirectoryVec[saveDirectoryVec.size()-1], saveDirectoryFinal );
    if( _parser.template get<bool> ("saving.removeOldRefinementSteps") ){
        for( int refinementStep = 0; refinementStep < saveDirectoryVec.size() - 1; ++refinementStep ) 
            pesopt::deleteDirectory( saveDirectoryVec[refinementStep] );
    }
    
  }
  
private :
    
  void saveInfo( const pesopt::NonlinearEnergyOp<DataTypeContainer> &totalEnergyOp,
                 const pesopt::NonlinearConstraintOps<DataTypeContainer> &constraintOps,
                 const VectorType &material, 
                 ParameterParserType & info
               ) const{

    totalEnergyOp.getEnergyInfo( material, info );
    constraintOps.getConstraintInfo( material, info );
    info.saveToFile( pesopt::strprintf ( "energyInfo.ini" ) );
  }
    
    
  template <typename OptimalDeformSolver>
  void saveResultsToFile( const ParameterParserType & parserRefinementStep,
                          const OptimalDeformSolver &OptDeformFinder,
                          const VectorType &material
                        ) const{
      
    const string saveDirectoryRefinementStep = parserRefinementStep.template get<string> ( "saving.saveDirectory" );
    
    //! save material vector
    pesopt::printVectorToFile<VectorType> ( material, pesopt::strprintf( "%s/materialCollabsed.txt", saveDirectoryRefinementStep.c_str() ).c_str(), 10 ); 
    
    //! save extended material vector
    VectorType materialExtended ( material );
    OptDeformFinder.getBoundaryHandlerMaterial().extendVector( materialExtended );
    pesopt::printVectorToFile<VectorType> ( materialExtended, pesopt::strprintf( "%s/materialExtended.txt", saveDirectoryRefinementStep.c_str() ).c_str(), 10 );
    
    
    //! save deformation vectors 
    const string saveDirectoryDeformation = parserRefinementStep.createSubDirectory( "Deformation" );
    for( int loadIdx=0; loadIdx< OptDeformFinder.getNumLoads(); ++loadIdx ){
        pesopt::printVectorToFile<VectorType> ( OptDeformFinder.getSolutionDisplacementCollabsed(loadIdx), 
                                                pesopt::strprintf( "%s/DisplacementPeriodic_Dir%d.txt", saveDirectoryDeformation.c_str(), loadIdx ).c_str(), 10 ); 
        pesopt::printVectorToFile<VectorType> ( OptDeformFinder.getAffineDisplacement(loadIdx), 
                                                pesopt::strprintf( "%s/DisplacementAffine_Dir%d.txt", saveDirectoryDeformation.c_str(), loadIdx ).c_str(), 10 ); 
    }
  }
    
  void saveResultsVTK( const ParameterParserType & parserRefinementStep ) const{
      
    MaterialOptimizationPeriodicCellVTKSaver<ConfiguratorType> vtkSaver( parserRefinementStep );
    vtkSaver.saveMaterial(  );
    vtkSaver.saveSingleCell( );
    vtkSaver.saveBlockOfCells(  parserRefinementStep.template get<int> ("saving.numBlocks") );
    if( parserRefinementStep.template get<bool> ( "saving.saveDeformations" ) )
        vtkSaver.saveDeformations( );
    if( parserRefinementStep.template get<bool> ("saving.saveInterface") ) 
        vtkSaver.saveInterface( );
    if( parserRefinementStep.template get<bool> ("saving.saveVonMisesOnInterface") )
        vtkSaver.saveVonMisesOnInterface( );
    if( parserRefinementStep.template get<bool>("saving.saveVonMises") ) 
        vtkSaver.saveVonMisesStresses( );

    
  }
 
  //NOTE can only plot to png if it was saved to vtk before, i.e. in function saveResultsVTK
  void plotResultsVTKtoPNG( const ParameterParserType & parserRefinementStep ) const{    
    
    //! plot vtk files to png
    MaterialOptimizationVTKtoPNGPlotter<DataTypeContainer> VTKtoPNGPlotter( parserRefinementStep );
    const int numLoads = parserRefinementStep.template get<int> ("AffineDisp.numLoads");
    
    if( parserRefinementStep.template get<bool> ("saving.plotPNGUndeformed") )
        VTKtoPNGPlotter.plotUndeformedToPNG( );  
    if( parserRefinementStep.template get<bool> ("saving.plotPNGBlockOfCells") )
        VTKtoPNGPlotter.plotBlockOfCellsToPNG( );  
    if( parserRefinementStep.template get<bool> ("saving.plotPNGDeformed") )
        VTKtoPNGPlotter.plotDeformationToPNG( numLoads ); 
    if( parserRefinementStep.template get<bool> ("saving.plotPNGInterface") ){
            VTKtoPNGPlotter.plotInterfaceToPNG( );
            if( parserRefinementStep.template get<bool> ("saving.plotPNGVonMisesOnInterface") )  
                VTKtoPNGPlotter.plotVonMisesStressesOnInterfaceToPNG( numLoads );
    }
    
  }

protected :
  
  template <typename OptimalDeformSolver>
  void saveAndPlotAllResults_RefinementStep( const string name, 
                                             const ParameterParserType & parserRefinementStep,
                                             const pesopt::NonlinearEnergyOp<DataTypeContainer> &totalEnergyOp,
                                             const pesopt::NonlinearConstraintOps<DataTypeContainer> &constraintOps,
                                             const OptimalDeformSolver &OptDeformFinder,
                                             const VectorType &material,
                                             const bool saveVTK,
                                             const bool plotVTKtoPNG
                                           ) const {
 
    pesopt::BoostParser parserMaterial ( parserRefinementStep );
    const string subDirectoryMaterial = parserMaterial.createSubDirectory( name.c_str() );
    parserMaterial.set ( "saving.saveDirectory", subDirectoryMaterial );
    
    pesopt::BoostParser energyInfo; energyInfo.set ( "saving.saveDirectory", subDirectoryMaterial );
    this->saveInfo( totalEnergyOp, constraintOps, material, energyInfo );
   
    this->template saveResultsToFile<OptimalDeformSolver>( parserMaterial, OptDeformFinder, material );
    
    //! save results as vtk
    if( saveVTK ) this->saveResultsVTK( parserMaterial );
    
    //!plot results to png
    if( plotVTKtoPNG ) this->plotResultsVTKtoPNG( parserMaterial );
                                                    
   }
  
      
};



#endif
