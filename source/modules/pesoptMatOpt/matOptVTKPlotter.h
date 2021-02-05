#ifndef __MATOPTVTKPLOTTER_H
#define __MATOPTVTKPLOTTER_H

#include <pesoptImageIO.h>
#include <pesopt_fe.h>
#include <pesopt_quocFE.h>
#include <pesopt_VTK.h>
#include "matOptDefines.h" 




template< typename ConfiguratorType>
class MaterialOptimizationVonMisesStressEvaluator{
    
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::RealVecChart RealVecChart;
  typedef typename ConfiguratorType::PointType PointType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::QuadRuleType QuadRuleType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
  typedef typename ConfiguratorType::DTContainer DataTypeContainer;
  
  PhaseFieldFunctions<RealType> _pfFcts;
  const ConfiguratorType &_conf;
  const MaterialProperties<RealType,ConfiguratorType::dimChartDomain> & _materialProperties; 
  RealType _mu, _lambda;
  const RealType _factorVoidMaterial;
  
public :
    
    MaterialOptimizationVonMisesStressEvaluator ( const ConfiguratorType &conf,
                                                  const MaterialProperties<RealType,ConfiguratorType::dimChartDomain> &materialProperites,
                                                  const RealType factorVoidMaterial
                                                ) :
      _conf ( conf ),
      _materialProperties ( materialProperites ), 
      _mu( _materialProperties.getMu() ), _lambda ( _materialProperties.getLambda() ),
      _factorVoidMaterial( factorVoidMaterial ) {    }
    
    
    void evaluateStress( const string fileName, const VectorType &material, 
                         const VectorType &dispPeriodic, const VectorType &dispAffine, 
                         std::vector<RealType> &stressVec, 
                         const bool isBlock ) const{
        
       PointType lengthSingleCell; 
       for( int i=0; i< lengthSingleCell.size(); ++i ) lengthSingleCell[i] = _conf.getMesh().getWidth(i);
        
       VTKSaver vtkHan;
       std::vector<PointType> pointVec;
       vtkHan.template getPoints<PointType>( fileName, pointVec );
       stressVec.resize ( pointVec.size() );
       
       //affine part
       const FEAffineFunctionEvaluator<DataTypeContainer, ConfiguratorType::dimChartDomain> dispAffineDFD ( dispAffine );
       const RealType div_Affine = dispAffineDFD.getDiv();
       const DerivativeVectorValuedType GradSym_Affine = dispAffineDFD.getSymGrad();
       
       FEScalarFunctionEvaluator<ConfiguratorType> materialDFD ( _conf, material );
       FEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( _conf, dispPeriodic );
       for( int p=0; p<pointVec.size(); ++p ){
           PointType pointVecSingleCell;
           if( isBlock ){
            for( int i=0; i< pointVecSingleCell.size(); ++i ){
                const RealType relativeCoord = pointVec[p][i] / lengthSingleCell[i];
                const int fac = static_cast<int> ( relativeCoord );
                pointVecSingleCell[i] = relativeCoord - fac;
            }
           }else{
             pointVecSingleCell = pointVec[p];   
           }
           int elementNumber; PointType LocalCoord;
           _conf.getLocalCoords ( pointVecSingleCell, elementNumber, LocalCoord );
           const RealType chi = _pfFcts.approxCharFct_material ( materialDFD.evaluate( _conf.getMesh().getElement(elementNumber), LocalCoord )  );
           const RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
           RealType div_Periodic = dispDFD.evaluateDivergence( _conf.getMesh().getElement(elementNumber), LocalCoord );
           DerivativeVectorValuedType GradSym_Periodic;
           dispDFD.evaluateSymmetrizedGradient( _conf.getMesh().getElement(elementNumber), LocalCoord, GradSym_Periodic );
           DerivativeVectorValuedType stress;
           stress.setZero();
           stress += 2. * _mu * ( GradSym_Periodic + GradSym_Affine );
           for( int i=0; i< _conf.dimChartDomain; ++i ) stress(i,i) += _lambda * ( div_Periodic + div_Affine );
           
           RealType vonMisesStress = 0.;
           for( int i=0; i< _conf.dimChartDomain; ++i ){
               vonMisesStress += pesopt::Sqr( stress(i,i) );
               for( int j=i+1; j< _conf.dimChartDomain; ++j ){
                     vonMisesStress -= stress(i,i) * stress(j,j);
                     vonMisesStress += 3 * pesopt::Sqr( stress(i,j) );
               }
           }
           // TODO use materialfactor???
           // stressVec[p] = materialfactor * std::sqrt( vonMisesStress );
           stressVec[p] = std::sqrt( vonMisesStress );
       }
    }
    
};




template< typename DataTypeContainer >
class MaterialOptimizationVTKtoPNGPlotter : 
public VTKtoPNGPlotter {
private:
  
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::PointType PointType;
  typedef pesopt::BoostParser ParameterParserType;
  
  const ParameterParserType &_parser;
  const string _saveDirectory;
  
public :
    
    MaterialOptimizationVTKtoPNGPlotter ( const ParameterParserType & parser ) :
    VTKtoPNGPlotter( parser.template get<string>( "saving.saveDirectory" ).c_str() ),
    _parser ( parser ),
    _saveDirectory( parser.template get<string>( "saving.saveDirectory" ).c_str() ) { }

    
   void plotUndeformedToPNG( ) const{
       pesopt::consoleOutputItem( "plot undeformed" );
       pesopt::BoostParser parserVTK ( _parser.template get<string> ( "saving.parserVTK" ) );
       ImageHandler imHandler;
       
       const bool useScalarData = true;
       this->template plotToPngWithParserInfo<pesopt::BoostParser>( 
                           pesopt::strprintf ( "%s/material.vtk", _saveDirectory.c_str() ),  
                           pesopt::strprintf ( "%s/_tmp_material.png", _saveDirectory.c_str() ).c_str(),
                           parserVTK, useScalarData, "material", VERTEX_DATA );
       
       imHandler.trimAndResizeImage( pesopt::strprintf ( "%s/_tmp_material.png", _saveDirectory.c_str() ).c_str(),
                                     pesopt::strprintf ( "%s/material.png", _saveDirectory.c_str() ).c_str(),
                                     parserVTK.template get<int> ("Window.resizeWidth"), parserVTK.template get<int> ("Window.resizeHeight")  );
       
       pesopt::deleteFile( "_tmp_material.png", _saveDirectory );
       
    }
    
    void plotBlockOfCellsToPNG( ) const{
       pesopt::consoleOutputItem( "plot block" );
       pesopt::BoostParser parserVTK ( _parser.template get<string> ( "saving.parserVTK" ) );
       ImageHandler imHandler;
       
       const bool useScalarData = true;
       this->template plotToPngWithParserInfo<pesopt::BoostParser>( 
                           pesopt::strprintf ( "%s/materialBlock.vtk", _saveDirectory.c_str() ),  
                           pesopt::strprintf ( "%s/_tmp_materialBlock.png", _saveDirectory.c_str() ).c_str(),
                           parserVTK,
                           useScalarData, "material", VERTEX_DATA );
       
       imHandler.trimAndResizeImage( pesopt::strprintf ( "%s/_tmp_materialBlock.png", _saveDirectory.c_str() ).c_str(),
                                     pesopt::strprintf ( "%s/materialBlock.png", _saveDirectory.c_str() ).c_str(),
                                     parserVTK.template get<int> ("Window.resizeWidth"), parserVTK.template get<int> ("Window.resizeHeight")  );
       
       pesopt::deleteFile( "_tmp_materialBlock.png", _saveDirectory );
       
    }
    
    
    void plotDeformationToPNG( const int numLoads ) const{
       pesopt::consoleOutputItem( "plot deformations" );
       const bool useScalarData = true;
       pesopt::BoostParser parserVTK ( _parser.template get<string> ( "saving.parserVTK" ) );
       for( int loadIdx=0; loadIdx< numLoads; ++loadIdx ){
        this->template plotToPngWithParserInfo<pesopt::BoostParser>( 
            pesopt::strprintf ( "%s/Deformation/Deformation_Dir%d.vtk", _saveDirectory.c_str(), loadIdx ), 
            pesopt::strprintf ( "%s/Deformation/Deformation_Dir%d.png", _saveDirectory.c_str(), loadIdx ).c_str(),
            parserVTK, useScalarData, "material", VERTEX_DATA);
       }
       
    }
    
    void plotInterfaceToPNG( ) const{
       pesopt::consoleOutputItem( "plot interface" );
       pesopt::BoostParser parserVTK ( _parser.template get<string> ( "saving.parserVTK" ) );
       ImageHandler imHandler;
       
       const int numSubdivLevels = _parser.template get<RealType> ( "saving.numSubdivLevelsInterface" );
       string subdivString = "";
       //TODO  maybe do subdivision here
       if( numSubdivLevels > 0 )
           subdivString = pesopt::strprintf( "_Subdiv%d", numSubdivLevels );
       
       this->template plotToPngWithParserInfo<pesopt::BoostParser>( 
                               pesopt::strprintf ( "%s/Interface/Interface%s.vtk", _saveDirectory.c_str(), subdivString.c_str() ), 
                               pesopt::strprintf ( "%s/Interface/_tmp_Interface.png", _saveDirectory.c_str() ).c_str(),
                               parserVTK );
       imHandler.trimAndResizeImage( pesopt::strprintf ( "%s/Interface/_tmp_Interface.png", _saveDirectory.c_str() ).c_str(),
                                     pesopt::strprintf ( "%s/Interface/Interface%s.png", _saveDirectory.c_str(), subdivString.c_str() ).c_str(),
                                     parserVTK.template get<int> ("Window.resizeWidth"), parserVTK.template get<int> ("Window.resizeHeight")  );
     
       pesopt::deleteFile( "_tmp_Interface.png", pesopt::strprintf("%s/Interface", _saveDirectory.c_str() ) );
       
       this->template plotToPngWithParserInfo<pesopt::BoostParser>( 
                               pesopt::strprintf ( "%s/Interface/BlockInterface%s.vtk", _saveDirectory.c_str(), subdivString.c_str() ), 
                               pesopt::strprintf ( "%s/Interface/_tmp_BlockInterface.png", _saveDirectory.c_str() ).c_str(),
                               parserVTK );
       imHandler.trimAndResizeImage( pesopt::strprintf ( "%s/Interface/_tmp_BlockInterface.png", _saveDirectory.c_str() ).c_str(),
                                     pesopt::strprintf ( "%s/Interface/BlockInterface%s.png", _saveDirectory.c_str(), subdivString.c_str() ).c_str(),
                                     parserVTK.template get<int> ("Window.resizeWidth"), parserVTK.template get<int> ("Window.resizeHeight") );
       
       pesopt::deleteFile( "_tmp_BlockInterface.png", pesopt::strprintf("%s/Interface", _saveDirectory.c_str() ) );
       
       //delete subdivision file 
       if( numSubdivLevels > 0 ){
           pesopt::deleteFile( pesopt::strprintf("Interface%s.vtk", subdivString.c_str()  ), pesopt::strprintf("%s/Interface", _saveDirectory.c_str() ) );
           pesopt::deleteFile( pesopt::strprintf("BlockInterface%s.vtk", subdivString.c_str()  ), pesopt::strprintf("%s/Interface", _saveDirectory.c_str() ) );
       }
    }
    
    
    void plotVonMisesStressesOnInterfaceToPNG( const int numLoads ) const{
       pesopt::consoleOutputItem( "plot von mises on interface" );
       pesopt::BoostParser parserVTKStress ( _parser.template get<string> ( "saving.parserVTKStress" ) );
       ImageHandler imHandler;
       for( int loadIdx=0; loadIdx<numLoads; ++loadIdx ){
          this->template  plotToPngWithParserInfo<pesopt::BoostParser>( 
            pesopt::strprintf( "%s/StressOnInterface/VonMises_Interface_Direction%d.vtp", _saveDirectory.c_str(), loadIdx ).c_str(), 
            pesopt::strprintf( "%s/StressOnInterface/_tmp_VonMises_Interface_Direction%d.png", _saveDirectory.c_str(), loadIdx ).c_str(),  parserVTKStress, true, "VonMisesStresses",  VERTEX_DATA );
          imHandler.trimAndResizeImage( 
            pesopt::strprintf( "%s/StressOnInterface/_tmp_VonMises_Interface_Direction%d.png", _saveDirectory.c_str(), loadIdx ).c_str(), 
            pesopt::strprintf( "%s/StressOnInterface/VonMises_Interface_Direction%d.png", _saveDirectory.c_str(), loadIdx ).c_str(), 
            parserVTKStress.template get<int> ("Window.resizeWidth"), parserVTKStress.template get<int> ("Window.resizeHeight")  );
          pesopt::deleteFile( pesopt::strprintf( "_tmp_VonMises_Interface_Direction%d.png", loadIdx ), pesopt::strprintf("%s/StressOnInterface", _saveDirectory.c_str() ) );
          
          this->template plotToPngWithParserInfo<pesopt::BoostParser>( 
                               pesopt::strprintf( "%s/StressOnInterface/VonMises_BlockInterface_Direction%d.vtp", _saveDirectory.c_str(), loadIdx ).c_str(), 
                               pesopt::strprintf( "%s/StressOnInterface/_tmp_VonMises_BlockInterface_Direction%d.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                               parserVTKStress, true, "VonMisesStresses",  VERTEX_DATA  );
          imHandler.trimAndResizeImage( 
            pesopt::strprintf( "%s/StressOnInterface/_tmp_VonMises_BlockInterface_Direction%d.png", _saveDirectory.c_str(), loadIdx ).c_str(),
            pesopt::strprintf( "%s/StressOnInterface/VonMises_BlockInterface_Direction%d.png", _saveDirectory.c_str(), loadIdx ).c_str(),
            parserVTKStress.template get<int> ("Window.resizeWidth"), parserVTKStress.template get<int> ("Window.resizeHeight")  );
          pesopt::deleteFile( pesopt::strprintf( "_tmp_VonMises_BlockInterface_Direction%d.png", loadIdx ), pesopt::strprintf("%s/StressOnInterface", _saveDirectory.c_str() ) );
        }
    }

};




// TODO replace by VTKDeformedMeshSaver
template< typename ConfiguratorType >
class VTKFEMeshSaver{
  
public:
  
  typedef typename ConfiguratorType::DTContainer            DataTypeContainer;
  typedef typename ConfiguratorType::RealType               RealType;
  typedef typename ConfiguratorType::InitType               MeshType;
  typedef typename ConfiguratorType::PointType              PointType;
  typedef typename ConfiguratorType::VectorType             VectorType;
  typedef pesopt::BoostParser ParameterParserType;
  
  const MeshType &_mesh;
  const int _numGlobalDofs;
  
public:
  
  VTKFEMeshSaver( const ConfiguratorType &conf ) : 
  _mesh( conf.getMesh() ), _numGlobalDofs ( conf.getNumGlobalDofs() )  { }
 
protected:
 
 void getDeformedMesh(  const VectorType &disp, MeshType& meshDeformed ) const {
    for( int i = 0; i < _mesh.getNumVertices(); ++i ){
        PointType coords = _mesh.getVertex( i );
        for( int comp = 0; comp < coords.size(); ++comp ) coords[comp] += disp[i + comp * _numGlobalDofs];
        meshDeformed.setVertex( i, coords );
    }
 }
    
public:
    
  void saveUndeformed ( const string outputFileNameVTK  ) const{
    VTKMeshSaver<MeshType> meshSaver ( _mesh );
    meshSaver.save( outputFileNameVTK, _mesh._VTKDATATYPEUNDEFORMED );
  }
 
  void saveDeformed ( const VectorType &disp, const string outputFileNameVTK ) const{
    MeshType meshDeformed ( _mesh );
    this->getDeformedMesh( disp, meshDeformed );
    VTKMeshSaver<MeshType> meshSaver ( meshDeformed );
    meshSaver.save( outputFileNameVTK, _mesh._VTKDATATYPEDEFORMED );
  }

};




template< typename ConfiguratorType >
class MaterialOptimizationVTKSaver
  :  public VTKFEMeshSaver<ConfiguratorType>{
  
public:
  
  typedef typename ConfiguratorType::DTContainer            DataTypeContainer;
  typedef typename ConfiguratorType::RealType               RealType;
  typedef typename ConfiguratorType::InitType               MeshType;
  typedef typename ConfiguratorType::PointType              PointType;
  typedef typename ConfiguratorType::VectorType             VectorType;
  typedef pesopt::BoostParser ParameterParserType;
  
  const MeshType &_mesh;
  const int _numGlobalDofs;
  
public:
  
  MaterialOptimizationVTKSaver( const ConfiguratorType &conf ) : 
  VTKFEMeshSaver<ConfiguratorType> ( conf ), 
  _mesh( conf.getMesh() ), _numGlobalDofs ( conf.getNumGlobalDofs() )  { }
    
public:

  void saveUndeformedWithMaterial ( const VectorType & material, 
                                    const string outputFileNameVTK,
                                    VTKDataSupp dataSupp = VERTEX_DATA ) const{
    VTKMeshSaver<MeshType> meshSaver ( _mesh );
    meshSaver.addScalarData ( material, "material", dataSupp );
    meshSaver.save( outputFileNameVTK, _mesh._VTKDATATYPEUNDEFORMED );
  }
 
  void saveDeformedWithMaterial ( const VectorType &disp, const VectorType & material, 
                                  const string outputFileNameVTK,
                                  VTKDataSupp dataSupp = VERTEX_DATA ) const{
    MeshType meshDeformed ( _mesh );
    this->getDeformedMesh( disp, meshDeformed );
    VTKMeshSaver<MeshType> meshSaver ( meshDeformed );
    meshSaver.addScalarData ( material, "material", dataSupp );
    meshSaver.save( outputFileNameVTK, _mesh._VTKDATATYPEDEFORMED );
  }

};



template< typename ConfiguratorType >
class MaterialOptimizationMacroVTKSaver {
    
private:
  
  typedef typename ConfiguratorType::RealType                               RealType;
  typedef typename ConfiguratorType::VectorType                             VectorType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType             DerivativeVectorValuedType;
  typedef typename ConfiguratorType::DTContainer                            DataTypeContainer;
  typedef typename DataTypeContainer::IntVecChart                           IntVecChart;
  typedef pesopt::BoostParser ParameterParserType;
  typedef typename DataTypeContainer::PointType                             PointType;
  typedef typename ConfiguratorType::SparseMatrixType                       SparseMatrixType;
  typedef typename ConfiguratorType::MaskType                               MaskType;
  typedef typename ConfiguratorType::InitType                               MeshType;
  
  const ParameterParserType &_parser;
  const string _saveDirectory;
  
  const MeshType _mesh;
  const ConfiguratorType _conf;
  
public :
    
    MaterialOptimizationMacroVTKSaver ( const ParameterParserType & parser ) :
    _parser ( parser ), _saveDirectory( parser.template get<string>( "saving.saveDirectory" ).c_str() ),
    _mesh ( generateMeshFromParser<MeshType>( parser ) ),
    _conf ( _mesh )   {  }
        
    void saveMaterial(  ) const{
        pesopt::consoleOutputItem( "save Material " );
            
        VectorType material ( _mesh.getNumVertices() );
        pesopt::loadVectorFromFile<VectorType>( material, pesopt::strprintf ( "%s/materialExtended.txt", _saveDirectory.c_str() ) ); 
            
        MaterialOptimizationVTKSaver<ConfiguratorType> meshSaver( _conf );
    
        meshSaver.saveUndeformedWithMaterial ( material, pesopt::strprintf( "%s/material.vtk", _saveDirectory.c_str() ).c_str() );
    }
        
    //=====================================================================================================================
    //save deformation
    //=====================================================================================================================
    void saveDeformations( ) const{        
        
        pesopt::consoleOutputItem( "save deformations" );
            
        MaterialOptimizationVTKSaver<ConfiguratorType> meshSaver( _conf );
            
        //TODO
        const int numLoads = 1;
        
        const string saveDirectoryDeformation = _parser.createSubDirectory( "Deformation" );
            
        for( int loadIdx=0; loadIdx < numLoads; ++loadIdx ){

            VectorType dispAtNodes ( _conf.dimChartDomain * _mesh.getNumVertices() );
            pesopt::loadVectorFromFile<VectorType>( dispAtNodes, pesopt::strprintf( "%s/Displacement_Dir%d.txt", saveDirectoryDeformation.c_str(), loadIdx ) ); 
            
            meshSaver.saveDeformed( dispAtNodes, pesopt::strprintf( "%s/Deformation_Dir%d.vtk", saveDirectoryDeformation.c_str(), loadIdx ).c_str() );
        }
    }
        
};



template< typename ConfiguratorType >
class MaterialOptimizationFineScaleVTKSaver {
    
private:
  
  typedef typename ConfiguratorType::RealType                               RealType;
  typedef typename ConfiguratorType::VectorType                             VectorType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType             DerivativeVectorValuedType;
  typedef typename ConfiguratorType::DTContainer                            DataTypeContainer;
  typedef typename DataTypeContainer::IntVecChart                           IntVecChart;              
  typedef pesopt::BoostParser ParameterParserType;
  typedef typename DataTypeContainer::PointType                             PointType;
  typedef typename ConfiguratorType::SparseMatrixType                       SparseMatrixType;
  typedef typename ConfiguratorType::MaskType                               MaskType;
  typedef typename ConfiguratorType::InitType                               MeshType;
  
  const ParameterParserType &_parser;
  const string _saveDirectory;
  
  const MeshType _mesh;
  const ConfiguratorType _conf;
  
public :
    
    MaterialOptimizationFineScaleVTKSaver ( const ParameterParserType & parser ) :
    _parser ( parser ), _saveDirectory( parser.template get<string>( "saving.saveDirectory" ).c_str() ),
    _mesh ( generateMeshFromParser<MeshType>( parser ) ),
    _conf ( _mesh )   {  }
        
    void saveMaterial(  ) const{
        pesopt::consoleOutputItem( "save Material " );
            
        VectorType material ( _mesh.getNumVertices() );
        pesopt::loadVectorFromFile<VectorType>( material, pesopt::strprintf ( "%s/materialExtended.txt", _saveDirectory.c_str() ) ); 
            
        MaterialOptimizationVTKSaver<ConfiguratorType> meshSaver( _conf );
    
        meshSaver.saveUndeformedWithMaterial ( material, pesopt::strprintf( "%s/material.vtk", _saveDirectory.c_str() ).c_str() );
    }
        
    void saveDeformations( ) const{        
        
        pesopt::consoleOutputItem( "save deformations" );
            
        MaterialOptimizationVTKSaver<ConfiguratorType> meshSaver( _conf );
            
        const int numLoads = _parser.template get<int> ( "AffineDisp.numLoads" );
        
        const string saveDirectoryDeformation = _parser.createSubDirectory( "Deformation" );
            
        for( int loadIdx=0; loadIdx < numLoads; ++loadIdx ){

            VectorType dispAtNodes ( _conf.dimChartDomain * _mesh.getNumVertices() );
            pesopt::loadVectorFromFile<VectorType>( dispAtNodes, pesopt::strprintf( "%s/Displacement_Dir%d.txt", saveDirectoryDeformation.c_str(), loadIdx ) ); 
            
            meshSaver.saveDeformed( dispAtNodes, pesopt::strprintf( "%s/Deformation_Dir%d.vtk", saveDirectoryDeformation.c_str(), loadIdx ).c_str() );
        }
    }
    
    
    void saveBoundaryMask(  ) const{
        pesopt::consoleOutputItem( "save BoundaryMask " );
            
        VectorType boundaryMaterial ( _mesh.getNumVertices() );
        pesopt::loadVectorFromFile<VectorType>( boundaryMaterial, pesopt::strprintf ( "%s/DirichletBoundaryMaterial.txt", _saveDirectory.c_str() ) ); 
        
        VectorType boundaryDisplacement ( _mesh.getNumVertices() );
        pesopt::loadVectorFromFile<VectorType>( boundaryDisplacement, pesopt::strprintf ( "%s/DirichletBoundaryDisplacement.txt", _saveDirectory.c_str() ) ); 
            
        MaterialOptimizationVTKSaver<ConfiguratorType> meshSaver( _conf );
        meshSaver.saveUndeformedWithMaterial ( boundaryMaterial, pesopt::strprintf( "%s/DirichletBoundaryMaterial.vtk", _saveDirectory.c_str() ).c_str() );
        meshSaver.saveUndeformedWithMaterial ( boundaryDisplacement, pesopt::strprintf( "%s/DirichletBoundaryDisplacement.vtk", _saveDirectory.c_str() ).c_str() );
    }
        
};


//TODO more general BoundaryConditionHandler
template< typename ConfiguratorType >
class MaterialOptimizationPeriodicCellVTKSaver {
    
private:
  
  typedef typename ConfiguratorType::RealType                               RealType;
  typedef typename ConfiguratorType::VectorType                             VectorType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType             DerivativeVectorValuedType;
  typedef typename ConfiguratorType::DTContainer                            DataTypeContainer;
  typedef typename DataTypeContainer::IntVecChart                           IntVecChart;
  typedef pesopt::BoostParser ParameterParserType;
  typedef typename DataTypeContainer::PointType                             PointType;
  typedef typename ConfiguratorType::SparseMatrixType                       SparseMatrixType;
  typedef typename ConfiguratorType::MaskType                               MaskType;
  typedef typename ConfiguratorType::InitType                               MeshType;
  
  const ParameterParserType &_parser;
  const string _saveDirectory;
  
  IntVecChart _numDofVec; 
  PointType _lengthVec;  
  const int _numDofsOffsetOutside = 3;
  
public :
    
    MaterialOptimizationPeriodicCellVTKSaver ( const ParameterParserType & parser ) :
    _parser ( parser ), _saveDirectory( parser.template get<string>( "saving.saveDirectory" ).c_str() ),
    _numDofsOffsetOutside( 3 ) {
        _parser.template getFixSizeVector<int,IntVecChart> ("InputMesh.NumDofVec", _numDofVec );
        _parser.template getFixSizeVector<RealType, PointType> ("InputMesh.LengthVec", _lengthVec );
    }
        
    void saveMaterial(  ) const{
        pesopt::consoleOutputItem( "save Material " );
        
        MeshType mesh ( _numDofVec, _lengthVec );
        ConfiguratorType conf ( mesh );
            
        VectorType material ( mesh.getNumVertices() );
        pesopt::loadVectorFromFile<VectorType>( material, pesopt::strprintf ( "%s/materialExtended.txt", _saveDirectory.c_str() ) ); 
            
        MaterialOptimizationVTKSaver<ConfiguratorType> meshSaver( conf );
    
        meshSaver.saveUndeformedWithMaterial ( material, pesopt::strprintf( "%s/material.vtk", _saveDirectory.c_str() ).c_str() );
    }
        
    //! ========================================================================================
    //! plot material (B or P) around single micro cell
    //! ========================================================================================
    void saveSingleCell( ) const{
            pesopt::consoleOutputItem(  "save single cell" );
            MeshType mesh ( _numDofVec, _lengthVec );
            ConfiguratorType conf ( mesh );
            
            VectorType material ( mesh.getNumVertices() );
            pesopt::loadVectorFromFile<VectorType>( material, pesopt::strprintf ( "%s/materialExtended.txt", _saveDirectory.c_str()  ) ); 
            FEScalarFunctionEvaluator<ConfiguratorType> discreteFctSingleCell ( conf, material );
        
            PointType offsetOutside; for( int i=0; i<offsetOutside.size(); ++i ) offsetOutside[i] = _numDofsOffsetOutside * mesh.getMeshSize( i );
            
            pesopt::consoleOutputItem(  "save single cell with material outside" );
            IntVecChart numDofVecOutside;
            PointType lengthVecOutside;
            for( int i=0; i<_numDofVec.size(); ++i ){
                numDofVecOutside[i] = _numDofVec[i] + 2 * _numDofsOffsetOutside;
                lengthVecOutside[i] = _lengthVec[i] + 2 * _numDofsOffsetOutside * mesh.getMeshSize( i ); 
            }
            MeshType meshOutside ( numDofVecOutside, lengthVecOutside );
            VectorType materialOutsideSoft ( meshOutside.getNumVertices() );
            //TODO read lower bound for phase field from parser
            for( int i=0; i<materialOutsideSoft.size(); ++i ) materialOutsideSoft[i] = 0.0;
            for( int nodeIdxOffset=0; nodeIdxOffset < meshOutside.getNumVertices(); nodeIdxOffset++ ){
                const PointType& GlobalCoordsOffset = meshOutside.getVertex ( nodeIdxOffset );
               
                bool inside = true;
                for( int i=0; i<GlobalCoordsOffset.size(); ++i ){
                  if( GlobalCoordsOffset[i] < offsetOutside[i] ) inside = false;
                  if( GlobalCoordsOffset[i] > _lengthVec[i] + offsetOutside[i] ) inside = false;
                }
                
                if( inside ) {
                    PointType GlobalCoordsSingleCell;
                    for( int i=0; i<GlobalCoordsSingleCell.size(); ++i){
                        GlobalCoordsSingleCell[i] = GlobalCoordsOffset[i] - _numDofsOffsetOutside * mesh.getMeshSize(i);
                    }
                    int elementNumberSingleCell; PointType LocalCoordSingleCell;
                    conf.getLocalCoords ( GlobalCoordsSingleCell, elementNumberSingleCell, LocalCoordSingleCell );
                    materialOutsideSoft[nodeIdxOffset] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
                }
            }
            
            VTKMeshSaver<MeshType> meshSaverOutsideSoft ( meshOutside );
            meshSaverOutsideSoft.addScalarData ( materialOutsideSoft, "material", VERTEX_DATA );
            meshSaverOutsideSoft.save ( pesopt::strprintf( "%s/materialOutsideSoft.vtk", _saveDirectory.c_str() ), offsetOutside, meshOutside._VTKDATATYPEUNDEFORMED );
    }
        
        
    //=====================================================================================================================
    //Plot block of mesh with material 
    //=====================================================================================================================
    void saveBlockOfCells( const int numBlocksPerDirection ) const{
            pesopt::consoleOutputItem( "save block of cells" );
            
            MeshType mesh ( _numDofVec, _lengthVec );
            ConfiguratorType conf ( mesh );
            
            VectorType material ( mesh.getNumVertices() );
            pesopt::loadVectorFromFile<VectorType>( material, pesopt::strprintf ( "%s/materialExtended.txt", _saveDirectory.c_str() ) ); 
            FEScalarFunctionEvaluator<ConfiguratorType> discreteFctSingleCell ( conf, material );
        
            PointType offsetOutside; for( int i=0; i<offsetOutside.size(); ++i ) offsetOutside[i] = _numDofsOffsetOutside * mesh.getMeshSize( i );
            
            //! ========================================================================================
            //! plot block
            //! ========================================================================================
            IntVecChart _numDofVecBlock;
            PointType _lengthVecBlock;
            for( int i=0; i<_numDofVec.size(); ++i ){
                int oldSize = _numDofVec[i];
                _numDofVecBlock[i] = static_cast<int> ( numBlocksPerDirection * (oldSize - 1) + 1 );
                _lengthVecBlock[i] = static_cast<RealType> ( numBlocksPerDirection ) * _lengthVec[i]; 
            }
            MeshType meshBlock ( _numDofVecBlock, _lengthVecBlock );
            VectorType materialBlock ( meshBlock.getNumVertices() );
            for( int nodeIdxBlock=0; nodeIdxBlock < meshBlock.getNumVertices(); nodeIdxBlock++ ){
                const PointType& GlobalCoordsBlock = meshBlock.getVertex ( nodeIdxBlock );
                PointType GlobalCoordsSingleCell;
                for( int i=0; i<GlobalCoordsSingleCell.size(); ++i){
                    int multDirection = static_cast<int> ( GlobalCoordsBlock[i] / _lengthVec[i] );
                    GlobalCoordsSingleCell[i] = GlobalCoordsBlock[i] - multDirection * _lengthVec[i];
                }
                int elementNumberSingleCell; PointType LocalCoordSingleCell;
                conf.getLocalCoords ( GlobalCoordsSingleCell, elementNumberSingleCell, LocalCoordSingleCell );
                materialBlock[nodeIdxBlock] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
            }
            VTKMeshSaver<MeshType> meshSaver ( meshBlock );
            meshSaver.addScalarData ( materialBlock, "material", VERTEX_DATA );
            meshSaver.save ( pesopt::strprintf( "%s/materialBlock.vtk", _saveDirectory.c_str() ), meshBlock._VTKDATATYPEUNDEFORMED );
            
            //! ========================================================================================
            //! plot material around block of micro cells
            //! ========================================================================================
            IntVecChart _numDofVecBlockOutside;
            PointType _lengthVecBlockOutside;
            for( int i=0; i<_numDofVec.size(); ++i ){
                _numDofVecBlockOutside[i] = _numDofVecBlock[i] + 2 * _numDofsOffsetOutside;
                _lengthVecBlockOutside[i] = _lengthVecBlock[i] + 2 * _numDofsOffsetOutside * mesh.getMeshSize( i ); 
            }
            MeshType meshBlockOutside ( _numDofVecBlockOutside, _lengthVecBlockOutside );
            VectorType materialBlockOutsideSoft ( meshBlockOutside.getNumVertices() );
            //TODO read lower bound for phase field from parser
            for( int i=0; i<materialBlockOutsideSoft.size(); ++i ) materialBlockOutsideSoft[i] = 0.0;
            for( int nodeIdxBlock=0; nodeIdxBlock < meshBlockOutside.getNumVertices(); nodeIdxBlock++ ){
                const PointType& GlobalCoordsBlockOffset = meshBlockOutside.getVertex ( nodeIdxBlock );
                
                bool inside = true;
                for( int i=0; i<GlobalCoordsBlockOffset.size(); ++i ){
                  if( GlobalCoordsBlockOffset[i] < offsetOutside[i] ) inside = false;
                  if( GlobalCoordsBlockOffset[i] > _lengthVecBlock[i] + offsetOutside[i] ) inside = false;
                }
                
                if( inside ) {
                    
                    PointType GlobalCoordsBlock;
                    for( int i=0; i<GlobalCoordsBlock.size(); ++i){
                        GlobalCoordsBlock[i] = GlobalCoordsBlockOffset[i] - _numDofsOffsetOutside * mesh.getMeshSize(i);
                    }
                    PointType GlobalCoordsSingleCell;
                    for( int i=0; i<GlobalCoordsSingleCell.size(); ++i){
                      int multDirection = static_cast<int> ( GlobalCoordsBlock[i] / _lengthVec[i] );
                      GlobalCoordsSingleCell[i] = GlobalCoordsBlock[i] - multDirection * _lengthVec[i];
                    }
                    int elementNumberSingleCell; PointType LocalCoordSingleCell;
                    conf.getLocalCoords ( GlobalCoordsSingleCell, elementNumberSingleCell, LocalCoordSingleCell );
                    materialBlockOutsideSoft[nodeIdxBlock] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
                }
            }
            
            VTKMeshSaver<MeshType> meshSaverBlockOutsideSoft ( meshBlockOutside );
            meshSaverBlockOutsideSoft.addScalarData ( materialBlockOutsideSoft, "material", VERTEX_DATA );
            meshSaverBlockOutsideSoft.save ( pesopt::strprintf( "%s/materialBlockOutsideSoft.vtk", _saveDirectory.c_str() ), offsetOutside, meshBlockOutside._VTKDATATYPEUNDEFORMED );
        }
        
        //=====================================================================================================================
        //save interface
        //=====================================================================================================================
        void saveDeformations( ) const{
        
            pesopt::consoleOutputItem( "save deformations" );
            
            MeshType mesh ( _numDofVec, _lengthVec );
            ConfiguratorType conf ( mesh );
            QuocPeriodicBoundaryConditionHandler<ConfiguratorType> bdryHandlerDisplacement( _parser, conf );
            
            MaterialOptimizationVTKSaver<ConfiguratorType> quocMeshSaver( conf );
            
            const int numLoads = _parser.template get<int> ( "AffineDisp.numLoads" );
            const int numAffineSymGradDofs = conf.numAffineSymGradDofs;
            
            const string saveDirectoryDeformation = _parser.createSubDirectory( "Deformation" );
            
            for( int loadIdx=0; loadIdx < numLoads; ++loadIdx ){

                VectorType dispAffineAtNodes ( conf.dimChartDomain * mesh.getNumVertices() );
                    
                VectorType displacementAffine ( numAffineSymGradDofs );
                pesopt::loadVectorFromFile<VectorType>( displacementAffine, pesopt::strprintf( "%s/DisplacementAffine_Dir%d.txt", saveDirectoryDeformation.c_str(), loadIdx ) ); 
                
                VectorType displacementPeriodic ( conf.dimChartDomain * conf.getNumGlobalDofs() ), displacementPeriodicExtended ( conf.dimChartDomain * conf.getNumGlobalDofs() );
                pesopt::loadVectorFromFile<VectorType>( displacementPeriodic, pesopt::strprintf( "%s/DisplacementPeriodic_Dir%d.txt", saveDirectoryDeformation.c_str(), loadIdx ) ); 
                bdryHandlerDisplacement.extendMultiVector( displacementPeriodic, displacementPeriodicExtended );
                
                FEAffineFunctionEvaluator<DataTypeContainer,ConfiguratorType::dimChartDomain> ( displacementAffine ).template getAffineDisplacementAtNodes<MeshType>( mesh, dispAffineAtNodes );
                quocMeshSaver.saveDeformed( displacementPeriodicExtended + dispAffineAtNodes, pesopt::strprintf( "%s/Deformation_Dir%d_AffinePeriodic.vtk", saveDirectoryDeformation.c_str(), loadIdx ).c_str() );
            }
        }
        
        
   //=====================================================================================================================
   //plot interface
   //=====================================================================================================================
   void saveInterface ( ) const {
        pesopt::consoleOutputItem( "save interface" );
        _parser.createSubDirectory("Interface");
            
        MeshType mesh ( _numDofVec, _lengthVec );
        ConfiguratorType conf ( mesh );
            
        VectorType material ( mesh.getNumVertices() );
        pesopt::loadVectorFromFile<VectorType>( material, pesopt::strprintf ( "%s/materialExtended.txt", _saveDirectory.c_str() ) ); 
        FEScalarFunctionEvaluator<ConfiguratorType> discreteFctSingleCell ( conf, material );
        
        PointType offsetOutside; for( int i=0; i<offsetOutside.size(); ++i ) offsetOutside[i] = _numDofsOffsetOutside * mesh.getMeshSize( i );
            
        //extract surface
        VTKInterfaceExtractor surfaceExtractor;
        const RealType threshold = _parser.template get<RealType> ( "saving.thresholdInterface" );
            
        surfaceExtractor.getInterfaceByMarchingCubes( pesopt::strprintf( "%s/materialOutsideSoft.vtk", _saveDirectory.c_str() ),
                                                          pesopt::strprintf( "%s/Interface/Interface.vtk", _saveDirectory.c_str() ), 
                                                          "material", threshold );

        surfaceExtractor.getInterfaceByMarchingCubes( pesopt::strprintf( "%s/materialBlockOutsideSoft.vtk", _saveDirectory.c_str() ),
                                                          pesopt::strprintf( "%s/Interface/BlockInterface.vtk", _saveDirectory.c_str() ), "material", threshold );
            
        //use loop subdivision for interfaces
        const int numSubdivLevels = _parser.template get<RealType> ( "saving.numSubdivLevelsInterface" );
        if( numSubdivLevels > 0 ){
            surfaceExtractor.loopSubdivision( pesopt::strprintf( "%s/Interface/Interface.vtk", _saveDirectory.c_str() ),
                                              pesopt::strprintf( "%s/Interface/Interface_Subdiv%d.vtk", _saveDirectory.c_str(), numSubdivLevels ), numSubdivLevels );
            surfaceExtractor.loopSubdivision( pesopt::strprintf( "%s/Interface/BlockInterface.vtk", _saveDirectory.c_str() ),
                                              pesopt::strprintf( "%s/Interface/BlockInterface_Subdiv%d.vtk", _saveDirectory.c_str(), numSubdivLevels ), numSubdivLevels );
        } 
            
    }
    
    
    
    //! ========================================================================================
    //! save von mises stresses on interface
    //! ========================================================================================
    void saveVonMisesOnInterface ( ) const {
            pesopt::consoleOutputItem( "save von mises stresses on interface" );
            _parser.createSubDirectory("StressOnInterface");
                
            MeshType mesh ( _numDofVec, _lengthVec );
            ConfiguratorType conf ( mesh );
            QuocPeriodicBoundaryConditionHandler<ConfiguratorType> bdryHandlerDisplacement( _parser, conf );
            
            VectorType material ( mesh.getNumVertices() );
            pesopt::loadVectorFromFile<VectorType>( material, pesopt::strprintf ( "%s/materialExtended.txt", _saveDirectory.c_str() ) ); 
            FEScalarFunctionEvaluator<ConfiguratorType> discreteFctSingleCell ( conf, material );
            
            
            //read material properties
            MaterialProperties<RealType,ConfiguratorType::dimChartDomain> HardMaterial ( "Hard", 1.0, 1.0, _parser.template get<RealType> ( "Material.ElastModulus" ), _parser.template get<RealType> ( "Material.PoissonRatio" ) );
            const RealType factorVoidMaterial = _parser.template get<RealType> ( "Material.factorVoidMaterial" );
                
            const int numLoads = _parser.template get<int> ( "AffineDisp.numLoads" );
            const int numAffineSymGradDofs = conf.numAffineSymGradDofs;

            for( int loadIdx=0; loadIdx<numLoads; ++loadIdx ){
                    
                VectorType displacementPeriodic ( conf.dimChartDomain * conf.getNumGlobalDofs() ), displacementPeriodicExtended ( conf.dimChartDomain * conf.getNumGlobalDofs() );
                pesopt::loadVectorFromFile<VectorType>( displacementPeriodic, pesopt::strprintf( "%s/Deformation/DisplacementPeriodic_Dir%d.txt", _saveDirectory.c_str(),  loadIdx ) ); 
                bdryHandlerDisplacement.extendMultiVector( displacementPeriodic, displacementPeriodicExtended );
                    
                VectorType displacementAffine ( numAffineSymGradDofs );
                pesopt::loadVectorFromFile<VectorType>( displacementAffine, pesopt::strprintf( "%s/Deformation/DisplacementAffine_Dir%d.txt", _saveDirectory.c_str(), loadIdx ) ); 
                    
                std::vector<RealType> stressVec, stressVecBlock;
                MaterialOptimizationVonMisesStressEvaluator<ConfiguratorType> vonMisesStressOp ( conf, HardMaterial, factorVoidMaterial );
                vonMisesStressOp.evaluateStress( pesopt::strprintf( "%s/Interface/Interface.vtk", _saveDirectory.c_str() ).c_str(), material, displacementPeriodicExtended, displacementAffine, stressVec, false  );
                vonMisesStressOp.evaluateStress( pesopt::strprintf( "%s/Interface/BlockInterface.vtk", _saveDirectory.c_str() ).c_str(), material, displacementPeriodicExtended, displacementAffine, stressVecBlock, true  );
                   
                vtkSaveScalarPointData<std::vector<RealType>>( pesopt::strprintf( "%s/Interface/Interface.vtk", _saveDirectory.c_str() ).c_str(), pesopt::strprintf( "%s/StressOnInterface/VonMises_Interface_Direction%d.vtp", _saveDirectory.c_str(), loadIdx ).c_str(), "VonMisesStresses", stressVec );
                vtkSaveScalarPointData<std::vector<RealType>>( pesopt::strprintf( "%s/Interface/BlockInterface.vtk", _saveDirectory.c_str() ).c_str(), pesopt::strprintf( "%s/StressOnInterface/VonMises_BlockInterface_Direction%d.vtp", _saveDirectory.c_str(), loadIdx ).c_str(), "VonMisesStresses", stressVecBlock );
            }
        }
    
        
        //=====================================================================================================================
        //plot von mises stresses
        //=====================================================================================================================
        
        void saveVonMisesStresses( ) const{
            pesopt::consoleOutputItem( "save von mises stresses" );
            _parser.createSubDirectory("Stress");
            
            MeshType mesh ( _numDofVec, _lengthVec );
            ConfiguratorType conf ( mesh );
            QuocPeriodicBoundaryConditionHandler<ConfiguratorType> bdryHandlerDisplacement( _parser, conf );
            
            VectorType material ( mesh.getNumVertices() );
            pesopt::loadVectorFromFile<VectorType>( material, pesopt::strprintf ( "%s/materialExtended.txt", _saveDirectory.c_str() ) ); 
            FEScalarFunctionEvaluator<ConfiguratorType> discreteFctSingleCell ( conf, material );
            
            //read material properties //TODO optional (since different for bone), maybe use name of material
            MaterialProperties<RealType,ConfiguratorType::dimChartDomain> HardMaterial ( "Hard", 1.0, 1.0, _parser.template get<RealType> ( "Material.ElastModulus" ),
                                                  _parser.template get<RealType> ( "Material.PoissonRatio" ) );
            const RealType factorVoidMaterial = _parser.template get<RealType> ( "Material.factorVoidMaterial" );
            
            const int numLoads = _parser.template get<int> ( "AffineDisp.numLoads" );
            const int numAffineSymGradDofs = conf.numAffineSymGradDofs;
            //for bone
            for( int loadIdx=0; loadIdx<numLoads; ++loadIdx ){
                
                VectorType displacementPeriodic ( conf.dimChartDomain * conf.getNumGlobalDofs() ), displacementPeriodicExtended ( conf.dimChartDomain * conf.getNumGlobalDofs() );
                pesopt::loadVectorFromFile<VectorType>( displacementPeriodic, pesopt::strprintf( "%s/Deformation/DisplacementPeriodic_Dir%d.txt", _saveDirectory.c_str(), loadIdx ) ); 
                bdryHandlerDisplacement.extendMultiVector( displacementPeriodic, displacementPeriodicExtended );
                
                VectorType displacementAffine ( numAffineSymGradDofs );
                pesopt::loadVectorFromFile<VectorType>( displacementAffine, pesopt::strprintf( "%s/Deformation/DisplacementAffine_Dir%d.txt", _saveDirectory.c_str(),  loadIdx ) ); 
                
                std::vector<RealType> stressVec3d, stressVecBlock3d;
                MaterialOptimizationVonMisesStressEvaluator<ConfiguratorType> vonMisesStressOp ( conf, HardMaterial, factorVoidMaterial );
                vonMisesStressOp.evaluateStress( pesopt::strprintf( "%s/material.vtk", _saveDirectory.c_str() ).c_str(), material, displacementPeriodicExtended, displacementAffine, stressVec3d, false  );
                vonMisesStressOp.evaluateStress( pesopt::strprintf( "%s/materialBlock.vtk", _saveDirectory.c_str() ).c_str(), material, displacementPeriodicExtended, displacementAffine, stressVecBlock3d, true );
                
                vtkSaveScalarPointData<std::vector<RealType>>( 
                    pesopt::strprintf( "%s/material.vtk", _saveDirectory.c_str() ).c_str(), 
                    pesopt::strprintf( "%s/Stress/VonMises_Direction%d.vtk", _saveDirectory.c_str(), loadIdx ).c_str(),
                    "VonMisesStresses",
                    stressVec3d );
                vtkSaveScalarPointData<std::vector<RealType>>(
                    pesopt::strprintf( "%s/materialBlock.vtk", _saveDirectory.c_str() ).c_str(), 
                    pesopt::strprintf( "%s/Stress/VonMises_Block_Direction%d.vtk", _saveDirectory.c_str(), loadIdx ).c_str(),
                    "VonMisesStresses", 
                    stressVecBlock3d );
            }
        }
        
};



#endif
