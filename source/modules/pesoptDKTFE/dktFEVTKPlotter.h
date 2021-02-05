#ifndef __DKTFEVTKPLOTTER_H
#define __DKTFEVTKPLOTTER_H

#include <pesopt_IO.h>
#include "dktFEMaterialOptimizationDefines.h"


template< typename DataTypeContainer >
class ShellVTKPlotter{
private:
  
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef pesopt::BoostParser ParameterParserType;
  typedef typename DataTypeContainer::Point3DType Point3DType;
  
  const string _saveDirectory;
  const string _extension;
  
public :
    
    ShellVTKPlotter ( const string saveDirectory, const string extension ) : _saveDirectory( saveDirectory ), _extension( extension ) { }
    
    
    void plotShellWithData( const string fileNameVTK, const ParameterParserType &parserVTK, 
                            const string nameData /*= "material"*/, 
                            const VTKDataSupp supp ) const{
       VTKtoPNGPlotter plotter( _saveDirectory );
       plotter.template plotToPngWithParserInfo<ParameterParserType>( 
          pesopt::strprintf ( "%s/%s.%s", _saveDirectory.c_str(), fileNameVTK.c_str(), _extension.c_str() ), 
          pesopt::strprintf ( "%s/%s_%s_old.png", _saveDirectory.c_str(), fileNameVTK.c_str(), nameData.c_str() ).c_str(), 
          parserVTK, true, nameData.c_str(), supp );
       
       ImageHandler imHandler; 
       imHandler.trimAndResizeImage( pesopt::strprintf ( "%s/%s_%s_old.png", _saveDirectory.c_str(), fileNameVTK.c_str(), nameData.c_str() ), 
                                     pesopt::strprintf ( "%s/%s_%s.png", _saveDirectory.c_str(), fileNameVTK.c_str(), nameData.c_str() ),
                                     parserVTK.template get<int> ("Window.resizeWidth"), parserVTK.template get<int> ("Window.resizeHeight") );
       const string NewDirectory = pesopt::createSubDirectory( _saveDirectory, "old" );
       pesopt::moveFile( pesopt::strprintf ( "%s_%s_old.png", fileNameVTK.c_str(), nameData.c_str() ), _saveDirectory.c_str(), NewDirectory );
    } 
    
     
    void plotStress( const ParameterParserType &parserVTKStress,
                     const string fileNameVTK,
                     const std::vector<string> &nameStressOnElementVec  ) const{
       VTKtoPNGPlotter plotterStress; 
       ImageHandler imHandler;
        
       //! plot stress on Elements
       for( int i=0; i<nameStressOnElementVec.size(); ++i ){
           plotterStress.template plotToPngWithParserInfo<ParameterParserType>( 
                pesopt::strprintf ( "%s/%s.%s", _saveDirectory.c_str(), fileNameVTK.c_str(), _extension.c_str() ),
                pesopt::strprintf ( "%s/%s_%s_old.png", _saveDirectory.c_str(), fileNameVTK.c_str(), nameStressOnElementVec[i].c_str() ).c_str(),
                parserVTKStress, true, nameStressOnElementVec[i].c_str(), FACE_DATA );
           imHandler.trimAndResizeImage( pesopt::strprintf ( "%s/%s_%s_old.png", _saveDirectory.c_str(), fileNameVTK.c_str(), nameStressOnElementVec[i].c_str() ),
                                         pesopt::strprintf ( "%s/%s_%s.png", _saveDirectory.c_str(), fileNameVTK.c_str(), nameStressOnElementVec[i].c_str() ),
                                         parserVTKStress.template get<int> ("Window.resizeWidth"), parserVTKStress.template get<int> ("Window.resizeHeight")  );
           const string NewDirectory = pesopt::createSubDirectory( _saveDirectory, "old" );
           pesopt::moveFile( pesopt::strprintf ( "%s_%s_old.png", fileNameVTK.c_str(), nameStressOnElementVec[i].c_str() ), _saveDirectory.c_str(), NewDirectory );
       }
   }
    
    void plotStressOnUndeformedShell( const ParameterParserType &parserVTK,
                                      const string fileNameVTK,
                                      const std::vector<string> &nameStressOnElementVec  ) const{
       VTKtoPNGPlotter plotterStress; 
       ImageHandler imHandler;
       ParameterParserType parserVTKStress ( parserVTK );
       parserVTKStress.set( "ScalarBar.useScalarData", parserVTK.template get<int> ( "ScalarBarStress.useScalarData" )  );
       parserVTKStress.set( "ScalarBar.useScalarBar", parserVTK.template get<int> ( "ScalarBarStress.useScalarBar" ) );
       parserVTKStress.set( "ScalarBar.ColorSpaceType", parserVTK.template get<int> ( "ScalarBarStress.ColorSpaceType" ) );
       parserVTKStress.set( "ScalarBar.ScaleType", parserVTK.template get<int> ( "ScalarBarStress.ScaleType" ) );
       
        const int scaleRGBPointsToRange = parserVTK.template get<int> ( "ScalarBarStress.scaleRGBPointsToRange" );
        parserVTKStress.set( "ScalarBar.scaleRGBPointsToRange", scaleRGBPointsToRange );
        const int numRGBPoints = parserVTK.template get<int> ( "ScalarBarStress.numRGBPoints" );
        parserVTKStress.set( "ScalarBar.numRGBPoints", numRGBPoints );
        for( int i=1; i<=numRGBPoints; ++i ){
            parserVTKStress.set( pesopt::strprintf("ScalarBar.RGBValue%d", i ).c_str(), parserVTK.template get<double>( pesopt::strprintf("ScalarBarStress.RGBValue%d", i ).c_str() ) );
            Point3DType RGBPoint; parserVTK.template getFixSizeVector<double,Point3DType> ( pesopt::strprintf("ScalarBarStress.RGBPoint%d", i ).c_str(), RGBPoint );
            parserVTKStress.template setFixSizeVector<Point3DType>( pesopt::strprintf("ScalarBar.RGBPoint%d", i ).c_str(), RGBPoint, false );
        }

    //        //! plot stress at Nodes
    //        const std::vector<string> nameStress = { "membraneStress", "bendingStress", "totalStress"  };
    //        for( int i=0; i<nameStress.size(); ++i ){
    //            plotterStress.template plotToPngWithParserInfo<ParameterParserType>( 
    //                 pesopt::strprintf ( "%s/Stress/allStressNodes.%s", _saveDirectory.c_str(), _extension.c_str() ),
    //                 pesopt::strprintf ( "%s/Stress/%sNodes_old.png", _saveDirectory.c_str(), nameStress[i].c_str() ).c_str(),
    //                 parserVTKStress, true, nameStress[i].c_str(), 1 );
    //            imHandler.trimAndResizeImage( pesopt::strprintf ( "%s/Stress/%sNodes_old.png", _saveDirectory.c_str(), nameStress[i].c_str() ),
    //                                          pesopt::strprintf ( "%s/Stress/%sNodes.png", _saveDirectory.c_str(), nameStress[i].c_str() ),
    //                                          parserVTKStress.template get<int> ("Window.resizeWidth"), parserVTKStress.template get<int> ("Window.resizeHeight")  );
    //        }

       //! plot stress on Elements
       for( int i=0; i<nameStressOnElementVec.size(); ++i ){
           plotterStress.template plotToPngWithParserInfo<ParameterParserType>( 
                pesopt::strprintf ( "%s/%s.%s", _saveDirectory.c_str(), fileNameVTK.c_str(), _extension.c_str() ),
                pesopt::strprintf ( "%s/%s_%s_old.png", _saveDirectory.c_str(), fileNameVTK.c_str(), nameStressOnElementVec[i].c_str() ).c_str(),
                parserVTKStress, true, nameStressOnElementVec[i].c_str(), FACE_DATA );
           imHandler.trimAndResizeImage( pesopt::strprintf ( "%s/%s_%s_old.png", _saveDirectory.c_str(), fileNameVTK.c_str(), nameStressOnElementVec[i].c_str() ),
                                         pesopt::strprintf ( "%s/%s_%s.png", _saveDirectory.c_str(), fileNameVTK.c_str(), nameStressOnElementVec[i].c_str() ),
                                         parserVTKStress.template get<int> ("Window.resizeWidth"), parserVTKStress.template get<int> ("Window.resizeHeight")  );
           const string NewDirectory = pesopt::createSubDirectory( _saveDirectory, "old" );
           pesopt::moveFile( pesopt::strprintf ( "%s_%s_old.png", fileNameVTK.c_str(), nameStressOnElementVec[i].c_str() ), _saveDirectory.c_str(), NewDirectory );
       }
   }
   
   void plotMarkedElementsOnUndeformedShell( const string fileNameVTK, const ParameterParserType &parserVTKUndeformed, const string nameVec  ) const{
       VTKtoPNGPlotter plotterUndeformed( _saveDirectory );
       plotterUndeformed.template plotToPngWithParserInfo<ParameterParserType>( 
        pesopt::strprintf ( "%s/%s.%s", _saveDirectory.c_str(), fileNameVTK.c_str(), _extension.c_str() ), 
        pesopt::strprintf ( "%s/%s_old.png", _saveDirectory.c_str(), fileNameVTK.c_str() ).c_str(),
        parserVTKUndeformed, true, nameVec.c_str(), FACE_DATA );
       
       ImageHandler imHandler; 
       imHandler.trimAndResizeImage( pesopt::strprintf ( "%s/%s_old.png", _saveDirectory.c_str(), fileNameVTK.c_str() ), pesopt::strprintf ( "%s/%s.png", _saveDirectory.c_str(), fileNameVTK.c_str() ),
                                     parserVTKUndeformed.template get<int> ("Window.resizeWidth"), parserVTKUndeformed.template get<int> ("Window.resizeHeight") );
       const string NewDirectory = pesopt::createSubDirectory( _saveDirectory, "old" );
       pesopt::moveFile( pesopt::strprintf ( "%s_old.png", fileNameVTK.c_str() ), _saveDirectory.c_str(), NewDirectory );
    }
   
   
  
};



#endif //__DKTFEVTKPLOTTER_H
