#ifndef __VTKTOPNGPLOTTER_H
#define __VTKTOPNGPLOTTER_H

#include <VTKSaver.h>

// #define __DEBUGVTKTOPNGPLOTTER
 
class VTKtoPNGPlotter {
    
public :
    
    VTKtoPNGPlotter( const string saveDirectory ) {}
    
    VTKtoPNGPlotter( ) {}

    template<typename ParameterParserType>
    void plotToPngWithParserInfo( const string &inputFileName, const string &outputFileName, const ParameterParserType &parser,
                                  const bool useScalarData = false,  const string &nameScalarData = "", const VTKDataSupp supp = VERTEX_DATA
                                ) const {
        
#ifdef __DEBUGVTKTOPNGPLOTTER
        cout << endl << "======================================" << endl
             << "plotToPngWithParserInfo" << endl
             << "inputFileName = " << inputFileName.c_str() << endl
             << "outputFileName = " << outputFileName.c_str() << endl;
#endif

        VTKSaver vtkReader;
        string file_extension; VTKDATATYPE dataType;
        vtkSmartPointer<vtkDataSet> grid; vtkReader.readVTKDataSet( inputFileName, grid, file_extension, dataType );
      
    //mapper
        auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper->SetInputData(grid);
        if( parser.template get<bool>( "Mapper.InterpolateScalarsBeforeMapping" ) ) mapper->InterpolateScalarsBeforeMappingOn();   
        else mapper->InterpolateScalarsBeforeMappingOff();   

              
        
    //coloring by point data
      double scalarRange[2];
      if( useScalarData ){
                   
        //TODO is this always the right option?
        mapper->SetColorModeToMapScalars();
        
#ifdef __DEBUGVTKTOPNGPLOTTER
        cout << "try to find array with name " << nameScalarData.c_str() << endl;
#endif
          
        switch( supp ){
            case VERTEX_DATA:{
            //TODO ???
            mapper->SetScalarModeToUsePointFieldData();
//             mapper->SetScalarModeToUsePointData();
            mapper->ScalarVisibilityOn();
            vtkSmartPointer<vtkDataArray> dataArray = grid->GetPointData()->GetArray(nameScalarData.c_str());
            if(dataArray){
#ifdef __DEBUGVTKTOPNGPLOTTER
                cout << "found data array " << nameScalarData.c_str() << endl;
                vtkIdType numVectors = dataArray->GetNumberOfTuples();
                const int numComponents = dataArray->GetNumberOfComponents();
                for (vtkIdType tupleIdx = 0; tupleIdx < numVectors; ++tupleIdx)
                    for( int compIdx = 0; compIdx < numComponents; ++compIdx )
                        cout << "dataArray_" << tupleIdx << "_" << compIdx << " = " << dataArray->GetComponent(tupleIdx, compIdx) << endl;
#endif
            }
            else throw invalid_argument( pesopt::strprintf ( "In File %s. Unknown Field for point scalar data %s. In File %s at line %d.", inputFileName.c_str(), nameScalarData.c_str(),  __FILE__, __LINE__ ).c_str() ); 
            mapper->SelectColorArray(nameScalarData.c_str());
            mapper->SetScalarRange(grid->GetPointData()->GetArray(nameScalarData.c_str())->GetRange());
            mapper->GetScalarRange (scalarRange);

            }
            break;
            
            case FACE_DATA:{
                mapper->SetScalarModeToUseCellFieldData();
                mapper->ScalarVisibilityOn();
                vtkSmartPointer<vtkDataArray> dataArray = grid->GetCellData()->GetArray(nameScalarData.c_str());
                if(dataArray){
    #ifdef __DEBUGVTKTOPNGPLOTTER
                    cout << "found data array " << nameScalarData.c_str() << endl;
    #endif
                }
                else throw invalid_argument( pesopt::strprintf ( "In File %s. Unknown Field for face scalar data %s. In File %s at line %d.", inputFileName.c_str(), nameScalarData.c_str(),  __FILE__, __LINE__ ).c_str() ); 
                mapper->SelectColorArray(nameScalarData.c_str());
                mapper->SetScalarRange(grid->GetCellData()->GetArray(nameScalarData.c_str())->GetRange());
                mapper->GetScalarRange (scalarRange);
            }
            break;

            default:
                throw std::invalid_argument( pesopt::strprintf ( "Unknown VTKDataSupp. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
         }  
      }
          
        
        
      //scalar bar
        auto scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
        auto ctf = vtkSmartPointer<vtkColorTransferFunction>::New(); 
        const int numRGBPoints = parser.template get<int> ( "ScalarBar.numRGBPoints" );
        for( int i=1; i<=numRGBPoints; ++i ){
                 double RGBValue = parser.template get<double> ( pesopt::strprintf("ScalarBar.RGBValue%d", i ).c_str() );
                 double RGBPoint[3]; parser.template getFixSizeVector<double,double[3]> ( pesopt::strprintf("ScalarBar.RGBPoint%d", i ).c_str(), RGBPoint );
                 if( parser.template get<bool> ( "ScalarBar.scaleRGBPointsToRange" ) )
                      ctf->AddRGBPoint( scalarRange[0] + RGBValue * (scalarRange[1] - scalarRange[0]), RGBPoint[0], RGBPoint[1], RGBPoint[2] );
                 else ctf->AddRGBPoint(RGBValue, RGBPoint[0], RGBPoint[1], RGBPoint[2]);
        }
             
        switch( parser.template get<int> ( "ScalarBar.ColorSpaceType" ) ){
                 case 1 :   ctf->SetColorSpaceToDiverging(); break;
                 case 2 :   ctf->SetColorSpaceToRGB(); break;
                 case 3 : { ctf->SetColorSpaceToHSV(); 
                            ctf->HSVWrapOff(); //ctf->HSVWrapOn();
                 }break;
                 case 4 :   ctf->SetColorSpaceToLab(); break;
                 default :  throw invalid_argument( pesopt::strprintf ( "Unknown ColorSpaceType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
        }
             
        switch( parser.template get<int> ( "ScalarBar.ScaleType" ) ){
                 case 1 : ctf->SetScaleToLinear(); break;
                 case 2 : ctf->SetScaleToLog10(); break;
                 default : throw invalid_argument( pesopt::strprintf ( "Unknown ScaleType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
        }
             
        mapper->SetLookupTable( ctf );
            
        scalarBar->SetLookupTable( vtkLookupTable::SafeDownCast(mapper->GetLookupTable()) );
        scalarBar->SetTitle("");
        scalarBar->SetNumberOfLabels( parser.template get<int> ( "ScalarBar.NumberOfLabels" )  ); 
        scalarBar->SetLookupTable( ctf );
        
        
      //actor
        auto actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetSpecular( parser.template get<double> ( "Actor.Specular" ) ); 
        actor->GetProperty()->SetSpecularPower( parser.template get<double> ( "Actor.SpecularPower" ) ); 
        double DiffuseColor[3]; parser.template getFixSizeVector<double,double[3]> ( "Actor.DiffuseColor", DiffuseColor );
        actor->GetProperty()->SetDiffuseColor( DiffuseColor ); 
        actor->GetProperty()->SetLineWidth( parser.template get<double> ( "Actor.LineWidth" ) );
        if( parser.template get<bool>( "Actor.useEdgeColor" ) ){
             double EdgeColor[3]; parser.template getFixSizeVector<double,double[3]> ( "Actor.EdgeColor", EdgeColor );
             actor->GetProperty()->SetEdgeColor( EdgeColor ); //(R,G,B)
             actor->GetProperty()->EdgeVisibilityOn();
        }
        
      //light
         auto light = vtkSmartPointer<vtkLightKit>::New();
         if( parser.template get<bool> ("Light.useLight") ){
            //light->SetPosition(lightPosition);         
            light->SetKeyLightWarmth(parser.template get<double>("Light.KeyWarmth")); light->SetKeyLightIntensity(parser.template get<double>("Light.KeyIntensity"));
            light->SetKeyLightAngle(parser.template get<double>("Light.KeyElevation"),parser.template get<double>("Light.KeyAzimuth"));
            light->SetFillLightWarmth(parser.template get<double>("Light.FillWarmth")); light->SetKeyToFillRatio(parser.template get<double>("Light.FillKeyRatio"));
            light->SetFillLightAngle(parser.template get<double>("Light.FillElevation"),parser.template get<double>("Light.FillAzimuth"));
            light->SetBackLightWarmth(parser.template get<double>("Light.BackWarmth")); light->SetKeyToBackRatio(parser.template get<double>("Light.BackKeyRatio"));
            light->SetBackLightAngle(parser.template get<double>("Light.BackElevation"),parser.template get<double>("Light.BackAzimuth"));
            light->SetHeadLightWarmth(parser.template get<double>("Light.HeadWarmth")); light->SetKeyToHeadRatio(parser.template get<double>("Light.HeadKeyRatio"));
         }
        
      //box at boundary
        vtkSmartPointer<vtkDataSet> gridBox;
        auto mapperBox = vtkSmartPointer<vtkDataSetMapper>::New();
        auto actorBox = vtkSmartPointer<vtkActor>::New();
        if( parser.template get<int> ("BoundaryBox.useBoundaryBox") ){
            string file_extensionBox; VTKDATATYPE dataTypeBox;
            vtkReader.readVTKDataSet( parser.template get<string>("BoundaryBox.fileName").c_str(), gridBox, file_extensionBox, dataTypeBox );
            mapperBox->SetInputData(gridBox);
            actorBox->SetMapper(mapperBox);
            if( parser.template get<int> ("BoundaryBox.useEdgeVisibility") ) actorBox->GetProperty()->EdgeVisibilityOn();
        }

      //Renderer
      //         renderer->SetInteractionModeTo2D();
        auto renderer = vtkSmartPointer<vtkRenderer>::New();
        renderer->AddActor(actor);
        if( parser.template get<int> ("BoundaryBox.useBoundaryBox") ) renderer->AddActor(actorBox);
        if( parser.template get<int> ("Light.useLight") ) light->AddLightsToRenderer(renderer);//for vtklight: renderer->AddLight(light);
        if( parser.template get<bool>("ScalarBar.useScalarBar")  ) renderer->AddActor2D(scalarBar);
        
        double Position[3]; parser.template getFixSizeVector<double,double[3]> ( "Camera.Position", Position );
        renderer->GetActiveCamera()->SetPosition(Position);
        double FocalPoint[3]; parser.template getFixSizeVector<double,double[3]> ( "Camera.FocalPoint", FocalPoint );
        renderer->GetActiveCamera()->SetFocalPoint(FocalPoint);
        double ViewUp[3]; parser.template getFixSizeVector<double,double[3]> ( "Camera.ViewUp", ViewUp );
        renderer->GetActiveCamera()->SetViewUp(ViewUp);
        renderer->GetActiveCamera()->SetParallelScale( parser.template get<double> ( "Camera.ParallelScale" ) );
        renderer->GetActiveCamera()->Roll( parser.template get<double> ( "Camera.Roll" ) );
        renderer->GetActiveCamera()->Azimuth( parser.template get<double> ( "Camera.Azimuth" ) );
        renderer->GetActiveCamera()->Elevation( parser.template get<double> ( "Camera.Elevation" ) );

        renderer->ResetCamera();
        renderer->GetActiveCamera()->Dolly( parser.template get<double> ( "Camera.Dolly" ) );
        renderer->ResetCameraClippingRange();
        double Background[3]; parser.template getFixSizeVector<double,double[3]> ( "Renderer.Background", Background );
        renderer->SetBackground(Background);
        renderer->UseFXAAOn(); //anti aliasing, renderWindow->SetAAFrames(10) is deprecated
        auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        
        auto renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);
        if( parser.template get<int>( "Window.FullScreenRendering" ) == 1 ) renderWindow->FullScreenOn();
        else renderWindow->SetSize( parser.template get<int>("Window.imageWidth"), parser.template get<int>("Window.imageHeight") );
        
        
        renderWindow->Render();
        // renderWindowInteractor->Start();
        renderWindow->OffScreenRenderingOn();
        renderWindow->Render();
        
      // Screenshot  
        auto windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
        windowToImageFilter->SetInput(renderWindow);         
        
        auto writer = vtkSmartPointer<vtkPNGWriter>::New();
        writer->SetFileName( outputFileName.c_str() );
        writer->SetCompressionLevel( parser.template get<int>("Writer.compressionLevel") ); //from 0-9: 0 no compression
        writer->SetInputConnection(windowToImageFilter->GetOutputPort());
        writer->Write();
        
        
        
        
        
      //Plot scalar bar separately
        if( parser.template get<bool>("ScalarBar.plotScalarBarSeparately")  ){
            
            size_t lastindex = outputFileName.find_last_of("."); 
            const string outputFileNameScalarBar = outputFileName.substr(0, lastindex) + "_ScalarBar.png"; 
            
            auto rendererScalarBar = vtkSmartPointer<vtkRenderer>::New();
            rendererScalarBar->AddActor2D(scalarBar);
            auto renderWindowScalarBar = vtkSmartPointer<vtkRenderWindow>::New();
            renderWindowScalarBar->AddRenderer(rendererScalarBar);
            auto renderWindowInteractorScalarBar = vtkSmartPointer<vtkRenderWindowInteractor>::New();
            renderWindowInteractorScalarBar->SetRenderWindow(renderWindowScalarBar);
            if( parser.template get<int>( "Window.FullScreenRendering" ) == 1 ) renderWindowScalarBar->FullScreenOn();
            else renderWindowScalarBar->SetSize( parser.template get<int>("Window.imageWidth"), parser.template get<int>("Window.imageHeight") );
            
            renderWindowScalarBar->Render();
            renderWindowScalarBar->OffScreenRenderingOn();
            renderWindowScalarBar->Render();
            
            auto windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
            windowToImageFilter->SetInput(renderWindowScalarBar);         
            
            auto writerScalarBar = vtkSmartPointer<vtkPNGWriter>::New();
            writerScalarBar->SetFileName( outputFileNameScalarBar.c_str() );
            writerScalarBar->SetCompressionLevel( parser.template get<int>("Writer.compressionLevel") ); //from 0-9: 0 no compression
            writerScalarBar->SetInputConnection(windowToImageFilter->GetOutputPort());
            writerScalarBar->Write();
        }
        
    }
    
    
    
    
    
     //scalar data type:  1 - point, 2 - cell //TODO replace by enum or string
    template<typename ParameterParserType>
    void plotMultipleToPngWithParserInfo( const string &inputFileName1, 
                                          const ParameterParserType &parser1, const bool useScalarData1,  const string &nameScalarData1, const VTKDataSupp supp1,
                                          const string &inputFileName2, 
                                          const ParameterParserType &parser2, const bool useScalarData2,  const string &nameScalarData2, const VTKDataSupp supp2,
                                          const ParameterParserType &parserLight, const ParameterParserType &parserRenderer,
                                          const string &outputFileName
                                        ) const {
        
#ifdef __DEBUGVTKTOPNGPLOTTER
        cout << endl << "======================================" << endl
             << "plotToPngWithParserInfo" << endl
             << "inputFileName1 = " << inputFileName1.c_str() << endl
             << "inputFileName2 = " << inputFileName2.c_str() << endl
             << "outputFileName = " << outputFileName.c_str() << endl;
#endif

        VTKSaver vtkReader;
        
        
   /////////////////////////////////////////////////     
   //             inputFile1     
   ////////////////////////////////////////////////     
        
        string file_extension1; VTKDATATYPE dataType1;
        vtkSmartPointer<vtkDataSet> grid1; vtkReader.readVTKDataSet( inputFileName1, grid1, file_extension1, dataType1 );
      
      //mapper
        auto mapper1 = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper1->SetInputData(grid1);
        if( parser1.template get<bool>( "Mapper.InterpolateScalarsBeforeMapping" ) ) mapper1->InterpolateScalarsBeforeMappingOn();   
        else mapper1->InterpolateScalarsBeforeMappingOff();   
        
      //scalar bar
//         auto scalarBar1 = vtkSmartPointer<vtkScalarBarActor>::New(); 
//         if( parser1.template get<bool>("ScalarBar.useScalarData") ){
            auto ctf1 = vtkSmartPointer<vtkColorTransferFunction>::New();
             const int numRGBPoints1 = parser1.template get<int> ( "ScalarBar.numRGBPoints" );
             for( int i=1; i<=numRGBPoints1; ++i ){
                 double RGBValue = parser1.template get<double> ( pesopt::strprintf("ScalarBar.RGBValue%d", i ).c_str() );
                 double RGBPoint[3]; parser1.template getFixSizeVector<double,double[3]> ( pesopt::strprintf("ScalarBar.RGBPoint%d", i ).c_str(), RGBPoint );
                 ctf1->AddRGBPoint(RGBValue, RGBPoint[0], RGBPoint[1], RGBPoint[2]);
             }
             switch( parser1.template get<int> ( "ScalarBar.ColorSpaceType" ) ){
                 case 1 :   ctf1->SetColorSpaceToDiverging(); break;
                 case 2 :   ctf1->SetColorSpaceToRGB(); break;
                 case 3 : { ctf1->SetColorSpaceToHSV(); 
                            ctf1->HSVWrapOff(); //ctf1->HSVWrapOn();
                 }break;
                 case 4 :   ctf1->SetColorSpaceToLab(); break;
                 default :  throw invalid_argument( pesopt::strprintf ( "Unknown ColorSpaceType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
             }
             switch( parser1.template get<int> ( "ScalarBar.ScaleType" ) ){
                 case 1 : ctf1->SetScaleToLinear(); break;
                 case 2 : ctf1->SetScaleToLog10(); break;
                 default : throw invalid_argument( pesopt::strprintf ( "Unknown ScaleType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
             }
            mapper1->SetLookupTable( ctf1 );

//             scalarBar1->SetLookupTable( vtkLookupTable::SafeDownCast(mapper1->GetLookupTable()) );
//             scalarBar1->SetTitle("");
//             scalarBar1->SetNumberOfLabels(4); 
//             scalarBar1->SetLookupTable( ctf1 );
//         }
        
      //coloring by point data
      if( useScalarData1 ){
        switch( supp1 ){
            case VERTEX_DATA:{
                mapper1->SetScalarModeToUsePointFieldData();
                mapper1->ScalarVisibilityOn();
                vtkSmartPointer<vtkDataArray> dataArray = grid1->GetPointData()->GetArray(nameScalarData1.c_str());
                if( !dataArray )  throw invalid_argument( pesopt::strprintf ( "Unknown Field for scalar data. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
                mapper1->SelectColorArray(nameScalarData1.c_str());
                mapper1->SetScalarRange(grid1->GetPointData()->GetArray(nameScalarData1.c_str())->GetRange());
            }
            break;
            
            case FACE_DATA:{
                mapper1->SetScalarModeToUseCellFieldData();
                mapper1->ScalarVisibilityOn();
                vtkSmartPointer<vtkDataArray> dataArray = grid1->GetCellData()->GetArray(nameScalarData1.c_str());
                if( !dataArray ) throw invalid_argument( pesopt::strprintf ( "Unknown Field for scalar data. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
                mapper1->SelectColorArray(nameScalarData1.c_str());
                mapper1->SetScalarRange(grid1->GetCellData()->GetArray(nameScalarData1.c_str())->GetRange());
            }
            break;

            default:
                throw std::invalid_argument( pesopt::strprintf ( "Unknown VTKDataSupp. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
         }
      } 
        
      //actor
        auto actor1 = vtkSmartPointer<vtkActor>::New();
        actor1->SetMapper(mapper1);
        actor1->GetProperty()->SetSpecular( parser1.template get<double> ( "Actor.Specular" ) ); 
        actor1->GetProperty()->SetSpecularPower( parser1.template get<double> ( "Actor.SpecularPower" ) ); 
        double DiffuseColor1[3]; parser1.template getFixSizeVector<double,double[3]> ( "Actor.DiffuseColor", DiffuseColor1 );
        actor1->GetProperty()->SetDiffuseColor( DiffuseColor1 ); 
        actor1->GetProperty()->SetLineWidth( parser1.template get<double> ( "Actor.LineWidth" ) );
        
        
        
   /////////////////////////////////////////////////     
   //             inputFile2   
   ////////////////////////////////////////////////     
        
        string file_extension2; VTKDATATYPE dataType2;
        vtkSmartPointer<vtkDataSet> grid2; vtkReader.readVTKDataSet( inputFileName2, grid2, file_extension2, dataType2 );
      
      //mapper
        auto mapper2 = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper2->SetInputData(grid2);
        if( parser2.template get<bool>( "Mapper.InterpolateScalarsBeforeMapping" ) ) mapper2->InterpolateScalarsBeforeMappingOn();   
        else mapper2->InterpolateScalarsBeforeMappingOff();   
        
      //scalar bar
//         auto scalarBar2 = vtkSmartPointer<vtkScalarBarActor>::New(); 
//         if( parser2.template get<bool>("ScalarBar.useScalarData") ){
            auto ctf2 = vtkSmartPointer<vtkColorTransferFunction>::New();
             const int numRGBPoints2 = parser2.template get<int> ( "ScalarBar.numRGBPoints" );
             for( int i=1; i<=numRGBPoints2; ++i ){
                 double RGBValue = parser2.template get<double> ( pesopt::strprintf("ScalarBar.RGBValue%d", i ).c_str() );
                 double RGBPoint[3]; parser2.template getFixSizeVector<double,double[3]> ( pesopt::strprintf("ScalarBar.RGBPoint%d", i ).c_str(), RGBPoint );
                 ctf2->AddRGBPoint(RGBValue, RGBPoint[0], RGBPoint[1], RGBPoint[2]);
             }
             switch( parser2.template get<int> ( "ScalarBar.ColorSpaceType" ) ){
                 case 1 :   ctf2->SetColorSpaceToDiverging(); break;
                 case 2 :   ctf2->SetColorSpaceToRGB(); break;
                 case 3 : { ctf2->SetColorSpaceToHSV(); 
                            ctf2->HSVWrapOff(); //ctf2->HSVWrapOn();
                 }break;
                 case 4 :   ctf2->SetColorSpaceToLab(); break;
                 default :  throw invalid_argument( pesopt::strprintf ( "Unknown ColorSpaceType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
             }
             switch( parser2.template get<int> ( "ScalarBar.ScaleType" ) ){
                 case 1 : ctf2->SetScaleToLinear(); break;
                 case 2 : ctf2->SetScaleToLog10(); break;
                 default : throw invalid_argument( pesopt::strprintf ( "Unknown ScaleType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
             }
            mapper2->SetLookupTable( ctf2 );

//             scalarBar2->SetLookupTable( vtkLookupTable::SafeDownCast(mapper2->GetLookupTable()) );
//             scalarBar2->SetTitle("");
//             scalarBar2->SetNumberOfLabels(4); 
//             scalarBar2->SetLookupTable( ctf2 );
//         }
        
      //coloring by point data
      if( useScalarData2 ){
        switch( supp2 ){
            case VERTEX_DATA:{
                mapper2->SetScalarModeToUsePointFieldData();
                mapper2->ScalarVisibilityOn();
                vtkSmartPointer<vtkDataArray> dataArray = grid2->GetPointData()->GetArray(nameScalarData2.c_str());
                if( !dataArray ) throw invalid_argument( pesopt::strprintf ( "Unknown Field for scalar data. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
                mapper2->SelectColorArray(nameScalarData2.c_str());
                mapper2->SetScalarRange(grid2->GetPointData()->GetArray(nameScalarData2.c_str())->GetRange());
            }
            break;
            
            case FACE_DATA:{
                mapper2->SetScalarModeToUseCellFieldData();
                mapper2->ScalarVisibilityOn();
                vtkSmartPointer<vtkDataArray> dataArray = grid2->GetCellData()->GetArray(nameScalarData2.c_str());
                if( !dataArray ) throw invalid_argument( pesopt::strprintf ( "Unknown Field for scalar data. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
                mapper2->SelectColorArray(nameScalarData2.c_str());
                mapper2->SetScalarRange(grid2->GetCellData()->GetArray(nameScalarData2.c_str())->GetRange());
            }
            break;

            default:
                throw std::invalid_argument( pesopt::strprintf ( "Unknown VTKDataSupp. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
         }
      } 
        
      //actor
        auto actor2 = vtkSmartPointer<vtkActor>::New();
        actor2->SetMapper(mapper2);
        actor2->GetProperty()->SetSpecular( parser2.template get<double> ( "Actor.Specular" ) ); 
        actor2->GetProperty()->SetSpecularPower( parser2.template get<double> ( "Actor.SpecularPower" ) ); 
        double DiffuseColor2[3]; parser2.template getFixSizeVector<double,double[3]> ( "Actor.DiffuseColor", DiffuseColor2 );
        actor2->GetProperty()->SetDiffuseColor( DiffuseColor2 ); 
        actor2->GetProperty()->SetLineWidth( parser2.template get<double> ( "Actor.LineWidth" ) );
        
        
      ///////////////////////////  
      //light
      ///////////////////////////  
         auto light = vtkSmartPointer<vtkLightKit>::New();
        //light->SetPosition(lightPosition);         
        light->SetKeyLightWarmth(parserLight.template get<double>("Light.KeyWarmth")); light->SetKeyLightIntensity(parserLight.template get<double>("Light.KeyIntensity"));
        light->SetKeyLightAngle(parserLight.template get<double>("Light.KeyElevation"),parserLight.template get<double>("Light.KeyAzimuth"));
        light->SetFillLightWarmth(parserLight.template get<double>("Light.FillWarmth")); light->SetKeyToFillRatio(parserLight.template get<double>("Light.FillKeyRatio"));
        light->SetFillLightAngle(parserLight.template get<double>("Light.FillElevation"),parserLight.template get<double>("Light.FillAzimuth"));
        light->SetBackLightWarmth(parserLight.template get<double>("Light.BackWarmth")); light->SetKeyToBackRatio(parserLight.template get<double>("Light.BackKeyRatio"));
        light->SetBackLightAngle(parserLight.template get<double>("Light.BackElevation"),parserLight.template get<double>("Light.BackAzimuth"));
        light->SetHeadLightWarmth(parserLight.template get<double>("Light.HeadWarmth")); light->SetKeyToHeadRatio(parserLight.template get<double>("Light.HeadKeyRatio"));
        
        
      //box at boundary //TODO
//         vtkSmartPointer<vtkDataSet> gridBox;
//         auto mapperBox = vtkSmartPointer<vtkDataSetMapper>::New();
//         auto actorBox = vtkSmartPointer<vtkActor>::New();
//         if( parser.template get<int> ("BoundaryBox.useBoundaryBox") ){
//             string file_extensionBox; VTKDATATYPE dataTypeBox;
//             vtkReader.readVTKDataSet( parser.template get<string>("BoundaryBox.fileName").c_str(), gridBox, file_extensionBox, dataTypeBox );
//             mapperBox->SetInputData(gridBox);
//             actorBox->SetMapper(mapperBox);
//         }
        
        
      //Renderer
        auto renderer = vtkSmartPointer<vtkRenderer>::New();
        renderer->AddActor(actor1);
        renderer->AddActor(actor2);
//         if( parserRenderer.template get<int> ("BoundaryBox.useBoundaryBox") ) renderer->AddActor(actorBox); //TODO
        light->AddLightsToRenderer(renderer);//for vtklight: renderer->AddLight(light);
//         if( parserRenderer.template get<bool>("ScalarBar.useScalarBar")  ) renderer->AddActor2D(scalarBar);
        
        double Position[3]; parserRenderer.template getFixSizeVector<double,double[3]> ( "Camera.Position", Position );
        renderer->GetActiveCamera()->SetPosition(Position);
        double FocalPoint[3]; parserRenderer.template getFixSizeVector<double,double[3]> ( "Camera.FocalPoint", FocalPoint );
        renderer->GetActiveCamera()->SetFocalPoint(FocalPoint);
        double ViewUp[3]; parserRenderer.template getFixSizeVector<double,double[3]> ( "Camera.ViewUp", ViewUp );
        renderer->GetActiveCamera()->SetViewUp(ViewUp);
        renderer->GetActiveCamera()->SetParallelScale( parserRenderer.template get<double> ( "Camera.ParallelScale" ) );
        renderer->GetActiveCamera()->Roll( parserRenderer.template get<double> ( "Camera.Roll" ) );
        renderer->GetActiveCamera()->Azimuth( parserRenderer.template get<double> ( "Camera.Azimuth" ) );
        renderer->GetActiveCamera()->Elevation( parserRenderer.template get<double> ( "Camera.Elevation" ) );

        renderer->ResetCamera();
        renderer->GetActiveCamera()->Dolly( parserRenderer.template get<double> ( "Camera.Dolly" ) );
        renderer->ResetCameraClippingRange();
        double Background[3]; parserRenderer.template getFixSizeVector<double,double[3]> ( "Renderer.Background", Background );
        renderer->SetBackground(Background);
        renderer->UseFXAAOn(); //anti aliasing, renderWindow->SetAAFrames(10) is deprecated
        auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        
        auto renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);
        if( parserRenderer.template get<int>( "Window.FullScreenRendering" ) == 1 ) renderWindow->FullScreenOn();
        else renderWindow->SetSize( parserRenderer.template get<int>("Window.imageWidth"), parserRenderer.template get<int>("Window.imageHeight") );
        
        
        renderWindow->Render();
        // renderWindowInteractor->Start();
        renderWindow->OffScreenRenderingOn();
        renderWindow->Render();
        
      // Screenshot  
        auto windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
        windowToImageFilter->SetInput(renderWindow);         
        
        auto writer = vtkSmartPointer<vtkPNGWriter>::New();
        writer->SetFileName( outputFileName.c_str() );
        writer->SetCompressionLevel( parserRenderer.template get<int>("Writer.compressionLevel") ); //from 0-9: 0 no compression
        writer->SetInputConnection(windowToImageFilter->GetOutputPort());
        writer->Write();
    }
    
    
#ifdef PESOPT_WITH_BOOST
    void plotToPngWithDefaultArgs( const string &inputFileName, const string &outputFileName, 
                                   const bool useScalarData = false, const string nameScalarData = "", 
                                   const VTKDataSupp supp = VERTEX_DATA ) const {
      
      pesopt::BoostParser parser;
      parser.set ( "Actor.Specular", 0. ); 
      parser.set ( "Actor.SpecularPower", 1. ); 
      parser.set ( "Actor.LineWidth", 1. );
      parser.set ( "Actor.useEdgeColor", 0 );
      std::vector<double> DiffuseColor = {0.5,0.5,0.5}; parser.setFixSizeVector<std::vector<double>> ( "Actor.DiffuseColor", DiffuseColor, false );
      parser.set("Light.useLight", 1 );
      parser.set("Light.KeyWarmth", 0.6 ); parser.set("Light.KeyIntensity", 0.75 );
      parser.set("Light.KeyElevation", 50 ); parser.set("Light.KeyAzimuth", 10 );
      parser.set("Light.FillWarmth", 0.4 ); parser.set("Light.FillKeyRatio", 3);
      parser.set("Light.FillElevation", -75); parser.set("Light.FillAzimuth", -10.);
      parser.set("Light.BackWarmth", 0.6); parser.set("Light.BackKeyRatio", 3.5);
      parser.set("Light.BackElevation", 0.); parser.set("Light.BackAzimuth", 110 );
      parser.set("Light.HeadWarmth", 0.5 ); parser.set("Light.HeadKeyRatio", 3.);
        
      const int numRGBPoints = 3;
      parser.set ( "ScalarBar.numRGBPoints", 3 );
      parser.set ( "ScalarBar.RGBValue1", 0. );
      std::vector<double> RGBPoint1 = {0.21960784313,0.29803921568,0.75294117647}; 
      parser.setFixSizeVector<std::vector<double>> ( "ScalarBar.RGBPoint1", RGBPoint1, false );
      parser.set ( "ScalarBar.RGBValue2", 0.5 );
      std::vector<double> RGBPoint2 = { 0.86666666666,0.86666666666,0.86666666666}; 
      parser.setFixSizeVector<std::vector<double>> ( "ScalarBar.RGBPoint2", RGBPoint2, false );
      parser.set ( "ScalarBar.RGBValue3", 1. );
      std::vector<double> RGBPoint3 = { 0.70588235294,0.0156862745,0.14901960784}; 
      parser.setFixSizeVector<std::vector<double>> ( "ScalarBar.RGBPoint3", RGBPoint3, false );
      
      parser.set ( "ScalarBar.ColorSpaceType", 1);
      parser.set ( "ScalarBar.ScaleType", 1);
      parser.set ( "ScalarBar.scaleRGBPointsToRange", 0 );
      parser.set ( "ScalarBar.plotScalarBarSeparately", 1 );
      parser.set ( "ScalarBar.NumberOfLabels", 0 ); 
          
      parser.set ("BoundaryBox.useBoundaryBox", 0 );
      parser.set ("ScalarBar.useScalarBar", false);
      std::vector<double> Position = {0.5,0.5,2.7320508075688776}; parser.setFixSizeVector<std::vector<double>> ( "Camera.Position", Position, false );
      std::vector<double> FocalPoint = { 0.5,0.5,0.0 }; parser.setFixSizeVector<std::vector<double>> ( "Camera.FocalPoint", FocalPoint, false );
      std::vector<double>ViewUp = {0.,0.,0.}; parser.setFixSizeVector<std::vector<double>> ( "Camera.ViewUp", ViewUp, false );
      parser.set ( "Camera.ParallelScale", 0.7071067811865476 );
      parser.set ( "Camera.Roll", 0. );
      parser.set ( "Camera.Azimuth", 0. );
      parser.set ( "Camera.Elevation", 0. );
      parser.set ( "Camera.Dolly", 1. );
      std::vector<double> Background = {1.0,1.0,1.0}; parser.setFixSizeVector<std::vector<double>> ( "Renderer.Background", Background, false );
      parser.set( "Mapper.InterpolateScalarsBeforeMapping", 1 );
      parser.set( "Window.FullScreenRendering", 1 );
      parser.set( "Writer.compressionLevel", 5 );
      
      this->template plotToPngWithParserInfo<pesopt::BoostParser> ( inputFileName, outputFileName, parser, useScalarData, nameScalarData, supp );
    }
#endif  //PESOPT_WITH_BOOST

};

#endif //__VTKtoPNGPlotter_H
