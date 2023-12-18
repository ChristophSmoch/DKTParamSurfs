#ifndef __VTKINTERFACETOSURFACE_H
#define __VTKINTERFACETOSURFACE_H

#include <VTKSaver.h>


class VTKInterfaceExtractor {
    
public :
    VTKInterfaceExtractor( ) {}
    
    //! \note  inputFile has to be vtkImageData (ie. either structuredpoint or uniformgrid)
    //         outputFile is polydata (ie. either .vtk with polydata or .vtp)
    bool getInterfaceByMarchingCubes( const string &inputFileName, const string &outputFileName, 
                                      const string &ActiveScalars, const float threshold = 0.5 ){
        
        VTKSaver vtkReader;
        string file_extension; VTKDATATYPE dataType;
        vtkSmartPointer<vtkDataSet> dataSet; vtkReader.readVTKDataSet( inputFileName, dataSet, file_extension, dataType );
        dataSet->GetPointData()->SetActiveScalars(ActiveScalars.c_str());
        
        if( (dataType == VTKSTRUCTUREDPOINTS) || (dataType == VTKSTRUCTUREDGRID ) || (dataType == VTKIMAGEDATA ) ){
          // using marching cubes algorithm
            auto mc = vtkSmartPointer<vtkMarchingCubes>::New();
            //mc->SetInputConnection(reader->GetOutputPort());
            mc->SetInputData(dataSet);
            mc->ComputeNormalsOn();
            //mc->ComputeGradientsOn();
            mc->ComputeScalarsOff();
            mc->SetValue(0, threshold);

            //save 
            auto writer = vtkSmartPointer<vtkPolyDataWriter>::New(); //.vtk
            writer->SetFileName(outputFileName.c_str());
            writer->SetFileTypeToASCII();         
            writer->SetInputConnection(mc->GetOutputPort());
            writer->Write();
//             VTKSaver vtkWriter;
//             string file_extension_output; //VTKDATATYPE dataType = VTKPOLYDATA;
//             /*vtkWriter.saveVTKDataSet( outputFileName, mc->GetOutput(), file_extension_output, dataType );*/ 
//             vtkWriter.writeVTKAlgorithmOutput( outputFileName, mc->GetOutputPort(), file_extension_output, VTKPOLYDATA ); 
            
            return true;
        }else{
            cout << "in getInterfaceByMarchingCubes: input is not of type vtkImageData" << endl;   
            return false;
        }

    }
    
    
    //         outputFile is polydata (ie. either .vtk with polydata or .vtp)
    void getInterfaceByContourFilter( const string &inputFileName, const string &outputFileName,
                                      const string &ActiveScalars, const float threshold = 0.5,
                                      const int extractionMode = 6, const bool colorRegions = false ){

        cout << "load" << endl;
        VTKSaver vtkReader;
        string file_extension; VTKDATATYPE dataType;
        vtkSmartPointer<vtkDataSet> dataSet; vtkReader.readVTKDataSet( inputFileName, dataSet, file_extension, dataType );
        dataSet->GetPointData()->SetActiveScalars(ActiveScalars.c_str());
        
      // using vtkContourFilter
        cout << "setup contourfilter" << endl;
        auto mc = vtkSmartPointer<vtkContourFilter>::New();
        //mc->SetInputConnection(reader->GetOutputPort());
        mc->SetInputData(dataSet);
        mc->ComputeNormalsOn();
        //mc->ComputeGradientsOn();
        mc->ComputeScalarsOff();
        mc->SetValue(0, threshold);
        
        auto connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
        connectivityFilter->SetInputConnection(mc->GetOutputPort());
        switch( extractionMode ){
            case 1 : connectivityFilter->SetExtractionModeToPointSeededRegions(); break;   
            case 2 : connectivityFilter->SetExtractionModeToCellSeededRegions(); break;
            case 3 : connectivityFilter->SetExtractionModeToLargestRegion(); break;
            case 4 : connectivityFilter->SetExtractionModeToSpecifiedRegions(); break;
            case 5 : connectivityFilter->SetExtractionModeToClosestPointRegion(); break;
            case 6 : connectivityFilter->SetExtractionModeToAllRegions(); break;     
            default: throw invalid_argument( pesopt::strprintf ( "Unknown extractionMode. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );  break;
        }
        
        if( colorRegions ) connectivityFilter->ColorRegionsOn();
        
        connectivityFilter->Update();

        //save 
        VTKSaver vtkWriter;
        string file_extension_output; //VTKDATATYPE dataType = VTKPOLYDATA;
        vtkWriter.saveVTKDataSet( outputFileName, connectivityFilter->GetOutput(), file_extension_output, VTKPOLYDATA  );
    }
    
      
      
      
//     void getTriangulatedInterface( const string &inputFileName, const string &outputFileName, const float threshold = 0.5 ){
//         VTKSaver vtkReader;
//         string file_extension; VTKDATATYPE dataType;
//         vtkSmartPointer<vtkDataSet> dataSet; vtkReader.readVTKDataSet( inputFileName, dataSet, file_extension, dataType );
//         
//         // using vtkContourFilter
//         auto mc = vtkSmartPointer<vtkContourFilter>::New();
// //         mc->SetInputConnection(reader->GetOutputPort());
//         mc->SetInputData( dataSet );
//         mc->ComputeNormalsOn();
//         //mc->ComputeGradientsOn();
//         mc->ComputeScalarsOff();
//         mc->SetValue(0, threshold);
//         
//         auto poly = vtkSmartPointer<vtkContourTriangulator>::New();
//         poly->SetInputConnection(mc->GetOutputPort());
//         
//         VTKSaver vtkWriter;
//         string file_extension_output;
//         vtkWriter.saveVTKDataSet( outputFileName, poly->GetOutput(), file_extension_output, dataType  );
//     }
    
    
    void getRegionByClipping( const string &inputFileName, const string &outputFileName, 
                              const string &ActiveScalars, const float threshold = 0.5, const bool useComplement = false ){
        VTKSaver vtkReader;
        string file_extension; VTKDATATYPE dataType;
        vtkSmartPointer<vtkDataSet> dataSet; vtkReader.readVTKDataSet( inputFileName, dataSet, file_extension, dataType );
        dataSet->GetPointData()->SetActiveScalars(ActiveScalars.c_str());
        
        auto clipper = vtkSmartPointer<vtkClipPolyData>::New();
        //old         clipper->SetInputConnection(reader->GetOutputPort());
        clipper->SetInputData(dataSet);
        clipper->SetValue(threshold);
        // clipper->GenerateClipScalarsOn();
        if( useComplement ) clipper->InsideOutOn();
        else                clipper->InsideOutOff();
        clipper->Update();
        
        //save 
        VTKSaver vtkWriter;
        string file_extension_output; // VTKDATATYPE dataType = VTKPOLYDATA;
        vtkWriter.saveVTKDataSet( outputFileName, clipper->GetOutput(), file_extension_output, dataType  );
        
    }
    
    
    
    
    //! \note inputfile has to be of type polydata 
    //        VTKSaver vtkReader;
    //    string file_extension; VTKDATATYPE dataType;
    //    vtkSmartPointer<vtkDataSet> dataSet; vtkReader.readVTKDataSet( inputFileName, dataSet, file_extension, dataType );
    void loopSubdivision( const string inputFileNameVTK, const string outputFileNameVTK, const int numberOfSubdivisions ) const{
        
        string file_extension_input  = pesopt::getFileNameExtension(inputFileNameVTK);

        vtkSmartPointer<vtkPolyData> originalMesh;
        if( file_extension_input == ".vtp" ){
            auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
            reader->SetFileName(inputFileNameVTK.c_str());
            reader->Update();
            // Subdivision filters only work on triangles
            auto triangles =  vtkSmartPointer<vtkTriangleFilter>::New();
            triangles->SetInputConnection(reader->GetOutputPort());
            triangles->Update();
            originalMesh = triangles->GetOutput();
        }
        if( file_extension_input == ".vtk" ){
               auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
               reader->SetFileName(inputFileNameVTK.c_str());
               reader->Update();
               
               // Subdivision filters only work on triangles
               auto triangles =  vtkSmartPointer<vtkTriangleFilter>::New();
               triangles->SetInputConnection(reader->GetOutputPort());
               triangles->Update();
               originalMesh = triangles->GetOutput();
        }
    
        auto subdivisionFilter = vtkSmartPointer<vtkLoopSubdivisionFilter>::New();
        subdivisionFilter->SetNumberOfSubdivisions(numberOfSubdivisions);
        subdivisionFilter->SetInputData(originalMesh);
        subdivisionFilter->Update();   
        
        VTKSaver vtkWriter;
        string file_extension_output; VTKDATATYPE dataType = VTKPOLYDATA;
        vtkWriter.saveVTKDataSet( outputFileNameVTK, subdivisionFilter->GetOutput(), file_extension_output, dataType  );
    }
    
    

};


#endif
