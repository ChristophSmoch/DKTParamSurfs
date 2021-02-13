#ifndef __VTKSAVER_H
#define __VTKSAVER_H

#include <pesopt_IO.h>
 
#include <includesVTK.h>

// #define __DEBUGVTKSaver

 
class VTKSaver {
    
public :
    VTKSaver( ) {}
    
    
    template<class TWriter> 
    void saveXMLVTKFile(const string &outputFileName, const vtkSmartPointer<vtkDataSet> &dataSet ) const{
       auto vtkwriter = vtkSmartPointer<TWriter>::New();
       vtkwriter->SetFileName( outputFileName.c_str() );
       vtkwriter->SetInputData(dataSet);
       vtkwriter->SetDataModeToAscii(); //vtkwriter->SetDataModeToBinary();
       vtkwriter->Update();
       vtkwriter->Write();
    }
    
    template<class TWriter> 
    void saveLegacyVTKFile(const string &outputFileName, const vtkSmartPointer<vtkDataSet> &dataSet ) const{
       auto vtkwriter = vtkSmartPointer<TWriter>::New();
       vtkwriter->SetFileName( outputFileName.c_str() );
       vtkwriter->SetInputData(dataSet);
       vtkwriter->SetFileTypeToASCII();
       vtkwriter->Update();
       vtkwriter->Write();
    }
    
//     template<class TWriter> 
//     void saveXMLVTKFile(const string &outputFileName, const vtkSmartPointer<vtkObjectBase> &algoOut ) const{
//        auto vtkwriter = vtkSmartPointer<TWriter>::New();
//        vtkwriter->SetFileName( outputFileName.c_str() );
//        vtkwriter->SetInputConnection(algoOut);
//        vtkwriter->SetDataModeToAscii();
//        vtkwriter->Update();
//        vtkwriter->Write();
//     }
//     
//     template<class TWriter> 
//     void saveLegacyVTKFile(const string &outputFileName, const vtkSmartPointer<vtkObjectBase> &algoOut ) const{
//        auto vtkwriter = vtkSmartPointer<TWriter>::New();
//        vtkwriter->SetFileName( outputFileName.c_str() );
//        vtkwriter->SetInputConnection(algoOut);
//        vtkwriter->SetFileTypeToASCII();
//        vtkwriter->Update();
//        vtkwriter->Write();
//     }
    
    
    void saveVTKDataSet( const string &outputFileName, 
                          const vtkSmartPointer<vtkDataSet> &dataSet, 
                          string &extension, 
                          const VTKDATATYPE &dataType  ) const{
#ifdef __DEBUGVTKSaver
        cout << endl 
             << "write dataset to file " << endl 
             << outputFileName.c_str() << endl
             << "extension = " << extension.c_str() << endl;
#endif
        extension = pesopt::getFileNameExtension( outputFileName );
        if      (extension == ".vtu")  saveXMLVTKFile<vtkXMLUnstructuredGridWriter> (outputFileName, dataSet );
        else if (extension == ".vtp")  saveXMLVTKFile<vtkXMLPolyDataWriter> (outputFileName, dataSet ); 
        else if (extension == ".vts")  saveXMLVTKFile<vtkXMLStructuredGridWriter> (outputFileName, dataSet );
        else if (extension == ".vtr")  saveXMLVTKFile<vtkXMLRectilinearGridWriter> (outputFileName, dataSet );
        else if (extension == ".vti")  saveXMLVTKFile<vtkXMLImageDataWriter> (outputFileName, dataSet );
        else if (extension == ".vtk"){
            if     ( dataType == VTKSTRUCTUREDGRID   )  saveLegacyVTKFile<vtkStructuredGridWriter>   ( outputFileName, dataSet );
            else if( dataType == VTKPOLYDATA         )  saveLegacyVTKFile<vtkPolyDataWriter>         ( outputFileName, dataSet );
            else if( dataType == VTKUNSTRUCTUREDGRID )  saveLegacyVTKFile<vtkUnstructuredGridWriter> ( outputFileName, dataSet );
            else if( dataType == VTKSTRUCTUREDPOINTS )  saveLegacyVTKFile<vtkStructuredPointsWriter> ( outputFileName, dataSet );
            else if( dataType == VTKRECTILINEARGRID  )  saveLegacyVTKFile<vtkRectilinearGridWriter>  ( outputFileName, dataSet );
            else throw std::invalid_argument( pesopt::strprintf ( "Unknown VTK-DataType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        }
        else throw std::invalid_argument( pesopt::strprintf ( "In function saveVTKDataSet, for outputfile %s, got unknown VTK-FileExtension %s. In File %s at line %d.", outputFileName.c_str(), extension.c_str(), __FILE__, __LINE__ ).c_str() );    
#ifdef __DEBUGVTKSaver
        cout << "finished writting dataset to file " << endl << outputFileName.c_str() << endl << endl;
#endif
     }
       
//     void writeVTKAlgorithmOutput( const string &outputFileName, const vtkSmartPointer<vtkObjectBase> &algoOut, string &extension, const VTKDATATYPE &dataType  ) const{
//         extension = pesopt::getFileNameExtension( outputFileName );
//         if      (extension == ".vtu")  saveXMLVTKFile<vtkXMLUnstructuredGridWriter> (outputFileName, algoOut );
//         else if (extension == ".vtp")  saveXMLVTKFile<vtkXMLPolyDataWriter> (outputFileName, algoOut ); 
//         else if (extension == ".vts")  saveXMLVTKFile<vtkXMLStructuredGridWriter> (outputFileName, algoOut );
//         else if (extension == ".vtr")  saveXMLVTKFile<vtkXMLRectilinearGridWriter> (outputFileName, algoOut );
//         else if (extension == ".vti")  saveXMLVTKFile<vtkXMLImageDataWriter> (outputFileName, algoOut );
//         else if (extension == ".vtk"){
//             if     ( dataType == VTKSTRUCTUREDGRID   )  saveLegacyVTKFile<vtkStructuredGridWriter>   ( outputFileName, algoOut );
//             else if( dataType == VTKPOLYDATA         )  saveLegacyVTKFile<vtkPolyDataWriter>         ( outputFileName, algoOut );
//             else if( dataType == VTKUNSTRUCTUREDGRID )  saveLegacyVTKFile<vtkUnstructuredGridWriter> ( outputFileName, algoOut );
//             else if( dataType == VTKSTRUCTUREDPOINTS )  saveLegacyVTKFile<vtkStructuredPointsWriter> ( outputFileName, algoOut );
//             else if( dataType == VTKRECTILINEARGRID  )  saveLegacyVTKFile<vtkRectilinearGridWriter>  ( outputFileName, algoOut );
//             else throw std::invalid_argument( pesopt::strprintf ( "Unknown VTK-DataType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
//         }
//         else throw std::invalid_argument( pesopt::strprintf ( "Unknown VTK-FileExtension. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );    
//      }
    
    
    
    template<class TReader> 
    void readXMLVTKFile(const string &inputFileName, vtkSmartPointer<vtkDataSet> &dataSet ) const{
        auto reader = vtkSmartPointer<TReader>::New();
        reader->SetFileName(inputFileName.c_str());
        reader->Update();
        reader->GetOutput()->Register(reader);
        dataSet = vtkDataSet::SafeDownCast(reader->GetOutput());
    }
    
    template<class TReader> 
    void readLegacyVTKFile(const string &inputFileName, vtkSmartPointer<vtkDataSet> &dataSet, VTKDATATYPE &dataType ) const{
        auto reader = vtkSmartPointer<TReader>::New();
        reader->SetFileName(inputFileName.c_str());
        reader->Update();
        reader->GetOutput()->Register(reader);
        dataSet = vtkDataSet::SafeDownCast(reader->GetOutput());
        
        if      ( reader->IsFileUnstructuredGrid() ) dataType = VTKUNSTRUCTUREDGRID;
        else if ( reader->IsFilePolyData() )         dataType = VTKPOLYDATA;
        else if ( reader->IsFileStructuredPoints() ) dataType = VTKSTRUCTUREDPOINTS;
        else if ( reader->IsFileRectilinearGrid()  ) dataType = VTKRECTILINEARGRID;
        else if ( reader->IsFileStructuredGrid()   ) dataType = VTKSTRUCTUREDGRID;
    }
    
    void readVTKDataSet( const string &inputFileName, vtkSmartPointer<vtkDataSet> &dataSet, string &extension, VTKDATATYPE &dataType  ) const{
#ifdef __DEBUGVTKSaver
        cout << endl << "read dataset from file " << endl << inputFileName.c_str() << endl;
#endif
        extension = pesopt::getFileNameExtension( inputFileName );
        if      (extension == ".vtu") { readXMLVTKFile<vtkXMLUnstructuredGridReader> (inputFileName, dataSet ); dataType = VTKUNSTRUCTUREDGRID;}
        else if (extension == ".vtp") { readXMLVTKFile<vtkXMLPolyDataReader> (inputFileName, dataSet ); dataType = VTKPOLYDATA;}
        else if (extension == ".vts") { readXMLVTKFile<vtkXMLStructuredGridReader> (inputFileName, dataSet ); dataType = VTKSTRUCTUREDGRID;}
        else if (extension == ".vtr") { readXMLVTKFile<vtkXMLRectilinearGridReader> (inputFileName, dataSet ); dataType = VTKRECTILINEARGRID; }
        else if (extension == ".vti") { readXMLVTKFile<vtkXMLImageDataReader> (inputFileName, dataSet ); dataType = VTKIMAGEDATA; }
        else if (extension == ".vtk")   readLegacyVTKFile<vtkDataSetReader> (inputFileName, dataSet, dataType );
        else throw std::invalid_argument( pesopt::strprintf ( "Unknown VTK-FileExtension for file %s. In File %s at line %d.", inputFileName.c_str(), __FILE__, __LINE__ ).c_str() ); 
#ifdef __DEBUGVTKSaver
        cout << "finished reading dataset from file " << endl << inputFileName.c_str() << endl << endl;
#endif
     }
    
    
    
    template<typename PointType>
    void getPoints( const string &inputFileName, std::vector<PointType> &pointVec ) const {
        // Get all data from the file
        string file_extension; VTKDATATYPE dataType;
        vtkSmartPointer<vtkDataSet> data; this->readVTKDataSet( inputFileName, data, file_extension, dataType );
        vtkIdType idNumPointsInFile = data->GetNumberOfPoints();
        PointType dummyPoint; const int dim = dummyPoint.size(); 
        for(int i = 0; i < idNumPointsInFile; i++){
            double point[dim];
            data->GetPoint(i, point);
            PointType p; for( int j=0;j<dim; ++j ) p[j] = point[j];
            pointVec.push_back( p );
        }
    }
    
    template<typename Indices3DType>
    void getCells( const string &inputFileName, std::vector<Indices3DType> &triangleVec ) const {
        string file_extension; VTKDATATYPE dataType;
        vtkSmartPointer<vtkDataSet> data; this->readVTKDataSet( inputFileName, data, file_extension, dataType );
        vtkIdType idNumCellsInFile = data->GetNumberOfCells();
        for(vtkIdType i = 0; i < idNumCellsInFile; i++){
            vtkCell* cell = data->GetCell(i);
            vtkIdList *listOfPointIdx = cell->GetPointIds();
            Indices3DType trianglePointIdx;
            for( int j=0; j < listOfPointIdx->GetNumberOfIds(); ++j ){
                trianglePointIdx[j] = listOfPointIdx->GetId(j);
            }
            triangleVec.push_back( trianglePointIdx );
        }
    }
    
    void getFreeSizeCells( const string &inputFileName, std::vector<std::vector<int>> &cellVec ) const {
        string file_extension; VTKDATATYPE dataType; 
        vtkSmartPointer<vtkDataSet> data; this->readVTKDataSet( inputFileName, data, file_extension, dataType );
        vtkIdType idNumCellsInFile = data->GetNumberOfCells();
        for(vtkIdType i = 0; i < idNumCellsInFile; i++){
            vtkCell* cell = data->GetCell(i);
            vtkIdList *listOfPointIdx = cell->GetPointIds();
            std::vector<int> cellPointIdx( listOfPointIdx->GetNumberOfIds() );
            for( int j=0; j < listOfPointIdx->GetNumberOfIds(); ++j ) cellPointIdx[j] = listOfPointIdx->GetId(j);
            cellVec.push_back( cellPointIdx );
        }
    }
    
    template <typename RealType>
    bool getScalarPointDataVec ( const string &inputFileName, const string &arrayName, std::vector<RealType> &dataVec ) const {
        string file_extension; VTKDATATYPE dataType;
        vtkSmartPointer<vtkDataSet> data; this->readVTKDataSet( inputFileName, data, file_extension, dataType );
        vtkIdType idNumPointsInFile = data->GetNumberOfPoints(); 
        dataVec.resize( idNumPointsInFile );
        vtkSmartPointer<vtkPointData> pd = data->GetPointData();
#ifdef __DEBUGVTKSaver
        if (pd){
            cout << " contains point data with " << pd->GetNumberOfArrays() << " arrays." << endl;
            for (int i = 0; i < pd->GetNumberOfArrays(); i++){
                cout << "\t Array " << i << " is named "  << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL") << endl;
                //  vtkSmartPointer<vtkDataArray> dataArray = data->GetPointData()->GetArray(i);
            }
        }
#endif
        vtkSmartPointer<vtkDataArray> dataArray = data->GetPointData()->GetArray(arrayName.c_str());
        if(dataArray){
#ifdef __DEBUGVTKSaver
            cout << "found data array " << arrayName.c_str() << endl;
#endif
            vtkIdType numTuples = dataArray->GetNumberOfTuples();
            for (vtkIdType tupleIdx = 0; tupleIdx < numTuples; ++tupleIdx){
                dataVec[tupleIdx] =  dataArray->GetComponent(tupleIdx, 0);
//                 cout << dataArray->GetComponent(tupleIdx, 0) << ", ";        
            }
            return true;
        }else{
#ifdef __DEBUGVTKSaver
            cout << "did not find data array " << arrayName.c_str() << endl;
#endif
            return false;
        }
  }
    
    //TODO does not work for legacy file format: only detects first array
    template<typename VecType>
    bool getPointDataVec ( const string &inputFileName, const string &arrayName, std::vector<VecType> &dataVec ) const {
        string file_extension; VTKDATATYPE dataType;
        vtkSmartPointer<vtkDataSet> data; this->readVTKDataSet( inputFileName, data, file_extension, dataType );
        vtkIdType idNumPointsInFile = data->GetNumberOfPoints(); 
        dataVec.resize( idNumPointsInFile );
        vtkSmartPointer<vtkPointData> pd = data->GetPointData();
#ifdef __DEBUGVTKSaver
        if (pd){
            cout << " contains point data with " << pd->GetNumberOfArrays() << " arrays." << endl;
            for (int i = 0; i < pd->GetNumberOfArrays(); i++){
                cout << "\t Array " << i << " is named "  << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL") << endl;
                // vtkSmartPointer<vtkDataArray> dataArray = data->GetPointData()->GetArray(i);
            }
        }
#endif
        vtkSmartPointer<vtkDataArray> dataArray = data->GetPointData()->GetArray(arrayName.c_str());
        if(dataArray){
#ifdef __DEBUGVTKSaver
            cout << "found data array " << arrayName.c_str() << endl;
#endif
            vtkIdType numVectors = dataArray->GetNumberOfTuples();
            const int numComponents = dataArray->GetNumberOfComponents();
#ifdef __DEBUGVTKSaver
            cout << "numVectors = " << numVectors << endl;
            cout << "numComponents = " << numComponents << endl;
#endif
            for (vtkIdType tupleIdx = 0; tupleIdx < numVectors; ++tupleIdx)
                for( int compIdx = 0; compIdx < numComponents; ++compIdx ){
                    dataVec[tupleIdx].resize(numComponents);
                    dataVec[tupleIdx][compIdx] = dataArray->GetComponent(tupleIdx, compIdx);         
                }
            return true;
        }else{
#ifdef __DEBUGVTKSaver
            cout << "did not find data array " << arrayName.c_str() << endl;
#endif
            return false;
        }
  }
  
  //TODO does not work for legacy file format: only detects first array
  template<typename VecType>
  bool getCellDataVec ( const string &inputFileName, const string &arrayName, std::vector<VecType> &dataVec ) const {
        string file_extension; VTKDATATYPE dataType; 
        vtkSmartPointer<vtkDataSet> data; this->readVTKDataSet( inputFileName, data, file_extension, dataType );
        vtkIdType idNumCellsInFile = data->GetNumberOfCells();
        dataVec.resize( idNumCellsInFile );
        vtkSmartPointer<vtkCellData> cd = data->GetCellData();
#ifdef __DEBUGVTKSaver
        if (cd) cout << " contains cell data with " << cd->GetNumberOfArrays() << " arrays." << endl;
#endif
        vtkSmartPointer<vtkDataArray> dataArray = data->GetCellData()->GetArray(arrayName.c_str());
        if(dataArray){
#ifdef __DEBUGVTKSaver
            cout << "found data array " << arrayName.c_str() << endl;
#endif
            vtkIdType numVectors = dataArray->GetNumberOfTuples();
            int numComponents = dataArray->GetNumberOfComponents();
            for (vtkIdType tupleIdx = 0; tupleIdx < numVectors; ++tupleIdx)
                for( int compIdx = 0; compIdx < numComponents; ++compIdx ){
                    dataVec[tupleIdx][compIdx] = dataArray->GetComponent(tupleIdx, compIdx);         
                }
             return true;
        }else{
#ifdef __DEBUGVTKSaver
            cout << "did not find data array " << arrayName.c_str() << endl;
#endif
            return false;
        }
  }

//   template<typename VecType>
//   void getFieldDataVec ( const string &inputFileName, const string &arrayName, std::vector<VecType> &dataVec ) const {
//         
//         // Get all data from the file
//         string file_extension; VTKDATATYPE dataType;
//         vtkSmartPointer<vtkDataSet> data; this->readVTKDataSet( inputFileName, data, file_extension, dataType );
//         
//         //TODO
// //         vtkIdType idNumCellsInFile = data->GetNumberOfCells();
// //         dataVec.resize( idNumCellsInFile );
// 
//         vtkSmartPointer<vtkPointData> fd = data->GetFieldData();
//         if (fd) cout << " contains field data with " << fd->GetNumberOfArrays() << " arrays." << endl;
//         vtkSmartPointer<vtkDataArray> dataArray = data->GetFieldData()->GetArray(arrayName.c_str());
//         
//         if(dataArray){
//             cout << "found data array " << arrayName.c_str() << endl;
//             vtkIdType numVectors = dataArray->GetNumberOfTuples();
//             int numComponents = dataArray->GetNumberOfComponents();
//             if( numVectors > 1 ){
//               for (vtkIdType tupleIdx = 0; tupleIdx < numVectors; ++tupleIdx)
//                 for( int compIdx = 0; compIdx < numComponents; ++compIdx ){
//                     //TODO
// //                     dataVec[tupleIdx][compIdx] = dataArray->GetComponent(tupleIdx, compIdx);         
//                 }
//               }
//         }else cout << "did not find data array " << arrayName.c_str() << endl;
//   }

};







template <typename VectorType>    
void vtkSaveScalarPointData( const string inputFileNameVTK, const string outputFileNameVTK, 
                             const string name, const VectorType &vec ) {
  
    VTKSaver saver;
    
    vtkSmartPointer<vtkDataSet> dataSet;
    string file_extension_input;
    VTKDATATYPE dataType;
    saver.readVTKDataSet( inputFileNameVTK, dataSet, file_extension_input, dataType );
        
    vtkSmartPointer<vtkDoubleArray> vecVTK = vtkSmartPointer<vtkDoubleArray>::New();
    vecVTK->SetName( name.c_str() );
    vecVTK->SetNumberOfTuples(dataSet->GetNumberOfPoints());
    for( int i=0; i<vec.size(); ++i) vecVTK->SetValue( i, vec[i] );
    dataSet->GetPointData()->AddArray(vecVTK);
        
    string file_extension_output;
    saver.saveVTKDataSet( outputFileNameVTK, dataSet, file_extension_output, dataType );

}



// template <typename VectorType>    
// void vtkSaveScalarPointData( const string inputFileNameVTK, const string outputFileNameVTK, 
//                     const string name, const VectorType &stressVec ) {
//     
//         string file_extension_input = pesopt::getFileNameExtension( inputFileNameVTK );
//         vtkSmartPointer<vtkDataSet> data;
//         if( file_extension_input == ".vtp" ){
//             auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
//             reader->SetFileName(inputFileNameVTK.c_str());
//             reader->Update();
//             data->DeepCopy( reader->GetOutput() );
//         }
//         if( file_extension_input == ".vtk" ){
//             auto reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
//             reader->SetFileName(inputFileNameVTK.c_str());
//             reader->Update();
//             if( reader->IsFileUnstructuredGrid() ){ data = vtkSmartPointer<vtkUnstructuredGrid>::New(); }
//             if( reader->IsFilePolyData() ){         data = vtkSmartPointer<vtkPolyData>::New(); }
//             if( reader->IsFileStructuredPoints() ){ data = vtkSmartPointer<vtkStructuredPoints>::New();  }
//             data->DeepCopy(reader->GetOutput());
//         }
//         
//         
//         vtkSmartPointer<vtkDoubleArray> stressVecVTK = vtkSmartPointer<vtkDoubleArray>::New();
//         stressVecVTK->SetName( name.c_str() );
//         stressVecVTK->SetNumberOfTuples(data->GetNumberOfPoints());
//         for( int i=0; i<stressVec.size(); ++i) stressVecVTK->SetValue( i, stressVec[i] );
//         data->GetPointData()->AddArray(stressVecVTK);
//         
//         string file_extension_output  = pesopt::getFileNameExtension(outputFileNameVTK);
//         if( file_extension_output == ".vtp" ){
//             auto vtkwriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//             vtkwriter->SetFileName( outputFileNameVTK.c_str() );
//             vtkwriter->SetInputData(data);
//             vtkwriter->SetDataModeToAscii(); //   vtkwriter->SetDataModeToBinary();
//             vtkwriter->Update();
//             vtkwriter->Write();
//         }
//         if( file_extension_output == ".vtk" ){
//             auto vtkwriter = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
//             vtkwriter->SetFileName( outputFileNameVTK.c_str() );
//             vtkwriter->SetInputData(data);
//             vtkwriter->SetFileTypeToASCII();
//             vtkwriter->Update();
//             vtkwriter->Write();
//         } 
// 
// }



void vtkConvertPNGToStructuredPoints ( const string &inputImage, const string &outputFile ){
  auto reader =  vtkSmartPointer<vtkPNGReader>::New();
  reader->SetFileName(inputImage.c_str());
  
  auto convertFilter = vtkSmartPointer<vtkImageToStructuredPoints>::New();
  convertFilter->SetInputConnection(reader->GetOutputPort());
  convertFilter->Update();  
  
  auto vtkwriter = vtkSmartPointer<vtkStructuredPointsWriter>::New();
  vtkwriter->SetFileName( outputFile.c_str() );
  vtkwriter->SetInputData(convertFilter->GetStructuredPointsOutput() );
  vtkwriter->SetFileTypeToASCII();
  vtkwriter->Update();
  vtkwriter->Write();
    
}

void vtkConvertPNGToStructuredPointsWithDim ( const string &inputImage, const string &outputFile, int *dimensions ){
  auto reader =  vtkSmartPointer<vtkPNGReader>::New();
  reader->SetFileName(inputImage.c_str());
  
  auto convertFilter = vtkSmartPointer<vtkImageToStructuredPoints>::New();
  convertFilter->SetInputConnection(reader->GetOutputPort());
  convertFilter->Update();  
  convertFilter->GetStructuredPointsOutput()->GetDimensions( dimensions );
  
  auto vtkwriter = vtkSmartPointer<vtkStructuredPointsWriter>::New();
  vtkwriter->SetFileName( outputFile.c_str() );
  vtkwriter->SetInputData(convertFilter->GetStructuredPointsOutput() );
  vtkwriter->SetFileTypeToASCII();
  vtkwriter->Update();
  vtkwriter->Write();
    
}


#endif
