#ifndef __VTKMESHSAVER_H
#define __VTKMESHSAVER_H

#include <pesopt_IO.h>
#include <includesVTK.h>
#include <VTKSaver.h>



//! This class is a container for a mesh plus data vectors on vertices (VERTEX_DATA) or faces (FACE_DATA)
//! vector-valued data can be saved as vectors or normals
//! We do not store any of the data vectors here, but only keep pointers to them.
template <class MeshType>
class VTKMeshSaver {
public:
  typedef typename MeshType::RealType RealType;
  typedef typename MeshType::PointType PointType;
  typedef typename MeshType::VectorType VectorType;
  

protected:
  struct ScalarData {
    string             _descr;
    const VectorType * _data;
  };

  struct VectorData {
    string             _descr;
    const VectorType * _data;
    const int          _numComponents;
  };

  std::vector<ScalarData> _scalarVertexData;
  std::vector<VectorData> _vectorVertexData;
  std::vector<VectorData> _normalVertexData;
  std::vector<ScalarData> _scalarFaceData;
  std::vector<VectorData> _vectorFaceData;
  std::vector<VectorData> _normalFaceData;
  
  const MeshType &_mesh;
  int _precision;   
  
public:
  
  VTKMeshSaver ( const MeshType & mesh ) : _mesh ( mesh ), _precision( 8 ) {}
  VTKMeshSaver ( const MeshType & mesh, int precision ) : _mesh ( mesh ), _precision( precision ) {}

  
  VTKMeshSaver & clearData ( ) {
      _scalarVertexData.clear();
      _vectorVertexData.clear();
      _normalVertexData.clear();
      _scalarFaceData.clear();
      _vectorFaceData.clear();
      _normalFaceData.clear();
      return *this;
  }
  
//   VTKMeshSaver & updateMesh ( const MeshType & mesh ) {
//       _mesh = mesh;
//       return *this;
//   }
  
  VTKMeshSaver & addScalarData ( const VectorType & data, string dataDescr, VTKDataSupp supp ) {
    
    ScalarData entry = { dataDescr, &data };

    switch ( supp ) {
    case VERTEX_DATA:
      if ( data.size() != static_cast<unsigned>( _mesh.getNumVertices () ) ){
          cout << "data.size()           = " << data.size() << endl;
          cout << "mesh.getNumVertices() = " << _mesh.getNumVertices() << endl;
          throw std::invalid_argument( pesopt::strprintf ( "Wrong size. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      }
      _scalarVertexData.push_back ( entry );
      break;

    case FACE_DATA:
      if ( data.size() != static_cast<unsigned>( _mesh.getNumElements () ) ){
          cout << "data.size()           = " << data.size() << endl;
          cout << "mesh.getNumElements() = " << _mesh.getNumElements() << endl;
          throw std::invalid_argument( pesopt::strprintf ( "Wrong size. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      }
      _scalarFaceData.push_back ( entry );
      break;
      
    default:
      throw std::invalid_argument( pesopt::strprintf ( "Unknown VTKDataSupp. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    return *this;
  }

  VTKMeshSaver & addVectorData ( const VectorType & data, const int numComponents, string dataDescr, VTKDataSupp supp ) {

    VectorData entry = { dataDescr, &data, numComponents };

    switch ( supp ) {
        case VERTEX_DATA:
            _vectorVertexData.push_back ( entry );
            break;
        case FACE_DATA:
            _vectorFaceData.push_back ( entry );
        break;
        default:  
            throw std::invalid_argument( pesopt::strprintf ( "Unknown VTKDataSupp. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    return *this;
  }
  
  VTKMeshSaver & addNormalData ( const VectorType & data, const int numComponents, string dataDescr, VTKDataSupp supp ) {

    VectorData entry = { dataDescr, &data, numComponents };

    switch ( supp ) {
        case VERTEX_DATA:
            _normalVertexData.push_back ( entry );
            break;
        case FACE_DATA:
            _normalFaceData.push_back ( entry );
        break;
        default:  
            throw std::invalid_argument( pesopt::strprintf ( "Unknown VTKDataSupp. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    return *this;
  }
  

  //! set precision in saving methods
  void setPrecisionTo( int prec ){ _precision = prec; }

  void saveVTKPolydata ( string outputFileNameVTK, const PointType &offsetPoints ) const {

    auto data = vtkSmartPointer<vtkPolyData>::New();
    const int numVertices = _mesh.getNumVertices();
    const int numCells = _mesh.getNumElements();
    
    //points
    auto points = vtkSmartPointer<vtkPoints>::New();
    for( int nodeIter=0; nodeIter<numVertices; ++nodeIter ){ 
        const PointType& point = _mesh.getVertex ( nodeIter );
        if( point.size() == 1 )
            points->InsertNextPoint ( point[0] - offsetPoints[0], 0., 0. );
        if( point.size() == 2 )
            points->InsertNextPoint ( point[0] - offsetPoints[0], point[1] - offsetPoints[1], 0. );
        if( point.size() == 3 )
            points->InsertNextPoint ( point[0] - offsetPoints[0], point[1] - offsetPoints[1], point[2] - offsetPoints[2] );
    }
    
    //cells
    auto cells = vtkSmartPointer<vtkCellArray>::New();   
    for( int elIdx=0; elIdx<numCells; ++elIdx ){
        const int numNodesOfElement = _mesh.getNumNodesOfElement(elIdx);
        vtkIdType nodesOfElement [numNodesOfElement];
        for( int i=0; i<numNodesOfElement; ++i ) nodesOfElement[i] = _mesh.getElementNodeIdx( elIdx, i );
        cells->InsertNextCell ( numNodesOfElement, nodesOfElement );
    }
    
    // Add the geometry and topology to the polydata
    data->SetPoints ( points );
    data->SetPolys ( cells );
    
      
    // scalar point data //TODO is saved as FIELD. How to save as POINTDATA?
    for (int i = 0; i < _scalarVertexData.size(); ++i) {
        auto pointDataVec = vtkSmartPointer<vtkDoubleArray>::New();
        pointDataVec->SetName( (_scalarVertexData[i]._descr).c_str() );
        pointDataVec->SetNumberOfTuples(numVertices);
        for( int nodeIdx=0; nodeIdx<(*_scalarVertexData[i]._data).size(); ++nodeIdx) 
            pointDataVec->SetValue( nodeIdx, (*_scalarVertexData[i]._data)[nodeIdx] );
        data->GetPointData()->AddArray(pointDataVec);
    }
    
    
    // vector point data
    for (int i = 0; i < _vectorVertexData.size(); ++i) {
        auto pointDataVec = vtkSmartPointer<vtkDoubleArray>::New();
        const int numComponents = _vectorVertexData[i]._numComponents;
        pointDataVec->SetNumberOfComponents( numComponents );
        pointDataVec->SetNumberOfTuples(numVertices);
        pointDataVec->SetName( (_vectorVertexData[i]._descr).c_str() );
        for( int nodeIdx=0; nodeIdx<numVertices; ++nodeIdx){
            double pN[numComponents];
            for( int comp=0; comp<numComponents; ++comp ){
              pN[comp] = (*_vectorVertexData[i]._data)(nodeIdx + comp * numVertices);
            }
            pointDataVec->SetTuple( nodeIdx, pN );
        }
        data->GetPointData()->AddArray(pointDataVec);
    }
    
    // point normals
    if( _normalVertexData.size() > 0 ){
        auto pointNormalsArray = vtkSmartPointer<vtkDoubleArray>::New();
        pointNormalsArray->SetNumberOfComponents(3); //3d normals (ie x,y,z)
        pointNormalsArray->SetNumberOfTuples(numVertices);
        for( int nodeIdx = 0; nodeIdx < numVertices; ++nodeIdx ){
            double pN[3] = {(*_normalVertexData[0]._data)(nodeIdx + 0 * numVertices), 
                            (*_normalVertexData[0]._data)(nodeIdx + 1 * numVertices), 
                            (*_normalVertexData[0]._data)(nodeIdx + 2 * numVertices)};
            pointNormalsArray->SetTuple(nodeIdx, pN) ;
        }
        data->GetPointData()->SetNormals(pointNormalsArray);
    }
    
    //scalar cell data
    for (int i = 0; i < _scalarFaceData.size(); ++i) {
        auto cellDataVec = vtkSmartPointer<vtkDoubleArray>::New();
        cellDataVec->SetName( (_scalarFaceData[i]._descr).c_str() );
        cellDataVec->SetNumberOfTuples(numCells);
        for( int elIdx=0; elIdx<numCells; ++elIdx ) 
            cellDataVec->SetValue( elIdx, (*_scalarFaceData[i]._data)[elIdx] );
        data->GetCellData()->AddArray(cellDataVec);
    }
    
        VTKSaver vtkWriter;
        string file_extension_output;
        vtkWriter.saveVTKDataSet( outputFileNameVTK, data, file_extension_output, VTKPOLYDATA  );
  }

  
  void saveVTKStructuredGrid ( string outputFileNameVTK, const PointType &offsetPoints ) const {

    auto data = vtkSmartPointer<vtkStructuredGrid>::New();
    const int numVertices = _mesh.getNumVertices();
    const int numCells = _mesh.getNumElements();
    
    data->SetDimensions( _mesh.getNumDofs(0), _mesh.getNumDofs(1), _mesh.getNumDofs(2) );
  
    //points
    auto points = vtkSmartPointer<vtkPoints>::New();
    for( int nodeIter=0; nodeIter<numVertices; ++nodeIter ){ 
        const PointType& point = _mesh.getVertex ( nodeIter );
        if( point.size() == 1 )
            points->InsertNextPoint ( point[0] - offsetPoints[0], 0., 0. );
        if( point.size() == 2 )
            points->InsertNextPoint ( point[0] - offsetPoints[0], point[1] - offsetPoints[1], 0. );
        if( point.size() == 3 )
            points->InsertNextPoint ( point[0] - offsetPoints[0], point[1] - offsetPoints[1], point[2] - offsetPoints[2] );
    }
    data->SetPoints ( points );
    
    // scalar point data //TODO is saved as FIELD. How to save as POINTDATA?
    for (int i = 0; i < _scalarVertexData.size(); ++i) {
        auto pointDataVec = vtkSmartPointer<vtkDoubleArray>::New();
        pointDataVec->SetName( (_scalarVertexData[i]._descr).c_str() );
        pointDataVec->SetNumberOfTuples(numVertices);
        for( int nodeIdx=0; nodeIdx<(*_scalarVertexData[i]._data).size(); ++nodeIdx) 
            pointDataVec->SetValue( nodeIdx, (*_scalarVertexData[i]._data)[nodeIdx] );
        data->GetPointData()->AddArray(pointDataVec);
    }

    // vector point data
    for (int i = 0; i < _vectorVertexData.size(); ++i) {
        auto pointDataVec = vtkSmartPointer<vtkDoubleArray>::New();
        const int numComponents = _vectorVertexData[i]._numComponents;
        pointDataVec->SetNumberOfComponents( numComponents );
        pointDataVec->SetNumberOfTuples(numVertices);
        pointDataVec->SetName( (_vectorVertexData[i]._descr).c_str() );
        for( int nodeIdx=0; nodeIdx<numVertices; ++nodeIdx){
            double pN[numComponents];
            for( int comp=0; comp<numComponents; ++comp ){
              pN[comp] = (*_vectorVertexData[i]._data)(nodeIdx + comp * numVertices);
            }
            pointDataVec->SetTuple( nodeIdx, pN );
        }
        data->GetPointData()->AddArray(pointDataVec);
    }

//     //scalar cell data
//     for (int i = 0; i < _scalarFaceData.size(); ++i) {
//         auto cellDataVec = vtkSmartPointer<vtkDoubleArray>::New();
//         cellDataVec->SetName( (_scalarFaceData[i]._descr).c_str() );
//         cellDataVec->SetNumberOfTuples(numCells);
//         for( int elIdx=0; elIdx<numCells; ++elIdx ) 
//             cellDataVec->SetValue( elIdx, (*_scalarFaceData[i]._data)[elIdx] );
//         data->GetCellData()->AddArray(cellDataVec);
//     }
    
        VTKSaver vtkWriter;
        string file_extension_output;
        vtkWriter.saveVTKDataSet( outputFileNameVTK, data, file_extension_output, VTKSTRUCTUREDGRID  );
    
  }

  void saveVTKStructuredPoints ( string outputFileNameVTK, const PointType &offsetPoints ) const {

    auto data = vtkSmartPointer<vtkStructuredPoints>::New();
    const int numVertices = _mesh.getNumVertices();
    const int numCells = _mesh.getNumElements();
    
    data->SetDimensions( _mesh.getNumDofs(0), _mesh.getNumDofs(1), _mesh.getNumDofs(2) );
    data->SetSpacing( _mesh.getMeshSize(0), _mesh.getMeshSize(1), _mesh.getMeshSize(2) );
    
    const int dim = offsetPoints.size();
    if( dim == 1 ) data->SetOrigin( -1. * offsetPoints[0], 0., 0. );
    if( dim == 2 ) data->SetOrigin( -1. * offsetPoints[0], -1. * offsetPoints[1], 0. );
    if( dim == 3 ) data->SetOrigin( -1. * offsetPoints[0], -1. * offsetPoints[1], -1. * offsetPoints[2] );

      
    // scalar point data //TODO is saved as FIELD. How to save as POINTDATA?
    for (int i = 0; i < _scalarVertexData.size(); ++i) {
        auto pointDataVec = vtkSmartPointer<vtkDoubleArray>::New();
        pointDataVec->SetName( (_scalarVertexData[i]._descr).c_str() );
        pointDataVec->SetNumberOfTuples(numVertices);
        for( int nodeIdx=0; nodeIdx<(*_scalarVertexData[i]._data).size(); ++nodeIdx) pointDataVec->SetValue( nodeIdx, (*_scalarVertexData[i]._data)[nodeIdx] );
        data->GetPointData()->AddArray(pointDataVec);
    }

    // vector point data
    for (int i = 0; i < _vectorVertexData.size(); ++i) {
        auto pointDataVec = vtkSmartPointer<vtkDoubleArray>::New();
        const int numComponents = _vectorVertexData[i]._numComponents;
        pointDataVec->SetNumberOfComponents( numComponents );
        pointDataVec->SetNumberOfTuples(numVertices);
        pointDataVec->SetName( (_vectorVertexData[i]._descr).c_str() );
        for( int nodeIdx=0; nodeIdx<numVertices; ++nodeIdx){
            double pN[numComponents];
            for( int comp=0; comp<numComponents; ++comp ){
              pN[comp] = (*_vectorVertexData[i]._data)(nodeIdx + comp * numVertices);
            }
            pointDataVec->SetTuple( nodeIdx, pN );
        }
        data->GetPointData()->AddArray(pointDataVec);
    }

    //scalar cell data
    for (int i = 0; i < _scalarFaceData.size(); ++i) {
        auto cellDataVec = vtkSmartPointer<vtkDoubleArray>::New();
        cellDataVec->SetName( (_scalarFaceData[i]._descr).c_str() );
        cellDataVec->SetNumberOfTuples(numCells);
        for( int elIdx=0; elIdx<numCells; ++elIdx ) 
            cellDataVec->SetValue( elIdx, (*_scalarFaceData[i]._data)[elIdx] );
        data->GetCellData()->AddArray(cellDataVec);
    }
    
        VTKSaver vtkWriter;
        string file_extension_output;
        vtkWriter.saveVTKDataSet( outputFileNameVTK, data, file_extension_output, VTKSTRUCTUREDPOINTS  );
    
  }
  
  
  void saveVTKUnstructuredGrid ( string outputFileNameVTK, const PointType &offsetPoints ) const {

    auto data = vtkSmartPointer<vtkUnstructuredGrid>::New();
    const int numVertices = _mesh.getNumVertices();
    const int numCells = _mesh.getNumElements();
    
    const int dim = offsetPoints.size();

    //points
    auto points = vtkSmartPointer<vtkPoints>::New();
    for( int nodeIter=0; nodeIter<numVertices; ++nodeIter ){ 
        const PointType& point = _mesh.getVertex ( nodeIter );
        if( dim == 1 ) points->InsertNextPoint ( point[0] - offsetPoints[0], 0., 0. );
        if( dim == 2 ) points->InsertNextPoint ( point[0] - offsetPoints[0], point[1] - offsetPoints[1], 0. );
        if( dim == 3 ) points->InsertNextPoint ( point[0] - offsetPoints[0], point[1] - offsetPoints[1], point[2] - offsetPoints[2] );
    }
    
    //cells
    auto cells = vtkSmartPointer<vtkCellArray>::New();   
    for( int elIdx=0; elIdx<numCells; ++elIdx ){
        const int numNodesOfElement = _mesh.getNumNodesOfElement(elIdx);
        vtkIdType nodesOfElement [numNodesOfElement];
        for( int i=0; i<numNodesOfElement; ++i ) nodesOfElement[i] = _mesh.getElementNodeIdx( elIdx, i );
        cells->InsertNextCell ( numNodesOfElement, nodesOfElement );
    }
    
    // Add the geometry and topology to the polydata
    data->SetPoints ( points );
    data->SetCells (  _mesh._VTKCELLTYPE, cells );
    
      
    // scalar point data //TODO is saved as FIELD. How to save as POINTDATA?
    for (int i = 0; i < _scalarVertexData.size(); ++i) {
        auto pointDataVec = vtkSmartPointer<vtkDoubleArray>::New();
        pointDataVec->SetName( (_scalarVertexData[i]._descr).c_str() );
        pointDataVec->SetNumberOfTuples(numVertices);
        for( int nodeIdx=0; nodeIdx<(*_scalarVertexData[i]._data).size(); ++nodeIdx) pointDataVec->SetValue( nodeIdx, (*_scalarVertexData[i]._data)[nodeIdx] );
        data->GetPointData()->AddArray(pointDataVec);
    }
    
    //scalar cell data
    for (int i = 0; i < _scalarFaceData.size(); ++i) {
        auto cellDataVec = vtkSmartPointer<vtkDoubleArray>::New();
        cellDataVec->SetName( (_scalarFaceData[i]._descr).c_str() );
        cellDataVec->SetNumberOfTuples(numCells);
        for( int elIdx=0; elIdx<numCells; ++elIdx ) 
            cellDataVec->SetValue( elIdx, (*_scalarFaceData[i]._data)[elIdx] );
        data->GetCellData()->AddArray(cellDataVec);
    }
    
    
    
    // vector point data
    for (int i = 0; i < _vectorVertexData.size(); ++i) {
        auto pointDataVec = vtkSmartPointer<vtkDoubleArray>::New();
        const int numComponents = _vectorVertexData[i]._numComponents;
        pointDataVec->SetNumberOfComponents( numComponents );
        pointDataVec->SetNumberOfTuples(numVertices);
        pointDataVec->SetName( (_vectorVertexData[i]._descr).c_str() );
        for( int nodeIdx=0; nodeIdx<numVertices; ++nodeIdx){
            double pN[numComponents];
            for( int comp=0; comp<numComponents; ++comp ){
              pN[comp] = (*_vectorVertexData[i]._data)(nodeIdx + comp * numVertices);
            }
            pointDataVec->SetTuple( nodeIdx, pN );
        }
        data->GetPointData()->AddArray(pointDataVec);
    }
    
    VTKSaver vtkWriter;
    string file_extension_output;
    vtkWriter.saveVTKDataSet( outputFileNameVTK, data, file_extension_output, VTKUNSTRUCTUREDGRID  );
    
  }


  
  void save ( const string outputFileNameVTK, const PointType &offsetPoints, const VTKDATATYPE dataType ) const {
      switch( dataType ){
          
          case VTKUNSTRUCTUREDGRID:{
              this->saveVTKUnstructuredGrid( outputFileNameVTK, offsetPoints );
          }break;
          
          case VTKPOLYDATA:{
              this->saveVTKPolydata( outputFileNameVTK, offsetPoints );
          }break;
          
          case VTKSTRUCTUREDGRID:{
              this->saveVTKStructuredGrid( outputFileNameVTK, offsetPoints );
          }break;
          
          case VTKSTRUCTUREDPOINTS:{
              this->saveVTKStructuredPoints( outputFileNameVTK, offsetPoints );
          }break;
          
          default : throw invalid_argument( pesopt::strprintf ( "function save in VTKMeshSaver is not implemented for given dataType = %d. In File %s at line %d.", dataType, __FILE__, __LINE__ ).c_str() ); 
          
      }
  }
  
  void save ( const string outputFileNameVTK, const VTKDATATYPE dataType ) const {
      PointType zeroOffset; zeroOffset.setZero();
      this->save( outputFileNameVTK, zeroOffset, dataType );
  }
  

};

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

template< typename MeshType >
class VTKDeformedMeshSaver{
  
  typedef typename MeshType::RealType RealType;
  typedef typename MeshType::RealVecChart RealVecChart;
  typedef typename MeshType::VectorType VectorType;
  
  const MeshType &_mesh;
  const int _numVertices;
  mutable VectorType _identityDeformVec;
  
public:
  
  VTKDeformedMeshSaver( const MeshType &mesh ) : 
  _mesh( mesh ), 
  _numVertices ( mesh.getNumVertices() ), 
  _identityDeformVec( static_cast<int>( MeshType::dimChartDomain ) * _numVertices ) {
    for ( int nodeIdx = 0; nodeIdx < _numVertices; ++nodeIdx ) {
      auto coords = _mesh.getVertex( nodeIdx );
      for( int comp = 0; comp < coords.size(); ++comp ) 
          _identityDeformVec(nodeIdx + comp * _numVertices ) = coords(comp);   
    }
  } 
 
 
  // TODO maybe use additional function xA?
  //  this is for example needed for shells
//   VTKDeformedMeshSaver( const MeshType &mesh,  const VectorType &xA ) : 
//   _mesh( mesh ), _numVertices ( mesh.getNumVertices() ),  _identityDeformVec( xA ) { }
 
protected:
 
 void getDeformedMesh( const VectorType & dispOrDeformVec, 
                       const string deformType,
                       MeshType &meshDeformed, 
                       VectorType &deformVec
                    ) const{
    
    if( deformType == "deformation" ) {
        deformVec = dispOrDeformVec;
    } else if( deformType == "displacement" ) { 
        deformVec = dispOrDeformVec + _identityDeformVec;
    } else if( deformType == "undeformed" ) {
        deformVec = _identityDeformVec; 
    }else {
       throw std::invalid_argument( pesopt::strprintf ( "Wrong channel for deformType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
    for( int nodeIdx = 0; nodeIdx < _numVertices; ++nodeIdx ){
        RealVecChart coords;
        for( int comp = 0; comp < 3; ++comp ) coords[comp] = deformVec[nodeIdx + comp * _numVertices];
        meshDeformed.setVertex( nodeIdx, coords );
    }

  }
 
  void getDeformedMesh( const VectorType & dispOrDeformVec, 
                        const string deformType,
                        MeshType &meshDeformed
                      ) const{
    VectorType deformVec ( dispOrDeformVec.size() );
    this->getDeformedMesh( dispOrDeformVec,  deformType,  meshDeformed,  deformVec );
 }
 
    
public:
    
  void save ( const VectorType & dispOrDeformVec, 
              const string deformType, 
              const string outputFileNameVTK ) const{
    MeshType meshDeformed ( _mesh );
    this->getDeformedMesh(  dispOrDeformVec,  deformType, meshDeformed );
    VTKMeshSaver<MeshType> meshSaver ( meshDeformed );
    if ( deformType ==  "undeformed" ) 
        meshSaver.save( outputFileNameVTK, _mesh._VTKDATATYPEUNDEFORMED );
    else
        meshSaver.save( outputFileNameVTK, _mesh._VTKDATATYPEDEFORMED );
  }
    
};



  

    
    
    
    

#endif
