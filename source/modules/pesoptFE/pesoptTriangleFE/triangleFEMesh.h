#ifndef __TRIANGLEMESH_H
#define __TRIANGLEMESH_H

#include <pesopt_IO.h>
#include <pesopt_VTK.h>
#include <feDefines.h>
#include "triangleFEDefines.h"
#include "triangleFEElement.h"
    
    
template< typename DataTypeContainer = DataTypeContainerTriangleFE,
          typename TriangleType = TriangleElement<DataTypeContainer>
//           typename EdgeType = TriangleEdge<DataTypeContainer> 
          >
class TriangleMesh {
public:
  typedef TriangleType                                      ElementType;
  typedef typename ElementType::ElementNodeIndicesType      ElementNodeIndicesType;
  typedef typename DataTypeContainer::RealType              RealType;
  typedef typename DataTypeContainer::RealVecChart          RealVecChart;
  typedef typename DataTypeContainer::IntVecChart           IntVecChart;
  typedef typename DataTypeContainer::PointType             PointType;
  typedef typename DataTypeContainer::VectorType            VectorType;
  typedef typename DataTypeContainer::MaskType              MaskType;
  typedef typename DataTypeContainer::TangentVecType        TangentVecType;  
  typedef pesopt::BoostParser ParameterParserType;
//   typedef typename EdgeType::EdgeNodeIndices                EdgeNodeIndices;

  static const int dimChartDomain = 2;
  static const int _VTKCELLTYPE = 5; //VTK_TRIANGLE
  static const VTKDATATYPE _VTKDATATYPEUNDEFORMED = VTKPOLYDATA;
  static const VTKDATATYPE _VTKDATATYPEDEFORMED =  VTKPOLYDATA;
  static const PESOPT_FEMeshType _FEMeshType = TRIANGLE;
  
protected:
  std::vector< PointType > _vertexIterator;
  std::vector< TriangleType > _triangIterator;
//   std::vector< EdgeType > _edgeIterator;
  
  mutable MaskType _boundaryMask;
  
  //neighbour elements of Elements
  mutable std::vector< ElementNodeIndicesType > _neighbourElementsOfElements; 
  
  mutable std::vector< std::vector<int>  > _vertexIteratorAdjacentTriangles;
  mutable std::vector< std::vector<int>  > _vertexIteratorAdjacentTrianglesLocalIdx;
public:
  
  //! Create empty TriangleMesh
//   TriangleMesh ( ) : 
//   _vertexIterator(),
//   _triangIterator()
// //   _edgeIterator() 
//   { }
  
  TriangleMesh ( const string& fileName ){
      this->loadFromVTK ( fileName );
      this->makeOrientationConsistent();
      this->generateBoundaryMask ( );   
      this->updateAllAdjacentElementsForNodes();
  }
  
  TriangleMesh ( std::vector< PointType > vertices, std::vector< ElementNodeIndicesType > triangles ){
      for( int i=0; i<vertices.size(); ++i ) this->pushBackVertex( vertices[i] );
      for( int i=0; i<triangles.size(); ++i ) this->pushBackTriang( triangles[i] );
      this->makeOrientationConsistent();
      //makeEdges();
      this->generateBoundaryMask ( );    
      this->updateAllAdjacentElementsForNodes();
  }
    
public:
  int getNumVertices ( ) const { return ( _vertexIterator.size() );}
  int getNumTriangs ( ) const { return ( static_cast<int> ( _triangIterator.size() ) );}
  int getNumElements ( ) const { return ( static_cast<int> ( _triangIterator.size() ) );}
  
//   int getNumEdges ( ) const { return _edgeIterator.size() ; }
  
  //! insert new vertex and return global index
  int pushBackVertex ( const PointType newVertex ) {
    _vertexIterator.push_back ( newVertex );
    return getNumVertices() - 1;
  }
  //! insert new triangle and return global index
  int pushBackTriang ( const ElementNodeIndicesType nodeIdx ) {
    int globalIdx = getNumTriangs();
    _triangIterator.push_back ( TriangleType( globalIdx, nodeIdx, _vertexIterator ) );
    return globalIdx;
  }
  //! insert new edge and return global index
//   int pushBackEdge ( const EdgeNodeIndices nodeIndices ) {
//     int globalIdx = getNumEdges();
//     _edgeIterator.push_back ( EdgeType( globalIdx, nodeIndices, _vertexIterator ) );
//     return globalIdx;
//   }
  //! insert new triangle and return global index
  void pushBackNeighbour ( const ElementNodeIndicesType &nodeIdx ) {
    _neighbourElementsOfElements.push_back( nodeIdx );
  }
  
  const PointType& getVertex ( const int num ) const { return _vertexIterator[num];}
  void setVertex ( const int num, const PointType Arg ) { _vertexIterator[num] = Arg;}

  void getVertexVector ( VectorType &x ) const {
      const int numVertices = this->getNumVertices();
      for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
        const PointType& coords ( this->getVertex(nodeIdx) );
        for( int comp=0; comp<coords.size(); ++comp )  x[ nodeIdx + numVertices * comp ] = coords[comp];
      }
  }

  
  int getNeighbour ( const int elementID, const int acrossLocalNode ) const { return _neighbourElementsOfElements[elementID][acrossLocalNode];}
  void setNeighbour ( const int elementID, const int acrossLocalNode, const int value ) const {  _neighbourElementsOfElements[elementID][acrossLocalNode] = value;}

  const TriangleType& getElement ( const int num ) const {return _triangIterator[num];}
//   TriangleType & getElement ( const int num ) {return _triangIterator[num];}
  void setElement ( const int num, const TriangleType Arg ) { _triangIterator[num] = Arg;}

  int getTriangNodeIdx  ( const int num, const int localNode ) const { return _triangIterator[num].getGlobalNodeIdx(localNode);}
  int getElementNodeIdx  ( const int num, const int localNode ) const { return this->getTriangNodeIdx(num,localNode);}
  int getNumNodesOfElement(const int elIdx ) const { return 3; }
  void setTriangNodeIdx ( const int num, const int localNode, const int newIdx ) {_triangIterator[num].setGlobalNodeIdx(localNode, newIdx );}
  void setElement ( const int num, const ElementNodeIndicesType nodeIdx ) {
    TriangleType triang ( num, nodeIdx, _vertexIterator );
    _triangIterator[num] = triang;
  }

  
  void updateAllAdjacentElementsForNodes () const {
    _vertexIteratorAdjacentTriangles.clear();
    _vertexIteratorAdjacentTrianglesLocalIdx.clear();
    _vertexIteratorAdjacentTriangles.resize( this->getNumVertices() );
    _vertexIteratorAdjacentTrianglesLocalIdx.resize( this->getNumVertices() );
    for ( int elementIndex = 0; elementIndex < getNumTriangs(); ++elementIndex  )
        for(int localNodeIdx=0; localNodeIdx<3; ++localNodeIdx){
           const int globalNodeIdx = _triangIterator[elementIndex].getGlobalNodeIdx(localNodeIdx);
           _vertexIteratorAdjacentTriangles[globalNodeIdx].push_back(_triangIterator[elementIndex].getGlobalElementIdx());
           _vertexIteratorAdjacentTrianglesLocalIdx[globalNodeIdx].push_back(localNodeIdx);
        }
  }
  
  //WARNING slow: iterates over all vertices
  void findClosestVertex( const PointType &coords, int &vertexIndex ) const{
      vertexIndex = 0;
      RealType distSqr = ( coords - _vertexIterator[0] ).squaredNorm();
      for( int i=0; i<this->getNumVertices(); ++i ){
          RealType tmp = ( coords - _vertexIterator[i] ).squaredNorm();
          if( tmp < distSqr ){
              vertexIndex = i; distSqr = tmp;
          }
      }
  }
  
  void getCommonElements( const int nodeIdx1, std::vector<int> &commonElements ) const{
      for(int i=0;i<_vertexIteratorAdjacentTriangles[nodeIdx1].size(); ++i ) commonElements.push_back(_vertexIteratorAdjacentTriangles[nodeIdx1][i]); 
  }
  
  void getCommonElements( const int nodeIdx1, const int nodeIdx2, std::vector<int> &commonElements ) const{
      for(int i=0;i<_vertexIteratorAdjacentTriangles[nodeIdx1].size(); ++i )
           for(int j=0;j<_vertexIteratorAdjacentTriangles[nodeIdx2].size(); ++j ){
              if(_vertexIteratorAdjacentTriangles[nodeIdx1][i] == _vertexIteratorAdjacentTriangles[nodeIdx2][j] ) 
                  commonElements.push_back(_vertexIteratorAdjacentTriangles[nodeIdx1][i]); 
           }
  }
  
  void getCommonElements( const int nodeIdx1, const int nodeIdx2, const int nodeIdx3, std::vector<int> &commonElements ) const{
      for(int i=0;i<_vertexIteratorAdjacentTriangles[nodeIdx1].size(); ++i )
           for(int j=0;j<_vertexIteratorAdjacentTriangles[nodeIdx2].size(); ++j )
               for(int k=0;k<_vertexIteratorAdjacentTriangles[nodeIdx3].size(); ++k ){
              if(    (_vertexIteratorAdjacentTriangles[nodeIdx1][i] == _vertexIteratorAdjacentTriangles[nodeIdx2][j] )
                  && (_vertexIteratorAdjacentTriangles[nodeIdx1][i] == _vertexIteratorAdjacentTriangles[nodeIdx3][k] ) ) commonElements.push_back(_vertexIteratorAdjacentTriangles[nodeIdx1][i]); 
           }
  }
  
  
  void getCommonElements( const std::vector<int> &nodeIdxVec, std::vector<int> &commonElements ) const{
      const int numNodes = nodeIdxVec.size();
      switch( numNodes ){
//         case 1: getCommonElements(nodeIdxVec[0],commonElements); break;
        case 2: getCommonElements(nodeIdxVec[0],nodeIdxVec[1],commonElements); break;
        case 3: getCommonElements(nodeIdxVec[0],nodeIdxVec[1],nodeIdxVec[2],commonElements); break;
          
        default : 
                throw std::invalid_argument ( pesopt::strprintf ( "wrong number of nodes to get common elements in file %s at line %d.", __FILE__, __LINE__ ).c_str() );
                break;
      }
  }
  
  int getLocalIdxWrtAdjacentElementIdx( const int nodeIdx, const int adjacentElementIdx ) const {
        return _vertexIteratorAdjacentTrianglesLocalIdx[nodeIdx][adjacentElementIdx];
  }
  
  int getLocalIdxWrtGlobalElementIdx( const int nodeIdx, const int elementIdx ) const {
      for( int adjacentElementIdx=0; adjacentElementIdx < _vertexIteratorAdjacentTriangles[nodeIdx].size(); ++adjacentElementIdx ){
       if( _vertexIteratorAdjacentTriangles[nodeIdx][adjacentElementIdx] == elementIdx) 
           return _vertexIteratorAdjacentTrianglesLocalIdx[nodeIdx][adjacentElementIdx];
      }
  }
  
  
  RealType getMinArea( ) const{
      RealType minarea = 1.e+100;
      for ( int elementIndex = 0; elementIndex < this->getNumTriangs(); ++elementIndex  ){
          RealType area = this->getElement( elementIndex ).getVolume();
          if( area < minarea ) minarea = area;
      }
      return minarea ;
  }

  RealType getInterfaceWith ( ) const {
    return std::sqrt( this->getMinArea( ) );
  }
  
  
  
  //! iterates over all triangles and collects all edges
//   void makeEdges() {
//     _edgeIterator.clear();
//     for ( int elementIdx = 0; elementIdx < this->getNumElements(); ++elementIdx )
//       for ( int locNodeIdx = 0; locNodeIdx < 3; locNodeIdx++ ) {
//         bool isNewEdge = true;
//         EdgeNodeIndices edge (  this->getTriangNodeIdx( elementIdx, ( locNodeIdx + 1 ) % 3 ), 
//                       this->getTriangNodeIdx( elementIdx, locNodeIdx % 3 ) );
//         EdgeNodeIndices edgeOpposite ( this->getTriangNodeIdx( elementIdx, locNodeIdx % 3 ), 
//                              this->getTriangNodeIdx( elementIdx, (locNodeIdx + 1 ) % 3 ) );
//         // iterate over all edges and check if already an edge
//         int oldEdgeIndex = -1;
//         for ( int k = 0; k < _edgeIterator.size(); k++ ) {
//           if      ( (edge(0) == _edgeIterator[k].getAdjacentNodeIdx(0)) && (edge(1) == _edgeIterator[k].getAdjacentNodeIdx(1)) ) { 
//               isNewEdge = false; oldEdgeIndex = k; break;} 
//           else if ( (edgeOpposite(0) == _edgeIterator[k].getAdjacentNodeIdx(0)) && (edgeOpposite(1) == _edgeIterator[k].getAdjacentNodeIdx(1)) ) { isNewEdge = false; oldEdgeIndex = k; break;}
//         }
//         if ( isNewEdge ) {
//             const int newEdgeIndex = this->pushBackEdge( edge );
//             _edgeIterator[newEdgeIndex].setOppositeNodeIdx( 0, this->getTriangNodeIdx( elementIdx, ( locNodeIdx + 2 ) % 3 ) );
//             _edgeIterator[newEdgeIndex].setTriangIdx( 0, elementIdx );
//         }else{
//             _edgeIterator[oldEdgeIndex].setOppositeNodeIdx( 1, this->getTriangNodeIdx( elementIdx, ( locNodeIdx + 2 ) % 3 ) );
//             _edgeIterator[oldEdgeIndex].setTriangIdx( 1, elementIdx );
//         }
//       }
//   }
  
  
  // refinement operations only update NodeIndex informations, hence have to update coords, edges
  void updateAllTriangles ( ) {
    for ( int elementIndex = 0; elementIndex < getNumTriangs(); ++elementIndex  ) 
        _triangIterator[elementIndex].update( _vertexIterator );
  }
  
  void printBasicInfo( ) const{
      cout << endl << "Mesh has " << endl 
           << this->getNumTriangs()  << " Elements" << endl
           << this->getNumVertices() << " Vertices" << endl;
//            << this->getNumEdges()    << " Edges" << endl;
  }

public:
    
        
 void loadFromVTK( const string& filename ) {
    VTKSaver vtkReader;
    int numVertices, numTriangles;
     
    std::vector<PointType> pointVec;
    vtkReader.getPoints<PointType>( filename, pointVec );
    numVertices = pointVec.size();    
    _vertexIterator.resize( numVertices );
    for( int i=0; i<numVertices; ++i ) this->setVertex ( i, pointVec[i] );

    std::vector<ElementNodeIndicesType> triangleVec;
    vtkReader.getCells<ElementNodeIndicesType>( filename, triangleVec );
    numTriangles = triangleVec.size();    
    _triangIterator.resize( numTriangles );
    for( int i = 0; i<numTriangles; ++i ) this->setElement ( i, triangleVec[i] ); 
    
  }    
    
    

  void makeNeighbour() const {
        int  *n_to_t, *num_of_n, *start_of_n;
        const int numNodes = getNumVertices(); const int numElements = getNumTriangs();
        
        if ( int(_neighbourElementsOfElements.size()) != numElements ) _neighbourElementsOfElements.resize ( numElements );

        n_to_t     = new int[numElements*3];
        num_of_n   = new int[numNodes];
        start_of_n = new int[numNodes+1];

        // iterate over all vertices and set num_of_n to zero, i.e. initialize
        for ( int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++ ) num_of_n[nodeIndex] = 0;

        // iterate over all triangles and all neighbours, i.e. 3, since every triangle has 3 neighbours for closed meshes
        for ( int elementIndex = 0; elementIndex < numElements; elementIndex++ )
            for ( int j = 0; j < 3; j++ ){
                num_of_n[getTriangNodeIdx ( elementIndex,j ) ]++;
                setNeighbour ( elementIndex, j, -1 );
            }
            
        // get Startindex for node
        int nodeIdx = 0; start_of_n[nodeIdx++] = 0;
        while ( nodeIdx < numNodes ) {
            start_of_n[nodeIdx] = start_of_n[nodeIdx-1] + num_of_n[nodeIdx-1];
            nodeIdx++;
        }
        start_of_n[numNodes] = 3 * numElements;

        // initialize reference
        for ( int elementIndex = 0; elementIndex < numElements*3; ++elementIndex ) n_to_t[elementIndex] = -1;

        // build reference
        for ( int elementIndex = 0; elementIndex < numElements; elementIndex++ )
            for ( int j = 0; j < 3; j++ ) {
            int k = start_of_n[getTriangNodeIdx ( elementIndex,j ) ];
            while ( n_to_t[k] > -1 ) k++;
            n_to_t[k] = elementIndex;
            }
        // find neighbour
        for ( int elementIndex = 0; elementIndex < numElements; elementIndex++ )
            for ( int j = 0; j < 3; j++ ) {
            int node1 = getTriangNodeIdx ( elementIndex, j        );
            int node2 = getTriangNodeIdx ( elementIndex, ( j + 1 ) % 3 );
            int node3 = getTriangNodeIdx ( elementIndex, ( j + 2 ) % 3 );

            for ( int k = start_of_n[node1]; k < start_of_n[node1+1]; k++ ) {
                int n = n_to_t[k]; // Tetraeder 
                if ( elementIndex < n ) // set neighborhood only once
                for ( int l = 0; l < 3; l++ ) {
                    if ( node3 == getTriangNodeIdx ( n, l ) ) {
                    setNeighbour ( elementIndex, ( j + 1 ) % 3, n );
                    for ( int v = 0; v < 3; v++ )
                        if ( v != l && getTriangNodeIdx ( n, v ) != node1 ) setNeighbour ( n, v, elementIndex );
                    }
                    else {
                    if ( node2 == getTriangNodeIdx ( n, l ) ) {
                        setNeighbour ( elementIndex, ( j + 2 ) % 3, n );
                        for ( int v = 0; v < 3; v++ )
                        if ( v != l && getTriangNodeIdx ( n, v ) != node1 ) setNeighbour ( n, v, elementIndex );
                    }
                    }
                }
            }
        }

        delete[] n_to_t; delete[] num_of_n; delete[] start_of_n;
  }



  void makeOrientationConsistent() {
        if ( int(_neighbourElementsOfElements.size()) != this->getNumTriangs() ) makeNeighbour();
        MaskType alreadyHandled( this->getNumTriangs() ); // true for all triangles T, whose neighboring triangles have already been oriented consistent with T
        MaskType inQueue( this->getNumTriangs() ); // true for all triangles who have already been dealt with or who are already waiting in queue
        std::queue<int> activeTriangles; // contains all triangles which are already oriented and whose neighbors will be dealt with next (may contain triangles twice for simplicity)
        activeTriangles.push( 0 );
        inQueue[0] = true;
        int currentTriangle; // the triangle whose neighbors are currently handled
        while( !activeTriangles.empty() ){
            currentTriangle = activeTriangles.front();
            activeTriangles.pop();
            // deal with all three neighbors of currentTriangle, i.e. orient them and add them to the list to deal with their neighbors
            for ( int i = 0; i < 3; i++ ){
                int neighbor = getNeighbour( currentTriangle, i );
                if ( neighbor >= 0 && neighbor < this->getNumTriangs() && !alreadyHandled[neighbor] ){
                    // compute the nodes "currentTriangle" and "neighbor" have in common
                    int node1 = getTriangNodeIdx ( currentTriangle, ( i + 1 ) % 3 );
                    int node2 = getTriangNodeIdx ( currentTriangle, ( i + 2 ) % 3 );
                    // check whether common nodes occur in reversed order in "neighbor", if not, change order
                    int j = 0; while ( getTriangNodeIdx ( neighbor, j ) != node2 ) j++;
                    if ( getTriangNodeIdx ( neighbor, ( j + 1 ) % 3 ) != node1 ){
                    // change order of nodes
                    int exchangeCache = getTriangNodeIdx ( neighbor, ( j + 1 ) % 3 );
                    setTriangNodeIdx( neighbor, ( j + 1 ) % 3, getTriangNodeIdx ( neighbor, ( j + 2 ) % 3 ) );
                    setTriangNodeIdx( neighbor, ( j + 2 ) % 3, exchangeCache );
                    // change order of corresponding neighbours
                    exchangeCache = getNeighbour( neighbor, ( j + 1 ) % 3 );
                    setNeighbour( neighbor, ( j + 1 ) % 3, getNeighbour( neighbor, ( j + 2 ) % 3 ) );
                    setNeighbour( neighbor, ( j + 2 ) % 3, exchangeCache );
                    }
                    if ( !inQueue[neighbor] ){
                        activeTriangles.push( neighbor );
                        inQueue[neighbor] = true;
                    }
                }
            }
            alreadyHandled[currentTriangle] = true;
        }
        
        for ( int elementIndex = 0; elementIndex < getNumTriangs(); ++elementIndex  ) 
            _triangIterator[elementIndex].update( _vertexIterator );
  }
    
    
  //! for compatibility
  int getNumDofs( const int i ) const{
        throw std::invalid_argument( pesopt::strprintf ( "getNumDofs for coordinate direction does not make sense for TriangleMesh. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        return 0;
  }
    
  RealType getMeshSize( const int i ) const{
         throw std::invalid_argument( pesopt::strprintf ( "getMeshSize for coordinate direction does not make sense for TriangleMesh. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        return 0.;
  }
    
    
  void generateBoundaryMask( ) const {
        _boundaryMask.resize( this->getNumVertices(), false );
        this->makeNeighbour();
        for( int elementIndex=0; elementIndex < this->getNumElements(); ++elementIndex )
            for( int localIdx = 0; localIdx < 3; ++localIdx ) 
                if ( this->getNeighbour( elementIndex , localIdx ) == -1 ){
                    _boundaryMask[_triangIterator[elementIndex].getGlobalNodeIdx( (localIdx+2) % 3 )] = true;
                    _boundaryMask[_triangIterator[elementIndex].getGlobalNodeIdx( (localIdx+1) % 3 )] = true;
                }
        //count numBoundaryNodes
//         for( int nodeIdx=0; nodeIdx < this->getNumVertices(); ++nodeIdx ){
//                if( mask[nodeIdx] ) ++numBoundaryNodes;    
//         }
  }
    
  const MaskType &getBoundaryMask( ) const { return _boundaryMask; }
    
    
  void saveToParser( ParameterParserType & parser, const string fileNameMesh ) const {
      parser.set ( "InputMesh.dimChartDomain", static_cast<int>(dimChartDomain) );
      parser.set ( "InputMesh.FEType", "Triangle" );
      parser.set ( "InputMesh.MeshType", "File" );
      parser.set ( "InputMesh.fileName", fileNameMesh ); 
  }

};


#endif
