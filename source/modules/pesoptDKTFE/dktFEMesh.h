#ifndef __DKTFEMESH_H
#define __DKTFEMESH_H

#include <pesopt_IO.h>
#include <pesopt_VTK.h>
#include <dktFEBoundary.h>
    
    
template< typename DataTypeContainer, typename TriangleType >
class TriangMeshWithTangentSpace {
public:
  typedef TriangleType                                  ElementType;
  typedef typename DataTypeContainer::RealType          RealType;
  typedef typename DataTypeContainer::Point3DType       Point3DType;
  typedef Point3DType                                   PointType;
  typedef typename DataTypeContainer::TangentVecType    TangentVecType;
  typedef typename DataTypeContainer::Indices3DType     Indices3DType;
  typedef Indices3DType                                 IndicesOfElementType;
  typedef typename DataTypeContainer::VectorType        VectorType;
  typedef typename DataTypeContainer::MaskType          MaskType;
  typedef typename DataTypeContainer::Vec2i             Vec2i;

  static const bool _hasTangentSpace = true;
  
  static const int _VTKCELLTYPE = 5; //VTK_TRIANGLE
  
protected:
  std::vector< Point3DType > _vertexIterator;
  std::vector< TriangleType > _triangIterator;
  
  mutable std::vector<TangentVecType> _tangentSpaceVec1;
  mutable std::vector<TangentVecType> _tangentSpaceVec2;
  mutable std::vector<TangentVecType> _normalSpaceVec;
  mutable std::vector<Indices3DType> _neighbour; //neighbour elements of Elements
  
  mutable std::vector< std::vector<int>  > _vertexIteratorAdjacentTriangles;
  mutable std::vector< std::vector<int>  > _vertexIteratorAdjacentTrianglesLocalIdx;
public:
  //! Create empty TriangMeshWithTangentSpace
  TriangMeshWithTangentSpace ( ) : _vertexIterator(), _triangIterator()  { }
  
  TriangMeshWithTangentSpace ( const string& fileName, const int tangentSpaceType ){
      this->loadFromVTK ( fileName, tangentSpaceType );
      //       this->makeOrientationConsistent();
      if( tangentSpaceType == 1 ) this->generateApproximativeTangentSpaceAtNodes();
      if( tangentSpaceType == 2 ){
        const int numVertices = this->getNumVertices ();
        MaskType boundaryMask; boundaryMask.resize( 3 * numVertices, false );
        int numBoundaryNodes;
        //ShellBoundaryType = 100 means boundary = topological boundary
        const bool ShellFETypeC1Dofs = false;
        generateDirichletBoundaryMaskUponShellBoundaryType<TriangMeshWithTangentSpace<DataTypeContainer,TriangleType>, ShellFETypeC1Dofs> ( static_cast<ShellBoundaryType>( 100 ), *this, boundaryMask, numBoundaryNodes );   
        this->generateApproximativeTangentSpaceAtNodes(boundaryMask);
      }
      updateAllProjectionCoefficients();
      updateAllAdjacentElementsForNodes();
  }
  
  
 TriangMeshWithTangentSpace ( std::vector< Point3DType > vertices, std::vector< Indices3DType > triangles ){
      for( int i=0; i<vertices.size(); ++i ) this->pushBackVertex( vertices[i] );
      for( int i=0; i<triangles.size(); ++i ) this->pushBackTriang( triangles[i] );
      makeOrientationConsistent();
      generateApproximativeTangentSpaceAtNodes();
      updateAllProjectionCoefficients();
      updateAllAdjacentElementsForNodes();
  }
    
  virtual ~TriangMeshWithTangentSpace ( ) {}  //why???
    
public:
  int getNumVertices ( ) const { return ( _vertexIterator.size() );}
  int getNumTriangs ( ) const { return ( static_cast<int> ( _triangIterator.size() ) );}
  int getNumElements ( ) const { return ( static_cast<int> ( _triangIterator.size() ) );}
  
  //! insert new vertex and return global index
  int pushBackVertex ( const Point3DType newVertex ) {
    _vertexIterator.push_back ( newVertex );
    return getNumVertices() - 1;
  }
  //! insert new triangle and return global index
  int pushBackTriang ( const Indices3DType nodeIdx ) {
    int globalIdx = getNumTriangs();
    _triangIterator.push_back ( TriangleType( globalIdx, nodeIdx, _vertexIterator ) );
    return globalIdx;
  }

  const Point3DType& getVertex ( const int num ) const { return _vertexIterator[num];}
  void setVertex ( const int num, const Point3DType Arg ) { _vertexIterator[num] = Arg;}

  int getNeighbour ( const int elementID, const int acrossLocalNode ) const { return _neighbour[elementID][acrossLocalNode];}
  void setNeighbour ( const int elementID, const int acrossLocalNode, const int value ) const {  _neighbour[elementID][acrossLocalNode] = value;}

  const TriangleType& getTriang ( const int num ) const {return _triangIterator[num];}
  TriangleType & getTriang ( const int num ) {return _triangIterator[num];}
  void setTriang ( const int num, const TriangleType Arg ) { _triangIterator[num] = Arg;}

  int getTriangNodeIdx  ( const int num, const int localNode ) const { return _triangIterator[num].getGlobalNodeIdx(localNode);}
  int getElementNodeIdx  ( const int num, const int localNode ) const { return this->getTriangNodeIdx(num,localNode);}
  int getNumNodesOfElement(const int elIdx ) const { return 3; }
  void setTriangNodeIdx ( const int num, const int localNode, const int newIdx ) {_triangIterator[num].setGlobalNodeIdx(localNode, newIdx );}
  void setTriang ( const int num, const Indices3DType nodeIdx ) {
    TriangleType triang ( num, nodeIdx, _vertexIterator );
    _triangIterator[num] = triang;
  }
  
  
  RealType getTangentVec1 ( int numComp, int VertexIndex ) const { return _tangentSpaceVec1[VertexIndex][numComp]; }
  RealType getTangentVec2 ( int numComp, int VertexIndex ) const { return static_cast<RealType>( _tangentSpaceVec2[VertexIndex][numComp] );  }
  const TangentVecType & getTangentVec1 (const int& vertexIndex ) const{ return _tangentSpaceVec1[vertexIndex];}
  const TangentVecType & getTangentVec2 (const int& vertexIndex ) const { return _tangentSpaceVec2[vertexIndex];}
  const TangentVecType & getNormalVec (const int& vertexIndex ) const{ return _normalSpaceVec[vertexIndex];}
  void setTangentVec1 (const int& vertexIndex, const TangentVecType& vec ) { _tangentSpaceVec1[vertexIndex] = vec;}
  void setTangentVec2 (const int& vertexIndex, const TangentVecType& vec ) { _tangentSpaceVec2[vertexIndex] = vec;}
  void setNormalVec (const int& vertexIndex, const TangentVecType& vec ) { _normalSpaceVec[vertexIndex] = vec;}
  
  
  
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
  
  //slow: iterates over all vertices
  void findClosestVertex( const Point3DType &coords, int &vertexIndex ) const{
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
      //TODO very slow!!! use that the element indices are ordered --> dont need the full loop
      for(int i=0;i<_vertexIteratorAdjacentTriangles[nodeIdx1].size(); ++i )
           for(int j=0;j<_vertexIteratorAdjacentTriangles[nodeIdx2].size(); ++j ){
              if(_vertexIteratorAdjacentTriangles[nodeIdx1][i] == _vertexIteratorAdjacentTriangles[nodeIdx2][j] ) commonElements.push_back(_vertexIteratorAdjacentTriangles[nodeIdx1][i]); 
           }
  }
  
  void getCommonElements( const int nodeIdx1, const int nodeIdx2, const int nodeIdx3, std::vector<int> &commonElements ) const{
      //TODO very slow!!! use that the element indices are ordered
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
      throw std::invalid_argument ( pesopt::strprintf ( "getLocalIdxWrtGlobalElementIdx fails in file %s at line %d.", __FILE__, __LINE__ ).c_str() );
      return -1;
  }
  
  
  
  // refinement operations only update NodeIndex informations, hence have to update coords and projection coefficients of elements
  void updateAllTriangles ( ) {
    for ( int elementIndex = 0; elementIndex < getNumTriangs(); ++elementIndex  ) 
        _triangIterator[elementIndex].updateNodesAndEdges( _vertexIterator );
    for ( int elementIndex = 0; elementIndex < getNumTriangs(); ++elementIndex  ) 
        _triangIterator[elementIndex].updateProjectionCoefficients( _tangentSpaceVec1, _tangentSpaceVec2 );   
  }
  
  void updateAllProjectionCoefficients ( ) {
    for ( int elementIndex = 0; elementIndex < getNumTriangs(); ++elementIndex  ) _triangIterator[elementIndex].updateProjectionCoefficients( _tangentSpaceVec1, _tangentSpaceVec2 );   
  }
  
  void printBasicInfo( ) const{
      cout << endl << "Mesh has " << endl 
           << this->getNumTriangs()  << " Elements" << endl
           << this->getNumVertices() << " Vertices" << endl;
  }
  
  void print ( const bool printNodes = true, const bool printElements = true ) const {
    if( printNodes ){
      cout << endl << "print nodes" << endl;
      for ( int nodeIndex = 0; nodeIndex < getNumVertices(); ++nodeIndex ){
       cout << endl << "nodeIndex = " << nodeIndex << endl;
       cout << "coords = " << _vertexIterator[nodeIndex].transpose() << endl;
       cout << "tv1    = " << _tangentSpaceVec1[nodeIndex].transpose() << endl;   
       cout << "tv2    = " << _tangentSpaceVec2[nodeIndex].transpose() << endl; 
       cout << "n      = " << _normalSpaceVec[nodeIndex].transpose() << endl; 
       cout << "v1xv2  = " << (_tangentSpaceVec1[nodeIndex].cross( _tangentSpaceVec2[nodeIndex] ) ).transpose() << endl;
      }
      cout << endl << "print elements adjacent to nodes" << endl;
      for ( int nodeIndex = 0; nodeIndex < getNumVertices(); ++nodeIndex ){
       cout << endl << "nodeIndex = " << nodeIndex << endl;
       for(int i=0;i<_vertexIteratorAdjacentTriangles[nodeIndex].size(); ++i)
           cout << "Element " << _vertexIteratorAdjacentTriangles[nodeIndex][i] 
                << ", localNodeIdx = " << _vertexIteratorAdjacentTrianglesLocalIdx[nodeIndex][i] << endl;
      }
    }
    if( printElements ){
        cout << endl << "print elements" << endl;
        for ( int elementIndex = 0; elementIndex < getNumTriangs(); ++elementIndex  ) _triangIterator[elementIndex].print(); 
        cout << endl << "finished print elements" << endl;
    }
  }

public:
    
        
 void loadFromVTK( const string& filename, const int tangentSpaceType ) {
    VTKSaver vtkReader;
    int numVertices, numTriangles;
     
    std::vector<PointType> pointVec;
    vtkReader.getPoints<PointType>( filename, pointVec );
    numVertices = pointVec.size();    
    _vertexIterator.resize( numVertices );
    _tangentSpaceVec1.resize ( numVertices ); _tangentSpaceVec2.resize ( numVertices ); _normalSpaceVec.resize( numVertices );
    for( int i=0; i<numVertices; ++i ) this->setVertex ( i, pointVec[i] );

    std::vector<Indices3DType> triangleVec;
    vtkReader.getCells<Indices3DType>( filename, triangleVec );
    numTriangles = triangleVec.size();    
    _triangIterator.resize( numTriangles );
    for( int i = 0; i<numTriangles; ++i ) this->setTriang ( i, triangleVec[i] ); 

    if( tangentSpaceType == 0 ){
        std::vector<TangentVecType> tangentSpaceVec1, tangentSpaceVec2;
        vtkReader.getPointDataVec<TangentVecType>( filename, "tangentVec1", tangentSpaceVec1 );
        vtkReader.getPointDataVec<TangentVecType>( filename, "tangentVec2", tangentSpaceVec2 );
        
        for( int i = 0; i<numVertices; ++i ){
            this->setTangentVec1 ( i, tangentSpaceVec1[i] );
            this->setTangentVec2 ( i, tangentSpaceVec2[i] );
            TangentVecType normal; normal = tangentSpaceVec1[i].cross( tangentSpaceVec2[i] );
            normal.normalize();
            this->setNormalVec ( i, normal );
        }
    }
  }    
    
    
  void generateApproximativeTangentSpaceAtNodes (){
      MaskType zeroBoundaryMask ( this->getNumVertices() );
      for( int nodeIdx=0; nodeIdx<zeroBoundaryMask.size(); ++nodeIdx ) zeroBoundaryMask[nodeIdx] = false;
      this->generateApproximativeTangentSpaceAtNodes( zeroBoundaryMask );
  }

  // for boundary node only averge both neighboring triangle normals s.t. the corresponding triangle has another boundary node
  void generateApproximativeTangentSpaceAtNodes ( const MaskType &boundaryMask ){
        const int numVertices = this->getNumVertices ();
        _tangentSpaceVec1.resize ( numVertices ); _tangentSpaceVec2.resize ( numVertices ); _normalSpaceVec.resize( numVertices );
        std::fill( _normalSpaceVec.begin(), _normalSpaceVec.end(), TangentVecType(0.,0.,0.) );
        for ( int elementIndex = 0; elementIndex < this->getNumTriangs(); ++elementIndex  ) {
            TangentVecType normalizedNormal; _triangIterator[elementIndex].getNormalizedNormalForFlattenedTriangle( normalizedNormal );
            RealType areaOfElement = _triangIterator[elementIndex].getAreaOfFlattenedTriangle ();
            for (int idx = 0; idx < 3; ++idx){ 
                int globIdx = _triangIterator[elementIndex].getGlobalNodeIdx(idx);
                if( boundaryMask[globIdx] == false ){ _normalSpaceVec[ globIdx ] += areaOfElement * normalizedNormal;
                }else {
                    // node is on boundary, so ask if element has another boundary node
                    if( boundaryMask[ _triangIterator[elementIndex].getGlobalNodeIdx( (idx + 1) % 3) ] || boundaryMask[ _triangIterator[elementIndex].getGlobalNodeIdx( (idx + 2) % 3) ] )
                        _normalSpaceVec[ globIdx ] += areaOfElement * normalizedNormal;
                }
            }
        }
        //check if normal is non zero and normalize
        for (int vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex){
            RealType normOfNormal = _normalSpaceVec[vertexIndex].norm();
            if( normOfNormal == 0. ){
                cout << pesopt::color::red << "normal is zero for vertex " << vertexIndex << pesopt::color::reset << endl;
                if( boundaryMask[vertexIndex] ) cout << "vertex on boundary" << endl;
                else                            cout << "vertex not on boundary" << endl;
                std::vector<int> commonElements;
                this->getCommonElements( vertexIndex, commonElements );
                for ( int elementIdx = 0; elementIdx < commonElements.size(); ++elementIdx  ) {
                    TangentVecType normalizedNormal; _triangIterator[commonElements[elementIdx]].getNormalizedNormalForFlattenedTriangle( normalizedNormal );
                    RealType areaOfElement = _triangIterator[commonElements[elementIdx]].getAreaOfFlattenedTriangle();
                    cout << "normal of El " << commonElements[elementIdx] << " is " << normalizedNormal.transpose() << endl;
                    cout << "area of El is " << areaOfElement << endl;
                 }
                _normalSpaceVec[vertexIndex] = TangentVecType(sqrt(1./3.),sqrt(1./3.),sqrt(1./3.));
            }
            else _normalSpaceVec[vertexIndex].normalize();
        }
        for (int vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex) getTangentSpaceFromNormal ( _tangentSpaceVec1[vertexIndex], _tangentSpaceVec2[vertexIndex], _normalSpaceVec[vertexIndex] );
    }


    //Compute orthonormal basis (v1, v2, n ) of (tangent space + normal space);   //since |n|^2 = nx^2 + ny^2 + nz^2 = 1  => at least one |ni| >= 1/sqrt(3) > 0.5 
    void getTangentSpaceFromNormal ( TangentVecType &tv1, TangentVecType &tv2, const TangentVecType& normal ) const {
        if ( (std::abs ( normal(2) ) > 0.5) ) tv1 = TangentVecType( normal(2), 0., (-1.) * normal(0) );
        else                                  tv1 = TangentVecType( normal(1), (-1.) * normal(0), 0. );
        tv1.normalize();
        tv2 = normal.cross( tv1 );
        tv2.normalize();
    }


    void makeNeighbour() const {
        int  *n_to_t, *num_of_n, *start_of_n;
        const int numNodes = getNumVertices(); const int numElements = getNumTriangs();
        
        if ( int(_neighbour.size()) != numElements ) _neighbour.resize ( numElements );

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
    
    
    void findTopologicalBoundaryMask( MaskType & mask ) const {
        const int numVertices =  this->getNumVertices();
        mask.resize( numVertices );
        for( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) mask[nodeIdx] = false;
        this->makeNeighbour();
        for( int elementIndex=0; elementIndex < this->getNumElements(); ++ elementIndex )
            for( int localIdx = 0; localIdx < 3; ++localIdx ) 
                if ( this->getNeighbour( elementIndex , localIdx ) == -1 ){
                    mask[_triangIterator[elementIndex].getGlobalNodeIdx( (localIdx+2) % 3 )] = true;
                    mask[_triangIterator[elementIndex].getGlobalNodeIdx( (localIdx+1) % 3 )] = true;
                }
    }



    void makeOrientationConsistent() {
        if ( int(_neighbour.size()) != this->getNumTriangs() ) makeNeighbour();
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
            _triangIterator[elementIndex].updateNodesAndEdges( _vertexIterator );
    }
    
    
    //! for compatibility
    int getNumDofs( const int i ) const{
        throw std::invalid_argument( pesopt::strprintf ( "getNumDofs for coordinate direction does not make sense for TriangMeshWithTangentSpace. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        return 0;
    }
    
    RealType getMeshSize( const int i ) const{
         throw std::invalid_argument( pesopt::strprintf ( "getMeshSize for coordinate direction does not make sense for TriangMeshWithTangentSpace. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        return 0.;
    }

};



#endif
