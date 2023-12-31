#ifndef __DKTFEADAPTIVEMESHWITHTANGENTSPACE_H
#define __DKTFEADAPTIVEMESHWITHTANGENTSPACE_H

#include <pesopt_IO.h>
#include <dktFEMesh.h>


template< typename DataTypeContainer, typename TriangleType >
class AdaptiveTriangMeshWithTangentSpace : public TriangMeshWithTangentSpace<DataTypeContainer, TriangleType> {
  
public:   
  typedef TriangleType ElementType;
  typedef typename DataTypeContainer::RealType       RealType;
  typedef typename DataTypeContainer::DomVecType     DomVecType;
  typedef typename DataTypeContainer::Matrix33       Matrix33;
  typedef typename DataTypeContainer::Point3DType    Point3DType;
  typedef typename DataTypeContainer::Indices3DType  Indices3DType;
  typedef typename DataTypeContainer::VectorType     VectorType;
  typedef typename DataTypeContainer::MaskType       MaskType;
  typedef typename DataTypeContainer::TangentVecType TangentVecType;
  typedef std::vector< Indices3DType >               VertexIndicesType;
  typedef std::vector< Point3DType >                 VertexCoordsType;
  
  typedef struct { int globalIndices[2];} ParentInformation;
  
  typedef short LocalIndex;
  typedef int GlobalIndex;
  static const int IndexNotSet = -1;
  
protected:  
// half-edge iterator (only for internal refinement!)
class DartIterator{
  
  typedef AdaptiveTriangMeshWithTangentSpace<DataTypeContainer, TriangleType> MeshType;
    
  protected:
    const MeshType & _mesh;
    GlobalIndex _triangle;
    LocalIndex  _node;
    LocalIndex  _edge;

  public:
    DartIterator( const MeshType & Mesh, GlobalIndex triangleIndex, LocalIndex localEdgeIndex, LocalIndex localNodeIndex ) : 
    _mesh( Mesh ), _triangle (triangleIndex), _node(localNodeIndex), _edge(localEdgeIndex){}    
    DartIterator( const MeshType & Mesh, GlobalIndex triangleIndex, LocalIndex localEdgeIndex ) : 
    _mesh( Mesh ), _triangle (triangleIndex), _node((localEdgeIndex + 1) % 3), _edge(localEdgeIndex){}

    void set ( const GlobalIndex triangle, const LocalIndex edge, const LocalIndex node )  { _triangle = triangle; _node = node; _edge = edge;}
    GlobalIndex getGlobalTriangleIndex() const { return _triangle; }
    GlobalIndex getGlobalNodeIndex()     const { return _mesh.getTriangNodeIdx( _triangle, _node );}
    LocalIndex  getLocalNodeIndex()      const { return _node; }
    LocalIndex  getLocalEdgeIndex()      const { return _edge;}
    LocalIndex  getNextNodeLocalIndex()  const { return 3 - (_node + _edge);}
    LocalIndex  getNextEdgeLocalIndex()  const { return 3 - (_node + _edge); }  
    GlobalIndex getNextNodeGlobalIndex() const { return _mesh.getTriangNodeIdx( _triangle, getNextNodeLocalIndex() );}
    GlobalIndex getNextTriangleIndex()   const { return _mesh.getNeighbour( _triangle, _edge ); }
 
    // returns local index of common node (first checking if "d" also refers to _triangle, but to a different node)
    template<typename DartIteratorType>
    LocalIndex getCommonNodeLocalIndex( const DartIteratorType & d ) const {
      if(_triangle != d.getGlobalTriangleIndex()) throw std::invalid_argument( pesopt::strprintf("T=%d, but this->tringle = %d. In File %s at line %d.", d.getGlobalTriangleIndex(), _triangle, __FILE__, __LINE__).c_str() );
      if ( _edge == d.getLocalEdgeIndex() ) throw std::invalid_argument( pesopt::strprintf("points to same edge=. In File %s at line %d.", __FILE__, __LINE__).c_str() );
      return 3 - ( getLocalEdgeIndex() + d.getLocalEdgeIndex() );
    }

    // returns global index of common node (also if "d" refers to a different triangle than _triangle )
    template<typename DartIteratorType>
    GlobalIndex getCommonNodeGlobalIndex( const DartIteratorType & d ) const {
      if ( getGlobalNodeIndex() == d.getGlobalNodeIndex() || getGlobalNodeIndex() == d.getNextNodeGlobalIndex() ) return getGlobalNodeIndex();
      if ( getNextNodeGlobalIndex() == d.getGlobalNodeIndex() || getNextNodeGlobalIndex() == d.getNextNodeGlobalIndex() ) return getNextNodeGlobalIndex();
      throw std::invalid_argument( pesopt::strprintf("Dart d and *this have no node in common. In File %s at line %d.", __FILE__, __LINE__).c_str() );
      return -1;
    }
    
    void print() const{ cout << "triangle = " << _triangle << ", node = " << _node << ", edge = " << _edge << endl; }

    // does _triangle have a neighbour across _edge?
    bool canFlipTriangle() const { return this->getNextTriangleIndex() != IndexNotSet ; }

    // moves to neighbouring triangle across _edge
    void flipTriangle() {
          if( !canFlipTriangle() ) throw std::invalid_argument( pesopt::strprintf("Cannot flip triangle. In File %s at line %d.", __FILE__, __LINE__).c_str() );
        
          GlobalIndex newTriangleIndex = getNextTriangleIndex();
          GlobalIndex globalNodeIndex = _mesh.getTriangNodeIdx( _triangle, _node );    
          GlobalIndex globalOtherNodesIndex = _mesh.getTriangNodeIdx( _triangle, 3 - (_node + _edge) );

          // find out which node of the new triangle is ours
          LocalIndex newNodeIndex = IndexNotSet;
          for (LocalIndex i = 0; i < 3; ++i){
            if ( _mesh.getTriangNodeIdx( newTriangleIndex, i) == globalNodeIndex ){ newNodeIndex = i; break;}
          }
            
          // find out which node of the new triangle lies on our edge
          LocalIndex newOtherNodesIndex = IndexNotSet;
          for (LocalIndex i = 1; i < 3; ++i){
            if ( _mesh.getTriangNodeIdx( newTriangleIndex,(newNodeIndex + i) % 3 ) == globalOtherNodesIndex ){ newOtherNodesIndex = (newNodeIndex + i) % 3; break; }
          }

          if ( min( newNodeIndex, newOtherNodesIndex ) < 0 ) 
              throw std::invalid_argument( pesopt::strprintf( "newTriang should neighbour of _triangle but has no common node .In File %s at line %d.", __FILE__, __LINE__ ).c_str() );

          LocalIndex newEdgeIndex = (-(newNodeIndex + newOtherNodesIndex)) % 3;
          if (newEdgeIndex < 0) newEdgeIndex += 3;

          // all information found, now set:
          _triangle = newTriangleIndex; _node = newNodeIndex; _edge = newEdgeIndex; 
    }
    
    // moves to other node along _edge (inside current triangle)
    void flipNode() { _node = getNextNodeLocalIndex();}
    // moves to other edge with _node (inside current triangle)
    void flipEdge() {_edge = getNextEdgeLocalIndex();}

};
  

protected:  
  int _refinementLevel;
  std::vector<bool> _markedForRefinement;
  std::map< int, ParentInformation > _interpolationMap;
  
public:
  
  AdaptiveTriangMeshWithTangentSpace ( ) : TriangMeshWithTangentSpace<DataTypeContainer,TriangleType>(), _refinementLevel(0) {}
    
  AdaptiveTriangMeshWithTangentSpace ( const string& fileName, const int tangentVecType ) : 
      TriangMeshWithTangentSpace<DataTypeContainer, TriangleType >( fileName, tangentVecType ), 
     _refinementLevel(0),
     _markedForRefinement( this->getNumTriangs(), false ) {
    this->makeNeighbour();
  }   
  
  AdaptiveTriangMeshWithTangentSpace ( std::vector< Point3DType > vertices, std::vector< Indices3DType > triangles ) :
      TriangMeshWithTangentSpace<DataTypeContainer, TriangleType >( vertices, triangles ), 
      _refinementLevel(0),
      _markedForRefinement( this->getNumTriangs(), false ){
     this->makeNeighbour();
  }

  const int getRefinementLevel() const{ return _refinementLevel; }
  const std::map< int, ParentInformation > & getInterpolationMap ( ) const { return _interpolationMap;}
  
  // mark
  void mark( int element ){ _markedForRefinement[ element ] = true;}
  void markAll(){ for( unsigned i = 0; i < _markedForRefinement.size(); i++ ) _markedForRefinement[i] = true; }
  void unmark( int element ){ _markedForRefinement[ element ] = false;}
  void unmarkAll(){ for( unsigned i = 0; i < _markedForRefinement.size(); i++ ) _markedForRefinement[i] = false; }
  bool isMarkedForRefinement( int element ) const {return _markedForRefinement[element];}
  
  // when adding a new face we have to expand the array _markedForRefinement correspondingly
  int pushBackTriang ( const Indices3DType newTriang ) {
    _markedForRefinement.push_back( false );
    return TriangMeshWithTangentSpace<DataTypeContainer, TriangleType >::pushBackTriang( newTriang );
  }
  
  //returns vector maskVec with 
  //maskVec[elementIndex] = 0 if not marked 
  //                      = 1 if marked
  void getMarkedTrianglesMask( VectorType &maskVec ) const{
      const int numElements = this->getNumTriangs();
      maskVec.resize( numElements ); maskVec.setZero();
      for ( int elementIndex = 0; elementIndex < numElements; ++elementIndex  )
        if ( isMarkedForRefinement( elementIndex ) ) maskVec[elementIndex] = 1.;
  }
  
  void refineMarkedTriangles( ) {
      _refinementLevel++;
      const int numElements = this->getNumTriangs();
      for ( int elementIndex = 0; elementIndex < numElements; ++elementIndex  )
        if ( isMarkedForRefinement( elementIndex ) ){ 
            refine( elementIndex ); 
        }
      cout << "makeOrientationConsistent" << endl; 
      this->makeOrientationConsistent();
      const int numVertices = this->getNumVertices();
      this->_tangentSpaceVec1.resize ( numVertices ); this->_tangentSpaceVec2.resize ( numVertices ); this->_normalSpaceVec.resize ( numVertices );
      
      this->updateAllAdjacentElementsForNodes();
  }
  
  //prolongate function on elements
  void prolongateScalarFct_P0( VectorType & function ) const {
      VectorType oldFunction = function;
      const int numElementsOld = function.size();
      function.resize( this->getNumElements() );
      for( int elementIndex = 0; elementIndex < numElementsOld; ++elementIndex ) function[elementIndex] = oldFunction[elementIndex];
      typename std::map< int, ParentInformation >::const_iterator iter;
      for( int elementIndex = numElementsOld; elementIndex < function.size(); ++elementIndex ){

      }  
  }
  
  // prolongate by linear interpolation
  void prolongateScalarFct_P1Linearly( VectorType & function ) const {
      VectorType oldFunction = function;
      const int numVerticesOld = function.size();
      function.resize( this->getNumVertices() );
      for( int i = 0; i < numVerticesOld; ++i ) function[i] = oldFunction[i];
      typename std::map< int, ParentInformation >::const_iterator iter;
      for( int i = numVerticesOld; i < function.size(); ++i ){
        iter = _interpolationMap.find( i );
        if( iter == _interpolationMap.end() ) throw std::invalid_argument ( pesopt::strprintf ( "unknown vertex in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        function[i] = 0.5 * ( function[iter->second.globalIndices[0]] + function[iter->second.globalIndices[1]] );
      }  
  }
  
  void prolongateVectorFct_P1Linearly( VectorType& function, const int numComp = 3 ) const {
      VectorType oldFunction = function;
      const int numVerticesOld = function.size() / numComp;
      const int numVerticesNew = this->getNumVertices();
      function.resize( numComp * numVerticesNew  );
      for( int i = 0; i < numVerticesOld; ++i ){
          for(int comp=0; comp < numComp; ++comp )
              function[i + comp * numVerticesNew] = oldFunction[i + comp * numVerticesOld];
      }
      typename std::map< int, ParentInformation >::const_iterator iter;
      for( int i = numVerticesOld; i < numVerticesNew; ++i ){
        iter = _interpolationMap.find( i );
        if( iter == _interpolationMap.end() ) throw std::invalid_argument ( pesopt::strprintf ( "unknown vertex in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        for(int comp=0; comp<numComp; ++comp ){ 
            function[i+comp*numVerticesNew] = 0.5 * ( function[iter->second.globalIndices[0] + comp*numVerticesNew] + function[iter->second.globalIndices[1] + comp*numVerticesNew] );
        }
      }  
  }
  
  void prolongateScalarFct_DKTLinearly( VectorType& function ) const {
      VectorType oldFunction = function;
      const int numVerticesOld = function.size() / 3;
      const int numVerticesNew = this->getNumVertices();
      function.resize( 3 * numVerticesNew  );
      for( int i = 0; i < numVerticesOld; ++i ){
          for(int k=0; k<3; ++k )
              function[i + k * numVerticesNew] = oldFunction[i+k*numVerticesOld];
      }
      typename std::map< int, ParentInformation >::const_iterator iter;
      for( int i = numVerticesOld; i < numVerticesNew; ++i ){
        iter = _interpolationMap.find( i );
        if( iter == _interpolationMap.end() ) throw std::invalid_argument ( pesopt::strprintf ( "unknown vertex in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        for(int k=0; k<3; ++k ){ 
            function[i+k*numVerticesNew] = 0.5 * ( function[iter->second.globalIndices[0] + k*numVerticesNew] + function[iter->second.globalIndices[1] + k*numVerticesNew] );
        }
      }  
  }
  
  void prolongateVectorFct_DKTLinearly( VectorType& function, const int numComp = 3 ) const {
      VectorType oldFunction = function;
      const int numVerticesOld = function.size() / ( 3 * numComp );
      const int numVerticesNew = this->getNumVertices();
      function.resize( 3 * numComp * numVerticesNew  );
      for( int i = 0; i < numVerticesOld; ++i ){
          for(int k=0; k < 3 * numComp; ++k )
              function[i + k * numVerticesNew] = oldFunction[i+k*numVerticesOld];
      }
      typename std::map< int, ParentInformation >::const_iterator iter;
      for( int i = numVerticesOld; i < numVerticesNew; ++i ){
        iter = _interpolationMap.find( i );
        if( iter == _interpolationMap.end() ) throw std::invalid_argument ( pesopt::strprintf ( "unknown vertex in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        for(int k=0; k<3 * numComp; ++k ){ 
            function[i+k*numVerticesNew] = 0.5 * ( function[iter->second.globalIndices[0] + k*numVerticesNew] + function[iter->second.globalIndices[1] + k*numVerticesNew] );
        }
      }  
  }
  
  //node is on boundary if both parents are boundary nodes
  void prolongateBoundaryMask( MaskType& mask, const bool onTopologicalBoundary ) const {
      MaskType oldMask = mask;
      const int numVerticesOld = oldMask.size();
      const int numVerticesNew = this->getNumVertices();
      mask.resize( numVerticesNew  );
      for( int i = 0; i < numVerticesOld; ++i ) mask[i] = oldMask[i];
      MaskType topologicalBoundaryMaskNew; this->findTopologicalBoundaryMask( topologicalBoundaryMaskNew );
      typename std::map< int, ParentInformation >::const_iterator iter;
      for( int i = numVerticesOld; i < numVerticesNew; ++i ){
        iter = _interpolationMap.find( i );
        if( iter == _interpolationMap.end() ) throw std::invalid_argument ( pesopt::strprintf ( "unknown vertex in File %s at line %d.", __FILE__, __LINE__ ).c_str() );           
        if( oldMask[iter->second.globalIndices[0]] && oldMask[iter->second.globalIndices[1]] ){
           if( onTopologicalBoundary ){
               if( topologicalBoundaryMaskNew[i] ) mask[i] = true;
           }else mask[i] = true;
        }
        else mask[i] = false;
      }  
  }

  void prolongateBoundaryMaskDKT( MaskType& mask, const bool onTopologicalBoundary ) const {
      MaskType oldMask = mask;
      const int numVerticesOld = oldMask.size() / 3;
      const int numVerticesNew = this->getNumVertices();
      mask.resize( 3 * numVerticesNew  );
      for( int i = 0; i < numVerticesOld; ++i ){
          for(int k=0; k<3; ++k ) mask[i + k * numVerticesNew] = oldMask[i+k*numVerticesOld];
      }
      MaskType topologicalBoundaryMaskNew; this->findTopologicalBoundaryMask( topologicalBoundaryMaskNew );
      typename std::map< int, ParentInformation >::const_iterator iter;
      for( int i = numVerticesOld; i < numVerticesNew; ++i ){
        iter = _interpolationMap.find( i );
        if( iter == _interpolationMap.end() ) throw std::invalid_argument ( pesopt::strprintf ( "unknown vertex in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        for(int k=0; k<3; ++k ){
            if( mask[iter->second.globalIndices[0]+k*numVerticesNew] && mask[iter->second.globalIndices[1]+k*numVerticesNew]){
                if( onTopologicalBoundary ){
                    if( topologicalBoundaryMaskNew[i] ){
                        mask[i+k*numVerticesNew] = true;
                    }
                } else{
                    mask[i+k*numVerticesNew] = true;
                }
            }
            else mask[i+k*numVerticesNew] = false;
        }
      }  
  }
  
  
  
public:
  // internal refinement functions
  GlobalIndex refine( GlobalIndex triangleToRefine ) {
    DartIterator d( *this, triangleToRefine, getLongestEdgeIndex( triangleToRefine ) );
    return refine( d );
  }
  
  GlobalIndex refine( const DartIterator& d) {
    // first cut the neighbour on edge l in two triangles:
    GlobalIndex P_new = addEdgeMidpoint( d ); 
    DartIterator d_prime = d;
    bool hasNeighbourOn_d_prime = d_prime.canFlipTriangle();
    GlobalIndex T_l_new = IndexNotSet;

    if ( hasNeighbourOn_d_prime ){
        d_prime.flipTriangle();
        T_l_new = refineOnlyThis(d_prime, P_new);
    }
    
    // create new triangle T_new and connect with neighbours
    const GlobalIndex T_new = refineOnlyThis(d, P_new);

    // take care of neighbours
    if (hasNeighbourOn_d_prime){
        // set d_prime to edge which lies in T_l_new, opposite to edge d2.edge of T_new
        d_prime.flipNode();
        d_prime.flipEdge();
        // flipTriangle() here will always work, as we have just created this triangle. Thus, no canFlipTriangle() call is needed.
        d_prime.flipTriangle();
        d_prime.flipEdge();
        // set neighbour GlobalIndex in T_l_new
        this->setNeighbour ( d_prime.getGlobalTriangleIndex(), d_prime.getLocalEdgeIndex(), T_new );
        // move d1 to T_new
        DartIterator d1 = d; d1.flipNode(); d1.flipEdge();
        d1.flipTriangle();
        d1.flipEdge();
        // set neighbour in T_new
        this->setNeighbour ( d1.getGlobalTriangleIndex(), d1.getLocalEdgeIndex(), T_l_new );
    }
    return T_new;
  }
    
  GlobalIndex refineOnlyThis( const DartIterator& d, GlobalIndex midpoint) {
    // in comments for this function, we write T for *this.
    LocalIndex l = getLongestEdgeIndex( d.getGlobalTriangleIndex() );
    if (l != d.getLocalEdgeIndex()) {
        DartIterator d_l(*this, d.getGlobalTriangleIndex(), l);
        LocalIndex commonNode = d.getCommonNodeLocalIndex(d_l);
        DartIterator d_prime( *this, d.getGlobalTriangleIndex(), l, commonNode);
        // bisect T along the longest edge
        refine(d_prime);
        // because d_prime.node also lies on d.edge, T is now modifies in such a way that d.edge has remained unchanged. Proceed by recursion until d.edge is the longest edge in T.
        return refineOnlyThis(d, midpoint);
    }
    // create DartIterators to all edges, starting from l
    DartIterator d1 = d;  d1.flipNode(); d1.flipEdge();
    DartIterator d2 = d1; d2.flipNode(); d2.flipEdge();
    // create new triangle T_new and connect with neighbours
    GlobalIndex T_newIndex = this->pushBackTriang ( Indices3DType( midpoint, d1.getGlobalNodeIndex(), d2.getGlobalNodeIndex() ) );
    this->_neighbour.push_back( Indices3DType ( d1.getNextTriangleIndex(), d.getGlobalTriangleIndex(), d.getNextTriangleIndex() ) );
    // set neighbour GlobalIndex in triangle lying opposite to d1.edge to T_new
    DartIterator d1_prime = d1;
    if (d1_prime.canFlipTriangle()){
        d1_prime.flipTriangle();
        this->setNeighbour ( d1_prime.getGlobalTriangleIndex(), d1_prime.getLocalEdgeIndex(), T_newIndex );
    }
    // connect ourselves (i. e., connect T) to T_new:
    this->setNeighbour( d1.getGlobalTriangleIndex(), d1.getLocalEdgeIndex(), T_newIndex );
    // do not have to refine this triangle again
    unmark( d1.getGlobalTriangleIndex() );
    // old triangle now has a new vertex, ie the midpoint
    this->getTriang( d1.getGlobalTriangleIndex() ).setGlobalNodeIdx( d1.getLocalNodeIndex(), midpoint);
    return T_newIndex;
  }
  
  // get longest edge index (starting search possibly with a preferred edge)
  LocalIndex getLongestEdgeIndex( GlobalIndex triangle ) const {
      // get global node indices of triangle
      const Indices3DType &triangIdx ( this->getTriang( triangle ).getGlobalNodeIdx() );
      // assume that edge 0 is longest edge
      LocalIndex localIdx = 0;
      Point3DType edge = this->getVertex( triangIdx[1] ); edge -= this->getVertex( triangIdx[2] ); RealType candidate0 = edge.squaredNorm(); 
      RealType maxLength = candidate0;
      // now check if edge 1 is longer
      edge = this->getVertex( triangIdx[2] ); edge -= this->getVertex( triangIdx[0] ); RealType candidate1 = edge.squaredNorm();
      if( candidate1 > maxLength ){ localIdx = 1; maxLength = candidate1; }
      // now check if edge 2 is longer
      edge = this->getVertex( triangIdx[0] ); edge -= this->getVertex( triangIdx[1] ); RealType candidate2 = edge.squaredNorm();
      if( candidate2 > maxLength ){ localIdx = 2; maxLength = candidate2; };
      return localIdx;  
   }
  
  // add edge midpoint on d.edge, update _interpolationMap and return global index of node
  GlobalIndex addEdgeMidpoint( const DartIterator& d ) {    
    DartIterator d1 = d; d1.flipNode();
    
    ParentInformation parents;
    parents.globalIndices[0] = d.getGlobalNodeIndex(); parents.globalIndices[1] = d1.getGlobalNodeIndex();
    
    // get edge midpoint between nodes n1 and n2
    Point3DType newVertex; newVertex.setZero();
    newVertex += this->getVertex ( d.getGlobalNodeIndex() );
    newVertex += this->getVertex ( d1.getGlobalNodeIndex() );
    newVertex /= 2.;
    
    // add new vertex
    int newNodeIdx = this->pushBackVertex( newVertex ); 
    _interpolationMap.insert( std::pair< int, ParentInformation >( newNodeIdx, parents ) );
    return newNodeIdx;
  }

};


#endif
