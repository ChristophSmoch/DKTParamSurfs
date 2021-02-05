#ifndef __TRIANGLEFEADAPTIVEREFINEMENT_H
#define __TRIANGLEFEADAPTIVEREFINEMENT_H


#include <feAdaptiveRefinement.h>

template<typename TriangleMeshType >
class FEAdaptiveMesh<TriangleMeshType, TRIANGLE> 
: public FEAdaptiveMeshBase<TriangleMeshType>{
    
  typedef typename TriangleMeshType::ElementType        ElementType;
  typedef typename ElementType::ElementNodeIndicesType  ElementNodeIndicesType;
  typedef typename TriangleMeshType::RealType           RealType;
  typedef typename TriangleMeshType::RealVecChart       RealVecChart;
  typedef typename TriangleMeshType::IntVecChart        IntVecChart;
  typedef typename TriangleMeshType::PointType          PointType;
  typedef typename TriangleMeshType::VectorType         VectorType;
  typedef typename TriangleMeshType::MaskType           MaskType;
    
  
  typedef struct { int globalIndices[2];} ParentInformation;
  
  typedef short LocalIndex;
  typedef int GlobalIndex;
  static const int IndexNotSet = -1;
  
protected:  
// half-edge iterator (only for internal refinement!)
class DartIterator{
    
  protected:
    const TriangleMeshType &_mesh;
    GlobalIndex _triangle;
    LocalIndex  _node;
    LocalIndex  _edge;

  public:
    DartIterator( const TriangleMeshType &Mesh, GlobalIndex triangleIndex, LocalIndex localEdgeIndex, LocalIndex localNodeIndex ) : 
      _mesh( Mesh ), _triangle (triangleIndex), _node(localNodeIndex), _edge(localEdgeIndex){}    
    DartIterator( const TriangleMeshType &Mesh, GlobalIndex triangleIndex, LocalIndex localEdgeIndex ) : 
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
      if(_triangle != d.getGlobalTriangleIndex()) 
          throw std::invalid_argument( pesopt::strprintf("T=%d, but tringle = %d. In File %s at line %d.", d.getGlobalTriangleIndex(), _triangle, __FILE__, __LINE__).c_str() );
      if ( _edge == d.getLocalEdgeIndex() )
          throw std::invalid_argument( pesopt::strprintf("points to same edge=. In File %s at line %d.", __FILE__, __LINE__).c_str() );
      return 3 - ( getLocalEdgeIndex() + d.getLocalEdgeIndex() );
    }

    // returns global index of common node (also if "d" refers to a different triangle than _triangle )
    template<typename DartIteratorType>
    GlobalIndex getCommonNodeGlobalIndex( const DartIteratorType & d ) const {
      if ( getGlobalNodeIndex() == d.getGlobalNodeIndex() || getGlobalNodeIndex() == d.getNextNodeGlobalIndex() ) return getGlobalNodeIndex();
      if ( getNextNodeGlobalIndex() == d.getGlobalNodeIndex() || getNextNodeGlobalIndex() == d.getNextNodeGlobalIndex() ) return getNextNodeGlobalIndex();
      throw std::invalid_argument( pesopt::strprintf("No common node ofDart d and mesh. In File %s at line %d.", __FILE__, __LINE__).c_str() );
      return -1;
    }

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
  mutable std::vector<bool> _markedForRefinement;
  mutable std::map< int, ParentInformation > _interpolationMap;
  
public:
    
  FEAdaptiveMesh( const TriangleMeshType &mesh ) :
    FEAdaptiveMeshBase<TriangleMeshType> ( mesh ),
    _markedForRefinement( mesh.getNumTriangs(), false ) { }
    
    
protected:
   const std::map< int, ParentInformation > & getInterpolationMap ( ) const { return _interpolationMap;}
  
  // mark
  void mark( int element ) const { _markedForRefinement[ element ] = true;}
  void markAll() const { for( unsigned i = 0; i < _markedForRefinement.size(); i++ ) _markedForRefinement[i] = true; }
  void unmark( int element ) const { _markedForRefinement[ element ] = false;}
  void unmarkAll() const { for( unsigned i = 0; i < _markedForRefinement.size(); i++ ) _markedForRefinement[i] = false; }
  const bool isMarkedForRefinement( int element ) const {return _markedForRefinement[element];}
     
public:
    
    
  void refineAll( ) const override{
      const MaskType markedElements ( this->_meshVector[this->_refinementLevel].getNumElements(), true );
      this->refineMarkedElements( markedElements );
  }
  
  void refineMarkedElements( const MaskType &markedElements ) const override {
      this->_refinementLevel++;
      TriangleMeshType refinedMesh ( this->_meshVector[this->_refinementLevel-1] );
      for(int elementIndex=0; elementIndex<markedElements.size(); ++elementIndex )
           mark( elementIndex );
      for ( int elementIndex = 0; elementIndex < markedElements.size(); ++elementIndex  )
        if ( isMarkedForRefinement( elementIndex ) ) refine( refinedMesh, elementIndex ); 

      refinedMesh.makeOrientationConsistent();
      refinedMesh.generateBoundaryMask();
      refinedMesh.updateAllAdjacentElementsForNodes();
      //       refinedMesh.makeEdges(); 
      this->_meshVector.push_back( refinedMesh );
  }
  
//   void refineByFactor( const RealType factor ) const override {
//   }
    
    
  template <typename ConfiguratorType>
  void prolongateScalarFunctionAtNodes ( VectorType &vec ) const {
     VectorType vecOld(vec);
     this->template prolongateScalarFunctionAtNodes<ConfiguratorType> ( vecOld, vec );
 }
 
 
  template <typename ConfiguratorType>
  void prolongateScalarFunctionAtNodes ( const VectorType &vecOld, VectorType &vecNew  ) const {
      const int numVerticesOld = vecOld.size();
      vecNew.resize( this->_meshVector[this->_refinementLevel].getNumVertices() );
      for( int i = 0; i < numVerticesOld; ++i ) vecNew[i] = vecOld[i];
      typename std::map< int, ParentInformation >::const_iterator iter;
      for( int i = numVerticesOld; i < vecNew.size(); ++i ){
        iter = _interpolationMap.find( i );
        if( iter == _interpolationMap.end() ) 
            throw std::invalid_argument ( pesopt::strprintf ( "unknown vertex in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        vecNew[i] = 0.5 * ( vecNew[iter->second.globalIndices[0]] + vecNew[iter->second.globalIndices[1]] );
      }  
      
 }
 
 
 
protected:
  // internal refinement functions
  GlobalIndex refine( TriangleMeshType &mesh, GlobalIndex triangleToRefine ) const {
    DartIterator d( mesh, triangleToRefine, getLongestEdgeIndex( mesh, triangleToRefine ) );
    return refine( mesh, d );
  }
  
  GlobalIndex refine( TriangleMeshType &mesh, const DartIterator& d) const{
    // first cut the neighbour on edge l in two triangles:
    GlobalIndex P_new = addEdgeMidpoint( mesh, d ); 
    DartIterator d_prime = d;
    bool hasNeighbourOn_d_prime = d_prime.canFlipTriangle();
    GlobalIndex T_l_new = IndexNotSet;

    if ( hasNeighbourOn_d_prime ){
        d_prime.flipTriangle();
        T_l_new = refineOnlyThis(mesh, d_prime, P_new);
    }
    
    // create new triangle T_new and connect with neighbours
    const GlobalIndex T_new = refineOnlyThis(mesh, d, P_new);

    // take care of neighbours
    if (hasNeighbourOn_d_prime){
        // set d_prime to edge which lies in T_l_new, opposite to edge d2.edge of T_new
        d_prime.flipNode();
        d_prime.flipEdge();
        // flipTriangle() will always work here, since we have just created the triangle. Thus, canFlipTriangle() is not needed.
        d_prime.flipTriangle();
        d_prime.flipEdge();
        // set neighbour GlobalIndex in T_l_new
        mesh.setNeighbour ( d_prime.getGlobalTriangleIndex(), d_prime.getLocalEdgeIndex(), T_new );
        // move d1 to T_new
        DartIterator d1 = d; d1.flipNode(); d1.flipEdge();
        d1.flipTriangle();
        d1.flipEdge();
        // set neighbour in T_new
        mesh.setNeighbour ( d1.getGlobalTriangleIndex(), d1.getLocalEdgeIndex(), T_l_new );
    }
    return T_new;
  }
    
  GlobalIndex refineOnlyThis( TriangleMeshType &mesh, const DartIterator& d, GlobalIndex midpoint) const{
    LocalIndex l = getLongestEdgeIndex( mesh, d.getGlobalTriangleIndex() );
    if (l != d.getLocalEdgeIndex()) {
        DartIterator d_l(mesh, d.getGlobalTriangleIndex(), l);
        LocalIndex commonNode = d.getCommonNodeLocalIndex(d_l);
        DartIterator d_prime( mesh, d.getGlobalTriangleIndex(), l, commonNode);
        // bisect T along the longest edge
        refine(mesh, d_prime);
        // because d_prime.node also lies on d.edge, T is now modifies in such a way that d.edge has remained unchanged. Proceed by recursion until d.edge is the longest edge in T.
        return refineOnlyThis(mesh, d, midpoint);
    }
    // create DartIterators to all edges, starting from l
    DartIterator d1 = d;  d1.flipNode(); d1.flipEdge();
    DartIterator d2 = d1; d2.flipNode(); d2.flipEdge();
    // create new triangle T_new and connect with neighbours
    GlobalIndex T_newIndex = mesh.pushBackTriang ( ElementNodeIndicesType( midpoint, d1.getGlobalNodeIndex(), d2.getGlobalNodeIndex() ) );
    _markedForRefinement.push_back( false );
    mesh.pushBackNeighbour( ElementNodeIndicesType ( d1.getNextTriangleIndex(), d.getGlobalTriangleIndex(), d.getNextTriangleIndex() ) );
    // set neighbour GlobalIndex in triangle lying opposite to d1.edge to T_new
    DartIterator d1_prime = d1;
    if (d1_prime.canFlipTriangle()){
        d1_prime.flipTriangle();
        mesh.setNeighbour ( d1_prime.getGlobalTriangleIndex(), d1_prime.getLocalEdgeIndex(), T_newIndex );
    }
    // connect ourselves (i. e., connect T) to T_new:
    mesh.setNeighbour( d1.getGlobalTriangleIndex(), d1.getLocalEdgeIndex(), T_newIndex );
    // do not have to refine this triangle again
    unmark( d1.getGlobalTriangleIndex() );
    // old triangle now has a new vertex, ie the midpoint
//     mesh.getElement( d1.getGlobalTriangleIndex() ).setGlobalNodeIdx( d1.getLocalNodeIndex(), midpoint);
    mesh.setTriangNodeIdx( d1.getGlobalTriangleIndex(), d1.getLocalNodeIndex(), midpoint);
   
    return T_newIndex;
  }
  
  // get longest edge index (starting search possibly with a preferred edge)
  LocalIndex getLongestEdgeIndex( TriangleMeshType &mesh, GlobalIndex triangle ) const {
      // get global node indices of triangle
      const ElementNodeIndicesType &triangIdx ( mesh.getElement( triangle ).getGlobalNodeIdx() );
      // assume that edge 0 is longest edge
      LocalIndex localIdx = 0;
      PointType edge = mesh.getVertex( triangIdx[1] ); edge -= mesh.getVertex( triangIdx[2] ); RealType candidate0 = edge.squaredNorm(); 
      RealType maxLength = candidate0;
      // now check if edge 1 is longer
      edge = mesh.getVertex( triangIdx[2] ); edge -= mesh.getVertex( triangIdx[0] ); RealType candidate1 = edge.squaredNorm();
      if( candidate1 > maxLength ){ localIdx = 1; maxLength = candidate1; }
      // now check if edge 2 is longer
      edge = mesh.getVertex( triangIdx[0] ); edge -= mesh.getVertex( triangIdx[1] ); RealType candidate2 = edge.squaredNorm();
      if( candidate2 > maxLength ){ localIdx = 2; maxLength = candidate2; };
      return localIdx;  
   }
  
  // add edge midpoint on d.edge, update _interpolationMap and return global index of node
  GlobalIndex addEdgeMidpoint( TriangleMeshType &mesh, const DartIterator& d ) const{    
    DartIterator d1 = d; d1.flipNode();
    
    ParentInformation parents;
    parents.globalIndices[0] = d.getGlobalNodeIndex(); parents.globalIndices[1] = d1.getGlobalNodeIndex();
    
    // get edge midpoint between nodes n1 and n2
    PointType newVertex; newVertex.setZero();
    newVertex += mesh.getVertex ( d.getGlobalNodeIndex() );
    newVertex += mesh.getVertex ( d1.getGlobalNodeIndex() );
    newVertex /= 2.;
    
    // add new vertex
    int newNodeIdx = mesh.pushBackVertex( newVertex ); 
    _interpolationMap.insert( std::pair< int, ParentInformation >( newNodeIdx, parents ) );
    return newNodeIdx;
  }
    
};


#endif
