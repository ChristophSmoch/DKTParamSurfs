#ifndef __DKTFEBASEFUNCTIONSETS_H
#define __DKTFEBASEFUNCTIONSETS_H

#include <pesopt_IO.h>
#include <pesopt_Eigen.h>


//! Inteface
template <typename DataTypeContainer, class QuadType >
class RefTriangleBaseFunctionSetInterface  {
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
    
public:
  RefTriangleBaseFunctionSetInterface( ) {}

  int numQuadPoints( ) const { return QuadType::numQuadPoints;}
  inline RealType getWeight ( int QuadPoint ) const { return _quadRule.getWeight ( QuadPoint );}
  inline const DomVecType& getRefCoord ( int QuadPoint ) const { return _quadRule.getRefCoord ( QuadPoint );}

protected:
  QuadType _quadRule;
};






template <typename DataTypeContainer, typename QuadType, typename TriangleType>
class RefTriangMeshBaseFunctionSetP0 : public RefTriangleBaseFunctionSetInterface< DataTypeContainer, QuadType >  {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;

  const TriangleType *_triangle;
  
public:
  RefTriangMeshBaseFunctionSetP0( ) : _triangle ( NULL ) { }

  void setTriangle ( const TriangleType &T ) { _triangle = &T;}
  
  enum { numBaseFuncs = 1 };
  
  RealType evaluateOnRefTriang ( int BaseFuncNum, const DomVecType &RefCoord ) const { return 1.;}
  inline const RealType evaluateOnRefTriang ( int BaseFuncNum, int QuadPoint ) const { return 1.;}
  
};






template <typename DataTypeContainer, typename QuadType, typename TriangleType, typename Imp>
class RefTriangMeshBaseFunctionSetPkInterface : public RefTriangleBaseFunctionSetInterface< DataTypeContainer, QuadType >  {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;

  const TriangleType *_triangle;

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
public:
  RefTriangMeshBaseFunctionSetPkInterface( ) : _triangle ( NULL ) { }

  void setTriangle ( const TriangleType &T ) { _triangle = &T;}
  
  //****************************************************
  // on unit triangle
  //****************************************************
   
  //! interface function, has to be provided in derived classes.
  RealType evaluateOnUnitTriang ( int BaseFuncNum, const DomVecType &RefCoord ) const { 
     throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
     return this->asImp().evaluateOnUnitTriang ( BaseFuncNum, RefCoord );
  }

  inline const RealType evaluateOnUnitTriang ( int BaseFuncNum, int QuadPoint ) const { 
      return this->asImp().evaluateOnUnitTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );
  }
  
  //! interface function, has to be provided in derived classes.
  void evaluateGradientOnUnitTriang ( int BaseFuncNum, const DomVecType &RefCoord, DomVecType &Gradient ) const {
     throw std::invalid_argument( pesopt::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
     this->asImp().evaluateGradientOnUnitTriang ( BaseFuncNum, RefCoord, Gradient );
  }
  
  inline DomVecType evaluateGradientOnUnitTriang ( int BaseFuncNum, int QuadPoint ) const {
     DomVecType Gradient;
     this->asImp().evaluateGradientOnUnitTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), Gradient );
     return Gradient;
  }
  //*****************************************************
  
  
  //****************************************************
  // on reference triangle
  //****************************************************
  RealType evaluateOnRefTriang ( int BaseFuncNum, const DomVecType &RefCoord ) const { return this->asImp().evaluateOnUnitTriang(BaseFuncNum,RefCoord);}
  inline const RealType evaluateOnRefTriang ( int BaseFuncNum, int QuadPoint ) const { return this->evaluateOnRefTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );}
  
  //TODO Gradient is now 3d Vector, before 2d vector!!!
  // Gradient = Dxref g^{-1} Grad_UnitTriang
  void evaluateGradientOnRefTriang ( int BaseFuncNum, const DomVecType &RefCoord, DomVecType &Gradient ) const {
    DomVecType GradUnitTriang; this->asImp().evaluateGradientOnUnitTriang( BaseFuncNum, RefCoord, GradUnitTriang );
    DomVecType ginvDx; ginvDx = _triangle->getInverseMetricRefTriang( ) * GradUnitTriang;  
    const DomVecType dir0 = _triangle->edgeRefTriang(0,1);
    const DomVecType dir1 = _triangle->edgeRefTriang(0,2); 
    //TODO for g^-1 = Dx^TDx: normalDirection = dir0 x dir1 / |..|
    Gradient = ginvDx[0] * dir0 + ginvDx[1] * dir1;
  }

  inline DomVecType evaluateGradientOnRefTriang ( int BaseFuncNum, int QuadPoint ) const {
    DomVecType Gradient;
    this->asImp().evaluateGradientOnRefTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), Gradient );
    return Gradient;
  }
  //*****************************************************
  
};







template <typename DataTypeContainer, typename QuadType, typename TriangleType>
class RefTriangMeshBaseFunctionSetP1 : 
public RefTriangMeshBaseFunctionSetPkInterface< DataTypeContainer, QuadType, TriangleType, RefTriangMeshBaseFunctionSetP1<DataTypeContainer, QuadType, TriangleType> >  {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
    
  static RealType _b1   ( const DomVecType &c ) { return 1. - c[0] - c[1]; }
  static RealType _b2   ( const DomVecType &c ) { return c[0]; }
  static RealType _b3   ( const DomVecType &c ) { return c[1]; }

  static RealType _d1_b1   ( const DomVecType &/*c*/ ) { return - 1.; }
  static RealType _d1_b2   ( const DomVecType &/*c*/ ) { return   1.; }
  static RealType _d1_b3   ( const DomVecType &/*c*/ ) { return   0.; }

  static RealType _d2_b1   ( const DomVecType &/*c*/ ) { return - 1.; }
  static RealType _d2_b2   ( const DomVecType &/*c*/ ) { return   0.; }
  static RealType _d2_b3   ( const DomVecType &/*c*/ ) { return   1.; }
  
  
  typedef RealType ( *BASIS_FUNC_TYPE ) ( const DomVecType &RefCoord );
  BASIS_FUNC_TYPE _basis[3];
  BASIS_FUNC_TYPE _dbasis[2][3];

public:
  RefTriangMeshBaseFunctionSetP1(  ) :
  RefTriangMeshBaseFunctionSetPkInterface< DataTypeContainer, QuadType, TriangleType, RefTriangMeshBaseFunctionSetP1<DataTypeContainer, QuadType, TriangleType> > () {
      _basis[0] = _b1; 
      _basis[1] = _b2; 
      _basis[2] = _b3;
      
      _dbasis[0][0] = _d1_b1;
      _dbasis[0][1] = _d1_b2;
      _dbasis[0][2] = _d1_b3;
    
      _dbasis[1][0] = _d2_b1;
      _dbasis[1][1] = _d2_b2;
      _dbasis[1][2] = _d2_b3;
  }
  
  enum { numBaseFuncs = 3 };
  
  RealType evaluateOnUnitTriang ( int BaseFuncNum, const DomVecType &RefCoord ) const { return _basis[BaseFuncNum] ( RefCoord );}
  
  void evaluateGradientOnUnitTriang ( int BaseFuncNum, const DomVecType &RefCoord, DomVecType &Gradient ) const {
    Gradient[0] = _dbasis[0][BaseFuncNum](RefCoord); Gradient[1] = _dbasis[1][BaseFuncNum](RefCoord);
  }
};







//! This class realizes the localToGlobal mapping when working with quadratic base functions, i.e. we have 6 (local) degrees of freedom (dofs) per element.
template< typename MeshType >
class P2Mapper {
  typedef Eigen::Vector2i KeyType;
  typedef std::pair<KeyType, int> PairType;
  typedef std::pair<int, KeyType> InvPairType;
      
  const MeshType& _mesh;
  int _globalDofs;
  
  std::map< KeyType, int > _map;        // maps triangle index and local node index to corresponding global node index
  std::map< KeyType, int > _edgeMap;    // maps global node indices of two nodes representing one edge to global index of dof on that edge 
  std::map< int, KeyType > _invEdgeMap; // maps one global node index to (i) itsef, if vertex; (ii) two adjacent global node indices, if on edge (neglecting order here!)
  
public:
  P2Mapper( const MeshType& mesh ) : _mesh( mesh ) { initMapper();}
  
  int localToGlobal( const int& T, const int& loc ) const { return typename std::map< KeyType, int >::const_iterator( _map.find( KeyType( T, loc ) ) )->second;}
  int getGlobalDofs () const{ return _globalDofs;}
  int getEdgeDOF( const int& v1, const int& v2 ) const {return typename std::map< KeyType, int >::const_iterator( _edgeMap.find( KeyType( v1, v2 ) ) )->second;}
  
  KeyType getGlobalPairing( const int& globIdx ) const {
    KeyType res( globIdx, globIdx );
    typename std::map< int, KeyType >::const_iterator it = _invEdgeMap.find( globIdx );
    if( it != _invEdgeMap.end() ) return it->second;
    return res;
  }
  
protected:
    
  void initMapper() {
    // first fill edge map to define global node indices for dofs on edges
    fillEdgeMaps();
    
    // add 6 local dofs for each element to the map
    for( int elementNumber = 0; elementNumber < _mesh.getNumTriangs(); ++elementNumber ){
      int node0 = _mesh.getElementNodeIdx(elementNumber,0);
      int node1 = _mesh.getElementNodeIdx(elementNumber,1);
      int node2 = _mesh.getElementNodeIdx(elementNumber,2);
      // ordering is (0,1,2,3,4,5) when going anti-clockwise
      _map.insert( PairType( KeyType(elementNumber,0), node0 ) );
      _map.insert( PairType( KeyType(elementNumber,1), _edgeMap[ KeyType( std::min(node0,node1), std::max(node0,node1) ) ] ) );
      _map.insert( PairType( KeyType(elementNumber,2), node1 ) );
      _map.insert( PairType( KeyType(elementNumber,3), _edgeMap[ KeyType( std::min(node1,node2), std::max(node1,node2) ) ] ) );
      _map.insert( PairType( KeyType(elementNumber,4), node2 ) );
      _map.insert( PairType( KeyType(elementNumber,5), _edgeMap[ KeyType( std::min(node2,node0),std::max(node2,node0) ) ] ) );
      
    }
  }
  
  void fillEdgeMaps () {     
    _globalDofs = _mesh.getNumVertices();
    // construct edge map and increase number of dofs
    for( int elementNumber = 0; elementNumber < _mesh.getNumTriangs(); ++elementNumber ){
      Eigen::Vector3i nodes;
      nodes[0] = _mesh.getElementNodeIdx(elementNumber,0);
      nodes[1] = _mesh.getElementNodeIdx(elementNumber,1);
      nodes[2] = _mesh.getElementNodeIdx(elementNumber,2);
      for( int i = 0; i < 3; i++ ){
        typename std::map< KeyType, int >::const_iterator it = _edgeMap.find( KeyType( std::min(nodes[i],nodes[(i+1)%3]), std::max(nodes[i],nodes[(i+1)%3]) ) );
        if( it == _edgeMap.end() ){
          _edgeMap.insert( PairType( KeyType( std::min(nodes[i],nodes[(i+1)%3]), std::max(nodes[i],nodes[(i+1)%3]) ), _globalDofs ) );
          _invEdgeMap.insert( InvPairType( _globalDofs++, KeyType( nodes[i], nodes[(i+1)%3] ) ) ); // order Vec does not matter as it is not the key!
        }
      }
    }
  }
};


//! BasefunctionSet for quadratic finite element functions on triangular meshes
template <typename DataTypeContainer, typename QuadType, typename TriangleType>
class RefTriangMeshBaseFunctionSetP2 :
public RefTriangMeshBaseFunctionSetPkInterface< DataTypeContainer, QuadType, TriangleType, RefTriangMeshBaseFunctionSetP2<DataTypeContainer, QuadType, TriangleType> >  {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
    
public:
  // ordering anti-clockwise
  static RealType _b1   ( const DomVecType &c ) { return 2. * ( 1. - c[0] - c[1] ) * ( 1./2. - c[0] - c[1] ); }
  static RealType _b2   ( const DomVecType &c ) { return 4. * c[0] * ( 1. - c[0] - c[1] ); }
  static RealType _b3   ( const DomVecType &c ) { return c[0] * ( 2. * c[0] - 1. ); }
  static RealType _b4   ( const DomVecType &c ) { return 4. * c[0] * c[1]; }
  static RealType _b5   ( const DomVecType &c ) { return c[1] * ( 2. * c[1] - 1. ); }
  static RealType _b6   ( const DomVecType &c ) { return 4. * c[1] * ( 1. - c[0] - c[1] ); }

  static RealType _d1_b1   ( const DomVecType &c ) { return - 4. * ( 3./4. - c[0] - c[1] ); }
  static RealType _d1_b2   ( const DomVecType &c ) { return 4. - 8.* c[0] - 4. * c[1]; }
  static RealType _d1_b3   ( const DomVecType &c ) { return 4. * c[0] - 1.; }
  static RealType _d1_b4   ( const DomVecType &c ) { return 4. * c[1]; }
  static RealType _d1_b5   ( const DomVecType &  ) { return 0.; }
  static RealType _d1_b6   ( const DomVecType &c ) { return -4. * c[1]; }

  static RealType _d2_b1   ( const DomVecType &c ) { return - 4. * ( 3./4. - c[0] - c[1] ); }
  static RealType _d2_b2   ( const DomVecType &c ) { return - 4. * c[0]; }
  static RealType _d2_b3   ( const DomVecType &  ) { return 0.; }
  static RealType _d2_b4   ( const DomVecType &c ) { return 4. * c[0]; }
  static RealType _d2_b5   ( const DomVecType &c ) { return 4. * c[1] - 1.; }
  static RealType _d2_b6   ( const DomVecType &c ) { return 4. - 4. * c[0] - 8. * c[1]; }
  
protected:
  typedef RealType ( *BASIS_FUNC_TYPE ) ( const DomVecType &RefCoord );
  BASIS_FUNC_TYPE _dbasis[2][6];
  BASIS_FUNC_TYPE _basis[6];

public:
  RefTriangMeshBaseFunctionSetP2(  ) : 
  RefTriangMeshBaseFunctionSetPkInterface< DataTypeContainer, QuadType, TriangleType, RefTriangMeshBaseFunctionSetP2<DataTypeContainer, QuadType, TriangleType> > () 
  {
  
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;

    _dbasis[0][0] = _d1_b1;
    _dbasis[0][1] = _d1_b2;
    _dbasis[0][2] = _d1_b3;
    _dbasis[0][3] = _d1_b4;
    _dbasis[0][4] = _d1_b5;
    _dbasis[0][5] = _d1_b6;

    _dbasis[1][0] = _d2_b1;
    _dbasis[1][1] = _d2_b2;
    _dbasis[1][2] = _d2_b3;
    _dbasis[1][3] = _d2_b4;
    _dbasis[1][4] = _d2_b5;
    _dbasis[1][5] = _d2_b6;
  }

  enum { numBaseFuncs = 6 };

  RealType evaluateOnUnitTriang ( int BaseFuncNum, const DomVecType &RefCoord ) const { return _basis[BaseFuncNum] ( RefCoord );}
  
  void evaluateGradientOnUnitTriang ( int BaseFuncNum, const DomVecType &RefCoord, DomVecType &Gradient ) const {
    Gradient[0] = _dbasis[0][BaseFuncNum](RefCoord); Gradient[1] = _dbasis[1][BaseFuncNum](RefCoord);
  }
  
};










//! This class realizes the localToGlobal mapping when working with reduced cubic base functions,
//! i.e. the dof in the center of the triangle has been removed
//! Hence, we have 9 (local) degrees of freedom (dofs) per element, which are all on the boundary of that element
//! and oredered in anti-clockwise direction
template< typename MeshType >
class P3redMapper {
  typedef Eigen::Vector2i KeyType;
  typedef std::pair<KeyType, int> PairType;
  typedef std::pair<int, KeyType> InvPairType;
      
  const MeshType& _mesh;
  int _globalDofs;
  
  std::map< KeyType, int > _map;        // maps triangle index and local node index to corresponding global node index
  std::map< KeyType, int > _edgeMap;    // maps global node indices of two nodes representing one edge to global index of dof on that edge 
  std::map< int, KeyType > _invEdgeMap; // maps one global node index to (i) itsef, if vertex; (ii) two adjacent global node indices, if on edge
  
public:
  P3redMapper( const MeshType& mesh ) : _mesh( mesh ) { initMapper();}
  int localToGlobal( const int& T, const int& loc ) const { return typename std::map< KeyType, int >::const_iterator( _map.find( KeyType( T, loc ) ) )->second;}
  int getGlobalDofs () const{ return _globalDofs;}
  int getEdgeDOF( const int& v1, const int& v2 ) const { return typename std::map< KeyType, int >::const_iterator( _edgeMap.find( KeyType( v1, v2 ) ) )->second;}
  
  KeyType getGlobalPairing( const int& globIdx ) const {
    KeyType res( globIdx, globIdx );
    typename std::map< int, KeyType >::const_iterator it = _invEdgeMap.find( globIdx );
    if( it != _invEdgeMap.end() ) return it->second;
    return res;
  }
  
protected:
  void initMapper() {
    // first fill edge map to define global node indices for dofs on edges
    fillEdgeMaps();
    
    // add 9 local dofs for each element to the map
    for( int elementNumber = 0; elementNumber < _mesh.getNumTriangs(); ++elementNumber ){
      Eigen::Vector3i nodes;
      nodes[0] = _mesh.getElementNodeIdx(elementNumber,0);
      nodes[1] = _mesh.getElementNodeIdx(elementNumber,1);
      nodes[2] = _mesh.getElementNodeIdx(elementNumber,2);

      // ordering is (0,1,2,3,4,5,6,7,8) when going anti-clockwise
      _map.insert( PairType( KeyType(elementNumber,0), nodes[0] ) );
      _map.insert( PairType( KeyType(elementNumber,1), _edgeMap[ KeyType( nodes[0], nodes[1] ) ] ) );
      _map.insert( PairType( KeyType(elementNumber,2), _edgeMap[ KeyType( nodes[1], nodes[0] ) ] ) );
      _map.insert( PairType( KeyType(elementNumber,3), nodes[1] ) );
      _map.insert( PairType( KeyType(elementNumber,4), _edgeMap[ KeyType( nodes[1], nodes[2] ) ] ) );
      _map.insert( PairType( KeyType(elementNumber,5), _edgeMap[ KeyType( nodes[2], nodes[1] ) ] ) );
      _map.insert( PairType( KeyType(elementNumber,6), nodes[2] ) );
      _map.insert( PairType( KeyType(elementNumber,7), _edgeMap[ KeyType( nodes[2], nodes[0] ) ] ) );
      _map.insert( PairType( KeyType(elementNumber,8), _edgeMap[ KeyType( nodes[0], nodes[2] ) ] ) );
      
    }
  }
  
  void fillEdgeMaps () {     
    _globalDofs = _mesh.getNumVertices();
    // construct edge map and increase number of dofs
    for( int elementNumber = 0; elementNumber < _mesh.getNumTriangs(); ++elementNumber ){
      Eigen::Vector3i nodes;
      nodes[0] = _mesh.getElementNodeIdx(elementNumber,0);
      nodes[1] = _mesh.getElementNodeIdx(elementNumber,1);
      nodes[2] = _mesh.getElementNodeIdx(elementNumber,2);
      for( int i = 0; i < 3; i++ ){
        typename std::map< KeyType, int >::const_iterator it = _edgeMap.find( KeyType( nodes[i], nodes[(i+1)%3] ) );
        if( it == _edgeMap.end() ){
          // first dof on edge (i,j) is closer to i, hence node i is first entry in key (j=i+1)
          _edgeMap.insert( PairType( KeyType( nodes[i], nodes[(i+1)%3] ), _globalDofs ) );
          _invEdgeMap.insert( InvPairType( _globalDofs++, KeyType( nodes[i], nodes[(i+1)%3] ) ) );
          // second dof on edge (i,j) is closer to j, hence node j is first entry in key (j=i+1)
          _edgeMap.insert( PairType( KeyType( nodes[(i+1)%3], nodes[i] ), _globalDofs ) );
          _invEdgeMap.insert( InvPairType( _globalDofs++, KeyType( nodes[(i+1)%3], nodes[i] ) ) );
        }
      }
    }
  }
};



//! BasefunctionSet for quadratic finite element functions on triangular meshes
template <typename DataTypeContainer, typename QuadType, typename TriangleType>
class RefTriangMeshBaseFunctionSetP3red : 
public RefTriangMeshBaseFunctionSetPkInterface< DataTypeContainer, QuadType, TriangleType, RefTriangMeshBaseFunctionSetP3red<DataTypeContainer, QuadType, TriangleType> >  {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
      
public:
  // ordering anti-clockwise
  static RealType _b1   ( const DomVecType &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 1./2. * c2 * ( 3. * c2 - 1. ) * ( 3. * c2 - 2. ) - 9./2. * c[0] * c[1] * c2;
  }
  static RealType _b2   ( const DomVecType &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 9./2. * c2 * c[0] * ( 3. * c2 - 1. ) + 27./4. * c[0] * c[1] * c2;
  }
  static RealType _b3   ( const DomVecType &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 9./2. * c2 * c[0] * ( 3. * c[0] - 1. ) + 27./4. * c[0] * c[1] * c2;
  }
  static RealType _b4   ( const DomVecType &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 1./2. * c[0] * ( 3. * c[0] - 1. ) * ( 3. * c[0] - 2. ) - 9./2. * c[0] * c[1] * c2;
  }
  static RealType _b5   ( const DomVecType &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 9./2. * c[0] * c[1] * ( 3. * c[0] - 1. ) + 27./4. * c[0] * c[1] * c2;
  }
  static RealType _b6   ( const DomVecType &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 9./2. * c[1] * c[0] * ( 3. * c[1] - 1. ) + 27./4. * c[0] * c[1] * c2;
  }
  static RealType _b7   ( const DomVecType &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 1./2. * c[1] * ( 3. * c[1] - 1. ) * ( 3. * c[1] - 2. ) - 9./2. * c[0] * c[1] * c2; 
  }
  static RealType _b8   ( const DomVecType &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 9./2. * c[1] * c2 * ( 3. * c[1] - 1. ) + 27./4. * c[0] * c[1] * c2; 
  }
  static RealType _b9   ( const DomVecType &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 9./2. * c2 * c[1] * ( 3. * c2 - 1. ) + 27./4. * c[0] * c[1] * c2;
  }

  //
  static RealType _d1_b1   ( const DomVecType &c ) { return -11./2. + 18. * c[0] - 27./2. * pesopt::Sqr(c[0]) + 27./2. * c[1] - 18. * c[0] * c[1] - 9. * pesopt::Sqr(c[1]); }
  static RealType _d1_b2   ( const DomVecType &c ) { return 9. - 45. * c[0] + 81./2. * pesopt::Sqr(c[0]) - 63./4. * c[1] + 81./2. * c[0] * c[1] + 27./4. * pesopt::Sqr(c[1]); }
  static RealType _d1_b3   ( const DomVecType &c ) { return - 9./2. + 36. * c[0] - 81./2. * pesopt::Sqr(c[0]) + 45./4. * c[1] - 81./2. * c[0] * c[1] - 27./4. * pesopt::Sqr(c[1]); }
  static RealType _d1_b4   ( const DomVecType &c ) { return 1. - 9. * c[0] + 27./2. * pesopt::Sqr(c[0]) - 9./2. * c[1] + 9. * c[0] * c[1] + 9./2. * pesopt::Sqr(c[1]); }
  static RealType _d1_b5   ( const DomVecType &c ) { return 9./4. * c[1] + 27./2. * c[0] * c[1] - 27./4. * pesopt::Sqr(c[1]); }
  static RealType _d1_b6   ( const DomVecType &c ) { return 9./4. * c[1] - 27./2. * c[0] * c[1] + 27./4. * pesopt::Sqr(c[1]); }
  static RealType _d1_b7   ( const DomVecType &c ) { return - 9./2. * c[1] + 9. * c[0] * c[1] + 9./2. * pesopt::Sqr(c[1]); }
  static RealType _d1_b8   ( const DomVecType &c ) { return 45./4. * c[1] - 27./2. * c[0] * c[1] - 81./4. * pesopt::Sqr(c[1]); }
  static RealType _d1_b9   ( const DomVecType &c ) { return - 63./4. * c[1] + 27./2. * c[0] * c[1] + 81./4. * pesopt::Sqr(c[1]); }
  //
  static RealType _d2_b1   ( const DomVecType &c ) { return -11./2. + 27./2. * c[0] - 9. * pesopt::Sqr(c[0]) + 18. * c[1] - 18. * c[0] * c[1] - 27./2. * pesopt::Sqr(c[1]); }
  static RealType _d2_b2   ( const DomVecType &c ) { return - 63./4. * c[0] + 81./4. * pesopt::Sqr(c[0]) + 27./2. * c[0] * c[1] ; }
  static RealType _d2_b3   ( const DomVecType &c ) { return 45./4. * c[0] - 81./4. * pesopt::Sqr(c[0]) - 27./2. * c[0] * c[1] ;}
  static RealType _d2_b4   ( const DomVecType &c ) { return -9./2. * c[0] + 9./2. * pesopt::Sqr(c[0]) + 9. * c[0] * c[1] ;}
  static RealType _d2_b5   ( const DomVecType &c ) { return 9./4. * c[0] + 27./4. * pesopt::Sqr(c[0]) - 27./2. * c[0] * c[1] ;}
  static RealType _d2_b6   ( const DomVecType &c ) { return 9./4. * c[0] - 27./4. * pesopt::Sqr(c[0]) + 27./2. * c[0] * c[1] ;}
  static RealType _d2_b7   ( const DomVecType &c ) { return 1. - 9./2. * c[0] + 9./2. * pesopt::Sqr(c[0]) - 9. * c[1] + 9. * c[0] * c[1] + 27./2. * pesopt::Sqr(c[1]);}
  static RealType _d2_b8   ( const DomVecType &c ) { return - 9./2. + 45./4. * c[0] - 27./4. * pesopt::Sqr(c[0]) + 36. * c[1] - 81./2. * c[0] * c[1] - 81./2. * pesopt::Sqr(c[1]);}
  static RealType _d2_b9   ( const DomVecType &c ) { return 9. - 63./4. * c[0] + 27./4. * pesopt::Sqr(c[0]) - 45. * c[1] + 81./2. * c[0] * c[1] + 81./2. * pesopt::Sqr(c[1]);}
  
  //
  static RealType _d11_b1   ( const DomVecType &c ) { return 18. - 27. * c[0] - 18. * c[1]; }
  static RealType _d11_b2   ( const DomVecType &c ) { return - 45. + 81. * c[0] + 81./2. * c[1]; }
  static RealType _d11_b3   ( const DomVecType &c ) { return 36. - 81. * c[0] - 81./2. * c[1]; }
  static RealType _d11_b4   ( const DomVecType &c ) { return - 9. + 27. * c[0] + 9. * c[1]; }
  static RealType _d11_b5   ( const DomVecType &c ) { return 27./2. * c[1] ; }
  static RealType _d11_b6   ( const DomVecType &c ) { return - 27./2. * c[1] ; }
  static RealType _d11_b7   ( const DomVecType &c ) { return 9. * c[1] ; }
  static RealType _d11_b8   ( const DomVecType &c ) { return - 27./2. * c[1] ; }
  static RealType _d11_b9   ( const DomVecType &c ) { return + 27./2. * c[1] ; }
  //
  static RealType _d12_b1   ( const DomVecType &c ) { return 27./2. - 18. * c[0] - 18. * c[1] ; }
  static RealType _d12_b2   ( const DomVecType &c ) { return - 63./4. + 81./2. * c[0] + 27./2. * c[1] ; }
  static RealType _d12_b3   ( const DomVecType &c ) { return 45./4. - 81./2. * c[0] - 27./2. * c[1] ;}
  static RealType _d12_b4   ( const DomVecType &c ) { return -9./2. + 9. * c[0] + 9. * c[1] ;}
  static RealType _d12_b5   ( const DomVecType &c ) { return 9./4. + 27./2. * c[0] - 27./2. * c[1] ;}
  static RealType _d12_b6   ( const DomVecType &c ) { return 9./4. - 27./2. * c[0] + 27./2. * c[1] ;}
  static RealType _d12_b7   ( const DomVecType &c ) { return - 9./2. + 9. * c[0]  + 9. * c[1];}
  static RealType _d12_b8   ( const DomVecType &c ) { return 45./4. - 27./2. * c[0] - 81./2. * c[1];}
  static RealType _d12_b9   ( const DomVecType &c ) { return - 63./4. + 27./2. * c[0] + 81./2. * c[1];}
  //
  static RealType _d22_b1   ( const DomVecType &c ) { return 18. - 18. * c[0] - 27. * c[1]; }
  static RealType _d22_b2   ( const DomVecType &c ) { return 27./2. * c[0]; }
  static RealType _d22_b3   ( const DomVecType &c ) { return - 27./2. * c[0];}
  static RealType _d22_b4   ( const DomVecType &c ) { return 9. * c[0];}
  static RealType _d22_b5   ( const DomVecType &c ) { return - 27./2. * c[0];}
  static RealType _d22_b6   ( const DomVecType &c ) { return 27./2. * c[0];}
  static RealType _d22_b7   ( const DomVecType &c ) { return - 9. + 9. * c[0] + 27. * c[1];}
  static RealType _d22_b8   ( const DomVecType &c ) { return 36. - 81./2. * c[0] - 81. * c[1];}
  static RealType _d22_b9   ( const DomVecType &c ) { return  - 45. + 81./2. * c[0] + 81. * c[1];}
protected:
  typedef RealType ( *BASIS_FUNC_TYPE ) ( const DomVecType &RefCoord );
  BASIS_FUNC_TYPE _basis[9];
  BASIS_FUNC_TYPE _dbasis[2][9];
  BASIS_FUNC_TYPE _d2basis[3][9];

public:
  RefTriangMeshBaseFunctionSetP3red ( ) :
  RefTriangMeshBaseFunctionSetPkInterface< DataTypeContainer, QuadType, TriangleType, RefTriangMeshBaseFunctionSetP3red<DataTypeContainer, QuadType, TriangleType> > () {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;
    _basis[6] = _b7;
    _basis[7] = _b8;
    _basis[8] = _b9;

    _dbasis[0][0] = _d1_b1;
    _dbasis[0][1] = _d1_b2;
    _dbasis[0][2] = _d1_b3;
    _dbasis[0][3] = _d1_b4;
    _dbasis[0][4] = _d1_b5;
    _dbasis[0][5] = _d1_b6;
    _dbasis[0][6] = _d1_b7;
    _dbasis[0][7] = _d1_b8;
    _dbasis[0][8] = _d1_b9;

    _dbasis[1][0] = _d2_b1;
    _dbasis[1][1] = _d2_b2;
    _dbasis[1][2] = _d2_b3;
    _dbasis[1][3] = _d2_b4;
    _dbasis[1][4] = _d2_b5;
    _dbasis[1][5] = _d2_b6;
    _dbasis[1][6] = _d2_b7;
    _dbasis[1][7] = _d2_b8;
    _dbasis[1][8] = _d2_b9;
    
    _d2basis[0][0] = _d11_b1;
    _d2basis[0][1] = _d11_b2;
    _d2basis[0][2] = _d11_b3;
    _d2basis[0][3] = _d11_b4;
    _d2basis[0][4] = _d11_b5;
    _d2basis[0][5] = _d11_b6;
    _d2basis[0][6] = _d11_b7;
    _d2basis[0][7] = _d11_b8;
    _d2basis[0][8] = _d11_b9;

    _d2basis[1][0] = _d22_b1;
    _d2basis[1][1] = _d22_b2;
    _d2basis[1][2] = _d22_b3;
    _d2basis[1][3] = _d22_b4;
    _d2basis[1][4] = _d22_b5;
    _d2basis[1][5] = _d22_b6;
    _d2basis[1][6] = _d22_b7;
    _d2basis[1][7] = _d22_b8;
    _d2basis[1][8] = _d22_b9;
    
    _d2basis[2][0] = _d12_b1;
    _d2basis[2][1] = _d12_b2;
    _d2basis[2][2] = _d12_b3;
    _d2basis[2][3] = _d12_b4;
    _d2basis[2][4] = _d12_b5;
    _d2basis[2][5] = _d12_b6;
    _d2basis[2][6] = _d12_b7;
    _d2basis[2][7] = _d12_b8;
    _d2basis[2][8] = _d12_b9;

  }

  enum { numBaseFuncs = 9 };
  
  RealType evaluateOnUnitTriang ( int BaseFuncNum, const DomVecType &RefCoord ) const { return _basis[BaseFuncNum] ( RefCoord );}
  
  void evaluateGradientOnUnitTriang ( int BaseFuncNum, const DomVecType &RefCoord, DomVecType &Gradient ) const {
    Gradient[0] = _dbasis[0][BaseFuncNum](RefCoord); Gradient[1] = _dbasis[1][BaseFuncNum](RefCoord);
  }
  
  void evaluateHessianSymOnUnitTriang ( int BaseFuncNum, const DomVecType &RefCoord, RealType &d11, RealType &d22, RealType &d12 ) const {
    d11 = _d2basis[0][BaseFuncNum](RefCoord); 
    d22 = _d2basis[1][BaseFuncNum](RefCoord);
    d12 = _d2basis[2][BaseFuncNum](RefCoord);
  }

};







template <typename DataTypeContainer, typename QuadType, typename TriangleType>
class RefTriangMeshBaseFunctionSetDKT : public RefTriangleBaseFunctionSetInterface< DataTypeContainer, QuadType > {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
  typedef typename DataTypeContainer::Matrix22 Matrix22;
  typedef RefTriangMeshBaseFunctionSetP3red<DataTypeContainer, QuadType, TriangleType> P3RedFct; 
 
  // ordering anti-clockwise
  static RealType _w1( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_b8(c) + 20.*P3RedFct::_b9(c) + 27.*P3RedFct::_b1(c) + 20.*P3RedFct::_b2(c) + 7.*P3RedFct::_b3(c) ) / 27.;     
  }
  
  static RealType _w2( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(2)[0]*P3RedFct::_b2(c) + 2.*t->getEdgeRefTriang(2)[0]*P3RedFct::_b3(c) - 2.*t->getEdgeRefTriang(1)[0]*P3RedFct::_b8(c) - 4.*t->getEdgeRefTriang(1)[0]*P3RedFct::_b9(c) ) / 27.;     
  }
  
  static RealType _w3( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(2)[1]*P3RedFct::_b2(c) + 2.*t->getEdgeRefTriang(2)[1]*P3RedFct::_b3(c) - 2.*t->getEdgeRefTriang(1)[1]*P3RedFct::_b8(c) - 4.*t->getEdgeRefTriang(1)[1]*P3RedFct::_b9(c) ) / 27.;     
  }
  
  static RealType _w4( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_b2(c) + 20.*P3RedFct::_b3(c) + 27.*P3RedFct::_b4(c) + 20.*P3RedFct::_b5(c) + 7.*P3RedFct::_b6(c) ) / 27.;     
  }  
    
  static RealType _w5( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(0)[0]*P3RedFct::_b5(c) + 2.*t->getEdgeRefTriang(0)[0]*P3RedFct::_b6(c) - 2.*t->getEdgeRefTriang(2)[0]*P3RedFct::_b2(c) - 4.*t->getEdgeRefTriang(2)[0]*P3RedFct::_b3(c) ) / 27.;     
  }
  
  static RealType _w6( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(0)[1]*P3RedFct::_b5(c) + 2.*t->getEdgeRefTriang(0)[1]*P3RedFct::_b6(c) - 2.*t->getEdgeRefTriang(2)[1]*P3RedFct::_b2(c) - 4.*t->getEdgeRefTriang(2)[1]*P3RedFct::_b3(c) ) / 27.;     
  }  
    
  static RealType _w7( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_b5(c) + 20.*P3RedFct::_b6(c) + 27.*P3RedFct::_b7(c) + 20.*P3RedFct::_b8(c) + 7.*P3RedFct::_b9(c) ) / 27.;     
  }  
    
  static RealType _w8( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(1)[0]*P3RedFct::_b8(c) + 2.*t->getEdgeRefTriang(1)[0]*P3RedFct::_b9(c) - 2.*t->getEdgeRefTriang(0)[0]*P3RedFct::_b5(c) - 4.*t->getEdgeRefTriang(0)[0]*P3RedFct::_b6(c) ) / 27.;     
  }
  
  static RealType _w9( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(1)[1]*P3RedFct::_b8(c) + 2.*t->getEdgeRefTriang(1)[1]*P3RedFct::_b9(c) - 2.*t->getEdgeRefTriang(0)[1]*P3RedFct::_b5(c) - 4.*t->getEdgeRefTriang(0)[1]*P3RedFct::_b6(c) ) / 27.;     
  }
  
  // derivative with respect to \xi
  static RealType _d1_w1( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d1_b8(c) + 20.*P3RedFct::_d1_b9(c) + 27.*P3RedFct::_d1_b1(c) + 20.*P3RedFct::_d1_b2(c) + 7.*P3RedFct::_d1_b3(c) ) / 27.;     
  }
  
  static RealType _d1_w2( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d1_b2(c) + 2.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d1_b3(c) - 2.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d1_b8(c) - 4.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d1_b9(c) ) / 27.;     
  }
  
  static RealType _d1_w3( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d1_b2(c) + 2.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d1_b3(c) - 2.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d1_b8(c) - 4.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d1_b9(c) ) / 27.;     
  }
  
  static RealType _d1_w4( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d1_b2(c) + 20.*P3RedFct::_d1_b3(c) + 27.*P3RedFct::_d1_b4(c) + 20.*P3RedFct::_d1_b5(c) + 7.*P3RedFct::_d1_b6(c) ) / 27.;     
  }  
    
  static RealType _d1_w5( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d1_b5(c) + 2.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d1_b6(c) - 2.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d1_b2(c) - 4.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d1_b3(c) ) / 27.;     
  }
  
  static RealType _d1_w6( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d1_b5(c) + 2.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d1_b6(c) - 2.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d1_b2(c) - 4.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d1_b3(c) ) / 27.;     
  }  
    
  static RealType _d1_w7( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d1_b5(c) + 20.*P3RedFct::_d1_b6(c) + 27.*P3RedFct::_d1_b7(c) + 20.*P3RedFct::_d1_b8(c) + 7.*P3RedFct::_d1_b9(c) ) / 27.;     
  }  
    
  static RealType _d1_w8( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d1_b8(c) + 2.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d1_b9(c) - 2.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d1_b5(c) - 4.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d1_b6(c) ) / 27.;     
  }
  
  static RealType _d1_w9( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d1_b8(c) + 2.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d1_b9(c) - 2.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d1_b5(c) - 4.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d1_b6(c) ) / 27.;     
  }
  
 // derivative with respect to \eta
  static RealType _d2_w1( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d2_b8(c) + 20.*P3RedFct::_d2_b9(c) + 27.*P3RedFct::_d2_b1(c) + 20.*P3RedFct::_d2_b2(c) + 7.*P3RedFct::_d2_b3(c) ) / 27.;     
  }
  
  static RealType _d2_w2( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d2_b2(c) + 2.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d2_b3(c) - 2.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d2_b8(c) - 4.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d2_b9(c) ) / 27.;     
  }
  
  static RealType _d2_w3( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d2_b2(c) + 2.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d2_b3(c) - 2.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d2_b8(c) - 4.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d2_b9(c) ) / 27.;     
  }
  
  static RealType _d2_w4( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d2_b2(c) + 20.*P3RedFct::_d2_b3(c) + 27.*P3RedFct::_d2_b4(c) + 20.*P3RedFct::_d2_b5(c) + 7.*P3RedFct::_d2_b6(c) ) / 27.;     
  }  
    
  static RealType _d2_w5( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d2_b5(c) + 2.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d2_b6(c) - 2.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d2_b2(c) - 4.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d2_b3(c) ) / 27.;     
  }
  
  static RealType _d2_w6( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d2_b5(c) + 2.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d2_b6(c) - 2.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d2_b2(c) - 4.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d2_b3(c) ) / 27.;     
  }  
    
  static RealType _d2_w7( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d2_b5(c) + 20.*P3RedFct::_d2_b6(c) + 27.*P3RedFct::_d2_b7(c) + 20.*P3RedFct::_d2_b8(c) + 7.*P3RedFct::_d2_b9(c) ) / 27.;     
  }  
    
  static RealType _d2_w8( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d2_b8(c) + 2.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d2_b9(c) - 2.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d2_b5(c) - 4.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d2_b6(c) ) / 27.;     
  }
  
  static RealType _d2_w9( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d2_b8(c) + 2.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d2_b9(c) - 2.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d2_b5(c) - 4.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d2_b6(c) ) / 27.;     
  }  
  
  
  // second derivatives
  static RealType _d11_w1( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d11_b8(c) + 20.*P3RedFct::_d11_b9(c) + 27.*P3RedFct::_d11_b1(c) + 20.*P3RedFct::_d11_b2(c) + 7.*P3RedFct::_d11_b3(c) ) / 27.;     
  }
  
  static RealType _d11_w2( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d11_b2(c) + 2.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d11_b3(c) - 2.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d11_b8(c) - 4.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d11_b9(c) ) / 27.;     
  }
  
  static RealType _d11_w3( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d11_b2(c) + 2.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d11_b3(c) - 2.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d11_b8(c) - 4.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d11_b9(c) ) / 27.;     
  }
  
  static RealType _d11_w4( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d11_b2(c) + 20.*P3RedFct::_d11_b3(c) + 27.*P3RedFct::_d11_b4(c) + 20.*P3RedFct::_d11_b5(c) + 7.*P3RedFct::_d11_b6(c) ) / 27.;     
  }  
    
  static RealType _d11_w5( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d11_b5(c) + 2.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d11_b6(c) - 2.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d11_b2(c) - 4.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d11_b3(c) ) / 27.;     
  }
  
  static RealType _d11_w6( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d11_b5(c) + 2.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d11_b6(c) - 2.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d11_b2(c) - 4.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d11_b3(c) ) / 27.;     
  }  
    
  static RealType _d11_w7( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d11_b5(c) + 20.*P3RedFct::_d11_b6(c) + 27.*P3RedFct::_d11_b7(c) + 20.*P3RedFct::_d11_b8(c) + 7.*P3RedFct::_d11_b9(c) ) / 27.;     
  }  
    
  static RealType _d11_w8( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d11_b8(c) + 2.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d11_b9(c) - 2.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d11_b5(c) - 4.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d11_b6(c) ) / 27.;     
  }
  
  static RealType _d11_w9( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d11_b8(c) + 2.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d11_b9(c) - 2.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d11_b5(c) - 4.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d11_b6(c) ) / 27.;     
  }
  
  
  static RealType _d22_w1( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d22_b8(c) + 20.*P3RedFct::_d22_b9(c) + 27.*P3RedFct::_d22_b1(c) + 20.*P3RedFct::_d22_b2(c) + 7.*P3RedFct::_d22_b3(c) ) / 27.;     
  }
  
  static RealType _d22_w2( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d22_b2(c) + 2.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d22_b3(c) - 2.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d22_b8(c) - 4.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d22_b9(c) ) / 27.;     
  }
  
  static RealType _d22_w3( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d22_b2(c) + 2.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d22_b3(c) - 2.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d22_b8(c) - 4.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d22_b9(c) ) / 27.;     
  }
  
  static RealType _d22_w4( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d22_b2(c) + 20.*P3RedFct::_d22_b3(c) + 27.*P3RedFct::_d22_b4(c) + 20.*P3RedFct::_d22_b5(c) + 7.*P3RedFct::_d22_b6(c) ) / 27.;     
  }  
    
  static RealType _d22_w5( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d22_b5(c) + 2.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d22_b6(c) - 2.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d22_b2(c) - 4.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d22_b3(c) ) / 27.;     
  }
  
  static RealType _d22_w6( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d22_b5(c) + 2.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d22_b6(c) - 2.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d22_b2(c) - 4.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d22_b3(c) ) / 27.;     
  }  
    
  static RealType _d22_w7( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d22_b5(c) + 20.*P3RedFct::_d22_b6(c) + 27.*P3RedFct::_d22_b7(c) + 20.*P3RedFct::_d22_b8(c) + 7.*P3RedFct::_d22_b9(c) ) / 27.;     
  }  
    
  static RealType _d22_w8( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d22_b8(c) + 2.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d22_b9(c) - 2.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d22_b5(c) - 4.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d22_b6(c) ) / 27.;     
  }
  
  static RealType _d22_w9( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d22_b8(c) + 2.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d22_b9(c) - 2.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d22_b5(c) - 4.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d22_b6(c) ) / 27.;     
  }  
  
  // derivative with respect to \xi
  static RealType _d12_w1( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d12_b8(c) + 20.*P3RedFct::_d12_b9(c) + 27.*P3RedFct::_d12_b1(c) + 20.*P3RedFct::_d12_b2(c) + 7.*P3RedFct::_d12_b3(c) ) / 27.;     
  }
  
  static RealType _d12_w2( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d12_b2(c) + 2.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d12_b3(c) - 2.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d12_b8(c) - 4.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d12_b9(c) ) / 27.;     
  }
  
  static RealType _d12_w3( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d12_b2(c) + 2.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d12_b3(c) - 2.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d12_b8(c) - 4.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d12_b9(c) ) / 27.;     
  }
  
  static RealType _d12_w4( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d12_b2(c) + 20.*P3RedFct::_d12_b3(c) + 27.*P3RedFct::_d12_b4(c) + 20.*P3RedFct::_d12_b5(c) + 7.*P3RedFct::_d12_b6(c) ) / 27.;     
  }  
    
  static RealType _d12_w5( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d12_b5(c) + 2.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d12_b6(c) - 2.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d12_b2(c) - 4.*t->getEdgeRefTriang(2)[0]*P3RedFct::_d12_b3(c) ) / 27.;     
  }
  
  static RealType _d12_w6( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d12_b5(c) + 2.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d12_b6(c) - 2.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d12_b2(c) - 4.*t->getEdgeRefTriang(2)[1]*P3RedFct::_d12_b3(c) ) / 27.;     
  }  
    
  static RealType _d12_w7( const TriangleType* /*t*/, const DomVecType &c ) { 
    return ( 7.*P3RedFct::_d12_b5(c) + 20.*P3RedFct::_d12_b6(c) + 27.*P3RedFct::_d12_b7(c) + 20.*P3RedFct::_d12_b8(c) + 7.*P3RedFct::_d12_b9(c) ) / 27.;     
  }  
    
  static RealType _d12_w8( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d12_b8(c) + 2.*t->getEdgeRefTriang(1)[0]*P3RedFct::_d12_b9(c) - 2.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d12_b5(c) - 4.*t->getEdgeRefTriang(0)[0]*P3RedFct::_d12_b6(c) ) / 27.;     
  }
  
  static RealType _d12_w9( const TriangleType* t, const DomVecType &c ) { 
    return ( 4.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d12_b8(c) + 2.*t->getEdgeRefTriang(1)[1]*P3RedFct::_d12_b9(c) - 2.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d12_b5(c) - 4.*t->getEdgeRefTriang(0)[1]*P3RedFct::_d12_b6(c) ) / 27.;     
  }
  
  
  //
  typedef RealType ( *BASIS_FUNC_TYPE ) ( const TriangleType* Triangle, const DomVecType &RefCoord );
  BASIS_FUNC_TYPE _basis[9];
  BASIS_FUNC_TYPE _dbasis[2][9];
  BASIS_FUNC_TYPE _d2basis[3][9];

  const TriangleType *_triangle;

public:
  RefTriangMeshBaseFunctionSetDKT(  ) : _triangle ( NULL ) {
    _basis[0] = _w1;
    _basis[1] = _w2;
    _basis[2] = _w3;
    _basis[3] = _w4;
    _basis[4] = _w5;
    _basis[5] = _w6;
    _basis[6] = _w7;
    _basis[7] = _w8;
    _basis[8] = _w9;

    _dbasis[0][0] = _d1_w1;
    _dbasis[0][1] = _d1_w2;
    _dbasis[0][2] = _d1_w3;
    _dbasis[0][3] = _d1_w4;
    _dbasis[0][4] = _d1_w5;
    _dbasis[0][5] = _d1_w6;
    _dbasis[0][6] = _d1_w7;
    _dbasis[0][7] = _d1_w8;
    _dbasis[0][8] = _d1_w9;

    _dbasis[1][0] = _d2_w1;
    _dbasis[1][1] = _d2_w2;
    _dbasis[1][2] = _d2_w3;
    _dbasis[1][3] = _d2_w4;
    _dbasis[1][4] = _d2_w5;
    _dbasis[1][5] = _d2_w6;
    _dbasis[1][6] = _d2_w7;
    _dbasis[1][7] = _d2_w8;
    _dbasis[1][8] = _d2_w9;
    
    
    _d2basis[0][0] = _d11_w1;
    _d2basis[0][1] = _d11_w2;
    _d2basis[0][2] = _d11_w3;
    _d2basis[0][3] = _d11_w4;
    _d2basis[0][4] = _d11_w5;
    _d2basis[0][5] = _d11_w6;
    _d2basis[0][6] = _d11_w7;
    _d2basis[0][7] = _d11_w8;
    _d2basis[0][8] = _d11_w9;

    _d2basis[1][0] = _d22_w1;
    _d2basis[1][1] = _d22_w2;
    _d2basis[1][2] = _d22_w3;
    _d2basis[1][3] = _d22_w4;
    _d2basis[1][4] = _d22_w5;
    _d2basis[1][5] = _d22_w6;
    _d2basis[1][6] = _d22_w7;
    _d2basis[1][7] = _d22_w8;
    _d2basis[1][8] = _d22_w9;
    
    _d2basis[2][0] = _d12_w1;
    _d2basis[2][1] = _d12_w2;
    _d2basis[2][2] = _d12_w3;
    _d2basis[2][3] = _d12_w4;
    _d2basis[2][4] = _d12_w5;
    _d2basis[2][5] = _d12_w6;
    _d2basis[2][6] = _d12_w7;
    _d2basis[2][7] = _d12_w8;
    _d2basis[2][8] = _d12_w9;
    
  }

  enum { numBaseFuncs = 9 };

  void setTriangle ( const TriangleType &T ) {_triangle = &T;}
  
  
  //****************************************************
  // on reference triangle
  //****************************************************
  RealType evaluateOnRefTriang ( int BaseFuncNum, const DomVecType &RefCoord ) const { return _basis[BaseFuncNum] ( _triangle, RefCoord ); }
  inline const RealType evaluateOnRefTriang ( int BaseFuncNum, int QuadPoint ) const { return this->evaluateOnRefTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );}
  
  void evaluateGradientOnRefTriang ( int BaseFuncNum, const DomVecType &RefCoord, DomVecType &Gradient ) const {
    DomVecType GradUnitTriang; 
    GradUnitTriang[0] = _dbasis[0][BaseFuncNum] ( _triangle, RefCoord ); 
    GradUnitTriang[1] = _dbasis[1][BaseFuncNum] ( _triangle, RefCoord );
    DomVecType ginvDx; ginvDx = _triangle->getInverseMetricRefTriang( ) * GradUnitTriang;  
    const DomVecType dir0 = _triangle->edgeRefTriang(0,1);
    const DomVecType dir1 = _triangle->edgeRefTriang(0,2); 
    Gradient = ginvDx[0] * dir0 + ginvDx[1] * dir1;
  }

  inline DomVecType evaluateGradientOnRefTriang ( int BaseFuncNum, int QuadPoint ) const {
    DomVecType Gradient;
    this->evaluateGradientOnRefTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), Gradient );
    return Gradient;
  }
  
  //Here: use that the map A:UnitTriang \to \RefTriang is linear (hence D^2A = 0)
  void evaluateHessianOnRefTriang ( int BaseFuncNum, const DomVecType &RefCoord, Matrix22 &Hessian ) const {
    Matrix22 HessianUnitTriang; 
    HessianUnitTriang(0,0) = _d2basis[0][BaseFuncNum] ( _triangle, RefCoord ); 
    HessianUnitTriang(1,1) = _d2basis[1][BaseFuncNum] ( _triangle, RefCoord );
    HessianUnitTriang(1,0) = _d2basis[2][BaseFuncNum] ( _triangle, RefCoord );
    HessianUnitTriang(0,1) = _d2basis[2][BaseFuncNum] ( _triangle, RefCoord );
    Matrix22 ginvSqrD2x; ginvSqrD2x = _triangle->getGradientMapInvRefTriang().transpose() * HessianUnitTriang * _triangle->getGradientMapInvRefTriang();
    Hessian = ginvSqrD2x;
  }

  inline Matrix22 evaluateHessianOnRefTriang ( int BaseFuncNum, int QuadPoint ) const {
    Matrix22 Hessian;
    this->evaluateHessianOnRefTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), Hessian );
    return Hessian;
  }
  
  //*****************************************************
  
  
};





/*! \brief BasefunctionSet for approximative gradient (i.e. \f$ \theta \approx \nabla w \f$) in Discrete Kirchhoff Triangle (DKT)
 *  The 9 DOFs per element T are denoted by \f$ W = [ w_0 w_{x0} w_{y0} w_1 w_{x1} w_{y1} w_2 w_{x2} w_{y2} ] \f$,
 *  where "w_i" refers to the function value of "w" at the vertex in "T" with local index "i"
 *  and "w_{zi}" to the corresponding derivative with respect to "z".
 *  
 *  The base functions are given by vector valued functions \f$ \theta \f$,
 *  i.e. \f$ \theta = (\theta^1, \theta^2): T \rightarrow \R^2, (\xi, \eta) \mapsto \theta(\xi, \eta) \f$.
 *  We shall use the notation \f$ w(\xi,\eta) = ( H_x W, H_y W )^T \f$.
 * 
 *  In applications with DKT the most relevant object is the stiffnes matrix \f$ L_ij = \int C \epsilon[\theta_i] : \epsilon[\theta_j] dx \f$,
 *  where \f$ \epsilon[\theta] \in \R^{2,2} \f$ denotes the symmetrized gradient and \f$ C \f$ is the elastic tensor, 
 *  i.e. \f$ \epsilon[\theta]_ij = 1/2 ( \partial_i\theta^j + \partial_j\theta^i) \f$.
 *  However, the evaluation of the gradient of \f$ \theta \f$ returns \f$ ( \partial_1\theta^1, \partial_2\theta^2, \partial_1\theta^2 + \partial_2\theta^1 ) \in \R^3 \f$.
 */
template <typename DataTypeContainer, typename QuadType, typename TriangleType>
class RefTriangMeshApproxGradientBaseFunctionSetDKT : public RefTriangleBaseFunctionSetInterface< DataTypeContainer, QuadType > {
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
  typedef typename DataTypeContainer::Point3DType Point3DType;
  typedef typename DataTypeContainer::Matrix22 Matrix22;
  typedef typename DataTypeContainer::Matrix32 Matrix32;  
  
  typedef RefTriangMeshBaseFunctionSetP2<DataTypeContainer, QuadType, TriangleType > P2Fct;
 
  const TriangleType *_triangle;
  
 public:
  // Numbering of local DOFs corresponds to order in W, i.e. 0 = w_0, 1 = w_{x0}, ..., 8 = w_{y2} 
  // Ordering of vertices anti-clockwise beginning at node zero, i.e. _triangle-> globNodeIdx( 0 ) 
  
  //P2Fct
  //static RealType _b1   ( const DomVecType &c ) { return 2. * ( 1. - c[0] - c[1] ) * ( 1./2. - c[0] - c[1] ); }
  //static RealType _b2   ( const DomVecType &c ) { return 4. * c[0] * ( 1. - c[0] - c[1] ); }
  //static RealType _b3   ( const DomVecType &c ) { return c[0] * ( 2. * c[0] - 1. ); }
  //static RealType _b4   ( const DomVecType &c ) { return 4. * c[0] * c[1]; }
  //static RealType _b5   ( const DomVecType &c ) { return c[1] * ( 2. * c[1] - 1. ); }
  //static RealType _b6   ( const DomVecType &c ) { return 4. * c[1] * ( 1. - c[0] - c[1] ); }  
  
     
  // Hx
  static RealType _tet11 ( const TriangleType* t, const DomVecType &c ) { 
      return 1.5 * ( t->edgeXNormSqr(1) * P2Fct::_b6(c) - t->edgeXNormSqr(2) * P2Fct::_b2(c) );
  }
  static RealType _tet12 ( const TriangleType* t, const DomVecType &c ) {
      return P2Fct::_b1(c) - 0.25 * ( t->edgeXXMin2YYNormSqr(2) * P2Fct::_b2(c) + t->edgeXXMin2YYNormSqr(1) * P2Fct::_b6(c) );
  }
  static RealType _tet13 ( const TriangleType* t, const DomVecType &c ) {
      return -0.75 * ( t->edgeXYNormSqr(2) * P2Fct::_b2(c) + t->edgeXYNormSqr(1) * P2Fct::_b6(c) );
  }
  static RealType _tet14 ( const TriangleType* t, const DomVecType &c ) {
      return 1.5 * ( t->edgeXNormSqr(2) * P2Fct::_b2(c) - t->edgeXNormSqr(0) * P2Fct::_b4(c) );
  } 
  static RealType _tet15 ( const TriangleType* t, const DomVecType &c ) { 
      return P2Fct::_b3(c) - 0.25 * ( t->edgeXXMin2YYNormSqr(2) * P2Fct::_b2(c) + t->edgeXXMin2YYNormSqr(0) * P2Fct::_b4(c) );
  }
  static RealType _tet16 ( const TriangleType* t, const DomVecType &c ) { 
      return -0.75 * ( t->edgeXYNormSqr(2) * P2Fct::_b2(c) + t->edgeXYNormSqr(0) * P2Fct::_b4(c) ); 
  }
  static RealType _tet17 ( const TriangleType* t, const DomVecType &c ) {
      return 1.5 * ( t->edgeXNormSqr(0) * P2Fct::_b4(c) - t->edgeXNormSqr(1) * P2Fct::_b6(c) );
  } 
  static RealType _tet18 ( const TriangleType* t, const DomVecType &c ) {
      return P2Fct::_b5(c) - 0.25 * ( t->edgeXXMin2YYNormSqr(0) * P2Fct::_b4(c) + t->edgeXXMin2YYNormSqr(1) * P2Fct::_b6(c) );
  }
  static RealType _tet19 ( const TriangleType* t, const DomVecType &c ) {
      return -0.75 * ( t->edgeXYNormSqr(0) * P2Fct::_b4(c) + t->edgeXYNormSqr(1) * P2Fct::_b6(c) );
  }

  //Hy
  static RealType _tet21 ( const TriangleType* t, const DomVecType &c ) { 
      return 1.5 * ( t->edgeYNormSqr(1) * P2Fct::_b6(c) - t->edgeYNormSqr(2) * P2Fct::_b2(c) );
  }
  static RealType _tet22 ( const TriangleType* t, const DomVecType &c ) {
      return -0.75 * ( t->edgeXYNormSqr(2) * P2Fct::_b2(c) + t->edgeXYNormSqr(1) * P2Fct::_b6(c) );
  }
  static RealType _tet23 ( const TriangleType* t, const DomVecType &c ) {
      return P2Fct::_b1(c) - 0.25 * ( t->edgeYYMin2XXNormSqr(2) * P2Fct::_b2(c) + t->edgeYYMin2XXNormSqr(1) * P2Fct::_b6(c) );
  }
  static RealType _tet24 ( const TriangleType* t, const DomVecType &c ) { 
      return 1.5 * ( t->edgeYNormSqr(2) * P2Fct::_b2(c) - t->edgeYNormSqr(0) * P2Fct::_b4(c) );
  }
  static RealType _tet25 ( const TriangleType* t, const DomVecType &c ) {
      return -0.75 * ( t->edgeXYNormSqr(2) * P2Fct::_b2(c) + t->edgeXYNormSqr(0) * P2Fct::_b4(c) );
  }
  static RealType _tet26 ( const TriangleType* t, const DomVecType &c ) {
      return P2Fct::_b3(c) - 0.25 * (t->edgeYYMin2XXNormSqr(2) * P2Fct::_b2(c) + t->edgeYYMin2XXNormSqr(0) * P2Fct::_b4(c) );
  }
  static RealType _tet27 ( const TriangleType* t, const DomVecType &c ) { 
      return 1.5 * ( t->edgeYNormSqr(0) * P2Fct::_b4(c) - t->edgeYNormSqr(1) * P2Fct::_b6(c) );
  }
  static RealType _tet28 ( const TriangleType* t, const DomVecType &c ) {
      return -0.75 * ( t->edgeXYNormSqr(0) * P2Fct::_b4(c) + t->edgeXYNormSqr(1) * P2Fct::_b6(c) );
  }
  static RealType _tet29 ( const TriangleType* t, const DomVecType &c ) { 
      return P2Fct::_b5(c) - 0.25 * (t->edgeYYMin2XXNormSqr(0) * P2Fct::_b4(c) + t->edgeYYMin2XXNormSqr(1) * P2Fct::_b6(c) );
  }
     
     
     
     
  // d Hx / d \xi
  static RealType _d1_tet11   ( const TriangleType* t, const DomVecType &c ) { 
      return - 6. * t->edgeXNormSqr(2) * ( 1. - 2. * c[0] ) + 6. * ( t->edgeXNormSqr(2) - t->edgeXNormSqr(1) ) * c[1];
  }
  static RealType _d1_tet12   ( const TriangleType* t, const DomVecType &c ) { 
      return -1. * t->edgeXXMin2YYNormSqr(2) * ( 1. - 2. * c[0] ) + 4. * ( c[0] + c[1] ) + ( t->edgeXXMin2YYNormSqr(1) + t->edgeXXMin2YYNormSqr(2) ) * c[1] - 3.;
  }
  static RealType _d1_tet13   ( const TriangleType* t, const DomVecType &c ) { 
      return -3. * t->edgeXYNormSqr(2) * ( 1. - 2. * c[0] ) + 3. * ( t->edgeXYNormSqr(1) + t->edgeXYNormSqr(2) ) * c[1];
  }
  static RealType _d1_tet14   ( const TriangleType* t, const DomVecType &c ) {
      return 6. * t->edgeXNormSqr(2) * ( 1. - 2. * c[0] ) - 6. * ( t->edgeXNormSqr(0) + t->edgeXNormSqr(2) ) * c[1] ;
  }  
  static RealType _d1_tet15   ( const TriangleType* t, const DomVecType &c ) {
      return -1. * t->edgeXXMin2YYNormSqr(2) * ( 1. - 2. * c[0] ) + ( t->edgeXXMin2YYNormSqr(2) - t->edgeXXMin2YYNormSqr(0) ) * c[1] + 4.*c[0] - 1.;
  }
  static RealType _d1_tet16   ( const TriangleType* t, const DomVecType &c ) { 
      return -3. * t->edgeXYNormSqr(2) * (1. - 2. * c[0] ) + 3. * ( t->edgeXYNormSqr(2) - t->edgeXYNormSqr(0) ) * c[1];
  }
  static RealType _d1_tet17   ( const TriangleType* t, const DomVecType &c ) { 
      return 6. * ( t->edgeXNormSqr(0) + t->edgeXNormSqr(1) ) * c[1];
  } 
  static RealType _d1_tet18   ( const TriangleType* t, const DomVecType &c ) { 
      return ( t->edgeXXMin2YYNormSqr(1) - t->edgeXXMin2YYNormSqr(0) ) * c[1];
  }
  static RealType _d1_tet19   ( const TriangleType* t, const DomVecType &c ) { 
      return 3. * ( t->edgeXYNormSqr(1) - t->edgeXYNormSqr(0) ) * c[1];
  }

  // Hx,\eta
  static RealType _d2_tet11   ( const TriangleType* t, const DomVecType &c ) { 
    return 6. * t->edgeXNormSqr(1) * ( 1. - 2. * c[1] ) + 6. * ( t->edgeXNormSqr(2) - t->edgeXNormSqr(1) ) * c[0];
  }
  static RealType _d2_tet12   ( const TriangleType* t, const DomVecType &c ) { 
    return -1. * t->edgeXXMin2YYNormSqr(1) * ( 1. - 2. * c[1] ) + 4. * ( c[0] + c[1] ) + ( t->edgeXXMin2YYNormSqr(2) + t->edgeXXMin2YYNormSqr(1) ) * c[0]  - 3.;
  }    
  static RealType _d2_tet13   ( const TriangleType* t, const DomVecType &c ) { 
    return -3. * t->edgeXYNormSqr(1) * ( 1. - 2 * c[1] ) + 3.*( t->edgeXYNormSqr(2) + t->edgeXYNormSqr(1) ) * c[0];
  }  
  static RealType _d2_tet14   ( const TriangleType* t, const DomVecType &c ) { 
    return -6. * ( t->edgeXNormSqr(2) + t->edgeXNormSqr(0) ) * c[0];
  }  
  static RealType _d2_tet15   ( const TriangleType* t, const DomVecType &c ) { 
    return ( t->edgeXXMin2YYNormSqr(2) - t->edgeXXMin2YYNormSqr(0) ) * c[0];
  }
  static RealType _d2_tet16   ( const TriangleType* t, const DomVecType &c ) { 
    return 3. * ( t->edgeXYNormSqr(2) - t->edgeXYNormSqr(0) ) * c[0];
  }
  static RealType _d2_tet17   ( const TriangleType* t, const DomVecType &c ) { 
    return -6. * t->edgeXNormSqr(1) * ( 1. - 2. * c[1] ) + 6. * ( t->edgeXNormSqr(0) + t->edgeXNormSqr(1) ) * c[0];
  }
  static RealType _d2_tet18   ( const TriangleType* t, const DomVecType &c ) {   
    return -1. * t->edgeXXMin2YYNormSqr(1) * ( 1. - 2. * c[1] ) + 4. * c[1]   + ( t->edgeXXMin2YYNormSqr(1) - t->edgeXXMin2YYNormSqr(0) ) * c[0] - 1.;
  }
  static RealType _d2_tet19   ( const TriangleType* t, const DomVecType &c ) { 
    return -3. * t->edgeXYNormSqr(1) * ( 1. - 2. * c[1] ) + 3.*( t->edgeXYNormSqr(1) - t->edgeXYNormSqr(0) ) * c[0];
  }
  
  
  // Hy,\xi
  static RealType _d1_tet21   ( const TriangleType* t, const DomVecType &c ) { 
    return  -6. * t->edgeYNormSqr(2) * (1. - 2. * c[0]) + 6. * (t->edgeYNormSqr(2) - t->edgeYNormSqr(1)) * c[1];
  }
  static RealType _d1_tet22   ( const TriangleType* t, const DomVecType &c ) { 
    return -3. * t->edgeXYNormSqr(2) * ( 1. - 2. * c[0] ) + 3. * ( t->edgeXYNormSqr(2) + t->edgeXYNormSqr(1) ) * c[1];
  }
  static RealType _d1_tet23   ( const TriangleType* t, const DomVecType &c ) {   
    return  -1. * t->edgeYYMin2XXNormSqr(2) * ( 1. - 2. * c[0] ) + 4. * ( c[0] + c[1] ) + ( t->edgeYYMin2XXNormSqr(2) + t->edgeYYMin2XXNormSqr(1) ) * c[1] - 3.;
  }
  static RealType _d1_tet24   ( const TriangleType* t, const DomVecType &c ) {   
    return 6. * t->edgeYNormSqr(2) * ( 1. - 2. * c[0] ) - 6. * ( t->edgeYNormSqr(2) + t->edgeYNormSqr(0) ) * c[1];
  }  
  static RealType _d1_tet25   ( const TriangleType* t, const DomVecType &c ) {  
    return -3. * t->edgeXYNormSqr(2) * ( 1. - 2. * c[0] ) + 3. * ( t->edgeXYNormSqr(2) - t->edgeXYNormSqr(0) ) * c[1];
  }
  static RealType _d1_tet26   ( const TriangleType* t, const DomVecType &c ) {  
    return  -1. * t->edgeYYMin2XXNormSqr(2) * ( 1. - 2. * c[0] ) + ( t->edgeYYMin2XXNormSqr(2) - t->edgeYYMin2XXNormSqr(0) ) *  c[1] + 4. * c[0] - 1.;
  }
  static RealType _d1_tet27   ( const TriangleType* t, const DomVecType &c ) {   
    return 6. * ( t->edgeYNormSqr(0) + t->edgeYNormSqr(1) ) *  c[1];
  } 
  static RealType _d1_tet28   ( const TriangleType* t, const DomVecType &c ) {   
    return 3. * ( t->edgeXYNormSqr(1) - t->edgeXYNormSqr(0) ) * c[1];
  }
  static RealType _d1_tet29   ( const TriangleType* t, const DomVecType &c ) { 
    return ( t->edgeYYMin2XXNormSqr(1) - t->edgeYYMin2XXNormSqr(0) ) * c[1];
  }
  
  // d Hy / d \eta
  static RealType _d2_tet21   ( const TriangleType* t, const DomVecType &c ) { 
    return  6. * t->edgeYNormSqr(1) * ( 1. - 2. * c[1] ) + 6. * ( t->edgeYNormSqr(2) - t->edgeYNormSqr(1) ) * c[0];
  }
  static RealType _d2_tet22   ( const TriangleType* t, const DomVecType &c ) { 
    return  -3. * t->edgeXYNormSqr(1) *  ( 1. - 2. * c[1] ) + 3. * ( t->edgeXYNormSqr(2) + t->edgeXYNormSqr(1) ) * c[0];
  }
  static RealType _d2_tet23   ( const TriangleType* t, const DomVecType &c ) {   
    return -1. * t->edgeYYMin2XXNormSqr(1) * ( 1. - 2. * c[1] ) + 4.*( c[0] + c[1] ) + ( t->edgeYYMin2XXNormSqr(2) + t->edgeYYMin2XXNormSqr(1) ) * c[0] - 3.;
  }
  static RealType _d2_tet24   ( const TriangleType* t, const DomVecType &c ) {   
    return -6. * ( t->edgeYNormSqr(2) + t->edgeYNormSqr(0) ) * c[0];
  }  
  static RealType _d2_tet25   ( const TriangleType* t, const DomVecType &c ) {  
    return 3. * ( t->edgeXYNormSqr(2) - t->edgeXYNormSqr(0) ) * c[0];
  }
  static RealType _d2_tet26   ( const TriangleType* t, const DomVecType &c ) {  
    return ( t->edgeYYMin2XXNormSqr(2) - t->edgeYYMin2XXNormSqr(0) ) * c[0];
  }
  static RealType _d2_tet27   ( const TriangleType* t, const DomVecType &c ) {   
    return  - 6. * t->edgeYNormSqr(1) * ( 1. - 2. * c[1] ) + 6. * ( t->edgeYNormSqr(0) + t->edgeYNormSqr(1) ) * c[0];
  }
  static RealType _d2_tet28   ( const TriangleType* t, const DomVecType &c ) {   
    return  -3. * t->edgeXYNormSqr(1) * ( 1. - 2. * c[1] ) + 3.*( t->edgeXYNormSqr(1) - t->edgeXYNormSqr(0) ) * c[0];
  }
  static RealType _d2_tet29   ( const TriangleType* t, const DomVecType &c ) { 
    return -1. * t->edgeYYMin2XXNormSqr(1) * ( 1. - 2. * c[1] ) + ( t->edgeYYMin2XXNormSqr(1) - t->edgeYYMin2XXNormSqr(0) ) * c[0]  + 4. * c[1] - 1.;
  }

protected:
  //
  typedef RealType ( *BASIS_FUNC_TYPE ) ( const TriangleType* Triangle, const DomVecType &RefCoord );
  BASIS_FUNC_TYPE _basis1[9];
  BASIS_FUNC_TYPE _basis2[9];
  BASIS_FUNC_TYPE _dbasis1[2][9];
  BASIS_FUNC_TYPE _dbasis2[2][9];


public:
  RefTriangMeshApproxGradientBaseFunctionSetDKT(  ) : _triangle ( NULL ) {
    // Hx
    _basis1[0] = _tet11;
    _basis1[1] = _tet12;
    _basis1[2] = _tet13;
    _basis1[3] = _tet14;
    _basis1[4] = _tet15;
    _basis1[5] = _tet16;
    _basis1[6] = _tet17;
    _basis1[7] = _tet18;
    _basis1[8] = _tet19;
    
    // Hy
    _basis2[0] = _tet21;
    _basis2[1] = _tet22;
    _basis2[2] = _tet23;
    _basis2[3] = _tet24;
    _basis2[4] = _tet25;
    _basis2[5] = _tet26;
    _basis2[6] = _tet27;
    _basis2[7] = _tet28;
    _basis2[8] = _tet29;
    
    // Hx,\xi
    _dbasis1[0][0] = _d1_tet11;
    _dbasis1[0][1] = _d1_tet12;
    _dbasis1[0][2] = _d1_tet13;
    _dbasis1[0][3] = _d1_tet14;
    _dbasis1[0][4] = _d1_tet15;
    _dbasis1[0][5] = _d1_tet16;
    _dbasis1[0][6] = _d1_tet17;
    _dbasis1[0][7] = _d1_tet18;
    _dbasis1[0][8] = _d1_tet19;
    
    // Hx,\eta
    _dbasis1[1][0] = _d2_tet11;
    _dbasis1[1][1] = _d2_tet12;
    _dbasis1[1][2] = _d2_tet13;
    _dbasis1[1][3] = _d2_tet14;
    _dbasis1[1][4] = _d2_tet15;
    _dbasis1[1][5] = _d2_tet16;
    _dbasis1[1][6] = _d2_tet17;
    _dbasis1[1][7] = _d2_tet18;
    _dbasis1[1][8] = _d2_tet19;
    
    // Hy,\xi
    _dbasis2[0][0] = _d1_tet21;
    _dbasis2[0][1] = _d1_tet22;
    _dbasis2[0][2] = _d1_tet23;
    _dbasis2[0][3] = _d1_tet24;
    _dbasis2[0][4] = _d1_tet25;
    _dbasis2[0][5] = _d1_tet26;
    _dbasis2[0][6] = _d1_tet27;
    _dbasis2[0][7] = _d1_tet28;
    _dbasis2[0][8] = _d1_tet29;
    
    // Hy, \eta
    _dbasis2[1][0] = _d2_tet21;
    _dbasis2[1][1] = _d2_tet22;
    _dbasis2[1][2] = _d2_tet23;
    _dbasis2[1][3] = _d2_tet24;
    _dbasis2[1][4] = _d2_tet25;
    _dbasis2[1][5] = _d2_tet26;
    _dbasis2[1][6] = _d2_tet27;
    _dbasis2[1][7] = _d2_tet28;
    _dbasis2[1][8] = _d2_tet29;
  }

  enum { numBaseFuncs = 9 };

  void setTriangle ( const TriangleType &T ) { _triangle = &T;}

  
  //****************************************************
  // on unit triangle (ATTENTION: still depends on triangle)
  //****************************************************
  
  //! Theta
  void evaluateApproxGradientOnUnitTriang ( int BaseFuncNum, const DomVecType &RefCoord, DomVecType &Gradient ) const {
    Gradient[0] = _basis1[BaseFuncNum](_triangle, RefCoord); Gradient[1] = _basis2[BaseFuncNum](_triangle, RefCoord);
  }

  inline DomVecType evaluateApproxGradientOnUnitTriang ( int BaseFuncNum, int QuadPoint ) const {
    DomVecType Gradient;
    this->evaluateApproxGradientOnUnitTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), Gradient );
    return Gradient;
  }
  //! D Theta
  void evaluateApproxHessianOnUnitTriang ( int BaseFuncNum, const DomVecType &RefCoord, Matrix22 &Hessian ) const {    
    Hessian(0,0) = _dbasis1[0][BaseFuncNum](_triangle, RefCoord);
    Hessian(1,0) = _dbasis1[1][BaseFuncNum](_triangle, RefCoord);
    Hessian(0,1) = _dbasis2[0][BaseFuncNum](_triangle, RefCoord);
    Hessian(1,1) = _dbasis2[1][BaseFuncNum](_triangle, RefCoord);
  }
  inline Matrix22 evaluateApproxHessianOnUnitTriang ( int BaseFuncNum, int QuadPoint ) const {    
    Matrix22 h;    
    evaluateApproxHessianOnUnitTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), h );
    return h;
  }
  //*****************************************************
  

  
  
  
  
  
  //****************************************************
  // on reference triangle
  //****************************************************
  
  void evaluateApproxGradientOnRefTriang ( int BaseFuncNum, const DomVecType &RefCoord, DomVecType &Gradient ) const {
      DomVecType GradUnitTriang; this->evaluateApproxGradientOnUnitTriang( BaseFuncNum, RefCoord, GradUnitTriang );
      Gradient[0] = GradUnitTriang[0];
      Gradient[1] = GradUnitTriang[1];
      //TODO
  }
  
  inline DomVecType evaluateApproxGradientOnRefTriang ( int BaseFuncNum, int QuadPoint ) const {
    DomVecType Gradient;
    this->evaluateApproxGradientOnRefTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), Gradient );
    return Gradient;
  }
  
  
   //DTheta  
   //TODO Matrix32 Hessian = Dx g^-1 HessianUnitTriang -> meaning in tangent vec directions??
//    void evaluateApproxHessianOnRefTriang ( int BaseFuncNum, const DomVecType &RefCoord, Matrix22 &Hessian ) const {    
//     const Matrix22& DxGinv = _triangle->getGradientTimesInverseMetricRefTriang();
//     Matrix22 HessianUnitTriang; this->evaluateApproxHessianOnUnitTriang( BaseFuncNum, RefCoord, HessianUnitTriang );
//     Hessian = DxGinv * HessianUnitTriang;
//   }
//   inline Matrix22 evaluateApproxHessianOnRefTriang ( int BaseFuncNum, int QuadPoint ) const {    
//     Matrix22 h;    
//     evaluateApproxHessianOnRefTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), h );
//     return h;
//   }
   
   
   
//    void evaluateApproxHessianAsVecOnRefTriang ( int BaseFuncNum, const DomVecType &RefCoord, Point3DType &Hessian ) const {    
//     const Matrix22& ginv = _triangle->getInverseMetricRefTriangDKT();
//     
//     RealType Hx1( _dbasis1[0][BaseFuncNum](_triangle, RefCoord) ), Hx2( _dbasis1[1][BaseFuncNum](_triangle, RefCoord) );
//     RealType Hy1( _dbasis2[0][BaseFuncNum](_triangle, RefCoord) ), Hy2( _dbasis2[1][BaseFuncNum](_triangle, RefCoord) );
//     
//     Hessian[0] = ginv(0,0) * Hx1 + ginv(1,0) * Hx2;
//     Hessian[1] = ginv(0,1) * Hy1 + ginv(1,1) * Hy2;
//     Hessian[2] = ginv(0,1) * Hx1 + ginv(1,1) * Hx2 + ginv(0,0) * Hy1 + ginv(1,0) * Hy2;
//   }
  void evaluateApproxHessianAsVecOnRefTriang ( int BaseFuncNum, const DomVecType &RefCoord, Point3DType &HessianVec ) const {    
    const Matrix22& DxGinv = _triangle->getGradientTimesInverseMetricRefTriang();
    
    RealType Hx1( _dbasis1[0][BaseFuncNum](_triangle, RefCoord) ), Hx2( _dbasis1[1][BaseFuncNum](_triangle, RefCoord) );
    RealType Hy1( _dbasis2[0][BaseFuncNum](_triangle, RefCoord) ), Hy2( _dbasis2[1][BaseFuncNum](_triangle, RefCoord) );
    
    HessianVec[0] = DxGinv(0,0) * Hx1 + DxGinv(0,1) * Hx2;
    HessianVec[1] = DxGinv(1,0) * Hy1 + DxGinv(1,1) * Hy2;
    HessianVec[2] = DxGinv(1,0) * Hx1 + DxGinv(1,1) * Hx2 + DxGinv(0,0) * Hy1 + DxGinv(0,1) * Hy2;
    //TODO
//        Matrix22 Hessian; this->evaluateApproxHessianOnRefTriang( BaseFuncNum, RefCoord, Hessian );
//        HessianVec(0) = Hessian(0,0);
//        HessianVec(1) = Hessian(1,1);
//        HessianVec(2) = Hessian(0,1) + Hessian(1,0);
  }
  inline Point3DType evaluateApproxHessianAsVecOnRefTriang ( int BaseFuncNum, int QuadPoint ) const {    
    Point3DType h;    
    evaluateApproxHessianAsVecOnRefTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), h );
    return h;
  }
   

  void evaluateApproxHessianOnRefTriang ( int BaseFuncNum, const DomVecType &RefCoord, Matrix22 &Hessian ) const { 
    const Matrix22& DxGinv = _triangle->getGradientTimesInverseMetricRefTriang();
    RealType Hx1( _dbasis1[0][BaseFuncNum](_triangle, RefCoord) ), Hx2( _dbasis1[1][BaseFuncNum](_triangle, RefCoord) );
    RealType Hy1( _dbasis2[0][BaseFuncNum](_triangle, RefCoord) ), Hy2( _dbasis2[1][BaseFuncNum](_triangle, RefCoord) );
    Hessian(0,0) = DxGinv(0,0) * Hx1 + DxGinv(0,1) * Hx2;
    Hessian(1,1) = DxGinv(1,0) * Hy1 + DxGinv(1,1) * Hy2;
    Hessian(0,1) = DxGinv(1,0) * Hx1 + DxGinv(1,1) * Hx2;
    Hessian(1,0) = DxGinv(0,0) * Hy1 + DxGinv(0,1) * Hy2;
  }
  inline Matrix22 evaluateApproxHessianOnRefTriang ( int BaseFuncNum, int QuadPoint ) const {    
    Matrix22 h;    
    evaluateApproxHessianOnRefTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), h );
    return h;
  }
   


  void evaluateApproxHessianSymOnRefTriang ( int BaseFuncNum, const DomVecType &RefCoord, Matrix22 &HessianSym ) const { 
      Point3DType HessianVec; evaluateApproxHessianAsVecOnRefTriang( BaseFuncNum, RefCoord, HessianVec );
      HessianSym(0,0) = HessianVec(0);
      HessianSym(1,1) = HessianVec(1);
      HessianSym(0,1) = 0.5 * HessianVec(2);
      HessianSym(1,0) = 0.5 * HessianVec(2);
  }
  inline Matrix22 evaluateApproxHessianSymOnRefTriang ( int BaseFuncNum, int QuadPoint ) const {    
    Matrix22 h;    
    evaluateApproxHessianSymOnRefTriang ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), h );
    return h;
  }
  
  //*****************************************************
  
};




#endif //__DKTFEBASEFUNCTIONSETS_H
