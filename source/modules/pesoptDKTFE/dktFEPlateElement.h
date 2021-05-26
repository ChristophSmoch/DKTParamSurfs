#ifndef __DKTFEPLATEELEMENT_H
#define __DKTFEPLATEELEMENT_H

#include <pesopt_IO.h>

// TODO need tangent space???
//TODO projection onto tangent space???
template<typename DataTypeContainer>
class DKTPlateElement {

public:
    
  typedef typename DataTypeContainer::RealType          RealType;
  typedef typename DataTypeContainer::Point3DType       Point3DType;
  typedef typename DataTypeContainer::TangentVecType    TangentVecType;
  typedef typename DataTypeContainer::Indices3DType     Indices3DType;
  typedef typename DataTypeContainer::Matrix22          Matrix22;
  typedef typename DataTypeContainer::Matrix32          Matrix32;
  typedef typename DataTypeContainer::Matrix33          Matrix33;
  typedef typename DataTypeContainer::DomVecType        DomVecType;
  
  typedef std::vector<Point3DType>                      VertexIterator;
  typedef std::vector<TangentVecType>                   TangentVecIterator;
  
public :   
      
    //****************************************
    //Flat triangle in R^3
    //****************************************  
    // global indices of element and nodes
    int _globIdx;
    Indices3DType _globNodeIdx;
    Point3DType _nodes[3];
    TangentVecType _edges[3]; // Notation:   e_i = N_{i-1} - N_{i+1} = N_{i+2} - N_{i+1}, i = 0,1,2
    Point3DType _midPoint;
    
    //****************************************
    //Flat Reference triangle n0, n1, n2
    //****************************************
    DomVecType _nodesRefTriang[3];
    DomVecType _edgesRefTriang[3]; // Notation:   e_i = N_{i-1} - N_{i+1} = N_{i+2} - N_{i+1}, i = 0,1,2
    RealType _edgesRefTriangLengthSqr[3];
    //given x: T_UnitTriang in R^2 -> T_RefTriang in R^3
    Matrix22 _ginvRefTriang; // (Dx^T Dx)^-1
    Matrix22 _GradMapRefTriang; // Dx
    Matrix22 _GradMapInvRefTriang; // Dx^-1
    Matrix22 _GradGinvRefTriang; // Dx g^-1 (TODO = Dx^{-T}
    RealType _areaRefTriang;
    //******************************************
    
    
//     //****************************************
//     //DKT Reference triangle (with edges projected on the tangent space 
//     //****************************************
//     TangentVecType _projectedEdgesRefTriang[3][2]; //
//     RealType _projectedEdgesRefTriangLengthSqr[3][2]; 
//     Matrix22 _projectedGinvRefTriang[3]; // (Dx^T Dx)^-1
// 
//     // corresponding coefficients: projVec = projCoeff1 v_1 + projCoeff2 v_2. Save as projcoeff(localNodeIndex)(Direction)
//     RealType _projCoeff1[3][2];
//     RealType _projCoeff2[3][2];
//     //******************************************

    
public:
    
  DKTPlateElement () : _globIdx(-1) {}
  
  DKTPlateElement( const int globalIdx, const Indices3DType globalNodeIndex, const VertexIterator &nodes  ) : 
     _globIdx(globalIdx), _globNodeIdx ( globalNodeIndex ) {
         // fill nodes and edges
        for ( int i = 0; i < 3; ++i ) _nodes[i] = nodes[ _globNodeIdx[i] ];
        for ( int i = 0; i < 3; ++i ) _edges[i] = getNode((i+2)%3) - getNode((i+1)%3);
        _midPoint = 1./3. * ( _nodes[0] + _nodes[1] + _nodes[2] );
        this->updateRefTriangNodes( nodes );
        for ( int i = 0; i < 3; ++i ) _edgesRefTriang[i] = getNodeRefTriang((i+2)%3) - getNodeRefTriang((i+1)%3);
        for ( int i = 0; i < 3; ++i ) _edgesRefTriangLengthSqr[i] = _edgesRefTriang[i].squaredNorm();
        this->setInverseMetricRefTriangle();
    }
  
  ~DKTPlateElement(){}
  
  
  //Let v1, v2 basis of tangent space, write w as linear combination of v1, v2 and (v1 x v2): w = a1 * v1 + a2 * v2 + a3 * (v1 x v2), return coefficient a1 and a2
//   void computeCoefficientsOfProjection (const TangentVecType& v1, const TangentVecType& v2, const TangentVecType& w, 
//                                         RealType& a1, RealType& a2 ){
//     RealType normV1Squared = v1.dot(v1); RealType normV2Squared = v2.dot(v2);
//     RealType scalarProductV1V2 = v1.dot(v2); RealType scalarProductV1W = v1.dot(w); RealType scalarProductV2W = v2.dot(w);
//     RealType scalarProductV1V2Squared = scalarProductV1V2 * scalarProductV1V2;
//     
//     a1 = normV2Squared * scalarProductV1W - scalarProductV1V2 * scalarProductV2W;
//     a1 /= (normV1Squared * normV2Squared - scalarProductV1V2Squared);
//     
//     a2 = scalarProductV2W * normV1Squared - scalarProductV1W * scalarProductV1V2;
//     a2 /= (normV1Squared * normV2Squared - scalarProductV1V2Squared);
//   }

  
  
  void updateRefTriangNodes( const VertexIterator &nodes ){
      for ( int i = 0; i < 3; ++i ) {
          if( std::abs( nodes[ _globNodeIdx[i] ](2) ) > 1.e-16 ){ 
               cout << "nodeIdx = " << _globNodeIdx[i] << endl
                    << "coords = " << nodes[_globNodeIdx[i] ].transpose() << endl;
               throw std::invalid_argument( pesopt::strprintf ( "Mesh is not contained in xy-Plane. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
           }
           _nodesRefTriang[i] << _nodes[i](0), _nodes[i](1);
       }
  }
  
  void updateNodesAndEdges( const VertexIterator &nodes ){
      for ( int i = 0; i < 3; ++i ) _nodes[i] = nodes[ _globNodeIdx[i] ];
      for ( int i = 0; i < 3; ++i ) _edges[i] = getNode((i+2)%3) - getNode((i+1)%3);
      this->updateRefTriangNodes( nodes );
      for ( int i = 0; i < 3; ++i ) _edgesRefTriang[i] = getNodeRefTriang((i+2)%3) - getNodeRefTriang((i+1)%3);
      for ( int i = 0; i < 3; ++i ) _edgesRefTriangLengthSqr[i] = _edgesRefTriang[i].squaredNorm();
      this->setInverseMetricRefTriangle();
  }
  
  
  //TODO
//   void updateProjectionCoefficients( const TangentVecIterator &tangentVec1, const TangentVecIterator &tangentVec2 ){ 
//     
//     // project edges to tangent space
//     RealType tmpCoeff1, tmpCoeff2;
//     TangentVecType vecToProject;
//     
//     // at node N_0 = (0,0) TODO -> replace by RefNode!!!
//     // vec corresponding to (1,0)
//     computeCoefficientsOfProjection ( tangentVec1 [ _globNodeIdx[0] ], tangentVec2 [ _globNodeIdx[0] ], _edges[2], _projCoeff1[0][0], _projCoeff2[0][0] );
//     // vec corresponding to (0,1)
//     vecToProject = -1.0 * _edges[1];
//     computeCoefficientsOfProjection ( tangentVec1 [ _globNodeIdx[0] ], tangentVec2 [ _globNodeIdx[0] ], vecToProject, _projCoeff1[0][1], _projCoeff2[0][1] );
//     
//     // at node N_1 = (1,0) TODO -> replace by RefNode!!!
//     // vec corresponding to (1,0)
//     vecToProject = -1.0 * _edges[2];
//     computeCoefficientsOfProjection ( tangentVec1 [ _globNodeIdx[1] ], tangentVec2 [ _globNodeIdx[1] ], vecToProject, tmpCoeff1 , tmpCoeff2  );
//     _projCoeff1 [1][0] = - 1.0 * tmpCoeff1;
//     _projCoeff2 [1][0] = - 1.0 * tmpCoeff2;
//     // vec corresponding to (0,1)
//     _projCoeff1 [1][1] = - 1.0 * tmpCoeff1;
//     _projCoeff2 [1][1] = - 1.0 * tmpCoeff2;
//     computeCoefficientsOfProjection ( tangentVec1 [ _globNodeIdx[1] ], tangentVec2 [ _globNodeIdx[1] ], _edges[0], tmpCoeff1 , tmpCoeff2 );
//     _projCoeff1 [1][1] +=  tmpCoeff1;
//     _projCoeff2 [1][1] +=  tmpCoeff2;
// 
//     // at node N_2 = (0,1) TODO -> replace by RefNode!!!
//     // vec corresponding to (0,1)
//     computeCoefficientsOfProjection ( tangentVec1 [ _globNodeIdx[2] ], tangentVec2 [ _globNodeIdx[2] ], _edges[1], tmpCoeff1 , tmpCoeff2  );
//     _projCoeff1 [2][1] = - 1.0 * tmpCoeff1;
//     _projCoeff2 [2][1] = - 1.0 * tmpCoeff2;
//     // vec corresponding to (1,0)
//     _projCoeff1 [2][0] = -1.0 * tmpCoeff1;
//     _projCoeff2 [2][0] = -1.0 * tmpCoeff2;
//     vecToProject = -1.0 * _edges[0];
//     computeCoefficientsOfProjection ( tangentVec1 [ _globNodeIdx[2] ], tangentVec2 [ _globNodeIdx[2] ], vecToProject, tmpCoeff1 , tmpCoeff2 );
//     _projCoeff1 [2][0] +=  tmpCoeff1;
//     _projCoeff2 [2][0] +=  tmpCoeff2;
//     
//   }

  void printNodes() const {
      for(int localIdx=0;localIdx<3;++localIdx){
      cout << "localNodeIndex = " << localIdx << " globalNodeIndex = " << getGlobalNodeIdx(localIdx) << endl 
           << "coords = " << getNode(localIdx).transpose() << endl
           << "coords ref triang = " << getNodeRefTriang(localIdx).transpose() << endl;
      }
  }
  
  void printEdges() const {
      for(int localIdx=0;localIdx<3;++localIdx){
      cout << "localEdge = " << localIdx << endl 
           << "edge = " << getEdge(localIdx).transpose() << endl
           << "edge ref triang = " << getEdgeRefTriang(localIdx).transpose() << endl;
      }
  }

  void print() const {
      cout << endl 
           << "---------------------------------------------------" << endl
           << "Element " << _globIdx << endl
           << "---------------------------------------------------" << endl;
      this->printNodes();
      this->printEdges();
//       cout << "projcoeff = " << endl;
//       for( int i=0; i<3; ++i )
//           for( int j=0; j<2; ++j )
//            cout << _projCoeff1[i][j] << "  ,    " << _projCoeff2[i][j] << endl;   
        
      cout << "area of FlattenedTriang = " << this->getAreaOfFlattenedTriangle() << endl;
      cout << "area of RefTriang = " << this->getAreaOfRefTriangle() << endl;
          
    for( int i = 0; i < 3; i++ ){
      cerr << "a_" << i+1 << ((i+1)%3)+1 << " = " << edgeYNormSqr ( i ) << endl;
      cerr << "b_" << i+1 << ((i+1)%3)+1 << " = " << edgeXYNormSqr ( i )<< endl;
      cerr << "c_" << i+1 << ((i+1)%3)+1 << " = " << edgeXXMin2YYNormSqr ( i )<< endl;
      cerr << "d_" << i+1 << ((i+1)%3)+1 << " = " << edgeXNormSqr ( i )<< endl;
      cerr << "e_" << i+1 << ((i+1)%3)+1 << " = " << edgeYYMin2XXNormSqr ( i )<< endl << endl;
    }
  }
  
  
  //The following methods compute the weighted/normalized normal of the triangle elements and its area under the assumption that the triangle is flat
  void getWeightedNormalForFlattenedTriangle ( TangentVecType &normal ) const {
    TangentVecType dx1 = getNode(1), dx2 = getNode(2);
    dx1 -= getNode(0);
    dx2 -= getNode(0);
    normal = TangentVecType( dx1.cross(dx2) );
  }

  void getNormalizedNormalForFlattenedTriangle ( TangentVecType &normal ) const {
    getWeightedNormalForFlattenedTriangle ( normal );
    normal.normalize();
  }

  inline RealType getAreaOfFlattenedTriangle ( ) const {
    TangentVecType n; getWeightedNormalForFlattenedTriangle ( n );
    return n.norm() / 2.;
  }
    
  // get and set functions
  int getGlobalElementIdx(  ) const { return _globIdx;}
  const Indices3DType & getGlobalNodeIdx( ) const{ return _globNodeIdx;}
  int getGlobalNodeIdx(int localIndex) const{ return _globNodeIdx[localIndex];}
  void setGlobalNodeIdx(int localIndex, int globalIndex) { _globNodeIdx[localIndex] = globalIndex;}
  
  const Point3DType& getNode ( int i ) const { return _nodes[i];}
  Point3DType& getNode ( int i ) { return _nodes[i];}
  void setNode ( int i, const Point3DType& node ) {_nodes[i] = node;}
  
  const TangentVecType& getEdge( int localEdgeIndex ) const { return _edges[localEdgeIndex];};
  const TangentVecType getEdge ( const int locNode1, const int locNode2 ) const { return (this->getNode(locNode2) - this->getNode(locNode1) );  }
  
  
  const Point3DType& getMidPoint (  ) const { return _midPoint;}
  
  const TangentVecType& operator[] ( int i ) const { return getNode(i);}
  TangentVecType& operator[] ( int i ) {return getNode(i);}
  
  //LocalCoord wrt. unit triangle (0,0) - (1,0) - (0,1)
  void getGlobalCoord ( const DomVecType &LocalCoord, TangentVecType &Coord ) const {
      TangentVecType dx1 = this->getEdge(0,1), dx2 = this->getEdge(0,2);
      Coord = LocalCoord[0] * dx1 + LocalCoord[1] * dx2;
      Coord += this->getNode(0);                            //Coord w.r.t. 0, 0
  }
  
  void getRefCoordFromGlobalCoord ( const TangentVecType &Coord,  DomVecType &LocalCoord ) const {
      Matrix33 mat;
      mat.col(0) = this->getEdge(0,1);
      mat.col(1) = this->getEdge(0,2);
      mat.col(2) = this->getEdge(0,1).cross( this->getEdge(0,2) );
      Matrix33 matInv = mat.inverse();
      TangentVecType tmp = matInv * ( Coord - this->getNode(0) );
      LocalCoord(0) = tmp(0); LocalCoord(1) = tmp(1);
  }
  
  
//   RealType getProjCoeff1 ( const int BaseFuncNum ) const {
//     int localNodeIndex = std::floor(BaseFuncNum/3);
//     int direction = BaseFuncNum%3;
//     direction -= 1;
//     return _projCoeff1[localNodeIndex][direction];
//   }
//   
//   RealType getProjCoeff2 ( const int BaseFuncNum ) const {
//     int localNodeIndex = std::floor(BaseFuncNum/3);
//     int direction = BaseFuncNum%3;
//     direction -= 1;
//     return _projCoeff2[localNodeIndex][direction];
//   }
  
  void getRefCoordsFromLocalIndex( const int localNodeIndex, DomVecType &RefCoords ) const{
      RefCoords.setZero();
      switch( localNodeIndex ){
          case 0 : break;
          case 1 : RefCoords[0] = 1.; break;
          case 2 : RefCoords[1] = 1.; break;
          default: throw std::invalid_argument( pesopt::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
  }
  
  void getRefCoordsFromGlobalNodeIdx( const int globalNodeIndex, DomVecType &RefCoords ) const{
     RefCoords.setZero();
     if( globalNodeIndex == _globNodeIdx[0] ) return;
     if( globalNodeIndex == _globNodeIdx[1] ) { RefCoords[0] = 1.; return;}
     if( globalNodeIndex == _globNodeIdx[2] ) { RefCoords[1] = 1.; return;}
     throw std::invalid_argument( pesopt::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  //================================================
  // Ref triangle 
  //================================================
  const DomVecType& getNodeRefTriang ( int i ) const { return _nodesRefTriang[i];}
  DomVecType& getNodeRefTriang ( int i ) { return _nodesRefTriang[i];}
  void setNodeRefTriang ( int i, const DomVecType& node ) {_nodesRefTriang[i] = node;}
  DomVecType edgeRefTriang ( const int locNode1, const int locNode2 ) const { return (this->getNodeRefTriang(locNode2) - this->getNodeRefTriang(locNode1) );  }
  const DomVecType& getEdgeRefTriang( int localEdgeIndex ) const { return _edgesRefTriang[localEdgeIndex];};
  RealType getAreaOfRefTriangle( ) const { return _areaRefTriang; }
  
  const Matrix22 &getInverseMetricRefTriang( ) const { return _ginvRefTriang; }
  const Matrix22 &getGradientTimesInverseMetricRefTriang( ) const { return _GradGinvRefTriang; }
//   const Matrix32 &getGradientMapRefTriang( ) const { return _GradMapRefTriang; }
  const Matrix22 &getGradientMapInvRefTriang( ) const { return _GradMapInvRefTriang; }
  
  // get d = ei_x / |ei|^2, i = local node index
  inline RealType edgeXNormSqr ( const int locNode ) const { return _edgesRefTriang[locNode][0] / _edgesRefTriangLengthSqr[locNode];}
  
  // get a = ei_y / |ei|^2, i = local node index
  inline RealType edgeYNormSqr ( const int locNode ) const { return _edgesRefTriang[locNode][1] / _edgesRefTriangLengthSqr[locNode];}  
  
  // get b = ei_x*ei_y / |ei|^2, i = local node index
  inline RealType edgeXYNormSqr ( const int locNode ) const { return _edgesRefTriang[locNode][0]*_edgesRefTriang[locNode][1] / _edgesRefTriangLengthSqr[locNode];}
  
  // get c = ( ei_x*ei_x - 2*ei_y*ei_y) / |ei|^2, i = local node index
  inline RealType edgeXXMin2YYNormSqr ( const int locNode ) const { return (_edgesRefTriang[locNode][0] * _edgesRefTriang[locNode][0] - 2.*_edgesRefTriang[locNode][1]*_edgesRefTriang[locNode][1]) / _edgesRefTriangLengthSqr[locNode];}
  
  // get e = ( ei_y*ei_y - 2*ei_x*ei_x) / |ei|^2, i = local node index
  inline RealType edgeYYMin2XXNormSqr ( const int locNode ) const { return (_edgesRefTriang[locNode][1]*_edgesRefTriang[locNode][1] - 2.*_edgesRefTriang[locNode][0]*_edgesRefTriang[locNode][0]) / _edgesRefTriangLengthSqr[locNode];}
  
  // get l^2 = |ei|^2, i = local node index
  inline RealType edgeLengthSqr( int locNode ) const { return _edgesRefTriangLengthSqr[locNode];} 

  
  void setInverseMetricRefTriangle() {
      //! Dx 
      DomVecType dx1 = this->edgeRefTriang(0,1), dx2 = this->edgeRefTriang(0,2);
      _GradMapRefTriang.col(0) = dx1; _GradMapRefTriang.col(1) = dx2;
      //! Dx^-1
      _GradMapInvRefTriang = _GradMapRefTriang.inverse();
      //! g = Dx^T Dx, here compute ginv
      RealType g11 = dx1.dot(dx1);
      RealType g12 = dx1.dot(dx2);
      RealType g22 = dx2.dot(dx2);
      RealType detg = g11 * g22 - pesopt::Sqr( g12 );
      if ( detg < 1e-15 ){
          pesopt::consoleOutput("Error");
          print();
          throw std::invalid_argument( pesopt::strprintf ( "det is 0. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      }
      _ginvRefTriang(0,0) =  g22/detg;
      _ginvRefTriang(1,1) =  g11/detg;
      _ginvRefTriang(1,0) = -g12/detg;
      _ginvRefTriang(0,1) = -g12/detg;
      _areaRefTriang = 0.5 * sqrt( detg );
      //! Dx ginv
      _GradGinvRefTriang = _GradMapRefTriang * _ginvRefTriang;
      
#ifdef _TESTSHELLELEMENTMETRIC
      cout << "========================================" << endl;
      RealType areaCheck = 0.5 * detgCheck;
      cout << "area       = " << _areaRefTriang << endl;
      cout << "area check = " << areaCheck << endl;
      cout << "g^-1 = (Dx^T Dx)^-1 = " << endl << _ginvRefTriang << endl;
      cout << "Dx^-1 = " << endl << _GradMapInvRefTriang << endl;
      cout << "Dx = " << endl << _GradMapRefTriang << endl;
      cout << "Dx g^1 = " << endl << _GradGinvRefTriang << endl;
#endif
      
   }
    
};



#endif                                                      // __DKTPLATEELEMENT_H
