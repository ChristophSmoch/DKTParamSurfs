#ifndef __TRIANGLEELEMENT_H
#define __TRIANGLEELEMENT_H

#include <pesopt_IO.h>
    
    
//****************************************
//  triangle in R^2
//
//
//          n_2
//         / \
//        /   \ /\
//       /     \ \
//  e_1 /       \e_0
//  /  /         \
// \/ /           \
//   /             \
//  n_0 ___________ n_1
//          e_2 ->
//
//
//****************************************  
template<typename DataTypeContainer>
class TriangleElement {

public:
    
  static const int numNodes = 3; 
   
  typedef typename DataTypeContainer::RealType        RealType;
  typedef typename DataTypeContainer::RealVecChart    RealVecChart;
  typedef typename DataTypeContainer::PointType       PointType;
  typedef typename DataTypeContainer::TangentVecType  TangentVecType;
  typedef typename DataTypeContainer::Matrix22        Matrix22;
  typedef typename DataTypeContainer::Matrix32        Matrix32;
  typedef typename DataTypeContainer::Matrix33        Matrix33;
   
  typedef std::vector<PointType>                      VertexIterator;
  typedef std::vector<TangentVecType>                 TangentVecIterator;
  typedef Eigen::Vector3i                             ElementNodeIndicesType;
  
protected :   
    
    //****************************************
    //triangle in R^2
    //****************************************  
    // global indices of element and nodes
    int _globIdx;
    ElementNodeIndicesType _globNodeIdx;
    PointType _nodes[3];
    TangentVecType _edges[3]; // Notation:   e_i = N_{i-1} - N_{i+1} = N_{i+2} - N_{i+1}, i = 0,1,2
    Matrix22 _ginv; // (Dx^T Dx)^-1
    RealType _area;

    
public:
   
  TriangleElement () : _globIdx(-1) {}
  
  TriangleElement( const int globalIdx, const ElementNodeIndicesType globalNodeIndex, const VertexIterator &nodes  ) : 
     _globIdx(globalIdx),
     _globNodeIdx ( globalNodeIndex ) {
        this->update( nodes );
  }
  
  void update( const VertexIterator &nodes ){
      for ( int i = 0; i < 3; ++i ) _nodes[i] = nodes[ _globNodeIdx[i] ];
      for ( int i = 0; i < 3; ++i ) _edges[i] = getNode((i+2)%3) - getNode((i+1)%3);
      this->setInverseMetric();
  }
    
  // get and set functions
  int getGlobalElementIdx(  ) const { return _globIdx;}
  const ElementNodeIndicesType & getGlobalNodeIdx( ) const{ return _globNodeIdx;}
  int getGlobalNodeIdx(int localIndex) const{ return _globNodeIdx[localIndex];}
  void setGlobalNodeIdx(int localIndex, int globalIndex) { _globNodeIdx[localIndex] = globalIndex;}
  
  const PointType& getNode ( int i ) const { return _nodes[i];}
  PointType& getNode ( int i ) { return _nodes[i];}
  void setNode ( int i, const PointType& node ) {_nodes[i] = node;}
  
//   const TangentVecType& getEdge( int localEdgeIndex ) const { return _edges[localEdgeIndex];};
  const TangentVecType getEdge ( const int locNode1, const int locNode2 ) const { return (this->getNode(locNode2) - this->getNode(locNode1) );  }
  
  
  void setInverseMetric() {
      //! Dx 
      RealVecChart dx1 = this->getEdge(0,1), dx2 = this->getEdge(0,2);
//       RealVecChart dx1 = this->getEdge(2), dx2 = this->getEdge(1);
      //! g = Dx^T Dx
      RealType g11 = dx1.dot(dx1);
      RealType g12 = dx1.dot(dx2);
      RealType g22 = dx2.dot(dx2);
      RealType detg = g11 * g22 - pesopt::Sqr( g12 );
      if ( detg < 1e-15 ){
          pesopt::consoleOutput("Error");
          throw std::invalid_argument( pesopt::strprintf ( "det is 0. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      }
      //! ginv
      _ginv(0,0) =  g22/detg;
      _ginv(1,1) =  g11/detg;
      _ginv(1,0) = -g12/detg;
      _ginv(0,1) = -g12/detg;
      
      _area = 0.5 * sqrt( detg );
      
      
//       Matrix22 A; 
//       A.col(0) = this->getEdge(0,1); 
//       A.col(1) = this->getEdge(0,2);
//       _area = 0.5 * A.determinant();
//       _ginv = A.inverse(); 
      
      
//       cout << "================================" << endl;
//       cout << "EL = " << this->getGlobalElementIdx() << endl;
//       cout << "node indices = " << this->getGlobalNodeIdx() << endl;
// //       cout << "ginv = " << _ginv << endl;
// //       cout << "area = " << _area << endl;
// //       cout << "================================" << endl;
//       
//       //!TEST
//       Matrix22 A; 
//       A.col(0) = this->getEdge(0,1); 
//       A.col(1) = this->getEdge(0,2);
//       Matrix22 alternative = A.inverse(); 
//       cout << "Dx ginv = " << endl << A * _ginv << endl;
//       cout << "alternative = " << endl << alternative << endl;
//       cout << "A = " << endl << A << endl;
//       cout << "0.5 * det(A) = " << 0.5 * A.determinant() << endl;
      
  }
  const Matrix22 &getInverseMetric( ) const { return _ginv; }
//   const RealType getArea ( ) const {
//       TangentVecType dx1 = getNode(1) - getNode(0), dx2 = getNode(2) - getNode(0);
//       return 0.5 * std::abs( dx1(0) * dx2(1) - dx1(1) * dx2(0) ); 
//   }
  //for compatibility
  const RealType getVolume() const {return _area;}
  
  
  //LocalCoord wrt. unit triangle (0,0) - (1,0) - (0,1)
  void getGlobalCoord ( const PointType &LocalCoord, PointType &Coord ) const {
      RealVecChart dx1 = this->getEdge(0,1), dx2 = this->getEdge(0,2);
      Coord = LocalCoord[0] * dx1 + LocalCoord[1] * dx2;
      Coord += this->getNode(0); //Coord Left Down
  }
  
  //LocalCoord wrt. unit triangle (0,0) - (1,0) - (0,1)
  void getLocalCoord ( const PointType &Coord, PointType &LocalCoord ) const {
      PointType tmp = Coord - this->getNode(0); 
      
      Matrix22 A; 
      A.col(0) = this->getEdge(0,1); 
      A.col(1) = this->getEdge(0,2);
      
      LocalCoord = A.inverse() * tmp;
  }

};





// template<typename DataTypeContainer>
// class TriangleEdge {
// 
// public:
//   typedef typename DataTypeContainer::PointType       PointType;
//   typedef typename DataTypeContainer::TangentVecType  TangentVecType;
//   typedef std::vector<PointType>                      VertexIterator;
//   typedef typename Eigen::Vector2i                    EdgeNodeIndices;
//   
// public :   
//     int _globEdgeIdx;
//     EdgeNodeIndices _globNodeIndices;
//     PointType _nodes[2];
//     EdgeNodeIndices _globOppositeNodeIndices;
//     TangentVecType _edge;
//     EdgeNodeIndices _neighbouringTriangIndices;
//    
//   TriangleEdge () : _globEdgeIdx(-1) {}
//   
//   TriangleEdge( const int globalEdgeIdx, const EdgeNodeIndices globalNodeIndices, const VertexIterator &nodes  ) : 
//      _globEdgeIdx(globalEdgeIdx),
//      _globNodeIndices ( globalNodeIndices ) {
//         _nodes[0] = nodes[ _globNodeIndices[0] ];
//         _nodes[1] = nodes[ _globNodeIndices[1] ];
//         _edge = _nodes[1] - _nodes[0];
//         _globOppositeNodeIndices[0] = -1; _globOppositeNodeIndices[1] = -1; //-1 means boundary edge
//         _neighbouringTriangIndices[0] = -1; _neighbouringTriangIndices[1] = -1; //-1 means boundary edge
//     }
//     
//     const TangentVecType &getEdge( ) const {return _edge; }
//     
//     const int getAdjacentNodeIdx( int localIdx ) const { return _globNodeIndices[localIdx]; }
//     const PointType &getAdjacentNode( int localIdx ) const { return _nodes[localIdx]; }
//     
//     const int getOppositeNodeIdx ( int localIdx ) const { return _globOppositeNodeIndices[localIdx]; }
//     void setOppositeNodeIdx ( int localIdx, int globalNodeIdx ) { _globOppositeNodeIndices[localIdx] = globalNodeIdx; }
//               
//     const int getTriangIdx ( int localIdx ) const { return _neighbouringTriangIndices[localIdx]; }
//     void setTriangIdx ( int localIdx, int globalTriangIdx ) { _neighbouringTriangIndices[localIdx] = globalTriangIdx; }
// };
    



#endif
