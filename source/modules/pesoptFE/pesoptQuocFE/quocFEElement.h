#ifndef __QUOCFEELEMENT_H
#define __QUOCFEELEMENT_H

#include <pesopt_IO.h>
#include "quocFEDefines.h"
    
    
//!
//  n_0 ________________ n_1
//
template<typename DataTypeContainer>
class QuocElement1D {

public : 
    
  static const int numNodes = 2;
    
  typedef typename DataTypeContainer::RealType          RealType;
  typedef typename DataTypeContainer::RealVecChart      RealVecChart;
  typedef typename DataTypeContainer::PointType         PointType;
  typedef std::vector<PointType>                        VertexIterator;
  typedef Eigen::Vector2i                               ElementNodeIndicesType;
 
protected:
    int _globIdx;
    ElementNodeIndicesType _globNodeIdx;
    PointType _nodes[2];
    RealType _hx;
    RealType _volume;
    
public:
  
  QuocElement1D () : _globIdx(-1) {}
  
  QuocElement1D( const int globalIdx, 
                 const ElementNodeIndicesType globalNodeIndex, 
                 const VertexIterator &nodes   ) : 
     _globIdx(globalIdx),
     _globNodeIdx ( globalNodeIndex ) {
        for ( int i = 0; i < 2; ++i ) _nodes[i] = nodes[ _globNodeIdx[i] ];
        PointType diff = _nodes[1] - _nodes[0];
        _hx = diff[0];
        _volume = _hx;
    }
    
  //
  const RealType getVolume() const {return _volume;}
    
  // get and set functions
  int getGlobalElementIdx(  ) const { return _globIdx;}
  const ElementNodeIndicesType & getGlobalNodeIdx( ) const{ return _globNodeIdx;}
  int getGlobalNodeIdx(int localIndex) const{ return _globNodeIdx[localIndex];}
  void setGlobalNodeIdx(int localIndex, int globalIndex) { _globNodeIdx[localIndex] = globalIndex;}
  
  const PointType& getNode ( int i ) const { return _nodes[i];}
  PointType& getNode ( int i ) { return _nodes[i];}
  void setNode ( int i, const PointType& node ) {_nodes[i] = node;}

  
  //LocalCoord wrt. [0,1]
  void getGlobalCoord ( const PointType &LocalCoord, PointType &Coord ) const {
      Coord = LocalCoord; 
      Coord *= _hx;
      Coord += this->getNode(0); //Coord Left Down
  }
  
  //LocalCoord wrt. [0,1]
  void getLocalCoord ( const PointType &Coord, PointType &LocalCoord ) const {
      LocalCoord = Coord - this->getNode(0); 
      LocalCoord /= _hx;
  }
  
};




// template<typename DataTypeContainer>
// class QuocBoundaryElement1D {
// public : 

// static const int numNodes = 1;

//   typedef typename DataTypeContainer::RealType              RealType;
//   typedef typename DataTypeContainer::RealVecChart            RealVecChart;
//   typedef typename DataTypeContainer::RealVecChartBoundary    RealVecChartBoundary;
//   typedef typename DataTypeContainer::PointType             PointType;
//   typedef typename ElementType::ElementNodeIndicesType      ElementNodeIndicesType;
//   typedef std::vector<PointType>                            VertexIterator;
//   typedef QuocElement1D<DataTypeContainer>                  ElementType;
// 
//   protected :   
//     // global indices of element and nodes
//     int _globBoundaryIdx;
//     const ElementType _element; //! \todo something wrong for const ElementType & _element;
//     const QuocBoundaryType _boundaryType;
//     RealVecChart _normal;
//     PointType _boundaryNodes[1];
//     PointType _boundaryNodesRefCoord[1];
//     const RealType _h;
//     int _nodeIndicesOf1DElement[1];
//     
// public:
//   
//   QuocBoundaryElement1D () : _globBoundaryIdx(-1) {}
//   
//   QuocBoundaryElement1D( const int globalBoundaryIdx, const ElementType & element, const QuocBoundaryType &boundaryType, const RealType h  ) : 
//      _globBoundaryIdx(globalBoundaryIdx), _element (element), _boundaryType ( boundaryType ), _h(h) {
//          
//          _normal.setZero();
//          switch( boundaryType ){
//              case LEFT:{
//                _boundaryNodesRefCoord[0] << 0.0;
//                _boundaryNodes[0] = element.getNode(0);
//                _normal[0] = -1.;
//                _nodeIndicesOf1DElement[0] = 0;
//              }break;
//                  
//              case RIGHT:{
//                _boundaryNodesRefCoord[0] << 1.0;
//                _boundaryNodes[0] = element.getNode(1);
//                _normal[0] = 1.;
//                _nodeIndicesOf1DElement[0] = 1;
//              }break;
//                  
//              default:
//                  break;
//          }
//     }
//   
//   ~QuocBoundaryElement1D(){}
// 
//   const ElementType & getElement() const {return _element;}
//   int getNodeIndexOfElement(const int index ) const{ return _nodeIndicesOf1DElement[index];}
//   const RealVecChart& getNormal() const {return _normal;}
//   const PointType& getBoundaryNode ( int i ) const { return _boundaryNodes[i];}
//   PointType& getBoundaryNode ( int i ) { return _boundaryNodes[i];}
//   const PointType& getBoundaryNodeRefCoord ( int i ) const { return _boundaryNodesRefCoord[i];}
//   PointType& getBoundaryNodeRefCoord ( int i ) { return _boundaryNodesRefCoord[i];}
//   const RealType getVolumeOfBoundaryElement() const {return _h;}
//   const QuocBoundaryType getBoundaryType() const {return _boundaryType;}
//   void getRefCoord( const RealVecChartBoundary &refCoordBoundary, RealVecChart & refCoord ) const {
//     refCoord = ( 1. - refCoordBoundary ) * this->getBoundaryNodeRefCoord(0) + refCoordBoundary * this->getBoundaryNodeRefCoord(1);   
//   }
// 
// };
    
    
    
    
    
//!
//     __________________
//  n_2                  n_3
//   |                    |
//   |                    |
//   |                    |
//   |                    |
//  n_0 ________________ n_1
//
template<typename DataTypeContainer>
class QuocElement2D {
public : 
    
  static const int numNodes = 4;
    
  typedef typename DataTypeContainer::RealType        RealType;
  typedef typename DataTypeContainer::RealVecChart    RealVecChart;
  typedef typename DataTypeContainer::PointType       PointType;
  typedef std::vector<PointType>                      VertexIterator;
  typedef Eigen::Vector4i                             ElementNodeIndicesType;
  
protected:
  int _globIdx;
  ElementNodeIndicesType _globNodeIdx;
  PointType _nodes[4];
  RealType _hx,_hy;
  RealType _volume;
    
public:
  
  QuocElement2D () : _globIdx(-1) {}
  
  QuocElement2D( const int globalIdx, 
                 const ElementNodeIndicesType globalNodeIndex, 
                 const VertexIterator &nodes  ) : 
     _globIdx(globalIdx),
     _globNodeIdx ( globalNodeIndex ) {
        for ( int i = 0; i < 4; ++i ) _nodes[i] = nodes[ _globNodeIdx[i] ];
        PointType diff = _nodes[3] - _nodes[0];
        _hx = diff[0]; _hy = diff[1];
        _volume = _hx * _hy;
    }
    
  //
  const RealType getVolume() const {return _volume;}
  
  // get and set functions
  int getGlobalElementIdx(  ) const { return _globIdx;}
  const ElementNodeIndicesType & getGlobalNodeIdx( ) const{ return _globNodeIdx;}
  int getGlobalNodeIdx(int localIndex) const{ return _globNodeIdx[localIndex];}
  void setGlobalNodeIdx(int localIndex, int globalIndex) { _globNodeIdx[localIndex] = globalIndex;}
  
  const PointType& getNode ( int i ) const { return _nodes[i];}
  PointType& getNode ( int i ) { return _nodes[i];}
  void setNode ( int i, const PointType& node ) {_nodes[i] = node;}
  
  //LocalCoord wrt. [0,1]^2
  void getGlobalCoord ( const PointType &LocalCoord, PointType &Coord ) const {
      Coord = LocalCoord; 
      Coord[0] *= _hx; Coord[1] *= _hy;
      Coord += this->getNode(0); //Coord Left Down
  }

  //LocalCoord wrt. [0,1]^2
  void getLocalCoord ( const PointType &Coord, PointType &LocalCoord ) const {
      LocalCoord = Coord - this->getNode(0); 
      LocalCoord[0] /= _hx; LocalCoord[1] /= _hy;
  }
  
};




template<typename DataTypeContainer>
class QuocBoundaryElement2D {

public : 
    
  static const int numNodes = 2;
    
  typedef typename DataTypeContainer::RealType              RealType;
  typedef typename DataTypeContainer::RealVecChart          RealVecChart;
  typedef typename DataTypeContainer::RealVecChartBoundary  RealVecChartBoundary;
  typedef typename DataTypeContainer::PointType             PointType;
  typedef std::vector<PointType>                            VertexIterator;
  typedef QuocElement2D<DataTypeContainer>                  ElementType;
  typedef typename ElementType::ElementNodeIndicesType      ElementNodeIndicesType;
  
  protected :   
    // global indices of element and nodes
    int _globBoundaryIdx;
    //! \todo something wrong for const ElementType & _element;
    //     const ElementType _element; 
    const ElementType &_element;
    const QuocBoundaryType _boundaryType;
    RealVecChart _normal;
    PointType _boundaryNodes[2];
    PointType _boundaryNodesRefCoord[2];
    const RealType _h;
    int _nodeIndicesOf2DElement[2];
    
public:
  
  QuocBoundaryElement2D () : _globBoundaryIdx(-1) {}
  
  QuocBoundaryElement2D( const int globalBoundaryIdx, 
                         const ElementType & element, 
                         const QuocBoundaryType &boundaryType, 
                         const RealType h  ) : 
     _globBoundaryIdx(globalBoundaryIdx), _element (element), _boundaryType ( boundaryType ), _h(h) {
         
         _normal.setZero();
         switch( boundaryType ){
             case LEFT:{
               _boundaryNodesRefCoord[0] << 0.0, 0.0; _boundaryNodesRefCoord[1] << 0.0, 1.0;
               _boundaryNodes[0] = element.getNode(0);_boundaryNodes[1] = element.getNode(2);
               _normal[0] = -1.;
               _nodeIndicesOf2DElement[0] = 0; _nodeIndicesOf2DElement[1] = 2;
             }break;
                 
             case RIGHT:{
               _boundaryNodesRefCoord[0] << 1.0, 0.0; _boundaryNodesRefCoord[1] << 1.0, 1.0;
               _boundaryNodes[0] = element.getNode(1);_boundaryNodes[1] = element.getNode(3);
               _normal[0] = 1.;
               _nodeIndicesOf2DElement[0] = 1; _nodeIndicesOf2DElement[1] = 3;
             }break;
                 
             case BOTTOM:{
               _boundaryNodesRefCoord[0] << 0.0, 0.0; _boundaryNodesRefCoord[1] << 1.0, 0.0;
               _boundaryNodes[0] = element.getNode(0); _boundaryNodes[1] = element.getNode(1);
               _normal[1] = -1.;
               _nodeIndicesOf2DElement[0] = 0; _nodeIndicesOf2DElement[1] = 1;
             }break;
                 
             case TOP:{
               _boundaryNodesRefCoord[0] << 0.0, 1.0; _boundaryNodesRefCoord[1] << 1.0, 1.0;
               _boundaryNodes[0] = element.getNode(2); _boundaryNodes[1] = element.getNode(3);
               _normal[1] = 1.;
               _nodeIndicesOf2DElement[0] = 2; _nodeIndicesOf2DElement[1] = 3;
             }break;
                 
             default:
                 break;
         }
    }

  const ElementType & getElement() const {return _element;}
  int getNodeIndexOfElement(const int index ) const{ return _nodeIndicesOf2DElement[index];}
  const RealVecChart& getNormal() const {return _normal;}
  const PointType& getBoundaryNode ( int i ) const { return _boundaryNodes[i];}
  PointType& getBoundaryNode ( int i ) { return _boundaryNodes[i];}
  const PointType& getBoundaryNodeRefCoord ( int i ) const { return _boundaryNodesRefCoord[i];}
  PointType& getBoundaryNodeRefCoord ( int i ) { return _boundaryNodesRefCoord[i];}
  const RealType getVolumeOfBoundaryElement() const {return _h;}
  const QuocBoundaryType getBoundaryType() const {return _boundaryType;}
  void getRefCoord( const RealVecChartBoundary &refCoordBoundary, RealVecChart & refCoord ) const {
    refCoord = ( 1. - refCoordBoundary ) * this->getBoundaryNodeRefCoord(0) + refCoordBoundary * this->getBoundaryNodeRefCoord(1);   
  }

};

// n0 (0,0,0)
// n1 (1,0,0)
// n2 (0,1,0)
// n3 (1,1,0)
// n4 (0,0,1)
// n5 (1,0,1)
// n6 (0,1,1)
// n7 (1,1,1)
//            __________________
//         n_6                  n_7
//         /                     /|
//        /|                    / |
//       / |                   /  |
//      /  | _________________/___|_
//     /   n_2               /   n_3
//    /____/______________  /     /
//  n_4   /               n_5    /
//   |   /                 |    /
//   |  /                  |   /
//   | /                   |  /
//   |/                    | /
//  n_0 ________________ n_1/
//
template<typename DataTypeContainer>
class QuocElement3D {

public : 
    
  static const int numNodes = 8;
    
  typedef typename DataTypeContainer::RealType                      RealType;
  typedef typename DataTypeContainer::RealVecChart                  RealVecChart;
  typedef typename DataTypeContainer::PointType                     PointType;
  typedef typename DataTypeContainer::MaskType                      MaskType;
  typedef std::vector<PointType>                                    VertexIterator;
  typedef Eigen::Matrix< int, 8, 1 >                                ElementNodeIndicesType;

  
protected:
    int _globIdx;
    ElementNodeIndicesType _globNodeIdx;
    PointType _nodes[8];
    RealType _hx,_hy,_hz;
    RealType _volume;
    
public:
  
  QuocElement3D () : _globIdx(-1) {}
  
  QuocElement3D( const int globalIdx, 
                 const ElementNodeIndicesType globalNodeIndex, 
                 const VertexIterator &nodes ) : 
     _globIdx(globalIdx),
     _globNodeIdx ( globalNodeIndex )
    {
        for ( int i = 0; i < 8; ++i ) _nodes[i] = nodes[ _globNodeIdx[i] ];
        PointType diff = _nodes[7] - _nodes[0];
        _hx = diff[0]; _hy = diff[1]; _hz = diff[2];
        _volume = _hx * _hy * _hz;
    }
  
  //
  const RealType getVolume( )  const { return _volume; }

  // get and set functions
  int getGlobalElementIdx(  ) const { return _globIdx;}
  const ElementNodeIndicesType & getGlobalNodeIdx( ) const{ return _globNodeIdx;}
  int getGlobalNodeIdx(int localIndex) const{ return _globNodeIdx[localIndex];}
  void setGlobalNodeIdx(int localIndex, int globalIndex) { _globNodeIdx[localIndex] = globalIndex;}
  
  const PointType& getNode ( int i ) const { return _nodes[i];}
  PointType& getNode ( int i ) { return _nodes[i];}
  void setNode ( int i, const PointType& node ) {_nodes[i] = node;}
  
  //LocalCoord wrt. [0,1]^3
  void getGlobalCoord ( const PointType &LocalCoord, PointType &Coord ) const {
      Coord = LocalCoord; 
      Coord[0] *= _hx; Coord[1] *= _hy; Coord[2] *= _hz;
      Coord += this->getNode(0); //Coord Left Down
  }
  
  //LocalCoord wrt. [0,1]^3
  void getLocalCoord ( const PointType &Coord, PointType &LocalCoord ) const {
      LocalCoord = Coord - this->getNode(0); 
      LocalCoord[0] /= _hx; LocalCoord[1] /= _hy; LocalCoord[2] /= _hz;
  }
  
};


template<typename DataTypeContainer>
class QuocBoundaryElement3D{
public : 
    
  static const int numNodes = 4;
    
  typedef typename DataTypeContainer::RealType                  RealType;
  typedef typename DataTypeContainer::RealVecChart              RealVecChart;
  typedef typename DataTypeContainer::RealVecChartBoundary      RealVecChartBoundary;
  typedef typename DataTypeContainer::PointType                 PointType;
  typedef std::vector<PointType>                                VertexIterator;
  typedef QuocElement3D<DataTypeContainer>                      ElementType;
  typedef typename ElementType::ElementNodeIndicesType          ElementNodeIndicesType;

  protected :   
    int _globBoundaryIdx;
    //! \todo something wrong for const ElementType & _element;
    // const ElementType _element; 
    const ElementType & _element;
    const QuocBoundaryType _boundaryType;
    RealVecChart _normal;
    PointType _boundaryNodes[4];
    PointType _boundaryNodesRefCoord[4];
    const RealType _h;
    int _nodeIndicesOf3DElement[4];
    
public:
  
  QuocBoundaryElement3D () : _globBoundaryIdx(-1) {}
  
  QuocBoundaryElement3D( const int globalBoundaryIdx, const ElementType &element, const QuocBoundaryType &boundaryType, const RealType h  ) : 
     _globBoundaryIdx(globalBoundaryIdx), _element (element), _boundaryType ( boundaryType ), _h(h) {
         
         _normal.setZero();
         switch( boundaryType ){
             case LEFT:{
               _nodeIndicesOf3DElement[0] = 0; _nodeIndicesOf3DElement[1] = 2; _nodeIndicesOf3DElement[2] = 4; _nodeIndicesOf3DElement[3] = 6;
               _boundaryNodesRefCoord[0] << 0.0, 0.0, 0.0; 
               _boundaryNodesRefCoord[1] << 0.0, 1.0, 0.0;
               _boundaryNodesRefCoord[2] << 0.0, 0.0, 1.0; 
               _boundaryNodesRefCoord[3] << 0.0, 1.0, 1.0;
               _boundaryNodes[0] = element.getNode(0);
               _boundaryNodes[1] = element.getNode(2);
               _boundaryNodes[2] = element.getNode(4);
               _boundaryNodes[3] = element.getNode(6);
               _normal[0] = -1.;
             }break;
                 
             case RIGHT:{
               _nodeIndicesOf3DElement[0] = 1; _nodeIndicesOf3DElement[1] = 3; _nodeIndicesOf3DElement[2] = 5; _nodeIndicesOf3DElement[3] = 7;
               _boundaryNodesRefCoord[0] << 1.0, 0.0, 0.0;
               _boundaryNodesRefCoord[1] << 1.0, 1.0, 0.0;
               _boundaryNodesRefCoord[2] << 1.0, 0.0, 1.0; 
               _boundaryNodesRefCoord[3] << 1.0, 1.0, 1.0;
               _boundaryNodes[0] = element.getNode(1);
               _boundaryNodes[1] = element.getNode(3);
               _boundaryNodes[2] = element.getNode(5);
               _boundaryNodes[3] = element.getNode(7);
               _normal[0] = 1.;
               
             }break;
                 
             case BOTTOM:{
               _nodeIndicesOf3DElement[0] = 0; _nodeIndicesOf3DElement[1] = 1; _nodeIndicesOf3DElement[2] = 2; _nodeIndicesOf3DElement[3] = 3;
               _boundaryNodesRefCoord[0] << 0.0, 0.0, 0.0;
               _boundaryNodesRefCoord[1] << 1.0, 0.0, 0.0;
               _boundaryNodesRefCoord[2] << 0.0, 1.0, 0.0; 
               _boundaryNodesRefCoord[3] << 1.0, 1.0, 0.0;
               _boundaryNodes[0] = element.getNode(0);
               _boundaryNodes[1] = element.getNode(1);
               _boundaryNodes[2] = element.getNode(2);
               _boundaryNodes[3] = element.getNode(3);
               _normal[2] = -1.;
             }break;
                 
             case TOP:{
               _nodeIndicesOf3DElement[0] = 4; _nodeIndicesOf3DElement[1] = 5; _nodeIndicesOf3DElement[2] = 6; _nodeIndicesOf3DElement[3] = 7;
               _boundaryNodesRefCoord[0] << 0.0, 0.0, 1.0;
               _boundaryNodesRefCoord[1] << 1.0, 0.0, 1.0;
               _boundaryNodesRefCoord[2] << 0.0, 1.0, 1.0; 
               _boundaryNodesRefCoord[3] << 1.0, 1.0, 1.0;
               _boundaryNodes[0] = element.getNode(4);
               _boundaryNodes[1] = element.getNode(5);
               _boundaryNodes[2] = element.getNode(6);
               _boundaryNodes[3] = element.getNode(7);
               _normal[2] = 1.;
             }break;
             
             case FRONT:{
               _nodeIndicesOf3DElement[0] = 0; _nodeIndicesOf3DElement[1] = 1; _nodeIndicesOf3DElement[2] = 4; _nodeIndicesOf3DElement[3] = 5;
               _boundaryNodesRefCoord[0] << 0.0, 0.0, 0.0;
               _boundaryNodesRefCoord[1] << 1.0, 0.0, 0.0;
               _boundaryNodesRefCoord[2] << 0.0, 0.0, 1.0; 
               _boundaryNodesRefCoord[3] << 1.0, 0.0, 1.0;
               _boundaryNodes[0] = element.getNode(0);
               _boundaryNodes[1] = element.getNode(1);
               _boundaryNodes[2] = element.getNode(4);
               _boundaryNodes[3] = element.getNode(5);
               _normal[1] = -1.;
             }break;
                 
             case BACK:{
               _nodeIndicesOf3DElement[0] = 2; _nodeIndicesOf3DElement[1] = 3; _nodeIndicesOf3DElement[2] = 6; _nodeIndicesOf3DElement[3] = 7;
               _boundaryNodesRefCoord[0] << 0.0, 1.0, 0.0;
               _boundaryNodesRefCoord[1] << 1.0, 1.0, 0.0;
               _boundaryNodesRefCoord[2] << 0.0, 1.0, 1.0; 
               _boundaryNodesRefCoord[3] << 1.0, 1.0, 1.0;
               _boundaryNodes[0] = element.getNode(2);
               _boundaryNodes[1] = element.getNode(3);
               _boundaryNodes[2] = element.getNode(6);
               _boundaryNodes[3] = element.getNode(7);
               _normal[1] = 1.;
             }break;
                 
             default:
                 break;
         }
    }
  
  const ElementType & getElement() const {return _element;}
  int getNodeIndexOfElement(const int index ) const{ return _nodeIndicesOf3DElement[index];}
  const RealVecChart& getNormal() const {return _normal;}
  const PointType& getBoundaryNode ( int i ) const { return _boundaryNodes[i];}
  PointType& getBoundaryNode ( int i ) { return _boundaryNodes[i];}
  const PointType& getBoundaryNodeRefCoord ( int i ) const { return _boundaryNodesRefCoord[i];}
  PointType& getBoundaryNodeRefCoord ( int i ) { return _boundaryNodesRefCoord[i];}
  const RealType getVolumeOfBoundaryElement() const {return _h;}
  const QuocBoundaryType getBoundaryType() const {return _boundaryType;}
  void getRefCoord( const RealVecChartBoundary &refCoordBoundary, RealVecChart & refCoord ) const {
    refCoord = (1. - refCoordBoundary[1]) * ( ( 1. - refCoordBoundary[0] ) * this->getBoundaryNodeRefCoord(0) 
               + refCoordBoundary[0] * this->getBoundaryNodeRefCoord(1) )
               + refCoordBoundary[1]  * ( ( 1. - refCoordBoundary[0] ) * this->getBoundaryNodeRefCoord(2) 
               + refCoordBoundary[0] * this->getBoundaryNodeRefCoord(3) ); 
  }

};


#endif
