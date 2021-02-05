#ifndef __QUOCFEMESHGENERATOR_H
#define __QUOCFEMESHGENERATOR_H

#include <feMeshGenerator.h>

template< >
class LineMesh<QuocMesh1D<>> {
protected:

  typedef QuocMesh1D<> QuocMeshType;
  typedef typename QuocMeshType::RealType                   RealType;
  typedef typename QuocMeshType::PointType                  PointType;
  typedef typename QuocMeshType::IntVecChart                IntVecChart;
  
public:

  LineMesh ( ) { }
 
  QuocMeshType generate( const int Nx, const RealType lx ) const {
      QuocMeshType mesh ( Nx, lx );
      return mesh;
  }
  
  QuocMeshType generate( const IntVecChart &numDofVec, const PointType &lengthVec ) const {
    return this->generate( numDofVec[0], lengthVec[0] );   
  }
  
};

template< >
class RectangleMesh<QuocMesh2D<>> {
protected:

  typedef QuocMesh2D<> QuocMeshType;
  typedef typename QuocMeshType::RealType                   RealType;
  typedef typename QuocMeshType::PointType                  PointType;
  typedef typename QuocMeshType::IntVecChart                IntVecChart;
  
public:

  RectangleMesh ( ) { }
 
  QuocMeshType generate( const int Nx, const int Ny, const RealType lx, const RealType ly ) const {
      QuocMeshType mesh ( Nx, Ny, lx, ly );
      return mesh;
  }
  
  QuocMeshType generate( const IntVecChart &numDofVec, const PointType &lengthVec ) const {
    return this->generate( numDofVec[0], numDofVec[1], lengthVec[0], lengthVec[1] );   
  }
  
};


template< >
class CuboidMesh<QuocMesh3D<>> {
protected:

  typedef QuocMesh3D<> QuocMeshType;
  typedef typename QuocMeshType::RealType                   RealType;
  typedef typename QuocMeshType::PointType                  PointType;
  typedef typename QuocMeshType::IntVecChart                IntVecChart;
  
public:

  CuboidMesh ( ) { }
 
  QuocMeshType generate( const int Nx, const int Ny, const int Nz, const RealType lx, const RealType ly, const RealType lz ) const {
      QuocMeshType mesh ( Nx, Ny, Nz, lx, ly, lz );
      return mesh;
  }
  
  QuocMeshType generate( const IntVecChart &numDofVec, const PointType &lengthVec ) const {
    return this->generate( numDofVec[0], numDofVec[1], numDofVec[2], lengthVec[0], lengthVec[1], lengthVec[2] );   
  }
  
};


// //stereographic projection P: flatMesh \subset \R^2 \to S^1
// template< typename MeshType >
// void stereographicProj ( const MeshType &flatMesh, MeshType &mesh ) {
//     
//    typedef typename MeshType::RealType           RealType;
//    typedef typename MeshType::PointType        PointType;
//     
//     mesh = flatMesh;
//     const int numVertices = flatMesh.getNumVertices();
//     
//     for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
//         const PointType& coords ( flatMesh.getVertex(nodeIdx) );
//         const RealType normSqr = coords[0] * coords[0] + coords[1] * coords[1];
//         PointType coordsOnSphere; 
//         coordsOnSphere(0) = 2. * coords(0) / (normSqr + 1.);
//         coordsOnSphere(1) = 2. * coords(1) / (normSqr + 1.);
//         coordsOnSphere(2) = (1. - normSqr) / (normSqr + 1.);
//         mesh.setVertex( nodeIdx, coordsOnSphere );
//         //!Tangentspace 
// //         if( ConfiguratorType::_ShellFEType == C1Dofs ){
// //                 RealType fac = 2. / pesopt::Sqr( normSqr + 1. );
// //                 TangentVecType firstTangentVecAtNode;  
// //                 firstTangentVecAtNode(0) = fac * ( - coords(0) * coords(0) + coords(1) * coords(1) + 1. ); 
// //                 firstTangentVecAtNode(1) = fac * -2. * coords(0) * coords(1);
// //                 firstTangentVecAtNode(2) = -2. * coords(0);
// //                 TangentVecType secondTangentVecAtNode; 
// //                 secondTangentVecAtNode(0) = fac * -2. * coords(0) * coords(1);
// //                 secondTangentVecAtNode(1) = fac * ( coords(0) * coords(0) - coords(1) * coords(1) + 1. );
// //                 secondTangentVecAtNode(2) = -2. * coords(1);
// //                 for( int comp=0; comp<3; ++comp ){
// //                   _xA[ nodeIdx +     numVertices + _numGlobalDofs * comp ] = firstTangentVecAtNode  [comp];
// //                   _xA[ nodeIdx + 2 * numVertices + _numGlobalDofs * comp ] = secondTangentVecAtNode [comp];
// //                 }
// //             }
//         }
// }
//     
//     
// 
// 
// template< typename MeshType >
// class TriangMeshGenerator_Circle {
// protected:
// 
//   typedef typename MeshType::RealType           RealType;
//   typedef typename MeshType::PointType        PointType;
//   typedef typename MeshType::TangentVecType     TangentVecType;
//   typedef typename MeshType::ElementNodeIndicesType      ElementNodeIndicesType;
//   typedef typename MeshType::VectorType         VectorType;
//   typedef typename MeshType::MaskType           MaskType;
//   typedef typename MeshType::ElementType        ElementType;
//   
// public:
// 
//   TriangMeshGenerator_Circle ( ) { }
// 
//   
//   MeshType generate( const int numDofs, const RealType radius ) const {
//       
//       const RealType pi = 4 * atan ( 1.0 );
//         
//       std::vector<PointType> vertices; vertices.reserve( numDofs + 1 );
//       vertices.push_back( PointType( 0., 0., 0. ) );
//       for( int nodeIdx = 0; nodeIdx < numDofs; ++nodeIdx ){
//          vertices.push_back( PointType( radius * cos( 2. * pi * static_cast<RealType>( nodeIdx ) / static_cast<RealType>( numDofs ) ), 
//                                           radius * sin( 2. * pi * static_cast<RealType>( nodeIdx ) / static_cast<RealType>( numDofs ) ),
//                                           0. ) );
//       }
//       
//       std::vector<ElementNodeIndicesType> triangles; triangles.reserve( numDofs );
//       for( int nodeIdx = 0; nodeIdx < numDofs; ++nodeIdx ){
//           if( nodeIdx < numDofs - 1 )
//              triangles.push_back( ElementNodeIndicesType( 0, nodeIdx + 1, nodeIdx + 2 ) );
//           if( nodeIdx == numDofs - 1 )
//              triangles.push_back( ElementNodeIndicesType( 0, nodeIdx + 1, 1 ) );
//       }
//       MeshType mesh ( vertices, triangles );
//       mesh.makeOrientationConsistent();
//       return mesh;
//   }
//   
// };


#endif
