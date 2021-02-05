#ifndef __TRIANGLEFEMESHGENERATOR_H
#define __TRIANGLEFEMESHGENERATOR_H

#include <feMeshGenerator.h>

template< >
class RectangleMesh<TriangleMesh<>> {
protected:

  typedef TriangleMesh<> TriangleMeshType;
  typedef typename TriangleMeshType::RealType                   RealType;
  typedef typename TriangleMeshType::PointType                  PointType;
  typedef typename TriangleMeshType::TangentVecType             TangentVecType;
  typedef typename TriangleMeshType::ElementNodeIndicesType     ElementNodeIndicesType;
  typedef typename TriangleMeshType::VectorType                 VectorType;
  typedef typename TriangleMeshType::MaskType                   MaskType;
  typedef typename TriangleMeshType::ElementType                ElementType;
  typedef typename TriangleMeshType::RealVecChart               RealVecChart;
  typedef typename TriangleMeshType::IntVecChart                IntVecChart;
public:

  RectangleMesh ( ) { }
 
  TriangleMeshType generate( const int Nx, const int Ny, const RealType lx, const RealType ly ) const {
      
      const RealType hx = lx / static_cast<RealType>( Nx - 1 );
      const RealType hy = ly / static_cast<RealType>( Ny - 1 );
      
      
      std::vector<PointType> vertices; vertices.reserve( Nx * Ny );
      for( int nodeIdx_y = 0; nodeIdx_y < Ny; ++nodeIdx_y )
          for( int nodeIdx_x = 0; nodeIdx_x < Nx; ++nodeIdx_x ){
              PointType node; node.setZero();
              node[0] = nodeIdx_x * hx;
              node[1] = nodeIdx_y * hy;
              vertices.push_back( node );
          }
      
      std::vector<ElementNodeIndicesType> triangles; triangles.reserve( 2 * (Nx-1) * (Ny-1) );
      int nodeIdx_x_y = -1;
      for( int nodeIdx_y= 0; nodeIdx_y < Ny - 1; ++nodeIdx_y )
          for( int nodeIdx_x = 0; nodeIdx_x < Nx; ++nodeIdx_x ){
              nodeIdx_x_y += 1;
              if( nodeIdx_x < Nx-1 ){
                  int nodeIdx_xP1_y = nodeIdx_x_y + 1;
                  int nodeIdx_x_yP1 = nodeIdx_x_y + Nx;
                  int nodeIdx_xP1_yP1 =  nodeIdx_xP1_y + Nx;
                  triangles.push_back( ElementNodeIndicesType( nodeIdx_x_y, nodeIdx_xP1_y, nodeIdx_xP1_yP1 ) );
                  triangles.push_back( ElementNodeIndicesType( nodeIdx_x_y, nodeIdx_x_yP1, nodeIdx_xP1_yP1 ) );
              }
          }
      
      TriangleMeshType mesh ( vertices, triangles );
      mesh.makeOrientationConsistent();
      return mesh;
  }
  
  TriangleMeshType generate( const IntVecChart &numDofVec, const RealVecChart &lengthVec ) const {
    return this->generate( numDofVec[0], numDofVec[1], lengthVec[0], lengthVec[1] );   
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
