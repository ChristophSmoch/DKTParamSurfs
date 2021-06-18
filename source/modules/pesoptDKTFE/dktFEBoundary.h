#ifndef __DKTFEBOUNDARY_H
#define __DKTFEBOUNDARY_H

#include <dktFEConfigurators.h>


enum ShellBoundaryType {
      ALLBOUNDARY = 100,
      NOBOUNDARY = 0,
  //Plate [0,1]^2
      PlateLeft = 1,
      PlateLeftTop = 2,
      PlateAll = 3,
      PlateHalfLeft = 4,
      PlateLeftAndRight = 5,
      PlateTopAndBottom = 6,
      PlateEpsilonLeft = 7,
      PlateRight = 8,
      PlateRightBottom = 9,
      PlateMiddleRightL = 1111,
  //Cylinder
      CylinderTopBottom = 11,
      CylinderBottom = 12,
  //Sphere
      SphereTopBottom = 21,
      SphereEquator = 22,
  //Circle
      CircleOutside = 31,
      CircleLeftAndRight = 32,
      CircleLeft = 33,
      CircleMidpoint = 34,
      CircleMidpointOneRing = 35,
  //PipeConnection
      PipeConnectionAll = 51,
  //Animals
      DragonTail = 211,
      DragonHead = 212
};


template< typename MeshType, bool ShellFETypeC1Dofs >
void generateDirichletBoundaryMaskUponShellBoundaryType ( const ShellBoundaryType ShellType, const MeshType &mesh,
                                                          typename MeshType::MaskType & mask, int & numBoundaryNodes,
                                                          const bool clampedBoundaryCondition = false ) {

//   typedef typename ConfiguratorType::InitType MeshType;
  typedef typename MeshType::RealType RealType;
  typedef typename MeshType::Point3DType    Point3DType;
  typedef typename MeshType::MaskType MaskType;

  numBoundaryNodes = 0;
  const int numVertices = mesh.getNumVertices();

  switch ( ShellType ){

      case NOBOUNDARY : {
      } break;

      case ALLBOUNDARY : {
            mesh.makeNeighbour();
            for( int ElementIndex = 0; ElementIndex < mesh.getNumTriangs(); ++ElementIndex ){
                for( int i = 0; i < 3; i++ ){
                    if ( mesh.getNeighbour( ElementIndex , i ) == -1 ){
                        mask[ mesh.getTriangNodeIdx( ElementIndex, (i+1)%3) ] = true;
                        mask[ mesh.getTriangNodeIdx( ElementIndex, (i+2)%3) ] = true;
                    }
                }
            }
            for( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ){
               if( mask[nodeIdx] ) ++numBoundaryNodes;
            }
      } break;


      //--------------------------------------------------------------------------------
      // Plates
      //---------------------------------------------------------------------------------
      case PlateLeft : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [0] == 0. ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;

      case PlateLeftTop : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [0] == 0. || coords[1] == 1. ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;

     case PlateAll : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [0] == 0. || coords [0] == 1. || coords [1] == 0. || coords[1] == 1. ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;

      case PlateHalfLeft : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if ( (coords[0] == 0.) && ( coords[1] > 0.5 ) ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;

      case PlateLeftAndRight : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if ( (coords [0] == 0.) || (coords[0] == 1.) ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;

     case PlateTopAndBottom : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if ( (coords [1] == 0.) || (coords[1] == 1.) ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;

      case PlateEpsilonLeft : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if ( (coords[0] == 0.) && ( coords[1] > 0.4 ) && ( coords[1] < 0.6 ) ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;

      case PlateRight : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [0] == 1. ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;

      case PlateRightBottom : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [0] == 1. || coords[1] == 0. ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;

      case PlateMiddleRightL : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [0] > 0.5 && coords[1] == 0.5 ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;


      //--------------------------------------------------------------------------------
      // Cylinder
      //---------------------------------------------------------------------------------
      case CylinderTopBottom : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [2] == 0. || coords [2] == 1.){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;

      case CylinderBottom : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if (coords [2] == 0. ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;

      //--------------------------------------------------------------------------------
      // Sphere
      //---------------------------------------------------------------------------------
      case SphereTopBottom : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if ( std::abs( coords[2] ) >= 1. - 1.e-15 ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;

      case SphereEquator : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if ( std::abs( coords[2] ) <= 1.e-15 ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;



      //--------------------------------------------------------------------------------
      // Circle
      //---------------------------------------------------------------------------------
      case CircleOutside : {
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              const RealType normSqr = coords(0) * coords(0) + coords(1) * coords(1);
              if ( normSqr >= 0.9 ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;


      case CircleLeftAndRight : {
           // first get topological boundary
            mesh.makeNeighbour();
            MaskType topologicalBdrMask ( mask.size() ); for( int i=0; i < mask.size(); ++i ) topologicalBdrMask[i] = false;
            for( int ElementIndex = 0; ElementIndex < mesh.getNumTriangs(); ++ElementIndex ){
                for( int i = 0; i < 3; i++ ){
                    if ( mesh.getNeighbour( ElementIndex , i ) == -1 ){
                        topologicalBdrMask[ mesh.getTriangNodeIdx( ElementIndex, (i+1)%3) ] = true;
                        topologicalBdrMask[ mesh.getTriangNodeIdx( ElementIndex, (i+2)%3) ] = true;
                    }
                }
            }
           // if coord < -0.9 or coord > 0.9 and on topological bdr
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if ( topologicalBdrMask[nodeIdx] && ( (coords[0] > 0.9) || (coords[0] < -0.9) ) ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;


     case CircleLeft : {
           // first get topological boundary
            mesh.makeNeighbour();
            MaskType topologicalBdrMask ( mask.size() ); for( int i=0; i < mask.size(); ++i ) topologicalBdrMask[i] = false;
            for( int ElementIndex = 0; ElementIndex < mesh.getNumTriangs(); ++ElementIndex ){
                for( int i = 0; i < 3; i++ ){
                    if ( mesh.getNeighbour( ElementIndex , i ) == -1 ){
                        topologicalBdrMask[ mesh.getTriangNodeIdx( ElementIndex, (i+1)%3) ] = true;
                        topologicalBdrMask[ mesh.getTriangNodeIdx( ElementIndex, (i+2)%3) ] = true;
                    }
                }
            }
           // if coord < -0.9 or coord > 0.9 and on topological bdr
            for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
              if ( topologicalBdrMask[nodeIdx] &&  (coords[0] < -0.9) ){
                ++numBoundaryNodes;
                mask[nodeIdx] = true;
              }
            }
      } break;


    case CircleMidpoint : {
        Point3DType midpoint ( 0.,0.,0.);
        int midpointVertexIndex;
        mesh.findClosestVertex( midpoint, midpointVertexIndex );
        mask[midpointVertexIndex] = true;
        ++numBoundaryNodes;
    } break;

    case CircleMidpointOneRing : {
        Point3DType midpoint ( 0.,0.,0.);
        int midpointVertexIndex;
        mesh.findClosestVertex( midpoint, midpointVertexIndex );
        mask[midpointVertexIndex] = true;
        ++numBoundaryNodes;

        std::vector<int> commonElements;
        mesh.getCommonElements( midpointVertexIndex, commonElements );
        for( int elIdx=0; elIdx < commonElements.size(); ++ elIdx ){
            for( int i=0; i<3; ++i ){
                int globalNodeIdx = mesh.getTriangNodeIdx( commonElements[elIdx], i );
                if( mask[globalNodeIdx] == false ){
                    mask[globalNodeIdx] = true;
                    ++numBoundaryNodes;
                }
            }
        }
    } break;

       //--------------------------------------------------------------------------------
      // other
      //---------------------------------------------------------------------------------
//      case PipeConnectionAll : {
//          for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
//              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
//              if (coords [2] == 0.5 || coords [2] == -0.5){
//                  ++numBoundaryNodes;
//                  mask[nodeIdx] = true;
//              }
//              if (coords [1] == 0.5 || coords [1] == -0.5){
//                  ++numBoundaryNodes;
//                  mask[nodeIdx] = true;
//              }
//          }
//       } break;
//
//       case DragonTail : {
//          for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
//              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
//              if (coords [0] <= 0.2){
//                  ++numBoundaryNodes;
//                  mask[nodeIdx] = true;
//              }
//          }
//       } break;
//
//
//        case DragonHead : {
//          for ( int nodeIdx=0; nodeIdx < numVertices; ++nodeIdx ) {
//              const Point3DType& coords ( mesh.getVertex(nodeIdx) );
//              if ( coords [0] >= 0.8){
//                  ++numBoundaryNodes;
//                  mask[nodeIdx] = true;
//              }
//          }
//       } break;

      default :
        throw std::invalid_argument( pesopt::strprintf ( "Wrong boundary condition. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        break;
    }

    //--------------------------------------------------------------------------------
    // clamped boundary condition
    //---------------------------------------------------------------------------------
    if( ShellFETypeC1Dofs &&  clampedBoundaryCondition  ){
        for( int i=0; i<numVertices; ++i ){
          if( mask[i] ){
            mask[ i + numVertices] = true;
            mask[ i + 2. * numVertices] = true;
          }
        }
    }

    //if non c1 dofs fix neighbouring nodes of all Dirichlet Boundary nodes
    if( (!ShellFETypeC1Dofs) && clampedBoundaryCondition  ){
        MaskType DirichletMask ( mask.size() ); DirichletMask = mask;
        for( int nodeIdx=0; nodeIdx<numVertices; ++nodeIdx ){
          if( DirichletMask[nodeIdx] ){
            std::vector<int> commonElements;
            mesh.getCommonElements( nodeIdx, commonElements );
            for( int i=0; i<commonElements.size(); ++i ){
                int elementIdx = commonElements[i];
                for( int localIndex = 0; localIndex <3; ++ localIndex ) mask[mesh.getTriang( elementIdx ).getGlobalNodeIdx(localIndex)] = true;
            }
          }
        }
    }


}


#endif // __DKTFEBOUNDARY_H
