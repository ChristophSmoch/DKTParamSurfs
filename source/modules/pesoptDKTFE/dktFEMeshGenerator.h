#ifndef __TRIANGSHELLMESHGENERATOR_H
#define __TRIANGSHELLMESHGENERATOR_H

#include <pesopt_IO.h>
#include <dktFEMeshInfo.h>
    
template< typename MeshType >
class TriangMeshGenerator {
protected:

  typedef typename MeshType::RealType           RealType;
  typedef typename MeshType::Point3DType        Point3DType;
  typedef typename MeshType::TangentVecType     TangentVecType;
  typedef typename MeshType::Indices3DType      Indices3DType;
  typedef typename MeshType::VectorType         VectorType;
  typedef typename MeshType::MaskType           MaskType;
  typedef typename MeshType::ElementType        ElementType;
  
public:

  TriangMeshGenerator ( ) { }
 
  MeshType generateRectangle( const int length_x, const int length_y ) const {
      
      const RealType meshsize = 1.;
      
      std::vector<Point3DType> vertices; vertices.reserve( (length_x + 1) * (length_y + 1) );
      for( int nodeIdx_y = 0; nodeIdx_y < length_y+1; ++nodeIdx_y )
          for( int nodeIdx_x = 0; nodeIdx_x < length_x+1; ++nodeIdx_x ){
              vertices.push_back( Point3DType( nodeIdx_x * meshsize, nodeIdx_y * meshsize, 0. ) );
          }
      
      std::vector<Indices3DType> triangles; triangles.reserve( 2 * length_x * length_y );
      int nodeIdx_x_y = -1;
      for( int nodeIdx_y= 0; nodeIdx_y < length_y; ++nodeIdx_y )
          for( int nodeIdx_x = 0; nodeIdx_x < length_x + 1; ++nodeIdx_x ){
              nodeIdx_x_y += 1;
              if( nodeIdx_x < length_x ){
                  int nodeIdx_xP1_y = nodeIdx_x_y + 1;
                  int nodeIdx_x_yP1 = nodeIdx_x_y + length_x + 1;
                  int nodeIdx_xP1_yP1 =  nodeIdx_xP1_y + length_x + 1;
                  triangles.push_back( Indices3DType( nodeIdx_x_y, nodeIdx_xP1_y, nodeIdx_xP1_yP1 ) );
                  triangles.push_back( Indices3DType( nodeIdx_x_y, nodeIdx_x_yP1, nodeIdx_xP1_yP1 ) );
              }
          }
      
      MeshType mesh ( vertices, triangles );
      mesh.makeOrientationConsistent();
      return mesh;
  }
  
  
    MeshType generateCircle( const int numDofs, const RealType radius ) const {
      
      const RealType pi = 4 * atan ( 1.0 );
        
      std::vector<Point3DType> vertices; vertices.reserve( numDofs + 1 );
      vertices.push_back( Point3DType( 0., 0., 0. ) );
      for( int nodeIdx = 0; nodeIdx < numDofs; ++nodeIdx ){
         vertices.push_back( Point3DType( radius * cos( 2. * pi * static_cast<RealType>( nodeIdx ) / static_cast<RealType>( numDofs ) ), 
                                          radius * sin( 2. * pi * static_cast<RealType>( nodeIdx ) / static_cast<RealType>( numDofs ) ),
                                          0. ) );
      }
      
      std::vector<Indices3DType> triangles; triangles.reserve( numDofs );
      for( int nodeIdx = 0; nodeIdx < numDofs; ++nodeIdx ){
          if( nodeIdx < numDofs - 1 )
             triangles.push_back( Indices3DType( 0, nodeIdx + 1, nodeIdx + 2 ) );
          if( nodeIdx == numDofs - 1 )
             triangles.push_back( Indices3DType( 0, nodeIdx + 1, 1 ) );
      }
      MeshType mesh ( vertices, triangles );
      mesh.makeOrientationConsistent();
      return mesh;
  }
  
};
    
    
    
template< typename MeshType >
class TriangMeshManipulator {
protected:

  typedef typename MeshType::RealType           RealType;
  typedef typename MeshType::Point3DType        Point3DType;
  typedef typename MeshType::TangentVecType     TangentVecType;
  typedef typename MeshType::Indices3DType      Indices3DType;
  typedef typename MeshType::VectorType         VectorType;
  typedef typename MeshType::MaskType           MaskType;
  
public:
  //! Create empty TriangMesh
  TriangMeshManipulator ( ) { }
 
  // (x,y,z) mapsto (x,y, cy + z)
  void shiftLinearlyInZDirection( MeshType &mesh, const RealType scalar ) const {
      for( int nodeIdx = 0; nodeIdx < mesh.getNumVertices(); ++ nodeIdx ){
          Point3DType newVertex = mesh.getVertex ( nodeIdx );
          newVertex[2] += scalar * newVertex[1];
          mesh.setVertex( nodeIdx, newVertex );
      }
      
      mesh.updateAllTriangles ( );
      mesh.updateAllProjectionCoefficients();
  }
  
  // (x,y,z) mapsto (x,y cos(phi) - z sin(phi), y sin(phi) + z cos(phi) )
  void rotateAroundX( MeshType &mesh, const RealType phi ) const {
      for( int nodeIdx = 0; nodeIdx < mesh.getNumVertices(); ++ nodeIdx ){
          const Point3DType oldVertex = mesh.getVertex ( nodeIdx );
          Point3DType newVertex = mesh.getVertex ( nodeIdx );
          newVertex[1] = oldVertex[1] * std::cos(phi) - oldVertex[2] * std::sin(phi);
          newVertex[2] = oldVertex[1] * std::sin(phi) + oldVertex[2] * std::cos(phi);
          mesh.setVertex( nodeIdx, newVertex );
      }
      
      mesh.updateAllTriangles ( );
      mesh.updateAllProjectionCoefficients();
  }
 
  
  void scale( MeshType &mesh, const Point3DType &newLowerbounds, const Point3DType &newUpperBounds ) const{
      DKTTriangMeshInfo<MeshType> meshInfo ( mesh );
      Point3DType lowerbounds, upperbounds;
      meshInfo.getBoundingBox( lowerbounds, upperbounds );
      Point3DType oldLength; oldLength = upperbounds - lowerbounds;
      Point3DType newLength; newLength = newUpperBounds - newLowerbounds;
      for( int nodeIdx = 0; nodeIdx < mesh.getNumVertices(); ++nodeIdx ){
        Point3DType oldVertex, newVertex;
        oldVertex = mesh.getVertex( nodeIdx );
        for( int coord=0; coord < 3; ++coord ){
             newVertex[coord] = newLength[coord] / oldLength[coord] * ( oldVertex[coord] - lowerbounds[coord] ) + newLowerbounds[coord];
        }
        mesh.setVertex( nodeIdx, newVertex );
      }
  }
  
  void scaleTo01( MeshType &mesh ) const{
   
      Point3DType newLowerbounds, newUpperBounds;
      for( int coord=0; coord<3; ++coord ){
         newLowerbounds[coord] = 0.;
         newUpperBounds[coord] = 1.;
      }
      this->scale( mesh, newLowerbounds, newUpperBounds );
  }
 
  void scaleTo01LengthPreserving( MeshType &mesh ) const{
   
      DKTTriangMeshInfo<MeshType> meshInfo ( mesh );
      Point3DType lowerbounds, upperbounds;
      meshInfo.getBoundingBox( lowerbounds, upperbounds );
      Point3DType oldLength; oldLength = upperbounds - lowerbounds;
      RealType maxLength = 0.; 
      for( int coord=0; coord<3; ++coord){
        if( oldLength[coord] > maxLength ) maxLength = oldLength[coord];   
      }
      
      Point3DType newLength;
      for( int coord=0; coord<3; ++coord) newLength[coord] = oldLength[coord] / maxLength;
      
      Point3DType newLowerbounds, newUpperBounds;
      for( int coord=0; coord<3; ++coord ){
         newLowerbounds[coord] = 0.5 - 0.5 * newLength[coord];
         newUpperBounds[coord] = 0.5 + 0.5 * newLength[coord];
      }
      
      this->scale( mesh, newLowerbounds, newUpperBounds );
  }
  
 //=========================================================================
 // Isometric deformations of plate
  
  // linear interpolation between init and final slope
  RealType getPhase_interpolateSlopeLinear( const RealType & t, const RealType initSlope, const RealType finalSlope ) const{
    const RealType init = atan( initSlope ), fin = atan ( finalSlope );
    return ( 1.0 - t ) * init + t * fin;
  }
  
  // computes int_0^x exp( i (K(t) + K(0) ) d t
  void getPoint_interpolateSlopeLinear( const RealType & t, RealType & gamma_x, RealType & gamma_z, const RealType initSlope, const RealType finalSlope ) const{
    
    gamma_x = 0.0;
    gamma_z = 0.0;
    
    int numSteps = 50; //TODO numQuadPoints
    RealType tau = t / static_cast<RealType> (numSteps);
    
    for(int i=0; i< numSteps; ++i){
      RealType coord_i = (static_cast<RealType> (i) + 0.5 ) * tau;
      RealType K_i = getPhase_interpolateSlopeLinear ( coord_i, initSlope, finalSlope );
      gamma_x += tau * cos( K_i );
      gamma_z += tau * sin( K_i );
    }
  }
  
  void IsometricDeformation_linearInterpolatedSlope( MeshType &mesh, const RealType initSlope, const RealType finalSlope  ) const {
//   void constructSimpleDisplacement( pesopt::MultiVector<RealType> &w, const RealType initSlope, const RealType finalSlope ) const{
    
    const int numVertices = mesh.getNumVertices();
    VectorType displacement ( numVertices );
    
    // construct displacement \phi(x,y,z) with 
    // \phi(x,y,z) = ( \int_0^x cos(K(t) + K_0) - x dt, 0, \int_0^x sin(K(t) + K_0) dt ) )
    // D \phi(x,y,z) = ( cos(K(x) + K_0) - 1    0
    //                        0                 0
    //                   sin(K(x) + K_0)        0
    for( int nodeIdx = 0; nodeIdx < numVertices; ++nodeIdx ){
        const Point3DType oldVertex = mesh.getVertex ( nodeIdx );
        Point3DType newVertex = mesh.getVertex ( nodeIdx );
        RealType K_index_x = getPhase_interpolateSlopeLinear( oldVertex[0], initSlope, finalSlope );
        RealType gamma_x, gamma_z;
        getPoint_interpolateSlopeLinear( oldVertex[0], gamma_x, gamma_z, initSlope, finalSlope );
        
        newVertex[0] = gamma_x;
        newVertex[1] = oldVertex[1];
        newVertex[2] = gamma_z;
        
        mesh.setVertex( nodeIdx, newVertex );
        
        
        if( MeshType::_hasTangentSpace ){
            
            const TangentVecType oldTangentVec1 = mesh.getTangentVec1 ( nodeIdx ), oldTangentVec2 = mesh.getTangentVec2( nodeIdx );
            TangentVecType newTangentVec1 ( oldTangentVec1 ), newTangentVec2 ( oldTangentVec2 );
            
            newTangentVec1[0] += cos( K_index_x ) - 1.0;
            newTangentVec1[2] += sin( K_index_x );
            
            mesh.setTangentVec1 ( nodeIdx, newTangentVec1 );
            mesh.setTangentVec2 ( nodeIdx, newTangentVec2 );
        }
    }
    
//     for(int i=0; i< numVertices; ++i){
//       if( _mask.get( i ) ){
//         for( int j=0; j<3; ++j){
//           w[j].set( i, 0.0 );   
//           w[j].set( i + 2 * numVertices, 0.0 ); 
//         }        
//         w[0].set( i + numVertices, cos ( atan(initSlope) ) - 1.0  ); 
//         w[1].set( i + numVertices, 0.0  ); 
//         w[2].set( i + numVertices, sin ( atan(initSlope) )  ); 
//       }
//     }

      mesh.updateAllTriangles ( );
      mesh.updateAllProjectionCoefficients();
    
  }
  
  
  
//   // use 1D-simple displacement f : \phi(x,y) = ( f^1(x), y, f^2(x) + a y ) 
//   void constructSimpleDisplacementWithNonZeroDirichletConditionLeft( pesopt::MultiVector<RealType> &w, const RealType initSlope, const RealType finalSlope, const RealType Y1  ) const{
//     
//     ConfiguratorType conf ( _mesh );
//     int numVertices = _mesh.getNumVertices();
//     
//     // construct displacement w = \phi - (x,y,0) with 
//     // w(x,y) = ( \int_0^x cos(K(t) + K_0) dt - x, 0, \int_0^x sin(K(t) + K_0) dt + Y1 y ) )
//     // D w(x,y) = ( cos(K(x) + K_0) - 1       0
//     //                   0                    0
//     //              sin(K(x) + K_0)           Y1
//     for( NodeIteratorType iter = conf.begin(); iter.notAtEnd(); ++iter ){
//         int index = iter.getIndex();
//         pesopt::Vec3<RealType> coords ( iter.getCoords() );
//         if( coords[0] == 0 ){
//           for( int j=0; j<2; ++j){
//             w[j].set( index, 0.0 );   
//             w[j].set( index + 2 * numVertices, 0.0 ); 
//           }    
//           w[2].set( index, Y1 * coords[1] );   
//           w[2].set( index + 2 * numVertices, Y1 ); 
//           w[0].set( index + numVertices, cos ( atan(initSlope) ) - 1.0  ); 
//           w[1].set( index + numVertices, 0.0  ); 
//           w[2].set( index + numVertices, sin ( atan(initSlope) )  ); 
//         } else {
//           RealType K_index_x = getPhase_interpolateSlopeLinear( coords[0], initSlope, finalSlope );
//           RealType gamma_x, gamma_z;
//           getPoint_interpolateSlopeLinear( coords[0], gamma_x, gamma_z, initSlope, finalSlope );
//           w[0].set( index , gamma_x - coords[0] );
//           w[0].set( index + numVertices, cos( K_index_x ) - 1.0 );
//           w[0].set( index + 2 * numVertices, 0.0 );
//           w[1].set( index , 0.0 );
//           w[1].set( index + numVertices, 0.0 );
//           w[1].set( index + 2 * numVertices, 0.0 );
//           w[2].set( index , gamma_z + Y1 * coords[1] );
//           w[2].set( index + numVertices, sin( K_index_x ) );
//           w[2].set( index + 2*numVertices, Y1 );
//         }
//         
//     }
//     
//   }
//   
//   
//   
//   
//   
//   
//   
//   void constructSimpleDisplacementFixedLeftAndRight( pesopt::MultiVector<RealType> &w, const RealType initSlope, const RealType finalSlope ) const{
//     
//     ConfiguratorType conf ( _mesh );
//     int numVertices = _mesh.getNumVertices();
//     
//     // construct displacement w = \phi - (x,y,0) with 
//     // w(x,y) = ( \int_0^x cos(K(t) + K_0) dt - x, 0, \int_0^x sin(K(t) + K_0) dt ) )
//     // D w(x,y) = ( cos(K(x) + K_0) - 1       0
//     //                   0                    0
//     //              sin(K(x) + K_0)           0
//     for( NodeIteratorType iter = conf.begin(); iter.notAtEnd(); ++iter ){
//         int index = iter.getIndex();
//         pesopt::Vec3<RealType> coords ( iter.getCoords() );
//         RealType K_index_x = getPhase_interpolateSlopeLinear( coords[0], initSlope, finalSlope );
//         RealType gamma_x, gamma_z;
//         getPoint_interpolateSlopeLinear( coords[0], gamma_x, gamma_z, initSlope, finalSlope );
//         w[0].set( index , gamma_x - coords[0] );
//         w[0].set( index + numVertices, cos( K_index_x ) - 1.0 );
//         w[0].set( index + 2 * numVertices, 0.0 );
//         w[1].set( index , 0.0 );
//         w[1].set( index + numVertices, 0.0 );
//         w[1].set( index + 2 * numVertices, 0.0 );
//         w[2].set( index , gamma_z );
//         w[2].set( index + numVertices, sin( K_index_x ) );
//         w[2].set( index + 2*numVertices, 0.0 );
//         
//         if( coords[0] == 0 ){
//             for( int j=0; j<3; ++j){
//             w[j].set( index, 0.0 );   
//             w[j].set( index + 2 * numVertices, 0.0 ); 
//             }        
//             w[0].set( index + numVertices, cos ( atan(initSlope) ) - 1.0  ); 
//             w[1].set( index + numVertices, 0.0  ); 
//             w[2].set( index + numVertices, sin ( atan(initSlope) )  ); 
//         }
//     }
//     
//   }
//   
//   
//   
//   
//   
//   
//   RealType getPhase_Twist( const RealType & coords_x, const RealType initSlope = 0.0 ) const {
// 
//     RealType init = atan( initSlope );
//     //RealType final = atan ( finalSlope );
//     RealType fin = 6 * atan ( 1.0 ); //=3pi/2
//     
//     RealType sizeFirstPart = 1. / 8.;
//     RealType endTwist = 1. / 2.;
//     
//     if( coords_x < sizeFirstPart ){
//       return init;
//     }
//     if( (coords_x >= sizeFirstPart) && (coords_x <= endTwist) ){
//       return  (coords_x - sizeFirstPart) * (fin - init) / (endTwist - sizeFirstPart ) + init;
//     } else {
//       return fin;
//     }
//     
//   }
//   
//   // computes int_0^x exp( i (K(t) + K(0) ) d t
//   void computeCurveFromPhaseForTwist( const RealType & coords_x, RealType & gamma_x, RealType & gamma_z, const RealType initSlope ) const{
//     
//     gamma_x = 0.0;
//     gamma_z = 0.0;
//     
//     int numSteps = 50; //TODO numQuadPoints
//     RealType tau = coords_x / static_cast<RealType> (numSteps);
//     
//     for(int i=0; i< numSteps; ++i){
//       RealType coord_i = ( static_cast<RealType> (i) + 0.5 ) * tau;
//       RealType K_i = getPhase_Twist ( coord_i, initSlope );
//       gamma_x += tau * cos( K_i );
//       gamma_z += tau * sin( K_i );
//     }
//   }
//   
//   void construcTwistDisplacement( pesopt::MultiVector<RealType> &w, const RealType initSlope = 0.0 ) const{
//     
//     ConfiguratorType conf ( _mesh );
//     
//     // construct displacement w = \phi - (x,y,0) with 
//     // w(x,y) = ( \int_0^x cos(K(t) + K_0) dt - x, 0, \int_0^x sin(K(t) + K_0) dt ) )
//     // D w(x,y) = ( cos(K(x) + K_0) - 1       0
//     //                   0                    0
//     //              sin(K(x) + K_0)           0
//     for( NodeIteratorType iter = conf.begin(); iter.notAtEnd(); ++iter ){
//         int index = iter.getIndex();
//         pesopt::Vec3<RealType> coords ( iter.getCoords() );
//         RealType K_index_x = getPhase_Twist( coords[0], initSlope );
//         RealType gamma_x, gamma_z;
//         computeCurveFromPhaseForTwist( coords[0], gamma_x, gamma_z, initSlope );
//         w[0].set( index , gamma_x - coords[0] );
//         w[0].set( index + _numVertices, cos( K_index_x ) - 1.0 );
//         w[0].set( index + 2 * _numVertices, 0.0 );
//         w[1].set( index , 0.0 );
//         w[1].set( index + _numVertices, 0.0 );
//         w[1].set( index + 2 * _numVertices, 0.0 );
//         w[2].set( index , gamma_z );
//         w[2].set( index + _numVertices, sin( K_index_x ) );
//         w[2].set( index + 2*_numVertices, 0.0 );
//     }
//     
//     for(int i=0; i< _numVertices; ++i){
//       if( _mask.get( i ) ){
//         for( int j=0; j<3; ++j){
//           w[j].set( i, 0.0 );   
//           w[j].set( i + 2 * _numVertices, 0.0 ); 
//         }        
//         w[0].set( i + _numVertices, cos ( atan(initSlope) ) - 1.0  ); 
//         w[1].set( i + _numVertices, 0.0  ); 
//         w[2].set( i + _numVertices, sin ( atan(initSlope) )  ); 
//       }
//     }
//     
//   }
//     
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////
  ////
//   void switchDisplacementType(  pesopt::MultiVector<RealType> &w ) const{
//       w.reallocate( 3 , _numGlobalDofs );
//       switch( _parser.getInt( "initDisplacementType" ) ){
//         // zero displacement 
//         case 0:
//           w[0].setZero();
//           w[1].setZero();
//           w[2].setZero();
//           break;
//         // simple displacement
//         case 1 :
//           constructSimpleDisplacement( w, _parser.getDouble( "initSlope" ), _parser.getDouble( "finalSlope" ) );
//           break;
//         // twist
//         case 2 :
//           construcTwistDisplacement( w, _parser.getDouble( "initSlope" ) );
//           break;
//         // simple displacement with left right fixed
//         case 3 :
//           constructSimpleDisplacementFixedLeftAndRight( w, _parser.getDouble( "initSlope" ), _parser.getDouble( "finalSlope" ) );
//           break;
//         case 4 :
//           constructSimpleDisplacementWithNonZeroDirichletConditionLeft ( w, _parser.getDouble( "initSlope" ), _parser.getDouble( "finalSlope" ), _parser.getDouble( "slopeInYDirection" ) );
//           break;
//         default:
//           throw pesopt::Exception ( "Wrong initDisplacementType !", __FILE__, __LINE__ );
//           break;
//       }
//   }
  
  
};

#endif
