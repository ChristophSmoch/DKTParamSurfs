#ifndef __DKTFEMESHINFO_H
#define __DKTFEMESHINFO_H

#include <dktFEMesh.h>
    
template< typename MeshType >
class DKTTriangMeshInfo {
protected:
  typedef typename MeshType::RealType RealType;
  typedef typename MeshType::PointType PointType;
  
  const MeshType& _mesh;
  
public:
  DKTTriangMeshInfo ( const MeshType& Mesh ) : _mesh( Mesh )  { }
  
  RealType getMaxAreaSqrt( ) const{
      RealType maxarea = 0.0;
      for ( int elementIndex = 0; elementIndex < _mesh.getNumTriangs(); ++elementIndex  ){
          RealType area = _mesh.getTriang( elementIndex ).getAreaOfFlattenedTriangle();
          if( area > maxarea ) maxarea = area;
      }
      return std::sqrt( maxarea );
  }

  RealType getMinAreaSqrt( ) const{
      RealType minarea = 1.e+100;
      for ( int elementIndex = 0; elementIndex < _mesh.getNumTriangs(); ++elementIndex  ){
          RealType area = _mesh.getTriang( elementIndex ).getAreaOfFlattenedTriangle();
          if( area < minarea ) minarea = area;
      }
      return std::sqrt( minarea );
  }

  RealType getInterfaceWith ( const int type, const RealType eps_factor ) const {
        RealType eps = 0.0;
        switch( type ){
        case 1 : 
            eps = this->getMaxAreaSqrt( );
            break;
        case 2 : 
            eps = this->getMinAreaSqrt( );
            break;
        default : 
            throw std::invalid_argument ( pesopt::strprintf ( "wrong method in file %s at line %d.", __FILE__, __LINE__ ).c_str() );
            break;
        }
        return eps * eps_factor;
  }
  
  
  void getBoundingBox( PointType &lowerBounds, PointType &upperBounds ) const{
      lowerBounds = _mesh.getVertex ( 0 );
      upperBounds = _mesh.getVertex ( 0 );
      for ( int nodeIndex = 0; nodeIndex < _mesh.getNumVertices(); ++nodeIndex  ){
          const PointType& vertex = _mesh.getVertex ( nodeIndex );
          for( int coord=0; coord < 3; ++coord ){
              if( vertex[coord] < lowerBounds[coord] ) lowerBounds[coord] = vertex[coord];
              if( vertex[coord] > upperBounds[coord] ) upperBounds[coord] = vertex[coord];
          }
      }
  }
  

};



#endif
