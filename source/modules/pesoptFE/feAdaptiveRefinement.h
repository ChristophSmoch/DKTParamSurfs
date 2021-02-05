#ifndef __FEADAPTIVEREFINEMENT_H
#define __FEADAPTIVEREFINEMENT_H

#include <pesopt_IO.h>
 
#include <feDefines.h>



template<typename MeshType>
class FEAdaptiveMeshBase {
    
protected:
    
  typedef typename MeshType::ElementType                    ElementType;
  typedef typename ElementType::ElementNodeIndicesType      ElementNodeIndicesType;
  typedef typename MeshType::RealType                       RealType;
  typedef typename MeshType::RealVecChart                   RealVecChart;
  typedef typename MeshType::IntVecChart                    IntVecChart;
  typedef typename MeshType::PointType                      PointType;
  typedef typename MeshType::VectorType                     VectorType;
  typedef typename MeshType::MaskType                       MaskType;
    
  const MeshType &_initMesh;
  mutable std::vector<MeshType> _meshVector;
  mutable int _refinementLevel;
    
public:
    
    
  FEAdaptiveMeshBase( const MeshType &mesh ) :
   _initMesh ( mesh ), _refinementLevel( 0 ) {
       _meshVector.push_back( _initMesh );
 }

  const MeshType& getCurrentMesh( ) const { return _meshVector[_refinementLevel]; }
    
  virtual void refineAll( ) const{
      throw std::logic_error ( pesopt::strprintf ( "FEAdaptiveMesh::refineAll is not defined in %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  virtual void refineMarkedElements( const MaskType &markedElements ) const{
     throw std::logic_error ( pesopt::strprintf ( "FEAdaptiveMesh::refineMarked is not defined in %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  virtual void refineByFactor( const RealType factor ) const {
    throw std::logic_error ( pesopt::strprintf ( "FEAdaptiveMesh::refineByFactor is not defined in %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
    
    
  template <typename ConfiguratorType>
  void prolongateScalarFunctionAtNodes ( VectorType &vec ) const {
     VectorType vecOld(vec);
     this->template prolongateScalarFunctionAtNodes<ConfiguratorType> ( vecOld, vec );
 }
 
 
  template <typename ConfiguratorType>
  void prolongateScalarFunctionAtNodes ( const VectorType &vecOld, VectorType &vecNew  ) const {
     throw std::logic_error ( pesopt::strprintf ( "FEAdaptiveMesh::prolongateScalarFunctionAtNodes is not defined in %s at line %d.", __FILE__, __LINE__ ).c_str() );
 }
 
 
  template <typename ConfiguratorType>
  void prolongateScalarFunctionOnElements ( VectorType &vec ) const {
     VectorType vecOld(vec);
     this->template prolongateScalarFunctionOnElements<ConfiguratorType> ( vecOld, vec );
 }
 
 
  template <typename ConfiguratorType>
  void prolongateScalarFunctionOnElements ( const VectorType &vecOld, VectorType &vecNew  ) const {
     throw std::logic_error ( pesopt::strprintf ( "FEAdaptiveMesh::prolongateScalarFunctionOnElements is not defined in %s at line %d.", __FILE__, __LINE__ ).c_str() );
 }
    
    
    
  template <typename ConfiguratorType>
  void prolongateVectorFunctionOnElements ( VectorType &vec ) const {
     VectorType vecOld(vec);
     this->template prolongateVectorFunctionOnElements<ConfiguratorType> ( vecOld, vec );
 }
 
 
  template <typename ConfiguratorType>
  void prolongateVectorFunctionOnElements ( const VectorType &vecOld, VectorType &vecNew  ) const {
     throw std::logic_error ( pesopt::strprintf ( "FEAdaptiveMesh::prolongateVectorFunctionOnElements is not defined in %s at line %d.", __FILE__, __LINE__ ).c_str() );
 }
    
};



template<typename MeshType, PESOPT_FEMeshType type = MeshType::_FEMeshType >
class FEAdaptiveMesh : public FEAdaptiveMeshBase<MeshType>{
    
public:
    
  FEAdaptiveMesh( const MeshType &mesh ) :
    FEAdaptiveMeshBase<MeshType> ( mesh ) { }

};


#endif
