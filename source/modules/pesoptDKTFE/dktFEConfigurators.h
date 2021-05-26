#ifndef __DKTFECONFIGURATORS_H
#define __DKTFECONFIGURATORS_H

#include <pesopt_IO.h>

#include <dktFEBaseFunctionSets.h>
  
    
enum ShellFEType {
  ElementValuedDofs,
  NodalValuedDofs,
  C1Dofs,
};



template <typename DataTypeContainer, typename MeshType, typename QuadType >
class RefTriangMeshConfiguratorP0 {
  
public:
  
  static const int maxNumLocalDofs = 1;
  static const ShellFEType _ShellFEType = ElementValuedDofs;

  typedef DataTypeContainer                                             DTContainer;
  typedef MeshType                                                      InitType;        
  typedef QuadType                                                      QuadRuleType;
  typedef typename MeshType::ElementType                                ElementType;
  typedef typename DataTypeContainer::RealType                          RealType;
  typedef typename DataTypeContainer::DomVecType                        DomVecType;
  typedef typename DataTypeContainer::TangentVecType                    TangentVecType;
  typedef typename DataTypeContainer::Point3DType                       Point3DType;
  typedef typename DataTypeContainer::Matrix22                          Matrix22;
  typedef typename DataTypeContainer::Matrix32                          Matrix32;
  typedef typename DataTypeContainer::Matrix33                          Matrix33;
  typedef typename DataTypeContainer::Tensor222Type                     Tensor222Type;
  typedef typename DataTypeContainer::Tensor322Type                     Tensor322Type;
  typedef typename DataTypeContainer::Tensor332Type                     Tensor332Type;
  typedef typename DataTypeContainer::Tensor333Type                     Tensor333Type;
  typedef typename DataTypeContainer::VectorType                        VectorType;
  typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
  typedef typename DataTypeContainer::TripletType                       TripletType;
  typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
  typedef typename DataTypeContainer::MaskType                          MaskType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;
  
  typedef RefTriangMeshBaseFunctionSetP0<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;

protected:
  
  const MeshType &_mesh;
  mutable BaseFuncSetType _baseFuncSet;
  
public:
  
  RefTriangMeshConfiguratorP0 ( const InitType &Mesh ) : _mesh ( Mesh ) {}

  const InitType& getInitializer( ) const { return this->_mesh; }

  inline int getNumLocalDofs ( const ElementType & ) const { return 1;}
  inline int getNumLocalDofs (  ) const { return 1;}

  int getNumGlobalDofs( ) const {return this->_mesh.getNumElements();}
  int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
  
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &T ) const {
      _baseFuncSet.setTriangle ( T );
      return _baseFuncSet;
  }
  
  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int /*localIndex*/ ) const {return T.getGlobalElementIdx( );}

};



template <typename DataTypeContainer, typename MeshType, typename QuadType >
class RefTriangMeshConfiguratorP1 {
public:
  
  static const int maxNumLocalDofs = 3;
  static const ShellFEType _ShellFEType = NodalValuedDofs;

  typedef DataTypeContainer                                             DTContainer;
  typedef MeshType                                                      InitType;      
  typedef QuadType                                                      QuadRuleType;
  typedef typename MeshType::ElementType                                ElementType;
  typedef typename DataTypeContainer::RealType                          RealType;
  typedef typename DataTypeContainer::DomVecType                        DomVecType;
  typedef typename DataTypeContainer::TangentVecType                    TangentVecType;
  typedef typename DataTypeContainer::Point3DType                       Point3DType;
  typedef typename DataTypeContainer::Matrix22                          Matrix22;
  typedef typename DataTypeContainer::Matrix32                          Matrix32;
  typedef typename DataTypeContainer::Matrix33                          Matrix33;
  typedef typename DataTypeContainer::Tensor222Type                     Tensor222Type;
  typedef typename DataTypeContainer::Tensor322Type                     Tensor322Type;
  typedef typename DataTypeContainer::Tensor332Type                     Tensor332Type;
  typedef typename DataTypeContainer::Tensor333Type                     Tensor333Type;
  typedef typename DataTypeContainer::VectorType                        VectorType;
  typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
  typedef typename DataTypeContainer::TripletType                       TripletType;
  typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
  typedef typename DataTypeContainer::MaskType                          MaskType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;
  
  typedef RefTriangMeshBaseFunctionSetP1<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  //typedef RefTriangMeshCachedBaseFunctionSetP1<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  
protected:
  
  const MeshType &_mesh;
  mutable BaseFuncSetType _baseFuncSet;
  
public:
  
  RefTriangMeshConfiguratorP1 ( const InitType &Mesh ) : _mesh ( Mesh ) {}

  const InitType& getInitializer( ) const { return this->_mesh; }

  inline int getNumLocalDofs ( const ElementType &/*T*/ ) const { return 3;}
  inline int getNumLocalDofs (  ) const { return 3;}

  int getNumGlobalDofs( ) const { return this->_mesh.getNumVertices();}
  int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &T ) const {
      _baseFuncSet.setTriangle ( T );
      return _baseFuncSet;
  }
  
  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int localIndex ) const {return T.getGlobalNodeIdx( localIndex );}

};



template <typename DataTypeContainer, typename MeshType, typename QuadType >
class RefTriangMeshConfiguratorDKT {
protected:

  const MeshType &_mesh;
  
public: 

  static const int maxNumLocalDofs = 9;
  static const ShellFEType _ShellFEType = C1Dofs;

  typedef DataTypeContainer                                             DTContainer;
  typedef MeshType                                                      InitType; 
  typedef QuadType                                                      QuadRuleType;
  typedef typename DataTypeContainer::RealType                          RealType;
  typedef typename DataTypeContainer::DomVecType                        DomVecType;
  typedef typename DataTypeContainer::TangentVecType                    TangentVecType;
  typedef typename DataTypeContainer::Point3DType                       Point3DType;
  typedef typename DataTypeContainer::Matrix22                          Matrix22;
  typedef typename DataTypeContainer::Matrix32                          Matrix32;
  typedef typename DataTypeContainer::Matrix33                          Matrix33;
  typedef typename DataTypeContainer::Tensor222Type                     Tensor222Type;
  typedef typename DataTypeContainer::Tensor322Type                     Tensor322Type;
  typedef typename DataTypeContainer::Tensor332Type                     Tensor332Type;
  typedef typename DataTypeContainer::Tensor333Type                     Tensor333Type;
  typedef typename DataTypeContainer::VectorType                        VectorType;
  typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
  typedef typename DataTypeContainer::TripletType                       TripletType;
  typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;
  typedef typename DataTypeContainer::MaskType                          MaskType;
  typedef typename MeshType::ElementType                                ElementType;
  
  typedef RefTriangMeshBaseFunctionSetDKT<DataTypeContainer, QuadType, ElementType>               BaseFuncSetType;
  typedef RefTriangMeshApproxGradientBaseFunctionSetDKT<DataTypeContainer, QuadType, ElementType> ApproxGradientBaseFuncSetType;
  
  RefTriangMeshConfiguratorDKT ( const InitType &Mesh ) : _mesh ( Mesh ) {}

  const InitType& getInitializer( ) const { return this->_mesh; }

  mutable BaseFuncSetType _baseFuncSet;
  mutable ApproxGradientBaseFuncSetType _approxGradBaseFuncSet;

  inline int getNumLocalDofs ( const ElementType & ) const { return 9;}
  inline int getNumLocalDofs (  ) const { return 9;}

  int getNumGlobalDofs( ) const {return 3 * this->_mesh.getNumVertices();}
  int maxNumQuadPoints( ) const {return QuadType::numQuadPoints;}
  
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &T ) const {
    _baseFuncSet.setTriangle ( T );
    return _baseFuncSet;
  }
  
  const ApproxGradientBaseFuncSetType& getApproxGradientBaseFunctionSet ( const ElementType &T ) const {
    _approxGradBaseFuncSet.setTriangle ( T );
    return _approxGradBaseFuncSet;
  }
  
  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int localIndex ) const {
    int aux = ( localIndex%3 == 0) ? 0 : (  ( (localIndex+2)%3 == 0) ? 1 : 2 );
    return aux * _mesh.getNumVertices() + T.getGlobalNodeIdx( std::floor(localIndex/3) );
  }

};



#endif //__DKTFECONFIGURATORS_H
