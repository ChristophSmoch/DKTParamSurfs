#ifndef __TRIANGLEFECONFIGURATORS_H
#define __TRIANGLEFECONFIGURATORS_H

#include <pesopt_IO.h>
#include "triangleFEDefines.h"
#include "triangleFEBaseFunctionSets.h"
#include "triangleFEMesh.h"
#include "triangleFEQuadrature.h"
  



template <typename DataTypeContainer = DataTypeContainerTriangleFE, 
          typename MeshType = TriangleMesh<>, 
          typename QuadType = CenterQuadrature<typename DataTypeContainer::RealType, typename DataTypeContainer::RealVecChart > >
class RefTriangMeshConfiguratorP0 {
  
public:
  
  static const int dimChartDomain = 2;
  static const int maxNumLocalDofs = 1;
  static const int maxNumLocalBoundaryDofs = 2;
  static const int numAffineSymGradDofs = 3;
  static const int sizeVoigtTensor = 6;
//   static const ShellFEType _ShellFEType = ElementValuedDofs;

  typedef DataTypeContainer                                             DTContainer;
  typedef MeshType                                                      InitType;        
  typedef QuadType                                                      QuadRuleType;
  typedef typename MeshType::ElementType                                ElementType;
  typedef typename DataTypeContainer::RealType                          RealType;
  typedef typename DataTypeContainer::IntVecChart                       IntVecChart;
  typedef typename DataTypeContainer::RealVecChart                      RealVecChart;
  typedef typename DataTypeContainer::TangentVecType                    TangentVecType;
  typedef typename DataTypeContainer::DerivativeVectorValuedType        DerivativeVectorValuedType;
  typedef typename DataTypeContainer::PointType                         PointType;
  typedef typename DataTypeContainer::Matrix22                          Matrix22;
  typedef typename DataTypeContainer::Matrix32                          Matrix32;
  typedef typename DataTypeContainer::Matrix33                          Matrix33;
  typedef typename DataTypeContainer::VectorType                        VectorType;
  typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
  typedef typename DataTypeContainer::TripletType                       TripletType;
  typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
  typedef typename DataTypeContainer::MaskType                          MaskType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;
  
  typedef RefTriangMeshBaseFunctionSetP0<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;

  
//   //! globally affine functions with constant gradient: \f$ u(x) = a \cdot x \f$
//   typedef Eigen::Matrix<RealType, dimChartDomain, dimChartDomain >                LocalMatrixTypeAffineGrad;
//   typedef Eigen::Matrix<RealType, dimChartDomain, maxNumLocalDofs >               LocalMatrixTypeMixedAffineGrad;
//   typedef GlobalAffineGradBaseFunctionSet2D<DataTypeContainer>                    GlobalAffineGradBaseFuncSet;
//   
//   //! globally affine functions with symmetric gradient
//   typedef Eigen::Matrix<RealType, numAffineSymGradDofs, numAffineSymGradDofs >    LocalMatrixTypeAffineSymGrad;
//   typedef Eigen::Matrix<RealType, numAffineSymGradDofs, maxNumLocalDofs >         LocalMatrixTypeMixedAffineSymGrad;
//   typedef GlobalAffineSymGradBaseFunctionSet2D<DataTypeContainer>                 GlobalAffineSymGradBaseFuncSet;
  
  
protected:
  
  const MeshType &_mesh;
  mutable BaseFuncSetType _baseFuncSet;
  
public:
  
  RefTriangMeshConfiguratorP0 ( const InitType &Mesh ) : _mesh ( Mesh ) {}

  const InitType& getMesh( ) const { return this->_mesh; }

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
  
  void getNodalValuesFromCoefficients( const VectorType &coeff,  VectorType nodalValues ) {
    throw std::invalid_argument( pesopt::strprintf ( "getNodalValuesFromCoefficients so far not implemented for RefTriangMeshConfiguratorP0. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );  
  }

};



template <typename DataTypeContainer = DataTypeContainerTriangleFE, 
          typename MeshType = TriangleMesh<>, 
          typename QuadType = TriQuadrature<typename DataTypeContainer::RealType, typename DataTypeContainer::RealVecChart > >
class RefTriangMeshConfiguratorP1 {
public:
  
  static const int dimChartDomain = 2;
  static const int maxNumLocalDofs = 3;
  static const int maxNumLocalBoundaryDofs = 2;
  static const int numAffineSymGradDofs = 3;
  static const int sizeVoigtTensor = 6;
  
  
//   static const ShellFEType _ShellFEType = NodalValuedDofs;

  typedef DataTypeContainer                                             DTContainer;
  typedef MeshType                                                      InitType;      
  typedef QuadType                                                      QuadRuleType;
  typedef typename MeshType::ElementType                                ElementType;
  typedef typename DataTypeContainer::RealType                          RealType;
  typedef typename DataTypeContainer::IntVecChart                       IntVecChart;
  typedef typename DataTypeContainer::RealVecChart                      RealVecChart;
  typedef typename DataTypeContainer::TangentVecType                    TangentVecType;
  typedef typename DataTypeContainer::PointType                         PointType;
  typedef typename DataTypeContainer::DerivativeVectorValuedType        DerivativeVectorValuedType;
  typedef typename DataTypeContainer::Matrix22                          Matrix22;
  typedef typename DataTypeContainer::Matrix32                          Matrix32;
  typedef typename DataTypeContainer::Matrix33                          Matrix33;
  typedef typename DataTypeContainer::VectorType                        VectorType;
  typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
  typedef typename DataTypeContainer::TripletType                       TripletType;
  typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
  typedef typename DataTypeContainer::MaskType                          MaskType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;
  typedef Eigen::Matrix<RealType, sizeVoigtTensor, 1 >                  VoigtTensorVec;
  
  typedef RefTriangMeshBaseFunctionSetP1<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  
    //! globally affine functions with constant gradient: \f$ u(x) = a \cdot x \f$
  typedef Eigen::Matrix<RealType, dimChartDomain, dimChartDomain >                LocalMatrixTypeAffineGrad;
  typedef Eigen::Matrix<RealType, dimChartDomain, maxNumLocalDofs >               LocalMatrixTypeMixedAffineGrad;
  typedef GlobalAffineGradBaseFunctionSet2D<DataTypeContainer>                    GlobalAffineGradBaseFuncSet;
  
  //! globally affine functions with symmetric gradient
  typedef Eigen::Matrix<RealType, numAffineSymGradDofs, numAffineSymGradDofs >    LocalMatrixTypeAffineSymGrad;
  typedef Eigen::Matrix<RealType, numAffineSymGradDofs, maxNumLocalDofs >         LocalMatrixTypeMixedAffineSymGrad;
  typedef GlobalAffineSymGradBaseFunctionSet2D<DataTypeContainer>                 GlobalAffineSymGradBaseFuncSet;
  
protected:
  
  const MeshType &_mesh;
  mutable BaseFuncSetType _baseFuncSet;
  
public:
  
  RefTriangMeshConfiguratorP1 ( const InitType &Mesh ) : _mesh ( Mesh ) {}

  const InitType& getMesh( ) const { return this->_mesh; }

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

  void getNodalValuesFromCoefficients( const VectorType &coeff,  VectorType nodalValues ) {
    nodalValues = coeff;   
  }
  
};


#endif
