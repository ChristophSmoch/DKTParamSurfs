#ifndef __FECONFIGURATOR_H
#define __FECONFIGURATOR_H

#include <pesopt_IO.h>

// template <typename DataTypeContainer, 
//           typename MeshType, 
//           typename QuadType
class FEConfigurator {
  
public:
  
//   static const int dimChartDomain = 2;
//   static const int maxNumLocalDofs = 1;
//   static const int maxNumLocalBoundaryDofs = 2;
//   static const int numAffineSymGradDofs = 3;
//   static const int sizeVoigtTensor = 6;

//   typedef DataTypeContainer                                             DTContainer;
//   typedef MeshType                                                      InitType;        
//   typedef QuadType                                                      QuadRuleType;
//   typedef typename MeshType::ElementType                                ElementType;
//   typedef typename DataTypeContainer::RealType                          RealType;
//   typedef typename DataTypeContainer::RealVecChart                      RealVecChart;
//   typedef typename DataTypeContainer::TangentVecType                    TangentVecType;
//   typedef typename DataTypeContainer::PointType                         PointType;
//   typedef typename DataTypeContainer::Matrix22                          Matrix22;
//   typedef typename DataTypeContainer::Matrix32                          Matrix32;
//   typedef typename DataTypeContainer::Matrix33                          Matrix33;
//   typedef typename DataTypeContainer::VectorType                        VectorType;
//   typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
//   typedef typename DataTypeContainer::TripletType                       TripletType;
//   typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
//   typedef typename DataTypeContainer::MaskType                          MaskType;
//   typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
//   typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;
//   typedef RefTriangMeshBaseFunctionSetP0<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;

// protected:
//   
//   const MeshType &_mesh;
//   mutable BaseFuncSetType _baseFuncSet;
  
public:
  
//   FEConfigurator ( const InitType &Mesh ) : _mesh ( Mesh ) {}

//   const InitType& getMesh( ) const { return this->_mesh; }
// 
//   inline int getNumLocalDofs ( const ElementType & ) const { return 1;}
//   inline int getNumLocalDofs (  ) const { return 1;}
// 
//   int getNumGlobalDofs( ) const {return this->_mesh.getNumElements();}
//   int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
//   
//   const BaseFuncSetType& getBaseFunctionSet ( const ElementType &T ) const {
//       _baseFuncSet.setTriangle ( T );
//       return _baseFuncSet;
//   }
//   
//   //! returns global index of the dof with number localIndex
//   inline int localToGlobal ( const ElementType &T, int /*localIndex*/ ) const {return T.getGlobalElementIdx( );}

};


#endif
