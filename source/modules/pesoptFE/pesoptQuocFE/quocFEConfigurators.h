#ifndef __QUOCFECONFIGURATORS_H
#define __QUOCFECONFIGURATORS_H

#include <pesopt_IO.h>
#include <quocFEBaseFunctionSets.h>
#include <quocFEMesh.h>
#include <quocFEQuadrature.h>


// template <typename DataTypeContainer,
//           typename MeshType,
//           typename QuadType = SimpsonQuadrature1D<typename DataTypeContainer::RealType, typename DataTypeContainer::RealVecChart> 
// //           typename QuadType = SimpsonQuadrature1D<typename DataTypeContainer::RealType> 
//           >
// class QuocConfigurator :
// public FEConfigurator {
//     QuocConfigurator () {}
//     
// }







template <typename DataTypeContainer = Quoc1DDataTypeContainer,
          typename MeshType = QuocMesh1D<DataTypeContainer>,
          typename QuadType = SimpsonQuadrature1D<typename DataTypeContainer::RealType, typename DataTypeContainer::RealVecChart> 
//           typename QuadType = SimpsonQuadrature1D<typename DataTypeContainer::RealType> 
          >
class QuocConfigurator1D {
public:
    static const int dimChartDomain = 1;
    static const int maxNumLocalDofs = 2;
//     static const int maxNumLocalBoundaryDofs = 2;
    static const int numAffineSymGradDofs = 1;
    static const int sizeVoigtTensor = 1;

  typedef DataTypeContainer                                             DTContainer;
  typedef QuadType                                                      QuadRuleType;
//   typedef BdrQuadType                                                   BoundaryQuadType;
  typedef MeshType                                                      InitType;               //!< that's the type that is needed by the constructor of the configurator
  typedef typename MeshType::ElementType                                ElementType;
//   typedef typename MeshType::BoundaryElementType                        BoundaryElementType;
  typedef typename DataTypeContainer::RealType                          RealType;
  typedef typename DataTypeContainer::RealVecChart                      RealVecChart;
  typedef typename DataTypeContainer::IntVecChart                        IntVecChart;
//   typedef typename DataTypeContainer::RealVecChartBoundary                RealVecChartBoundary;
  typedef typename DataTypeContainer::PointType                         PointType;
  typedef typename DataTypeContainer::Matrix22                          Matrix22;
  typedef typename DataTypeContainer::Matrix32                          Matrix32;
  typedef typename DataTypeContainer::Matrix33                          Matrix33;
  typedef typename DataTypeContainer::DerivativeVectorValuedType        DerivativeVectorValuedType;
  typedef typename DataTypeContainer::VectorType                        VectorType;
  typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
  typedef typename DataTypeContainer::TripletType                       TripletType;
  typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
  typedef typename DataTypeContainer::MaskType                          MaskType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;
  typedef Eigen::Matrix<RealType, sizeVoigtTensor, 1 >                  VoigtTensorVec;
  
//   typedef QuocBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  typedef QuocCachedBaseFunctionSet1D<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  
  //! globally affine functions with symmetric gradient: u(x) = A x
  typedef Eigen::Matrix<RealType, numAffineSymGradDofs, numAffineSymGradDofs >    LocalMatrixTypeAffineSymGrad;
  typedef Eigen::Matrix<RealType, numAffineSymGradDofs, maxNumLocalDofs >         LocalMatrixTypeMixedAffineSymGrad; //old: LocalMatrixTypeMixed
  typedef GlobalAffineSymGradBaseFunctionSet2D<DataTypeContainer> GlobalAffineSymGradBaseFuncSet;
  
  //! globally affine functions with constant gradient: \f$ u(x) = a \cdot x \f$
  typedef Eigen::Matrix<RealType, dimChartDomain, dimChartDomain >               LocalMatrixTypeAffineGrad;
  typedef Eigen::Matrix<RealType, dimChartDomain, maxNumLocalDofs >         LocalMatrixTypeMixedAffineGrad;
  typedef GlobalAffineGradBaseFunctionSet2D<DataTypeContainer> GlobalAffineGradBaseFuncSet;
  
protected:  
  const MeshType &_mesh;
  const QuadType _quad;
  const BaseFuncSetType _baseFuncSet;
  
public:
  
  QuocConfigurator1D ( const InitType &Mesh ) :
  _mesh ( Mesh ), 
  _baseFuncSet( _mesh.getMeshSize(0) )
  {}

  const InitType& getMesh( ) const { return this->_mesh; }

  inline int getNumLocalDofs ( const ElementType & ) const { return maxNumLocalDofs;}
  inline int getNumLocalDofs (  ) const { return maxNumLocalDofs;}

  int getNumGlobalDofs( ) const {return this->_mesh.getNumVertices();}
  int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &/*T*/ ) const {return _baseFuncSet;}
  const BaseFuncSetType& getBaseFunctionSet ( ) const {return _baseFuncSet;}
  
  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &El, int localIndex ) const {return El.getGlobalNodeIdx( localIndex );}
  
  void getLocalCoords ( const PointType &GlobalCoord, int &elementNumber, PointType &LocalCoord ) const {
      PointType BasePoint;
      std::vector<int> BasePointFac ( 1 );
      for ( int comp = 0; comp < 1; ++comp ) {
          const RealType sc = GlobalCoord[comp] / _mesh.getMeshSize(comp);
          BasePointFac[comp] = static_cast<int> ( sc );
          BasePoint[comp] = BasePointFac[comp] * _mesh.getMeshSize(comp);
          LocalCoord[comp] = sc - BasePointFac[comp];
          //special case: "right" boundary
          if( std::abs( _mesh.getWidth( comp ) - GlobalCoord[comp] ) < 2.e-16 ){
              BasePoint[comp] = _mesh.getWidth(comp) - _mesh.getMeshSize(comp);
              BasePointFac[comp]  = _mesh.getNumDofs( comp ) - 2.;
              LocalCoord[comp] = 1.;
          }
      }
      elementNumber = _mesh.getGlobalElementIndex( BasePointFac[0] );
  }
  
//   void getGlobalCoordsForBoundaryElement ( const BoundaryElementType &El, const PointType &LocalCoord, PointType &Coord ) const {
//       Coord = LocalCoord; 
//       Coord[0] *= _mesh.getMeshSize(0); Coord[1] *= _mesh.getMeshSize(1);
//       Coord += El.getElement().getNode(0); //Coord Left Down
//   }

  void getNodalValuesFromCoefficients ( const VectorType& coefficients , VectorType& nodalValues ) const {
      nodalValues = coefficients;
  }

};
    
    
    
    
    
    
template <typename DataTypeContainer = Quoc2DDataTypeContainer,
          typename MeshType = QuocMesh2D<DataTypeContainer>,
          typename QuadType = SimpsonQuadrature2D<typename DataTypeContainer::RealType, typename DataTypeContainer::RealVecChart, typename DataTypeContainer::RealVecChart1D>,
          typename BdrQuadType = SimpsonQuadrature1D<typename DataTypeContainer::RealType, typename DataTypeContainer::RealType> >
class QuocConfigurator2D {
public:
    static const int dimChartDomain = 2;
    static const int maxNumLocalDofs = 4;
    static const int maxNumLocalBoundaryDofs = 2;
    static const int numAffineSymGradDofs = 3;
    static const int sizeVoigtTensor = 6;
    
  typedef DataTypeContainer                                             DTContainer;
  typedef QuadType                                                      QuadRuleType;
  typedef BdrQuadType                                                   BoundaryQuadType;
  typedef MeshType                                                      InitType;
  typedef typename MeshType::ElementType                                ElementType;
  typedef typename MeshType::BoundaryElementType                        BoundaryElementType;
  typedef typename DataTypeContainer::RealType                          RealType;
  typedef typename DataTypeContainer::RealVecChart                      RealVecChart;
  typedef typename DataTypeContainer::IntVecChart                       IntVecChart;
  typedef typename DataTypeContainer::RealVecChartBoundary              RealVecChartBoundary;
  typedef typename DataTypeContainer::PointType                         PointType;
  typedef typename DataTypeContainer::Matrix22                          Matrix22;
  typedef typename DataTypeContainer::Matrix32                          Matrix32;
  typedef typename DataTypeContainer::Matrix33                          Matrix33;
  typedef typename DataTypeContainer::DerivativeVectorValuedType        DerivativeVectorValuedType;
  //   typedef typename DataTypeContainer::Tensor222Type                     Tensor222Type;
  //   typedef typename DataTypeContainer::Tensor322Type                     Tensor322Type;
  typedef typename DataTypeContainer::VectorType                        VectorType;
  typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
  typedef typename DataTypeContainer::TripletType                       TripletType;
  typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
  typedef typename DataTypeContainer::MaskType                          MaskType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;
  typedef Eigen::Matrix<RealType, sizeVoigtTensor, 1 >                  VoigtTensorVec;
  
//   typedef QuocBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  typedef QuocCachedBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType>   BaseFuncSetType;
  
  
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
  const QuadType _quad;
  const BaseFuncSetType _baseFuncSet;
  
public:
  
  QuocConfigurator2D ( const InitType &Mesh ) : 
  _mesh ( Mesh ),
  _baseFuncSet( _mesh.getMeshSize(0), _mesh.getMeshSize(1) ) {}

  const InitType& getMesh( ) const { return this->_mesh; }

  inline int getNumLocalDofs ( const ElementType & ) const { return maxNumLocalDofs;}
  inline int getNumLocalDofs (  ) const { return maxNumLocalDofs;}
  inline int getNumLocalBoundaryDofs (  ) const { return maxNumLocalBoundaryDofs;}

  int getNumGlobalDofs( ) const {return this->_mesh.getNumVertices();}
  int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
  int maxNumBoundaryQuadPoints() const {return BoundaryQuadType::numQuadPoints;}
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &/*T*/ ) const {return _baseFuncSet;}
  const BaseFuncSetType& getBaseFunctionSet ( ) const {return _baseFuncSet;}
  
  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &El, int localIndex ) const {return El.getGlobalNodeIdx( localIndex );}
  
  void getLocalCoords ( const PointType &GlobalCoord, int &elementNumber, PointType &LocalCoord ) const {
      PointType BasePoint;
      std::vector<int> BasePointFac ( 2 );
      for ( int comp = 0; comp < 2; ++comp ) {
          const RealType sc = GlobalCoord[comp] / _mesh.getMeshSize(comp);
          BasePointFac[comp] = static_cast<int> ( sc );
          BasePoint[comp] = BasePointFac[comp] * _mesh.getMeshSize(comp);
          LocalCoord[comp] = sc - BasePointFac[comp];
          //special case: "right" boundary
          if( std::abs( _mesh.getWidth( comp ) - GlobalCoord[comp] ) < 2.e-16 ){
              BasePoint[comp] = _mesh.getWidth(comp) - _mesh.getMeshSize(comp);
              BasePointFac[comp]  = _mesh.getNumDofs( comp ) - 2.;
              LocalCoord[comp] = 1.;
          }
      }
      
      elementNumber = _mesh.getGlobalElementIndex( BasePointFac[0], BasePointFac[1] );
      
      // The position is not in the domain.
      // //     if ( ( GlobalCoord[c] < 0. ) || ( El[c] < 0 ) || ( El[c] >= ( (Grid.getSize())[c] - 1 ) ) )  return false;
  }
  
  void getGlobalCoordsForBoundaryElement ( const BoundaryElementType &El, const PointType &LocalCoord, PointType &Coord ) const {
      Coord = LocalCoord; 
      Coord[0] *= _mesh.getMeshSize(0); Coord[1] *= _mesh.getMeshSize(1);
      Coord += El.getElement().getNode(0); //Coord Left Down
  }
  
  void getNodalValuesFromCoefficients ( const VectorType& coefficients , VectorType& nodalValues ) const {
      nodalValues = coefficients;
  }

};






template <typename DataTypeContainer = Quoc3DDataTypeContainer, 
          typename MeshType = QuocMesh3D<DataTypeContainer>,
          typename QuadType = SimpsonQuadrature3D<typename DataTypeContainer::RealType, typename DataTypeContainer::RealVecChart, typename DataTypeContainer::RealVecChart1D>,
          typename BdrQuadType = SimpsonQuadrature2D<typename DataTypeContainer::RealType, typename DataTypeContainer::RealVecChartBoundary, typename DataTypeContainer::RealVecChart1D> >
class QuocConfigurator3D {
public:
  static const int dimChartDomain = 3;
  static const int maxNumLocalDofs = 8;
  static const int maxNumLocalBoundaryDofs = 4;
  static const int numAffineSymGradDofs = 6;
  static const int sizeVoigtTensor = 21;

  typedef DataTypeContainer                                             DTContainer;
  typedef QuadType                                                      QuadRuleType;
  typedef BdrQuadType                                                   BoundaryQuadType;
  typedef MeshType                                                      InitType;        
  typedef typename MeshType::ElementType                                ElementType;
  typedef typename MeshType::BoundaryElementType                        BoundaryElementType;
  typedef typename DataTypeContainer::RealType                          RealType;
  typedef typename DataTypeContainer::RealVecChart                      RealVecChart;
  typedef typename DataTypeContainer::IntVecChart                       IntVecChart;
  typedef typename DataTypeContainer::RealVecChartBoundary              RealVecChartBoundary;
  typedef typename DataTypeContainer::PointType                         PointType;
  typedef typename DataTypeContainer::Matrix22                          Matrix22;
  typedef typename DataTypeContainer::Matrix32                          Matrix32;
  typedef typename DataTypeContainer::Matrix33                          Matrix33;
  typedef typename DataTypeContainer::DerivativeVectorValuedType        DerivativeVectorValuedType;
  typedef typename DataTypeContainer::VectorType                        VectorType;
  typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
  typedef typename DataTypeContainer::TripletType                       TripletType;
  typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
  typedef typename DataTypeContainer::MaskType                          MaskType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;
  typedef Eigen::Matrix<RealType, sizeVoigtTensor, 1 >                  VoigtTensorVec;
  
//   typedef QuocBaseFunctionSet3D<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  typedef QuocCachedBaseFunctionSet3D<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  
  //! globally affine functions with symmetric gradient
  typedef Eigen::Matrix<RealType, numAffineSymGradDofs, numAffineSymGradDofs >    LocalMatrixTypeAffineSymGrad;
  typedef Eigen::Matrix<RealType, numAffineSymGradDofs, maxNumLocalDofs >         LocalMatrixTypeMixedAffineSymGrad;
  typedef GlobalAffineSymGradBaseFunctionSet3D<DataTypeContainer>             GlobalAffineSymGradBaseFuncSet;
  
  //! globally affine functions with constant gradient: \f$ u(x) = a \cdot x \f$
  typedef Eigen::Matrix<RealType, dimChartDomain, dimChartDomain >               LocalMatrixTypeAffineGrad;
  typedef Eigen::Matrix<RealType, dimChartDomain, maxNumLocalDofs >         LocalMatrixTypeMixedAffineGrad;
  typedef GlobalAffineGradBaseFunctionSet3D<DataTypeContainer>     GlobalAffineGradBaseFuncSet;

protected:  
  const MeshType &_mesh;
  const QuadType _quad;
  mutable BaseFuncSetType _baseFuncSet;
  
public:
  
  QuocConfigurator3D ( const InitType &Mesh ) :
  _mesh ( Mesh ), 
  _baseFuncSet( _mesh.getMeshSize(0), _mesh.getMeshSize(1), _mesh.getMeshSize(2) )
  {}

  const InitType& getMesh( ) const { return this->_mesh; }

  inline int getNumLocalDofs ( const ElementType & ) const { return maxNumLocalDofs;}
  inline int getNumLocalDofs ( ) const { return maxNumLocalDofs;}
  inline int getNumLocalBoundaryDofs (  ) const { return maxNumLocalBoundaryDofs;}
  
  int getNumGlobalDofs( ) const {return this->_mesh.getNumVertices();}
  int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
  int maxNumBoundaryQuadPoints() const {return BoundaryQuadType::numQuadPoints;}
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &/*T*/ ) const {return _baseFuncSet;}
  const BaseFuncSetType& getBaseFunctionSet ( ) const {return _baseFuncSet;}
  
  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int localIndex ) const {return T.getGlobalNodeIdx( localIndex );}
  
  //TODO option either left, right, inside 
  // divide globalCoord = BasePoint + h * LocalCoord
  // BasePoint = ceil( globalCoord / h ) * h
  // LocalCoord = globalCoord / length - ceil( globalCoord / h )
  // BasePointFac = ceil( globalCoord / h )
  //TODO also for 2d case
  void getLocalCoords ( const PointType &GlobalCoord, int &elementNumber, PointType &LocalCoord ) const {
      PointType BasePoint;
      std::vector<int> BasePointFac ( 3 );
      
      for ( int comp = 0; comp < 3; ++comp ) {
          
         bool leftBoundary = false, rightBoundary = false, outside = false;
         if( GlobalCoord[comp] < 2.e-16 ) { leftBoundary = true; outside = true; }
         if( GlobalCoord[comp] > (_mesh.getWidth( comp ) - 2.e-16) ) { rightBoundary = true; outside = true;}
         
         if( outside ){
            if( leftBoundary ){
                BasePoint[comp] = 0.;
                BasePointFac[comp]  = 0.;
                LocalCoord[comp] = 0.;   
            }
            if( rightBoundary ){
                BasePoint[comp] = _mesh.getWidth(comp) - _mesh.getMeshSize(comp);
                BasePointFac[comp]  = _mesh.getNumDofs( comp ) - 2.;
                LocalCoord[comp] = 1.;   
            }
          }else{
              const RealType relativeCoord = GlobalCoord[comp] / _mesh.getMeshSize(comp);
              BasePointFac[comp] = static_cast<int> ( relativeCoord );
              BasePoint[comp] = BasePointFac[comp] * _mesh.getMeshSize(comp);
              LocalCoord[comp] = relativeCoord - BasePointFac[comp];
          }
      }
      
      elementNumber = _mesh.getGlobalElementIndex( BasePointFac[0], BasePointFac[1], BasePointFac[2] );
      
  }
  
  void getGlobalCoordsForBoundaryElement ( const BoundaryElementType &El, const PointType &LocalCoord, PointType &Coord ) const {
      Coord = LocalCoord; 
      Coord[0] *= _mesh.getMeshSize(0); Coord[1] *= _mesh.getMeshSize(1); Coord[2] *= _mesh.getMeshSize(2);
      Coord += El.getElement().getNode(0); //Coord Left Down
  }
  
  void getNodalValuesFromCoefficients ( const VectorType& coefficients , VectorType& nodalValues ) const {
      nodalValues = coefficients;
  }
  
};


#endif //__QUOCMESHCONFIGURATORS_H
