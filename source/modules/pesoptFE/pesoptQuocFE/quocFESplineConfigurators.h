#ifndef __QUOCFESPLINECONFIGURATORS_H
#define __QUOCFESPLINECONFIGURATORS_H

# include <linearSystemSolver.h>
# include "quocFESplineBaseFunctionSets.h"


static const unsigned short int ElementsInSupportOf1DSpline = 4;

static unsigned short int determine1DSplineElementNumber ( const int Index,  const int numDofs ) {
    if ( Index <= 0 ) return 0;
    if ( Index >= numDofs - 2 ) return 2;
    return 1;
}

static int localToGlobal1DSpline ( const int BasePoint, const int localIndex,  const int numDofs) {
      int globalIndex = BasePoint;
      if ( BasePoint > 0 ) {
          if ( BasePoint < numDofs - 2 ) globalIndex -= 1;
          else                       globalIndex -= 2;
      }
      return globalIndex + localIndex;
}



template <typename DataTypeContainer, 
//           = Quoc1DDataTypeContainer,
          typename MeshType, 
//             = QuocMesh1D<DataTypeContainer>,
          typename QuadType 
//           typename QuadType = SimpsonQuadrature1D<typename DataTypeContainer::RealType> 
          >
class QuocConfiguratorBSpline1D {
    
public:
    static const int dimChartDomain = 1;
    static const int maxNumLocalDofs = 4;
//     static const int maxNumLocalBoundaryDofs = 2;
    static const int numAffineSymGradDofs = 1;
    static const int sizeVoigtTensor = 1;

  typedef DataTypeContainer                                             DTContainer;
  typedef QuadType                                                      QuadRuleType;
//   typedef BdrQuadType                                                   BoundaryQuadType;
  typedef MeshType                                                      InitType;           
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
  
  typedef QuocSplineBaseFunctionSet1D<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;

protected:  
  const MeshType &_mesh;
  const int _Nx;
  const QuadType _quad;
  std::vector < BaseFuncSetType > _baseFuncSet;

public:
  explicit QuocConfiguratorBSpline1D ( const InitType& mesh )
  : _mesh ( mesh ), _Nx ( _mesh.getNumDofs( 0 ) ), 
    _baseFuncSet ( ) {
        for ( unsigned short int i = 0; i < BaseFuncSetType::numberOfDifferentBaseFunctionSets; ++i )
            _baseFuncSet.emplace_back ( _mesh.getMeshSize(0), i );
     }
 
 
  const InitType& getMesh( ) const { return this->_mesh; }

  inline int getNumLocalDofs ( const ElementType & ) const { return maxNumLocalDofs;}
  inline int getNumLocalDofs (  ) const { return maxNumLocalDofs;}

  int getNumGlobalDofs( ) const {return this->_mesh.getNumVertices();}
  int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
  
  
  
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
  
  //! returns global index of the dof with number localIndex
  // Warning: some elements appear twice due to localToGlobal
  // In most cases, the right "elementNumber" is obtained by "getConsecutiveElementNumber"
  inline int localToGlobal ( const ElementType &El, const int localIndex ) const {
    return localToGlobal1DSpline ( El.getGlobalNodeIdx( 0 ), localIndex, _Nx );
  }
 
  unsigned short int determineSplineElementNumber ( const ElementType& El ) const {
    return determine1DSplineElementNumber ( El.getGlobalNodeIdx( 0 ), _Nx );
  }
  
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType& El ) const {
    return _baseFuncSet[determineSplineElementNumber ( El )];
  }

  void convertNodalValuesToSplineCoefficients ( const VectorType &nodalValues, 
                                                VectorType &splineCoefficients ) const{
    std::vector<TripletType> tripletList;
    
    //Dirichlet nodes
    tripletList.push_back( TripletType( 0,  0,  1.0 ) );
    tripletList.push_back( TripletType( _Nx-1, _Nx-1, 1. ) );
            
    //
    for ( int i = 1; i < _Nx-1; ++i ) {
      tripletList.push_back( TripletType( i, i-1, 1./6. ) );
      tripletList.push_back( TripletType( i, i, 4./6. ) );
      tripletList.push_back( TripletType( i, i+1, 1./6. ) );
    }
         
    SparseMatrixType SystemMatrix ( _Nx, _Nx );
    SystemMatrix.setFromTriplets( tripletList.begin(), tripletList.end() );
    
//     IterativeLinearSystemSolver<DataTypeContainer,iterativeLSMethod> iterativeLinearSystemSolver(  10. * Eigen::NumTraits<RealType>::epsilon(),  1. );
//     iterativeLinearSystemSolver.solve( SystemMatrix,  splineCoefficients,  nodalValues );

#ifdef PESOPT_WITH_SUITESPARSE
    DirectLinearSystemSolver<DataTypeContainer,EigenUmfPackLU > directLinearSystemSolver;
#else 
   DirectLinearSystemSolver<DataTypeContainer,EigenSparseLU > directLinearSystemSolver;
#endif
    directLinearSystemSolver.prepareSolver( SystemMatrix );
    directLinearSystemSolver.solve( splineCoefficients,  nodalValues );

  }
  
  
  void getNodalValuesFromCoefficients ( const VectorType& splineCoefficients , VectorType& nodalValues ) const {
    nodalValues.setZero ( );
    for ( int nodeIdx = 0; nodeIdx < _mesh.getNumVertices(); ++nodeIdx ) {
        auto GlobalCoord = _mesh.getVertex(nodeIdx);
        RealVecChart LocalCoord;
        int elementNumber;
        this->getLocalCoords ( GlobalCoord, elementNumber, LocalCoord );
        auto El = _mesh.getElement(elementNumber);
        const BaseFuncSetType &bfs = this->getBaseFunctionSet ( El );
        for ( unsigned short int b = 0; b < this->getNumLocalDofs ( ); ++b )
            nodalValues[nodeIdx] += splineCoefficients[ this->localToGlobal ( El, b ) ] * bfs.evaluate ( b, LocalCoord );
    }
    
  }
  
};




template <typename DataTypeContainer, 
//           = Quoc2DDataTypeContainer,
          typename MeshType, 
//             = QuocMesh2D<DataTypeContainer>,
          typename QuadType 
//             = SimpsonQuadrature2D<typename DataTypeContainer::RealType, typename DataTypeContainer::RealVecChart, typename DataTypeContainer::RealVecChart1D>,
//           typename BdrQuadType = SimpsonQuadrature1D<typename DataTypeContainer::RealType, typename DataTypeContainer::RealType> 
          >
class QuocConfiguratorBSpline2D {
    
public:
    static const int dimChartDomain = 2;
    static const int maxNumLocalDofs = 16;
    static const int maxNumLocalBoundaryDofs = 4;
    static const int numAffineSymGradDofs = 3;
    static const int sizeVoigtTensor = 6;


  typedef DataTypeContainer                                             DTContainer;
  typedef QuadType                                                      QuadRuleType;
//   typedef BdrQuadType                                                   BoundaryQuadType;
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
  
  typedef QuocSplineBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  
  static const int numberOfDifferentBaseFunctionSets1D = BaseFuncSetType::BaseFunc1DType::numberOfDifferentBaseFunctionSets;

protected:  
  const MeshType &_mesh;
  const int _Nx,  _Ny;
  const QuadType _quad;
  std::vector < BaseFuncSetType > _baseFuncSet;

public:
  explicit QuocConfiguratorBSpline2D ( const InitType& mesh )
  : _mesh ( mesh ), 
    _Nx ( _mesh.getNumDofs( 0 ) ), _Ny ( _mesh.getNumDofs( 1 ) ), 
    _baseFuncSet ( ) {
//         for ( unsigned short int i = 0; i < numberOfDifferentBaseFunctionSets1D; ++i )
//             for ( unsigned short int j = 0; j < numberOfDifferentBaseFunctionSets1D; ++j )
//                _baseFuncSet.emplace_back ( _mesh.getMeshSize(0), _mesh.getMeshSize(1), i,  j);
        for ( unsigned short int i = 0; i < numberOfDifferentBaseFunctionSets1D; ++i )
            for ( unsigned short int j = 0; j < numberOfDifferentBaseFunctionSets1D; ++j )
               _baseFuncSet.emplace_back ( _mesh.getMeshSize(0), _mesh.getMeshSize(1), j,  i);
     }
     
     
  const InitType& getMesh( ) const { return this->_mesh; }

  inline int getNumLocalDofs ( const ElementType & ) const { return maxNumLocalDofs;}
  inline int getNumLocalDofs (  ) const { return maxNumLocalDofs;}
  inline int getNumLocalBoundaryDofs (  ) const { return maxNumLocalBoundaryDofs;}

  int getNumGlobalDofs( ) const {return this->_mesh.getNumVertices();}
  int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
//   int maxNumBoundaryQuadPoints() const {return BoundaryQuadType::numQuadPoints;}

  
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
  
  //! returns global index of the dof with number localIndex
  // Warning: some elements appear twice due to localToGlobal
  // In most cases, the right "elementNumber" is obtained by "getConsecutiveElementNumber"
  inline int localToGlobal ( const ElementType &El, const int localIndex ) const {
    const int localIndexX = localIndex % ElementsInSupportOf1DSpline;
    const int localIndexY = localIndex / ElementsInSupportOf1DSpline;
    auto NodeIndices = _mesh.getNodeIndices ( El.getGlobalNodeIdx( 0 ) );
    return localToGlobal1DSpline ( NodeIndices[0], localIndexX,  _Nx) + _Nx * localToGlobal1DSpline( NodeIndices[1], localIndexY,  _Ny );
  }
  
  

  unsigned short int determineSplineElementNumber ( const ElementType& El ) const {
    auto NodeIndices = _mesh.getNodeIndices ( El.getGlobalNodeIdx( 0 ) );
    return determine1DSplineElementNumber( NodeIndices[1],  _Ny) * numberOfDifferentBaseFunctionSets1D + determine1DSplineElementNumber ( NodeIndices[0],  _Nx );
  }
  
  
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType& El ) const {
    return _baseFuncSet[determineSplineElementNumber ( El )];
  }
  
  
  void convertNodalValuesToSplineCoefficients ( const VectorType& nodalValues, VectorType& splineCoefficients ) const {

    std::vector<TripletType> tripletList;

    //! Dirichlet Nodes   
    tripletList.push_back( TripletType( _mesh.getGlobalNodeIndex (0,0),         _mesh.getGlobalNodeIndex (0,0),         1. ));
    tripletList.push_back( TripletType( _mesh.getGlobalNodeIndex (0,_Ny-1),     _mesh.getGlobalNodeIndex (0,_Ny-1),     1. ));
    tripletList.push_back( TripletType( _mesh.getGlobalNodeIndex (_Nx-1,0),     _mesh.getGlobalNodeIndex (_Nx-1,0),     1. ));
    tripletList.push_back( TripletType( _mesh.getGlobalNodeIndex (_Nx-1,_Ny-1), _mesh.getGlobalNodeIndex (_Nx-1,_Ny-1), 1. ));
    
    for ( int i = 1; i < _Nx-1; ++i )
      for ( int j = 1; j < _Ny-1; ++j ) {
          const int currentIndex = _mesh.getGlobalNodeIndex (i,j);

          tripletList.push_back( TripletType( currentIndex, currentIndex, 16./36. ) );

          tripletList.push_back( TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j) , 4./36. ) );
          tripletList.push_back( TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j) , 4./36. ));
          tripletList.push_back( TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,j-1) , 4./36. ));
          tripletList.push_back( TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,j+1) , 4./36. ));

          tripletList.push_back( TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j-1) , 1./36. ));
          tripletList.push_back( TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j-1) , 1./36. ));
          tripletList.push_back( TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j+1) , 1./36. ));
          tripletList.push_back( TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j+1) , 1./36. ));
      }

    for ( int i = 1; i < _Nx-1; ++i ){
        const int currentIndex = _mesh.getGlobalNodeIndex (i,0);
        tripletList.push_back( TripletType( currentIndex, currentIndex , 24./36. ));
        tripletList.push_back( TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,0) , 6./36. ));
        tripletList.push_back( TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,0) , 6./36. ));

        const int currentIndexEnd = _mesh.getGlobalNodeIndex (i,_Ny-1);
        tripletList.push_back( TripletType( currentIndexEnd, currentIndexEnd , 24./36. ));
        tripletList.push_back( TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i-1,_Ny-1) , 6./36. ));
        tripletList.push_back( TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i+1,_Ny-1) , 6./36. ));
    }

    for ( int j = 1; j < _Ny-1; ++j ){
        const int currentIndex = _mesh.getGlobalNodeIndex (0,j);
        tripletList.push_back( TripletType( currentIndex, currentIndex , 24./36. ));
        tripletList.push_back( TripletType( currentIndex, _mesh.getGlobalNodeIndex (0,j-1) , 6./36. ));
        tripletList.push_back( TripletType( currentIndex, _mesh.getGlobalNodeIndex (0,j+1) , 6./36. ));

        const int currentIndexEnd = _mesh.getGlobalNodeIndex (_Nx-1,j);
        tripletList.push_back( TripletType( currentIndexEnd, currentIndexEnd , 24./36. ));
        tripletList.push_back( TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (_Nx-1,j-1) , 6./36. ));
        tripletList.push_back( TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (_Nx-1,j+1) , 6./36. ));
    }
      
      
    SparseMatrixType SystemMatrix ( this->getNumGlobalDofs( ), this->getNumGlobalDofs( ) );
    SystemMatrix.setFromTriplets( tripletList.begin(), tripletList.end() );
    
    
//     IterativeLinearSystemSolver<DataTypeContainer,iterativeLSMethod> iterativeLinearSystemSolver (  10. * Eigen::NumTraits<RealType>::epsilon(),  1. );
//     iterativeLinearSystemSolver.solve( SystemMatrix,  splineCoefficients,  nodalValues );

    
#ifdef PESOPT_WITH_SUITESPARSE
    DirectLinearSystemSolver<DataTypeContainer,EigenUmfPackLU > directLinearSystemSolver;
#else 
   DirectLinearSystemSolver<DataTypeContainer,EigenSparseLU > directLinearSystemSolver;
#endif
    directLinearSystemSolver.prepareSolver( SystemMatrix );
    directLinearSystemSolver.solve( splineCoefficients,  nodalValues );
  }


  void getNodalValuesFromCoefficients ( const VectorType& splineCoefficients , VectorType& nodalValues ) const {
    nodalValues.setZero ( );
    for ( int nodeIdx = 0; nodeIdx < _mesh.getNumVertices(); ++nodeIdx ) {
        auto GlobalCoord = _mesh.getVertex(nodeIdx);
        RealVecChart LocalCoord;
        int elementNumber;
        this->getLocalCoords ( GlobalCoord, elementNumber, LocalCoord );
        auto El = _mesh.getElement(elementNumber);
        const BaseFuncSetType &bfs = this->getBaseFunctionSet ( El );
        for ( unsigned short int b = 0; b < this->getNumLocalDofs ( ); ++b )
            nodalValues[nodeIdx] += splineCoefficients[ this->localToGlobal ( El, b ) ] * bfs.evaluate ( b, LocalCoord );
    }
    
  }
  
};



template <typename DataTypeContainer, 
//          = Quoc3DDataTypeContainer,
          typename MeshType, 
//             = QuocMesh3D<DataTypeContainer>,
          typename QuadType 
//             = SimpsonQuadrature3D<typename DataTypeContainer::RealType, typename DataTypeContainer::RealVecChart, typename DataTypeContainer::RealVecChart1D>,
//           typename BdrQuadType = SimpsonQuadrature2D<typename DataTypeContainer::RealType, typename DataTypeContainer::RealVecChartBoundary, typename DataTypeContainer::RealVecChart1D> 
          >
class QuocConfiguratorBSpline3D {
    
public:
    static const int dimChartDomain = 3;
    static const int maxNumLocalDofs = 64;
    static const int maxNumLocalBoundaryDofs = 16;
    static const int numAffineSymGradDofs = 6;
    static const int sizeVoigtTensor = 21;


  typedef DataTypeContainer                                             DTContainer;
  typedef QuadType                                                      QuadRuleType;
//   typedef BdrQuadType                                                   BoundaryQuadType;
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
  
  typedef QuocSplineBaseFunctionSet3D<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  
  static const int numberOfDifferentBaseFunctionSets1D = BaseFuncSetType::BaseFunc1DType::numberOfDifferentBaseFunctionSets;

protected:  
  const MeshType &_mesh;
  const int _Nx, _Ny, _Nz;
  const QuadType _quad;
  std::vector < BaseFuncSetType > _baseFuncSet;

public:
  explicit QuocConfiguratorBSpline3D ( const InitType& mesh )
  : _mesh ( mesh ), 
    _Nx ( _mesh.getNumDofs(0) ), _Ny(  _mesh.getNumDofs(1) ), _Nz( _mesh.getNumDofs(2) ), 
    _baseFuncSet ( ) {
//         for ( unsigned short int i = 0; i < numberOfDifferentBaseFunctionSets1D; ++i )
//             for ( unsigned short int j = 0; j < numberOfDifferentBaseFunctionSets1D; ++j )
//                 for ( unsigned short int k = 0; k < numberOfDifferentBaseFunctionSets1D; ++k )
//                     _baseFuncSet.emplace_back ( _mesh.getMeshSize(0), _mesh.getMeshSize(1), _mesh.getMeshSize(2), i, j, k);
        for ( unsigned short int i = 0; i < numberOfDifferentBaseFunctionSets1D; ++i )
            for ( unsigned short int j = 0; j < numberOfDifferentBaseFunctionSets1D; ++j )
                for ( unsigned short int k = 0; k < numberOfDifferentBaseFunctionSets1D; ++k )
                    _baseFuncSet.emplace_back ( _mesh.getMeshSize(0), _mesh.getMeshSize(1), _mesh.getMeshSize(2), k, j, i);
     }
     
     
  const InitType& getMesh( ) const { return this->_mesh; }

  inline int getNumLocalDofs ( const ElementType & ) const { return maxNumLocalDofs;}
  inline int getNumLocalDofs (  ) const { return maxNumLocalDofs;}
  inline int getNumLocalBoundaryDofs (  ) const { return maxNumLocalBoundaryDofs;}

  int getNumGlobalDofs( ) const {return this->_mesh.getNumVertices();}
  int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
//   int maxNumBoundaryQuadPoints() const {return BoundaryQuadType::numQuadPoints;}

  
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
  
  
  
    //! returns global index of the dof with number localIndex
  // Warning: some elements appear twice due to localToGlobal
  // In most cases, the right "elementNumber" is obtained by "getConsecutiveElementNumber"
  inline int localToGlobal ( const ElementType &El, const int localIndex ) const {
    const int localIndexX = localIndex % ElementsInSupportOf1DSpline;
    const int localIndexY = ( localIndex % pesopt::Sqr ( ElementsInSupportOf1DSpline ) ) / ElementsInSupportOf1DSpline;
    const int localIndexZ = localIndex / pesopt::Sqr ( ElementsInSupportOf1DSpline );
    auto NodeIndices = _mesh.getNodeIndices ( El.getGlobalNodeIdx( 0 ) );
    return localToGlobal1DSpline( NodeIndices[0], localIndexX, _Nx ) + _Nx * localToGlobal1DSpline( NodeIndices[1], localIndexY,  _Ny) + _Nx * _Ny * localToGlobal1DSpline( NodeIndices[2], localIndexZ,  _Nz );
  }

  unsigned short int determineSplineElementNumber ( const ElementType& El ) const {
    auto NodeIndices = _mesh.getNodeIndices ( El.getGlobalNodeIdx( 0 ) );
    return determine1DSplineElementNumber ( NodeIndices[2],  _Nz ) * pesopt::Sqr ( numberOfDifferentBaseFunctionSets1D ) + determine1DSplineElementNumber ( NodeIndices[1],  _Ny ) * numberOfDifferentBaseFunctionSets1D + determine1DSplineElementNumber( NodeIndices[0], _Nx );
  }
  
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType& El ) const {
    return _baseFuncSet[determineSplineElementNumber ( El )];
  }
  

  void convertNodalValuesToSplineCoefficients ( const VectorType& nodalValues, VectorType& splineCoefficients ) const {

    std::vector<TripletType> tripletList;
    
    const int Nx = _mesh.getNumDofs(0); 
    const int Ny = _mesh.getNumDofs(1);
    const int Nz = _mesh.getNumDofs(2);

    //! Dirichlet Nodes   
    tripletList.push_back ( TripletType( _mesh.getGlobalNodeIndex ( 0, 0, 0 ), _mesh.getGlobalNodeIndex ( 0, 0, 0 ), 1. ) );
    tripletList.push_back ( TripletType( _mesh.getGlobalNodeIndex ( 0, Ny-1, 0 ), _mesh.getGlobalNodeIndex ( 0, Ny-1, 0 ), 1. ) );
    tripletList.push_back ( TripletType( _mesh.getGlobalNodeIndex ( Nx-1, 0, 0 ), _mesh.getGlobalNodeIndex ( Nx-1, 0, 0 ), 1. ) );
    tripletList.push_back ( TripletType( _mesh.getGlobalNodeIndex ( Nx-1, Ny-1, 0 ), _mesh.getGlobalNodeIndex ( Nx-1, Ny-1, 0 ), 1. ) );
    tripletList.push_back ( TripletType( _mesh.getGlobalNodeIndex ( 0, 0, Nz-1 ), _mesh.getGlobalNodeIndex ( 0, 0, Nz-1 ), 1. ) );
    tripletList.push_back ( TripletType( _mesh.getGlobalNodeIndex ( 0, Ny-1, Nz-1 ), _mesh.getGlobalNodeIndex ( 0, Ny-1, Nz-1 ), 1. ) );
    tripletList.push_back ( TripletType( _mesh.getGlobalNodeIndex ( Nx-1, 0, Nz-1 ), _mesh.getGlobalNodeIndex ( Nx-1, 0, Nz-1 ), 1. ) );
    tripletList.push_back ( TripletType( _mesh.getGlobalNodeIndex ( Nx-1, Ny-1, Nz-1 ), _mesh.getGlobalNodeIndex ( Nx-1, Ny-1, Nz-1 ), 1. ) );
    
    // !
    for ( int i = 1; i < Nx - 1; ++i ) {
      for ( int j = 1; j < Ny - 1; ++j ) {
        for ( int k = 1; k < Nz - 1; ++k ) {

          const int currentIndex = _mesh.getGlobalNodeIndex  ( i, j, k );

          tripletList.push_back (  TripletType( currentIndex, currentIndex, 64./216. ) );

          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j,k), 16./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j,k), 16./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,j-1,k), 16./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,j+1,k), 16./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,j,k-1), 16./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,j,k+1), 16./216. ) );

          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j-1,k), 4./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j-1,k), 4./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j+1,k), 4./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j+1,k), 4./216. ) );

          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j,k-1), 4./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j,k-1), 4./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j,k+1), 4./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j,k+1), 4./216. ) );

          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,j-1,k-1), 4./216.) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,j+1,k-1), 4./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,j-1,k+1), 4./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,j+1,k+1), 4./216. ) );

          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j-1,k-1), 1./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j-1,k-1), 1./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j+1,k-1), 1./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j+1,k-1), 1./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j-1,k+1), 1./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j-1,k+1), 1./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j+1,k+1), 1./216. ) );
          tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j+1,k+1), 1./216. ) );
        }
      }
    }
    
    for ( int i = 1; i < Nx - 1; ++i ) {
      for ( int j = 1; j < Ny - 1; ++j ) {
        const int currentIndex = _mesh.getGlobalNodeIndex (i,j,0);
        tripletList.push_back (  TripletType( currentIndex, currentIndex, 16./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j,0), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j,0), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,j-1,0), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,j+1,0), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j-1,0), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j-1,0), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,j+1,0), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,j+1,0), 1./36. ) );

        const int currentIndexEnd = _mesh.getGlobalNodeIndex (i,j,Nz-1);
        tripletList.push_back (  TripletType( currentIndexEnd, currentIndexEnd, 16./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i-1,j,Nz-1), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i+1,j,Nz-1), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i,j-1,Nz-1), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i,j+1,Nz-1), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i-1,j-1,Nz-1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i+1,j-1,Nz-1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i-1,j+1,Nz-1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i+1,j+1,Nz-1), 1./36. ) );

      }
    }
    for ( int i = 1; i < Nx - 1; ++i ) {
      for ( int k = 1; k < Nz - 1; ++k ) {
        const int currentIndex = _mesh.getGlobalNodeIndex (i,0,k);
        tripletList.push_back (  TripletType( currentIndex, currentIndex, 16./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,0,k), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,0,k), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,0,k-1), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i,0,k+1), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,0,k-1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,0,k-1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i-1,0,k+1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (i+1,0,k+1), 1./36. ) );

        const int currentIndexEnd = _mesh.getGlobalNodeIndex (i,Ny-1,k);
        tripletList.push_back (  TripletType( currentIndexEnd, currentIndexEnd, 16./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i-1,Ny-1,k), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i+1,Ny-1,k), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i,Ny-1,k-1), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i,Ny-1,k+1), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i-1,Ny-1,k-1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i+1,Ny-1,k-1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i-1,Ny-1,k+1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (i+1,Ny-1,k+1), 1./36. ) );
      }
    }
    for ( int j = 1; j < Ny - 1; ++j ) {
      for ( int k = 1; k < Nz - 1; ++k ) {
        const int currentIndex = _mesh.getGlobalNodeIndex (0,j,k);
        tripletList.push_back (  TripletType( currentIndex, currentIndex, 16./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (0,j-1,k), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (0,j+1,k), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (0,j,k-1), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (0,j,k+1), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (0,j-1,k-1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (0,j+1,k-1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (0,j-1,k+1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndex, _mesh.getGlobalNodeIndex (0,j+1,k+1), 1./36. ) );

        const int currentIndexEnd = _mesh.getGlobalNodeIndex (Nx-1,j,k);
        tripletList.push_back (  TripletType( currentIndexEnd, currentIndexEnd, 16./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (Nx-1,j-1,k), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (Nx-1,j+1,k), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (Nx-1,j,k-1), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (Nx-1,j,k+1), 4./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (Nx-1,j-1,k-1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (Nx-1,j+1,k-1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (Nx-1,j-1,k+1), 1./36. ) );
        tripletList.push_back (  TripletType( currentIndexEnd, _mesh.getGlobalNodeIndex (Nx-1,j+1,k+1), 1./36. ) );
      }
    }

    for ( int i = 1; i < Nx - 1; ++i ) {
        const int index00 = _mesh.getGlobalNodeIndex (i,0,0);
        tripletList.push_back (  TripletType( index00, _mesh.getGlobalNodeIndex (i,0,0), 4./6. ) );
        tripletList.push_back (  TripletType( index00, _mesh.getGlobalNodeIndex (i-1,0,0), 1./6. ) );
        tripletList.push_back (  TripletType( index00, _mesh.getGlobalNodeIndex (i+1,0,0), 1./6. ) );
          
        const int index01 = _mesh.getGlobalNodeIndex (i,0,Nz-1);
        tripletList.push_back (  TripletType( index01, _mesh.getGlobalNodeIndex (i,0,Nz-1), 4./6. ) );
        tripletList.push_back (  TripletType( index01, _mesh.getGlobalNodeIndex (i-1,0,Nz-1), 1./6. ) );
        tripletList.push_back (  TripletType( index01, _mesh.getGlobalNodeIndex (i+1,0,Nz-1), 1./6. ) );
          
        const int index10 = _mesh.getGlobalNodeIndex (i,Ny-1,0);
        tripletList.push_back (  TripletType( index10, _mesh.getGlobalNodeIndex (i,Ny-1,0), 4./6. ) );  
        tripletList.push_back (  TripletType( index10, _mesh.getGlobalNodeIndex (i-1,Ny-1,0), 1./6. ) );
        tripletList.push_back (  TripletType( index10, _mesh.getGlobalNodeIndex (i+1,Ny-1,0), 1./6. ) );
          
        const int index11 = _mesh.getGlobalNodeIndex (i,Ny-1,Nz-1);
        tripletList.push_back (  TripletType( index11, _mesh.getGlobalNodeIndex (i,Ny-1,Nz-1), 4./6. ) );
        tripletList.push_back (  TripletType( index11, _mesh.getGlobalNodeIndex (i-1,Ny-1,Nz-1), 1./6. ) );
        tripletList.push_back (  TripletType( index11, _mesh.getGlobalNodeIndex (i+1,Ny-1,Nz-1), 1./6. ) );
    }
    for ( int j = 1; j < Ny - 1; ++j ) {
        const int index00 = _mesh.getGlobalNodeIndex (0,j,0);
        tripletList.push_back (  TripletType( index00, _mesh.getGlobalNodeIndex (0,j,0), 4./6. ) );
        tripletList.push_back (  TripletType( index00, _mesh.getGlobalNodeIndex (0,j-1,0), 1./6. ) );
        tripletList.push_back (  TripletType( index00, _mesh.getGlobalNodeIndex (0,j+1,0), 1./6. ) );
          
        const int index01 = _mesh.getGlobalNodeIndex (0,j,Nz-1);
        tripletList.push_back (  TripletType( index01, _mesh.getGlobalNodeIndex (0,j,Nz-1), 4./6. ) );
        tripletList.push_back (  TripletType( index01, _mesh.getGlobalNodeIndex (0,j-1,Nz-1), 1./6. ) );
        tripletList.push_back (  TripletType( index01, _mesh.getGlobalNodeIndex (0,j+1,Nz-1), 1./6. ) );
          
        const int index10 = _mesh.getGlobalNodeIndex (Nx-1,j,0);
        tripletList.push_back (  TripletType( index10, _mesh.getGlobalNodeIndex (Nx-1,j,0), 4./6. ) );
        tripletList.push_back (  TripletType( index10, _mesh.getGlobalNodeIndex (Nx-1,j-1,0), 1./6. ) );
        tripletList.push_back (  TripletType( index10, _mesh.getGlobalNodeIndex (Nx-1,j+1,0), 1./6. ) );
          
        const int index11 = _mesh.getGlobalNodeIndex (Nx-1,j,Nz-1);
        tripletList.push_back (  TripletType( index11, _mesh.getGlobalNodeIndex (Nx-1,j,Nz-1), 4./6. ) );
        tripletList.push_back (  TripletType( index11, _mesh.getGlobalNodeIndex (Nx-1,j-1,Nz-1), 1./6. ) );
        tripletList.push_back (  TripletType( index11, _mesh.getGlobalNodeIndex (Nx-1,j+1,Nz-1), 1./6. ) );
    }
    for ( int k = 1; k < Nz - 1; ++k ) {
        const int index00 = _mesh.getGlobalNodeIndex (0,0,k);
        tripletList.push_back (  TripletType( index00, _mesh.getGlobalNodeIndex (0,0,k), 4./6. ) );
        tripletList.push_back (  TripletType( index00, _mesh.getGlobalNodeIndex (0,0,k-1), 1./6. ) );
        tripletList.push_back (  TripletType( index00, _mesh.getGlobalNodeIndex (0,0,k+1), 1./6. ) );
          
        const int index01 = _mesh.getGlobalNodeIndex (0,Ny-1,k);
        tripletList.push_back (  TripletType( index01, _mesh.getGlobalNodeIndex (0,Ny-1,k), 4./6. ) );
        tripletList.push_back (  TripletType( index01, _mesh.getGlobalNodeIndex (0,Ny-1,k-1), 1./6. ) );
        tripletList.push_back (  TripletType( index01, _mesh.getGlobalNodeIndex (0,Ny-1,k+1), 1./6. ) );
          
        const int index10 = _mesh.getGlobalNodeIndex (Nx-1,0,k);
        tripletList.push_back (  TripletType( index10, _mesh.getGlobalNodeIndex (Nx-1,0,k), 4./6. ) );
        tripletList.push_back (  TripletType( index10, _mesh.getGlobalNodeIndex (Nx-1,0,k-1), 1./6. ) );
        tripletList.push_back (  TripletType( index10, _mesh.getGlobalNodeIndex (Nx-1,0,k+1), 1./6. ) );
          
        const int index11 = _mesh.getGlobalNodeIndex (Nx-1,Ny-1,k);
        tripletList.push_back (  TripletType( index11, _mesh.getGlobalNodeIndex (Nx-1,Ny-1,k), 4./6. ) );
        tripletList.push_back (  TripletType( index11, _mesh.getGlobalNodeIndex (Nx-1,Ny-1,k-1), 1./6. ) );
        tripletList.push_back (  TripletType( index11, _mesh.getGlobalNodeIndex (Nx-1,Ny-1,k+1), 1./6. ) );
    }
    
    SparseMatrixType SystemMatrix ( this->getNumGlobalDofs( ), this->getNumGlobalDofs( ) );
    SystemMatrix.setFromTriplets( tripletList.begin(), tripletList.end() );
    
//     IterativeLinearSystemSolver<DataTypeContainer,iterativeLSMethod> iterativeLinearSystemSolver (  10. * Eigen::NumTraits<RealType>::epsilon(),  1. );
//     iterativeLinearSystemSolver.solve( SystemMatrix,  splineCoefficients,  nodalValues );
    
#ifdef PESOPT_WITH_SUITESPARSE
    DirectLinearSystemSolver<DataTypeContainer,EigenUmfPackLU > directLinearSystemSolver;
#else 
   DirectLinearSystemSolver<DataTypeContainer,EigenSparseLU > directLinearSystemSolver;
#endif
    directLinearSystemSolver.prepareSolver( SystemMatrix );
    directLinearSystemSolver.solve( splineCoefficients,  nodalValues );

  }
  



  void getNodalValuesFromCoefficients ( const VectorType& splineCoefficients , VectorType& nodalValues ) const {
    nodalValues.setZero ( );
    for ( int nodeIdx = 0; nodeIdx < _mesh.getNumVertices(); ++nodeIdx ) {
        auto GlobalCoord = _mesh.getVertex(nodeIdx);
        RealVecChart LocalCoord;
        int elementNumber;
        this->getLocalCoords ( GlobalCoord, elementNumber, LocalCoord );
        auto El = _mesh.getElement(elementNumber);
        const BaseFuncSetType &bfs = this->getBaseFunctionSet ( El );
        for ( unsigned short int b = 0; b < this->getNumLocalDofs ( ); ++b )
            nodalValues[nodeIdx] += splineCoefficients[ this->localToGlobal ( El, b ) ] * bfs.evaluate ( b, LocalCoord );
    }
    
  }


};

#endif
