#ifndef __SHELLCURVATURE_H
#define __SHELLCURVATURE_H

#include <pesopt_DKT.h>

// computes int \sqrt(g) |det(g^-1 h )|
template<typename ConfiguratorType>
class GaussCurvatureL1 :
public RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xStorage;
  public:
    GaussCurvatureL1 ( const ConfiguratorType &conf,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xStorage ) :
      RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xStorage ( xStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gInv = _xStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& h = _xStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType GaussCurv = ( gInv * h ).determinant();
        return _xStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * std::abs( GaussCurv );
    }
};



// computes int \sqrt(g) |det(g_coarse^-1 h_coarse ) - det(g_fine^-1 h_fine )|
template<typename ConfiguratorType>
class GaussCurvatureL1DiffConf :
public RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1DiffConf<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::DomVecType DomVecType;
       typedef typename ConfiguratorType::TangentVecType TangentVecType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
       typedef typename ConfiguratorType::QuadRuleType QuadType;
 public:
      QuadType _quadRule;
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_confCoarse,  &_confFine;
       const DKTFEVectorFunctionEvaluator<ConfiguratorType> _xDFEvaluatorCoarse;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xStorageFine;

       vtkSmartPointer<vtkPolyData> _coarseMeshData;
       vtkSmartPointer<vtkCellLocator> _coarseMeshCellLocator;

  public:
    GaussCurvatureL1DiffConf (
        const ConfiguratorType &confCoarse,
        const VectorType &FEDofsCoarse,
        const ConfiguratorType &confFine,
        const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xStorageFine ) :
      RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1DiffConf<ConfiguratorType>> (confFine),
      _confCoarse( confCoarse ), _confFine ( confFine ),
      _xDFEvaluatorCoarse ( confCoarse,  FEDofsCoarse, 3 ),
      _xStorageFine( xStorageFine ),
      _coarseMeshData( vtkSmartPointer<vtkPolyData>::New() ),
      _coarseMeshCellLocator ( vtkSmartPointer<vtkCellLocator>::New() )
      {

        // save mesh data on fine scale to vtk
        const int numVerticesCoarse = _confCoarse.getInitializer().getNumVertices();
        const int numCellsCoarse = _confCoarse.getInitializer().getNumElements();
        auto points = vtkSmartPointer<vtkPoints>::New();
        for( int nodeIter=0; nodeIter<numVerticesCoarse; ++nodeIter ){
            auto point = _confCoarse.getInitializer().getVertex ( nodeIter );
            if( point.size() == 1 ) points->InsertNextPoint ( point[0], 0., 0. );
            if( point.size() == 2 ) points->InsertNextPoint ( point[0], point[1], 0. );
            if( point.size() == 3 ) points->InsertNextPoint ( point[0], point[1], point[2] );
        }
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        for( int elIdx=0; elIdx<numCellsCoarse; ++elIdx ){
            const int numNodesOfElement = _confCoarse.getInitializer().getNumNodesOfElement(elIdx);
            vtkIdType nodesOfElement [numNodesOfElement];
            for( int i=0; i<numNodesOfElement; ++i ) nodesOfElement[i] = _confCoarse.getInitializer().getElementNodeIdx( elIdx, i );
            cells->InsertNextCell ( numNodesOfElement, nodesOfElement );
        }
        _coarseMeshData->SetPoints ( points );
        _coarseMeshData->SetPolys ( cells );

        // Create the cell locator
        _coarseMeshCellLocator->SetDataSet(_coarseMeshData);
        _coarseMeshCellLocator->BuildLocator();
    }

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        // evaluate on fine element
        auto gInvFine = _xStorageFine.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        auto hFine = _xStorageFine.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType GaussCurvFine = ( gInvFine * hFine ).determinant();
        // find coarseElIdx, coarseRefCoord
        Point3DType globalCoordsChartFineTmp;
        auto refCoordFine = _quadRule.getRefCoord ( QuadPoint );
        El.getGlobalCoord ( refCoordFine, globalCoordsChartFineTmp );
        double globalCoordsChartFine[3] = {globalCoordsChartFineTmp[0], globalCoordsChartFineTmp[1], globalCoordsChartFineTmp[2] };
        double closestPoint[3];
        double closestPointDistSqr;
        vtkIdType closestPointCellIdx;
        int closestPointSubIdx;
        _coarseMeshCellLocator->FindClosestPoint(globalCoordsChartFine, closestPoint, closestPointCellIdx, closestPointSubIdx, closestPointDistSqr);
        int coarseElIdx = closestPointCellIdx;
        Point3DType closestPointTmp;
        closestPointTmp(0) = closestPoint[0];
        closestPointTmp(1) = closestPoint[1];
        closestPointTmp(2) = closestPoint[2];
        DomVecType coarseRefCoord;
        _confCoarse.getInitializer().getTriang(coarseElIdx).getRefCoordFromGlobalCoord( closestPointTmp , coarseRefCoord);
        if ( closestPointDistSqr > 1.e-8 )
            std::cout << "Error: Squared distance to closest point: " << closestPointDistSqr << std::endl;
        // evaluate on coarse element
        Matrix32 DxCoarse;
        _xDFEvaluatorCoarse.evaluateGradient ( _confCoarse.getInitializer().getTriang(coarseElIdx), coarseRefCoord, DxCoarse );
        Tensor322Type DDxCoarse;
        _xDFEvaluatorCoarse.evaluateApproxHessianSym ( _confCoarse.getInitializer().getTriang(coarseElIdx), coarseRefCoord, DDxCoarse );
        PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> pointwiseEvaluator;
        TangentVecType normalCoarse;
        pointwiseEvaluator.evaluateNormal ( DxCoarse, normalCoarse  );
        Matrix22 gCoarse;
        pointwiseEvaluator.evaluateFirstFundamentalForm ( DxCoarse, gCoarse );
        Matrix22 gInvCoarse;
        gInvCoarse = gCoarse.inverse();
        Matrix22 hCoarse;
        pointwiseEvaluator.evaluateSecondFundamentalForm ( DDxCoarse, normalCoarse, hCoarse );
        const RealType GaussCurvCoarse = ( gInvCoarse * hCoarse ).determinant();
        // summarize
        return _xStorageFine.getArea(El.getGlobalElementIdx(),QuadPoint) * std::abs( GaussCurvFine - GaussCurvCoarse );
    }
};


// computes int \sqrt(g) |det(g_coarse^-1 h_coarse ) - det(g_fine^-1 h_fine )|
template<typename ConfiguratorType>
class GaussCurvatureL1DiffConfWithoutLeft :
public RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1DiffConfWithoutLeft<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::DomVecType DomVecType;
       typedef typename ConfiguratorType::TangentVecType TangentVecType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
       typedef typename ConfiguratorType::QuadRuleType QuadType;
 public:
      QuadType _quadRule;
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_confCoarse,  &_confFine;
       const DKTFEVectorFunctionEvaluator<ConfiguratorType> _xDFEvaluatorCoarse;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xStorageFine;

       vtkSmartPointer<vtkPolyData> _coarseMeshData;
       vtkSmartPointer<vtkCellLocator> _coarseMeshCellLocator;

  public:
    GaussCurvatureL1DiffConfWithoutLeft (
        const ConfiguratorType &confCoarse,
        const VectorType &FEDofsCoarse,
        const ConfiguratorType &confFine,
        const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xStorageFine ) :
      RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1DiffConfWithoutLeft<ConfiguratorType>> (confFine),
      _confCoarse( confCoarse ), _confFine ( confFine ),
      _xDFEvaluatorCoarse ( confCoarse,  FEDofsCoarse, 3 ),
      _xStorageFine( xStorageFine ),
      _coarseMeshData( vtkSmartPointer<vtkPolyData>::New() ),
      _coarseMeshCellLocator ( vtkSmartPointer<vtkCellLocator>::New() )
      {

        // save mesh data on fine scale to vtk
        const int numVerticesCoarse = _confCoarse.getInitializer().getNumVertices();
        const int numCellsCoarse = _confCoarse.getInitializer().getNumElements();
        auto points = vtkSmartPointer<vtkPoints>::New();
        for( int nodeIter=0; nodeIter<numVerticesCoarse; ++nodeIter ){
            auto point = _confCoarse.getInitializer().getVertex ( nodeIter );
            if( point.size() == 1 ) points->InsertNextPoint ( point[0], 0., 0. );
            if( point.size() == 2 ) points->InsertNextPoint ( point[0], point[1], 0. );
            if( point.size() == 3 ) points->InsertNextPoint ( point[0], point[1], point[2] );
        }
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        for( int elIdx=0; elIdx<numCellsCoarse; ++elIdx ){
            const int numNodesOfElement = _confCoarse.getInitializer().getNumNodesOfElement(elIdx);
            vtkIdType nodesOfElement [numNodesOfElement];
            for( int i=0; i<numNodesOfElement; ++i ) nodesOfElement[i] = _confCoarse.getInitializer().getElementNodeIdx( elIdx, i );
            cells->InsertNextCell ( numNodesOfElement, nodesOfElement );
        }
        _coarseMeshData->SetPoints ( points );
        _coarseMeshData->SetPolys ( cells );

        // Create the cell locator
        _coarseMeshCellLocator->SetDataSet(_coarseMeshData);
        _coarseMeshCellLocator->BuildLocator();
    }

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{

      bool isRight = true;
      for (int vertexIdx = 0; vertexIdx < 3; ++vertexIdx){
        // const Point3DType& vertexCoord = El.getVertex(vertexIdx);
        if (El.getNode(vertexIdx)[0] < 0.25){
          isRight = false;
        }
      }
      if (isRight){
        // evaluate on fine element
        auto gInvFine = _xStorageFine.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        auto hFine = _xStorageFine.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType GaussCurvFine = ( gInvFine * hFine ).determinant();
        // find coarseElIdx, coarseRefCoord
        Point3DType globalCoordsChartFineTmp;
        auto refCoordFine = _quadRule.getRefCoord ( QuadPoint );
        El.getGlobalCoord ( refCoordFine, globalCoordsChartFineTmp );
        double globalCoordsChartFine[3] = {globalCoordsChartFineTmp[0], globalCoordsChartFineTmp[1], globalCoordsChartFineTmp[2] };
        double closestPoint[3];
        double closestPointDistSqr;
        vtkIdType closestPointCellIdx;
        int closestPointSubIdx;
        _coarseMeshCellLocator->FindClosestPoint(globalCoordsChartFine, closestPoint, closestPointCellIdx, closestPointSubIdx, closestPointDistSqr);
        int coarseElIdx = closestPointCellIdx;
        Point3DType closestPointTmp;
        closestPointTmp(0) = closestPoint[0];
        closestPointTmp(1) = closestPoint[1];
        closestPointTmp(2) = closestPoint[2];
        DomVecType coarseRefCoord;
        _confCoarse.getInitializer().getTriang(coarseElIdx).getRefCoordFromGlobalCoord( closestPointTmp , coarseRefCoord);
        if ( closestPointDistSqr > 1.e-8 )
            std::cout << "Error: Squared distance to closest point: " << closestPointDistSqr << std::endl;
        // evaluate on coarse element
        Matrix32 DxCoarse;
        _xDFEvaluatorCoarse.evaluateGradient ( _confCoarse.getInitializer().getTriang(coarseElIdx), coarseRefCoord, DxCoarse );
        Tensor322Type DDxCoarse;
        _xDFEvaluatorCoarse.evaluateApproxHessianSym ( _confCoarse.getInitializer().getTriang(coarseElIdx), coarseRefCoord, DDxCoarse );
        PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> pointwiseEvaluator;
        TangentVecType normalCoarse;
        pointwiseEvaluator.evaluateNormal ( DxCoarse, normalCoarse  );
        Matrix22 gCoarse;
        pointwiseEvaluator.evaluateFirstFundamentalForm ( DxCoarse, gCoarse );
        Matrix22 gInvCoarse;
        gInvCoarse = gCoarse.inverse();
        Matrix22 hCoarse;
        pointwiseEvaluator.evaluateSecondFundamentalForm ( DDxCoarse, normalCoarse, hCoarse );
        const RealType GaussCurvCoarse = ( gInvCoarse * hCoarse ).determinant();
        // summarize
        return _xStorageFine.getArea(El.getGlobalElementIdx(),QuadPoint) * std::abs( GaussCurvFine - GaussCurvCoarse );
      }
      else {
        return 0.0;
      }
    }
};


// computes int \sqrt(g) |det(g_coarse^-1 h_coarse )|
template<typename ConfiguratorType>
class GaussCurvatureL1Conf :
public RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1Conf<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::DomVecType DomVecType;
       typedef typename ConfiguratorType::TangentVecType TangentVecType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
       typedef typename ConfiguratorType::QuadRuleType QuadType;
 public:
      QuadType _quadRule;
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_confCoarse,  &_confFine;
       const DKTFEVectorFunctionEvaluator<ConfiguratorType> _xDFEvaluatorCoarse;

       vtkSmartPointer<vtkPolyData> _coarseMeshData;
       vtkSmartPointer<vtkCellLocator> _coarseMeshCellLocator;

  public:
    GaussCurvatureL1Conf (
        const ConfiguratorType &confCoarse,
        const VectorType &FEDofsCoarse,
        const ConfiguratorType &confFine) :
      RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1Conf<ConfiguratorType>> (confFine),
      _confCoarse( confCoarse ), _confFine ( confFine ),
      _xDFEvaluatorCoarse ( confCoarse,  FEDofsCoarse, 3 ),
      _coarseMeshData( vtkSmartPointer<vtkPolyData>::New() ),
      _coarseMeshCellLocator ( vtkSmartPointer<vtkCellLocator>::New() )
      {

        // save mesh data on fine scale to vtk
        const int numVerticesCoarse = _confCoarse.getInitializer().getNumVertices();
        const int numCellsCoarse = _confCoarse.getInitializer().getNumElements();
        auto points = vtkSmartPointer<vtkPoints>::New();
        for( int nodeIter=0; nodeIter<numVerticesCoarse; ++nodeIter ){
            auto point = _confCoarse.getInitializer().getVertex ( nodeIter );
            if( point.size() == 1 ) points->InsertNextPoint ( point[0], 0., 0. );
            if( point.size() == 2 ) points->InsertNextPoint ( point[0], point[1], 0. );
            if( point.size() == 3 ) points->InsertNextPoint ( point[0], point[1], point[2] );
        }
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        for( int elIdx=0; elIdx<numCellsCoarse; ++elIdx ){
            const int numNodesOfElement = _confCoarse.getInitializer().getNumNodesOfElement(elIdx);
            vtkIdType nodesOfElement [numNodesOfElement];
            for( int i=0; i<numNodesOfElement; ++i ) nodesOfElement[i] = _confCoarse.getInitializer().getElementNodeIdx( elIdx, i );
            cells->InsertNextCell ( numNodesOfElement, nodesOfElement );
        }
        _coarseMeshData->SetPoints ( points );
        _coarseMeshData->SetPolys ( cells );

        // Create the cell locator
        _coarseMeshCellLocator->SetDataSet(_coarseMeshData);
        _coarseMeshCellLocator->BuildLocator();
    }

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{

        // find coarseElIdx, coarseRefCoord
        Point3DType globalCoordsChartFineTmp;
        auto refCoordFine = _quadRule.getRefCoord ( QuadPoint );
        El.getGlobalCoord ( refCoordFine, globalCoordsChartFineTmp );
        double globalCoordsChartFine[3] = {globalCoordsChartFineTmp[0], globalCoordsChartFineTmp[1], globalCoordsChartFineTmp[2] };
        double closestPoint[3];
        double closestPointDistSqr;
        vtkIdType closestPointCellIdx;
        int closestPointSubIdx;
        _coarseMeshCellLocator->FindClosestPoint(globalCoordsChartFine, closestPoint, closestPointCellIdx, closestPointSubIdx, closestPointDistSqr);
        int coarseElIdx = closestPointCellIdx;
        Point3DType closestPointTmp;
        closestPointTmp(0) = closestPoint[0];
        closestPointTmp(1) = closestPoint[1];
        closestPointTmp(2) = closestPoint[2];
        DomVecType coarseRefCoord;
        _confCoarse.getInitializer().getTriang(coarseElIdx).getRefCoordFromGlobalCoord( closestPointTmp , coarseRefCoord);
        if ( closestPointDistSqr > 1.e-8 )
            std::cout << "Error: Squared distance to closest point: " << closestPointDistSqr << std::endl;
        // evaluate on coarse element
        Matrix32 DxCoarse;
        _xDFEvaluatorCoarse.evaluateGradient ( _confCoarse.getInitializer().getTriang(coarseElIdx), coarseRefCoord, DxCoarse );
        Tensor322Type DDxCoarse;
        _xDFEvaluatorCoarse.evaluateApproxHessianSym ( _confCoarse.getInitializer().getTriang(coarseElIdx), coarseRefCoord, DDxCoarse );
        PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> pointwiseEvaluator;
        TangentVecType normalCoarse;
        pointwiseEvaluator.evaluateNormal ( DxCoarse, normalCoarse  );
        Matrix22 gCoarse;
        pointwiseEvaluator.evaluateFirstFundamentalForm ( DxCoarse, gCoarse );
        Matrix22 gInvCoarse;
        gInvCoarse = gCoarse.inverse();
        Matrix22 hCoarse;
        pointwiseEvaluator.evaluateSecondFundamentalForm ( DDxCoarse, normalCoarse, hCoarse );
        const RealType GaussCurvCoarse = ( gInvCoarse * hCoarse ).determinant();
        const RealType AreaCoarse = sqrt(gCoarse.determinant());
        // summarize
        return AreaCoarse * std::abs( GaussCurvCoarse );
    }
};




// computes int \sqrt(g) |det(g_coarse^-1 h_coarse )| only on [0.25,1.0]x[0.0,1.0]
template<typename ConfiguratorType>
class GaussCurvatureL1ConfWithoutLeft :
public RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1ConfWithoutLeft<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::DomVecType DomVecType;
       typedef typename ConfiguratorType::TangentVecType TangentVecType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
       typedef typename ConfiguratorType::QuadRuleType QuadType;
 public:
      QuadType _quadRule;
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_confCoarse,  &_confFine;
       const DKTFEVectorFunctionEvaluator<ConfiguratorType> _xDFEvaluatorCoarse;

       vtkSmartPointer<vtkPolyData> _coarseMeshData;
       vtkSmartPointer<vtkCellLocator> _coarseMeshCellLocator;

  public:
    GaussCurvatureL1ConfWithoutLeft (
        const ConfiguratorType &confCoarse,
        const VectorType &FEDofsCoarse,
        const ConfiguratorType &confFine) :
      RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1ConfWithoutLeft<ConfiguratorType>> (confFine),
      _confCoarse( confCoarse ), _confFine ( confFine ),
      _xDFEvaluatorCoarse ( confCoarse,  FEDofsCoarse, 3 ),
      _coarseMeshData( vtkSmartPointer<vtkPolyData>::New() ),
      _coarseMeshCellLocator ( vtkSmartPointer<vtkCellLocator>::New() )
      {

        // save mesh data on fine scale to vtk
        const int numVerticesCoarse = _confCoarse.getInitializer().getNumVertices();
        const int numCellsCoarse = _confCoarse.getInitializer().getNumElements();
        auto points = vtkSmartPointer<vtkPoints>::New();
        for( int nodeIter=0; nodeIter<numVerticesCoarse; ++nodeIter ){
            auto point = _confCoarse.getInitializer().getVertex ( nodeIter );
            if( point.size() == 1 ) points->InsertNextPoint ( point[0], 0., 0. );
            if( point.size() == 2 ) points->InsertNextPoint ( point[0], point[1], 0. );
            if( point.size() == 3 ) points->InsertNextPoint ( point[0], point[1], point[2] );
        }
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        for( int elIdx=0; elIdx<numCellsCoarse; ++elIdx ){
            const int numNodesOfElement = _confCoarse.getInitializer().getNumNodesOfElement(elIdx);
            vtkIdType nodesOfElement [numNodesOfElement];
            for( int i=0; i<numNodesOfElement; ++i ) nodesOfElement[i] = _confCoarse.getInitializer().getElementNodeIdx( elIdx, i );
            cells->InsertNextCell ( numNodesOfElement, nodesOfElement );
        }
        _coarseMeshData->SetPoints ( points );
        _coarseMeshData->SetPolys ( cells );

        // Create the cell locator
        _coarseMeshCellLocator->SetDataSet(_coarseMeshData);
        _coarseMeshCellLocator->BuildLocator();
    }

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{

      bool isRight = true;
      for (int vertexIdx = 0; vertexIdx < 3; ++vertexIdx){
        // const Point3DType& vertexCoord = El.getVertex(vertexIdx);
        if (El.getNode(vertexIdx)[0] < 0.25){
          isRight = false;
        }
      }
      if (isRight){
        // find coarseElIdx, coarseRefCoord
        Point3DType globalCoordsChartFineTmp;
        auto refCoordFine = _quadRule.getRefCoord ( QuadPoint );
        El.getGlobalCoord ( refCoordFine, globalCoordsChartFineTmp );
        double globalCoordsChartFine[3] = {globalCoordsChartFineTmp[0], globalCoordsChartFineTmp[1], globalCoordsChartFineTmp[2] };
        double closestPoint[3];
        double closestPointDistSqr;
        vtkIdType closestPointCellIdx;
        int closestPointSubIdx;
        _coarseMeshCellLocator->FindClosestPoint(globalCoordsChartFine, closestPoint, closestPointCellIdx, closestPointSubIdx, closestPointDistSqr);
        int coarseElIdx = closestPointCellIdx;
        Point3DType closestPointTmp;
        closestPointTmp(0) = closestPoint[0];
        closestPointTmp(1) = closestPoint[1];
        closestPointTmp(2) = closestPoint[2];
        DomVecType coarseRefCoord;
        _confCoarse.getInitializer().getTriang(coarseElIdx).getRefCoordFromGlobalCoord( closestPointTmp , coarseRefCoord);
        if ( closestPointDistSqr > 1.e-8 )
            std::cout << "Error: Squared distance to closest point: " << closestPointDistSqr << std::endl;
        // evaluate on coarse element
        Matrix32 DxCoarse;
        _xDFEvaluatorCoarse.evaluateGradient ( _confCoarse.getInitializer().getTriang(coarseElIdx), coarseRefCoord, DxCoarse );
        Tensor322Type DDxCoarse;
        _xDFEvaluatorCoarse.evaluateApproxHessianSym ( _confCoarse.getInitializer().getTriang(coarseElIdx), coarseRefCoord, DDxCoarse );
        PointWiseVectorFunctionEvaluatorShellFE<ConfiguratorType> pointwiseEvaluator;
        TangentVecType normalCoarse;
        pointwiseEvaluator.evaluateNormal ( DxCoarse, normalCoarse  );
        Matrix22 gCoarse;
        pointwiseEvaluator.evaluateFirstFundamentalForm ( DxCoarse, gCoarse );
        Matrix22 gInvCoarse;
        gInvCoarse = gCoarse.inverse();
        Matrix22 hCoarse;
        pointwiseEvaluator.evaluateSecondFundamentalForm ( DDxCoarse, normalCoarse, hCoarse );
        const RealType GaussCurvCoarse = ( gInvCoarse * hCoarse ).determinant();
        const RealType AreaCoarse = sqrt(gCoarse.determinant());
        // summarize
        return AreaCoarse * std::abs( GaussCurvCoarse );
      }
      else {
        return 0.0;
      }
    }
};




// computes int \sqrt(g) det(g^-1 h )
template<typename ConfiguratorType>
class GaussCurvatureInt :
public RefTriangleIntegrator<ConfiguratorType, GaussCurvatureInt<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xStorage;
  public:
    GaussCurvatureInt ( const ConfiguratorType &conf,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xStorage ) :
      RefTriangleIntegrator<ConfiguratorType, GaussCurvatureInt<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xStorage ( xStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gInv = _xStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& h = _xStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType GaussCurv = ( gInv * h ).determinant();
        return _xStorage.getArea(El.getGlobalElementIdx(),QuadPoint) *  GaussCurv ;
    }
};


// computes int \sqrt(g_A) |det(g_B^-1 h_B) - det(g_A^-1 h_A )|
template<typename ConfiguratorType>
class GaussCurvatureL1Diff :
public RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1Diff<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xBStorage;
  public:
    GaussCurvatureL1Diff ( const ConfiguratorType &conf,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xBStorage ) :
      RefTriangleIntegrator<ConfiguratorType, GaussCurvatureL1Diff<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xAStorage ( xAStorage ),
      _xBStorage ( xBStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gAInv = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& hA = _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType GaussCurvA = ( gAInv * hA ).determinant();
        const Matrix22& gBInv = _xBStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& hB = _xBStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType GaussCurvB = ( gBInv * hB ).determinant();
        return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * std::abs( GaussCurvB - GaussCurvA );
        // return _xBStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * std::abs( GaussCurvB );
    }
};


// computes int \sqrt(g) |tr(g^-1 h )|
template<typename ConfiguratorType>
class MeanCurvatureL1 :
public RefTriangleIntegrator<ConfiguratorType, MeanCurvatureL1<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xStorage;
  public:
    MeanCurvatureL1 ( const ConfiguratorType &conf,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xStorage ) :
      RefTriangleIntegrator<ConfiguratorType, MeanCurvatureL1<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xStorage ( xStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gInv = _xStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& h = _xStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType MeanCurv = ( gInv * h ).trace();
        return _xStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * std::abs( MeanCurv );
    }
};


// computes int \sqrt(g_A) |g_A^-1 (h_B -h_A )|^2
template<typename ConfiguratorType>
class RelativeShapeOperatorL2 :
public RefTriangleIntegrator<ConfiguratorType, RelativeShapeOperatorL2<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage, &_xBStorage;
  public:
    RelativeShapeOperatorL2 ( const ConfiguratorType &conf,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xBStorage ) :
      RefTriangleIntegrator<ConfiguratorType, RelativeShapeOperatorL2<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xAStorage ( xAStorage ),
      _xBStorage ( xBStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gAInv = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& hA = _xAStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& hB = _xBStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType aux = ( gAInv * (hB - hA) ).squaredNorm();
        return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * aux;
    }
};


// computes int \sqrt(g_A) |h_A |^2
template<typename ConfiguratorType>
class SecondFFL2 :
public RefTriangleIntegrator<ConfiguratorType, SecondFFL2<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xStorage;
  public:
    SecondFFL2 ( const ConfiguratorType &conf,
                           const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xStorage) :
      RefTriangleIntegrator<ConfiguratorType, SecondFFL2<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xStorage ( xStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gInv = _xStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Matrix22& h = _xStorage.getSecondFF( El.getGlobalElementIdx(), QuadPoint );
        const RealType aux =  h.squaredNorm();
        return _xStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * aux;
    }
};



// computes int \sqrt(g_A) |D^2 (x_B - x_A)|^2
template<typename ConfiguratorType>
class SecondDerivativeEnergy :
public RefTriangleIntegrator<ConfiguratorType, SecondDerivativeEnergy<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xAStorage, &_xBStorage;
  public:
    SecondDerivativeEnergy ( const ConfiguratorType &conf,
                             const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xAStorage,
                             const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xBStorage ) :
      RefTriangleIntegrator<ConfiguratorType, SecondDerivativeEnergy<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xAStorage ( xAStorage ),
      _xBStorage ( xBStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gAInv = _xAStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Tensor322Type& D2xA = _xAStorage.getHessian( El.getGlobalElementIdx(), QuadPoint );
        const Tensor322Type& D2xB = _xBStorage.getHessian( El.getGlobalElementIdx(), QuadPoint );
        Tensor322Type D2u;
        for( int i=0; i<3; ++i )
            for( int j=0; j<3; ++j )
                for( int k = 0; k <2; ++k )
                    D2u.set(i,j,k, D2xB.get(i,j,k) - D2xA.get(i,j,k) );

        const RealType aux =  D2u.normSqr();
        return _xAStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * aux;
    }
};


// computes int \sqrt(g_A) |D^2 (x_B - x_A)|^2
template<typename ConfiguratorType>
class SecondDerivativeEnergyConf :
public RefTriangleIntegrator<ConfiguratorType, SecondDerivativeEnergyConf<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::DomVecType DomVecType;
       typedef typename ConfiguratorType::TangentVecType TangentVecType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
       typedef typename ConfiguratorType::QuadRuleType QuadType;
 public:
      QuadType _quadRule;
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_confCoarse,  &_confFine;
       const DKTFEVectorFunctionEvaluator<ConfiguratorType> _xDFEvaluatorCoarse;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xStorageFine;

       vtkSmartPointer<vtkPolyData> _coarseMeshData;
       vtkSmartPointer<vtkCellLocator> _coarseMeshCellLocator;

  public:
    SecondDerivativeEnergyConf (
        const ConfiguratorType &confCoarse,
        const VectorType &FEDofsCoarse,
        const ConfiguratorType &confFine,
        const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xStorageFine ) :
      RefTriangleIntegrator<ConfiguratorType, SecondDerivativeEnergyConf<ConfiguratorType>> (confFine),
      _confCoarse( confCoarse ), _confFine ( confFine ),
      _xDFEvaluatorCoarse ( confCoarse,  FEDofsCoarse, 3 ),
      _xStorageFine( xStorageFine ),
      _coarseMeshData( vtkSmartPointer<vtkPolyData>::New() ),
      _coarseMeshCellLocator ( vtkSmartPointer<vtkCellLocator>::New() )
      {

        // save mesh data on fine scale to vtk
        const int numVerticesCoarse = _confCoarse.getInitializer().getNumVertices();
        const int numCellsCoarse = _confCoarse.getInitializer().getNumElements();
        auto points = vtkSmartPointer<vtkPoints>::New();
        for( int nodeIter=0; nodeIter<numVerticesCoarse; ++nodeIter ){
            auto point = _confCoarse.getInitializer().getVertex ( nodeIter );
            if( point.size() == 1 ) points->InsertNextPoint ( point[0], 0., 0. );
            if( point.size() == 2 ) points->InsertNextPoint ( point[0], point[1], 0. );
            if( point.size() == 3 ) points->InsertNextPoint ( point[0], point[1], point[2] );
        }
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        for( int elIdx=0; elIdx<numCellsCoarse; ++elIdx ){
            const int numNodesOfElement = _confCoarse.getInitializer().getNumNodesOfElement(elIdx);
            vtkIdType nodesOfElement [numNodesOfElement];
            for( int i=0; i<numNodesOfElement; ++i ) nodesOfElement[i] = _confCoarse.getInitializer().getElementNodeIdx( elIdx, i );
            cells->InsertNextCell ( numNodesOfElement, nodesOfElement );
        }
        _coarseMeshData->SetPoints ( points );
        _coarseMeshData->SetPolys ( cells );

        // Create the cell locator
        _coarseMeshCellLocator->SetDataSet(_coarseMeshData);
        _coarseMeshCellLocator->BuildLocator();
    }

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        // evaluate on fine element
        auto DDxFine = _xStorageFine.getHessian(El.getGlobalElementIdx(), QuadPoint);
        // find coarseElIdx, coarseRefCoord
        Point3DType globalCoordsChartFineTmp;
        auto refCoordFine = _quadRule.getRefCoord ( QuadPoint );
        El.getGlobalCoord ( refCoordFine, globalCoordsChartFineTmp );
        double globalCoordsChartFine[3] = {globalCoordsChartFineTmp[0], globalCoordsChartFineTmp[1], globalCoordsChartFineTmp[2] };
        double closestPoint[3];
        double closestPointDistSqr;
        vtkIdType closestPointCellIdx;
        int closestPointSubIdx;
        _coarseMeshCellLocator->FindClosestPoint(globalCoordsChartFine, closestPoint, closestPointCellIdx, closestPointSubIdx, closestPointDistSqr);
        int coarseElIdx = closestPointCellIdx;
        Point3DType closestPointTmp;
        closestPointTmp(0) = closestPoint[0];
        closestPointTmp(1) = closestPoint[1];
        closestPointTmp(2) = closestPoint[2];
        DomVecType coarseRefCoord;
        _confCoarse.getInitializer().getTriang(coarseElIdx).getRefCoordFromGlobalCoord( closestPointTmp , coarseRefCoord);
        if ( closestPointDistSqr > 1.e-8 )
            std::cout << "Error: Squared distance to closest point: " << closestPointDistSqr << std::endl;
        // evaluate on coarse element
        Tensor322Type DDxCoarse;
        _xDFEvaluatorCoarse.evaluateApproxHessianSym ( _confCoarse.getInitializer().getTriang(coarseElIdx), coarseRefCoord, DDxCoarse );
        // summarize
        Tensor322Type D2u;
        for( int i=0; i<3; ++i )
            for( int j=0; j<3; ++j )
                for( int k = 0; k <2; ++k )
                    D2u.set(i,j,k, DDxFine.get(i,j,k) - DDxCoarse.get(i,j,k) );

        const RealType aux =  D2u.normSqr();
        return _xStorageFine.getArea(El.getGlobalElementIdx(),QuadPoint) * aux;
    }
};



// computes int \sqrt(g_A) |D^2 (x_B - x_A)|^2 only on [0.25,1.0]x[0.0,1.0]
template<typename ConfiguratorType>
class SecondDerivativeEnergyConfWithoutLeft :
public RefTriangleIntegrator<ConfiguratorType, SecondDerivativeEnergyConfWithoutLeft<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::DomVecType DomVecType;
       typedef typename ConfiguratorType::TangentVecType TangentVecType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
       typedef typename ConfiguratorType::QuadRuleType QuadType;
 public:
      QuadType _quadRule;
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_confCoarse,  &_confFine;
       const DKTFEVectorFunctionEvaluator<ConfiguratorType> _xDFEvaluatorCoarse;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xStorageFine;

       vtkSmartPointer<vtkPolyData> _coarseMeshData;
       vtkSmartPointer<vtkCellLocator> _coarseMeshCellLocator;

  public:
    SecondDerivativeEnergyConfWithoutLeft (
        const ConfiguratorType &confCoarse,
        const VectorType &FEDofsCoarse,
        const ConfiguratorType &confFine,
        const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xStorageFine ) :
      RefTriangleIntegrator<ConfiguratorType, SecondDerivativeEnergyConfWithoutLeft<ConfiguratorType>> (confFine),
      _confCoarse( confCoarse ), _confFine ( confFine ),
      _xDFEvaluatorCoarse ( confCoarse,  FEDofsCoarse, 3 ),
      _xStorageFine( xStorageFine ),
      _coarseMeshData( vtkSmartPointer<vtkPolyData>::New() ),
      _coarseMeshCellLocator ( vtkSmartPointer<vtkCellLocator>::New() )
      {

        // save mesh data on fine scale to vtk
        const int numVerticesCoarse = _confCoarse.getInitializer().getNumVertices();
        const int numCellsCoarse = _confCoarse.getInitializer().getNumElements();
        auto points = vtkSmartPointer<vtkPoints>::New();
        for( int nodeIter=0; nodeIter<numVerticesCoarse; ++nodeIter ){
            auto point = _confCoarse.getInitializer().getVertex ( nodeIter );
            if( point.size() == 1 ) points->InsertNextPoint ( point[0], 0., 0. );
            if( point.size() == 2 ) points->InsertNextPoint ( point[0], point[1], 0. );
            if( point.size() == 3 ) points->InsertNextPoint ( point[0], point[1], point[2] );
        }
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        for( int elIdx=0; elIdx<numCellsCoarse; ++elIdx ){
            const int numNodesOfElement = _confCoarse.getInitializer().getNumNodesOfElement(elIdx);
            vtkIdType nodesOfElement [numNodesOfElement];
            for( int i=0; i<numNodesOfElement; ++i ) nodesOfElement[i] = _confCoarse.getInitializer().getElementNodeIdx( elIdx, i );
            cells->InsertNextCell ( numNodesOfElement, nodesOfElement );
        }
        _coarseMeshData->SetPoints ( points );
        _coarseMeshData->SetPolys ( cells );

        // Create the cell locator
        _coarseMeshCellLocator->SetDataSet(_coarseMeshData);
        _coarseMeshCellLocator->BuildLocator();
    }

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        // evaluate on fine element
        bool isRight = true;
        for (int vertexIdx = 0; vertexIdx < 3; ++vertexIdx){
          // const Point3DType& vertexCoord = El.getVertex(vertexIdx);
          if (El.getNode(vertexIdx)[0] < 0.25){
            isRight = false;
          }
        }
        if (isRight){
          auto DDxFine = _xStorageFine.getHessian(El.getGlobalElementIdx(), QuadPoint);
          // find coarseElIdx, coarseRefCoord
          Point3DType globalCoordsChartFineTmp;
          auto refCoordFine = _quadRule.getRefCoord ( QuadPoint );
          El.getGlobalCoord ( refCoordFine, globalCoordsChartFineTmp );
          double globalCoordsChartFine[3] = {globalCoordsChartFineTmp[0], globalCoordsChartFineTmp[1], globalCoordsChartFineTmp[2] };
          double closestPoint[3];
          double closestPointDistSqr;
          vtkIdType closestPointCellIdx;
          int closestPointSubIdx;
          _coarseMeshCellLocator->FindClosestPoint(globalCoordsChartFine, closestPoint, closestPointCellIdx, closestPointSubIdx, closestPointDistSqr);
          int coarseElIdx = closestPointCellIdx;
          Point3DType closestPointTmp;
          closestPointTmp(0) = closestPoint[0];
          closestPointTmp(1) = closestPoint[1];
          closestPointTmp(2) = closestPoint[2];
          DomVecType coarseRefCoord;
          _confCoarse.getInitializer().getTriang(coarseElIdx).getRefCoordFromGlobalCoord( closestPointTmp , coarseRefCoord);
          if ( closestPointDistSqr > 1.e-8 )
              std::cout << "Error: Squared distance to closest point: " << closestPointDistSqr << std::endl;
          // evaluate on coarse element
          Tensor322Type DDxCoarse;
          _xDFEvaluatorCoarse.evaluateApproxHessianSym ( _confCoarse.getInitializer().getTriang(coarseElIdx), coarseRefCoord, DDxCoarse );
          // summarize
          Tensor322Type D2u;
          for( int i=0; i<3; ++i )
              for( int j=0; j<3; ++j )
                  for( int k = 0; k <2; ++k )
                      D2u.set(i,j,k, DDxFine.get(i,j,k) - DDxCoarse.get(i,j,k) );

          const RealType aux =  D2u.normSqr();
          return _xStorageFine.getArea(El.getGlobalElementIdx(),QuadPoint) * aux;
        }
        else {
          return 0.0;
        }
    }
};



// computes int \sqrt(g_A) |D^2  x_A|^2
template<typename ConfiguratorType>
class SecondDerivativeL2:
public RefTriangleIntegrator<ConfiguratorType, SecondDerivativeL2<ConfiguratorType> >
{
  protected:
       typedef typename ConfiguratorType::RealType RealType;
       typedef typename ConfiguratorType::Point3DType Point3DType;
       typedef typename ConfiguratorType::Matrix22  Matrix22;
       typedef typename ConfiguratorType::Matrix32  Matrix32;
       typedef typename ConfiguratorType::Matrix33  Matrix33;
       typedef typename ConfiguratorType::Tensor322Type Tensor322Type;
       typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
       typedef typename ConfiguratorType::ElementType ElementType;
       typedef typename ConfiguratorType::VectorType VectorType;
 public:
      static const DiscreteFunctionCacheType _DiscreteFunctionCacheType = FirstAndSecondOrder;
 protected:
       const ConfiguratorType &_conf;
       const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &_xStorage;
  public:
    SecondDerivativeL2 ( const ConfiguratorType &conf,
                             const DiscreteVectorFunctionStorage<ConfiguratorType,_DiscreteFunctionCacheType> &xStorage ) :
      RefTriangleIntegrator<ConfiguratorType, SecondDerivativeL2<ConfiguratorType>> (conf),
      _conf ( conf ),
      _xStorage ( xStorage ) {}

    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const{
        const Matrix22& gInv = _xStorage.getFirstFFInv( El.getGlobalElementIdx(), QuadPoint );
        const Tensor322Type& D2x = _xStorage.getHessian( El.getGlobalElementIdx(), QuadPoint );
        const RealType aux =  D2x.normSqr();
        return _xStorage.getArea(El.getGlobalElementIdx(),QuadPoint) * aux;
    }
};


#endif
