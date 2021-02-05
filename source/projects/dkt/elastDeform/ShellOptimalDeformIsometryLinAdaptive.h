#ifndef __SHELLBUCKLINGOPTIMALDEFORMISOMETRYADADAPTIVE_H
#define __SHELLBUCKLINGOPTIMALDEFORMISOMETRYADADAPTIVE_H

#include <pesopt_VTK.h>
#include <energyDefines.h>

#include "ShellOptimalDeformSolverInterfaceAdaptive.h"
#include "ShellOptimalDeformSolverIsometryLin.h"
#include "ShellCurvature.h"


template<typename MatOptConfigurator, typename LinElastEnergyType, bool IsometryConstraint>
class optimalDeformSolverIsometryAdaptive 
: public optimalDeformSolverInterfaceAdaptive<MatOptConfigurator>{
  
  typedef typename MatOptConfigurator::ConfiguratorType         ConfiguratorType;
  typedef typename MatOptConfigurator::ConfiguratorTypePf       ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType                   RealType;
  typedef typename ConfiguratorType::MaskType                   MaskType;
  typedef typename ConfiguratorType::InitType                   MeshType;
  typedef typename ConfiguratorTypePf::DTContainer              DataTypeContainer;
  typedef typename ConfiguratorType::VectorType                 VectorType;
  typedef pesopt::BoostParser   ParameterParserType;
  typedef typename DataTypeContainer::Matrix22                  Matrix22;
  typedef typename DataTypeContainer::Matrix32                  Matrix32;
  typedef typename DataTypeContainer::DomVecType                DomVecType;
  typedef typename DataTypeContainer::TangentVecType            TangentVecType;
  typedef typename LinElastEnergyType::EvaluationHelper         EvaluationHelper;
  
public:
  optimalDeformSolverIsometryAdaptive ( const ParameterParserType &Parser ) :
  optimalDeformSolverInterfaceAdaptive<MatOptConfigurator> ( Parser ) { } 

  void computeOptimalDeformation_OnCurrentMesh( 
           const MeshType & mesh,                    
           const VectorType& material,  
           const VectorType &initDisp, VectorType &solDisp,
           MaskType &DirichletMask,
           VectorType &membraneStressVec, 
           VectorType &bendingStressVec, 
           VectorType &totalStressVec,
           VectorType &GaussCurvVec,
           DeformationOptimizationShellEnergyInfo<RealType> &energyInfo,
           IsometryInfo<RealType> &isometryInfo,
           const string saveDirectory,
           const int refinementStep = 0,
           const RealType factorForce = 1.
    ) const {
       
    ParameterParserType parser; parser = this->_parser;
    parser.set ( "saving.saveDirectory", saveDirectory );
    parser.saveToFile( pesopt::strprintf ( "ParameterParser.ini" ) );
    mesh.printBasicInfo();
    auto startTime = std::chrono::high_resolution_clock::now(); 
    
    //initialize configurators on current level
    ConfiguratorType conf ( mesh ); ConfiguratorTypePf confpf ( mesh );
    MatOptConfigurator matOptConf ( parser, conf, confpf );
      
    ShellHandler<ConfiguratorType,LinElastEnergyType::_DiscreteFunctionCacheType> shellHandler ( parser, conf );
    if( parser.template get<bool> ( "ConstraintProblemAdaptiveProlongationTypes.constructProlongatedBoundaryByShellHandler" ) )  DirichletMask = shellHandler.getDirichletMask();
    else shellHandler.setDirichletMask( DirichletMask ); 
    
    //! construct rhs 
    VectorType rhs_force ( 3 * conf.getNumGlobalDofs() );
    this->constructForceRHS( parser, matOptConf, shellHandler.getChartToUndeformedShell_Cache(), shellHandler.getDirichletMask(), rhs_force, factorForce );
    
    //! Compute Optimal Deformation
    optimalDeformSolver<MatOptConfigurator, true, LinElastEnergyType, IsometryConstraint> OptDeformFinder ( parser, matOptConf, material, shellHandler, rhs_force, initDisp );
    
    solDisp.resize( OptDeformFinder.getSolutionDisplacement().size() ); solDisp = OptDeformFinder.getSolutionDisplacement();  
    OptDeformFinder.getEnergyInfo( energyInfo );
    OptDeformFinder.getLastStressOnElements( membraneStressVec, bendingStressVec, totalStressVec );
    
    //! print elapsed time
    std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << endl << "duration of refinementStep = " << diff.count() << " sec" << endl;
    pesopt::consoleOutput( pesopt::strprintf( "finished computation optimal design for refinementStep %d", refinementStep ).c_str() );
    
    //! save results
    VectorType areaOnElementsVec ( mesh.getNumElements() );
    evaluateAreaOfElements<ConfiguratorType,LinElastEnergyType::_DiscreteFunctionCacheType>( conf, shellHandler, areaOnElementsVec ); 
    VectorType averagedMembraneStressVec ( membraneStressVec ), averagedBendingStressVec( bendingStressVec ), averagedTotalStressVec( totalStressVec );
    for( int elIdx=0; elIdx < mesh.getNumElements(); ++elIdx ){
        averagedMembraneStressVec[elIdx] = membraneStressVec[elIdx] / areaOnElementsVec[elIdx];
        averagedBendingStressVec[elIdx] = bendingStressVec[elIdx] / areaOnElementsVec[elIdx];
        averagedTotalStressVec[elIdx] = totalStressVec[elIdx] / areaOnElementsVec[elIdx];  
    }
    const std::vector<VectorType> StressOnElementVec = {membraneStressVec, bendingStressVec, totalStressVec,
                                                        averagedMembraneStressVec, averagedBendingStressVec, averagedTotalStressVec,
                                                        areaOnElementsVec };
    const std::vector<string> nameStressOnElementVec = {"membraneStressOnElements", "bendingStressOnElements", "totalStressOnElements",
                                                        "averagedMembraneStressOnElements", "averagedBendingStressOnElements", "averagedTotalStressOnElements", 
                                                        "areaOnElements"};
    this->template saveAllResults<LinElastEnergyType::_DiscreteFunctionCacheType,IsometryConstraint> ( saveDirectory, matOptConf, shellHandler, material, solDisp, StressOnElementVec, nameStressOnElementVec, totalStressVec, areaOnElementsVec, energyInfo );
    
    
    //BEGIN TEST CURVATURE
    
    const int numVertices = mesh.getNumVertices();
    const int numGlobalDofsDeform = conf.getNumGlobalDofs();
    VectorType xB ( solDisp.size() ); xB = solDisp + shellHandler.getChartToUndeformedShell();
    DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xBStorage ( matOptConf._conf, xB, 3 );
    
    //Isometry Constraint 
    IsometryConstraintIntegrated<MatOptConfigurator> isoOpIntegrated ( matOptConf, shellHandler.getChartToUndeformedShell_Cache() );
    RealType isometryErrorL2Sqr = 0.; isoOpIntegrated.apply( solDisp, isometryErrorL2Sqr );
    VectorType isoErrorL2SqrOnElements ( totalStressVec.size() ); isoOpIntegrated.applyOnElements( solDisp, isoErrorL2SqrOnElements );
    VectorType avIsoErrorL2SqrOnElements ( totalStressVec.size() );
    for( int elIdx=0; elIdx < mesh.getNumElements(); ++elIdx ) avIsoErrorL2SqrOnElements[elIdx] = isoErrorL2SqrOnElements[elIdx] / areaOnElementsVec[elIdx];
    isometryInfo.setIsometryErrorL2( sqrt(isometryErrorL2Sqr) );
    
    IsometryConstraintL1<MatOptConfigurator> isoOpL1 ( matOptConf, shellHandler.getChartToUndeformedShell_Cache() );
    RealType isometryErrorL1 = 0.; isoOpL1.apply( solDisp, isometryErrorL1 );
    VectorType isoErrorL1OnElements ( totalStressVec.size() ); isoOpL1.applyOnElements( solDisp, isoErrorL1OnElements );
    VectorType avIsoErrorL1OnElements ( totalStressVec.size() );
    for( int elIdx=0; elIdx < mesh.getNumElements(); ++elIdx ) avIsoErrorL1OnElements[elIdx] = isoErrorL1OnElements[elIdx] / areaOnElementsVec[elIdx];
    isometryInfo.setIsometryErrorL1( isometryErrorL1 );
    
    //Gauss and Mean Curvature, Relative Shape Operator
    VectorType GaussCurvVector ( totalStressVec.size() ), MeanCurvVector ( totalStressVec.size() ), RelShapeOpVector ( totalStressVec.size() );
    VectorType avGaussCurvVector ( totalStressVec.size() ), avMeanCurvVector ( totalStressVec.size() ), avRelShapeOpVector ( totalStressVec.size() );
    GaussCurvatureL1<ConfiguratorType> ( matOptConf._conf, xBStorage ).assembleOnElements( GaussCurvVector ); 
    MeanCurvatureL1<ConfiguratorType> ( matOptConf._conf, xBStorage ).assembleOnElements( MeanCurvVector ); 
    RelativeShapeOperatorL2<ConfiguratorType> ( matOptConf._conf, shellHandler.getChartToUndeformedShell_Cache(), xBStorage ).assembleOnElements( RelShapeOpVector ); 
    for( int elIdx=0; elIdx < mesh.getNumElements(); ++elIdx ){
        avGaussCurvVector[elIdx] = GaussCurvVector[elIdx] / areaOnElementsVec[elIdx];
        avMeanCurvVector[elIdx] = MeanCurvVector[elIdx] / areaOnElementsVec[elIdx];
        avRelShapeOpVector[elIdx] = RelShapeOpVector[elIdx] / areaOnElementsVec[elIdx];
    }
    RealType integralGauss = 0.;
    GaussCurvatureL1<ConfiguratorType> ( matOptConf._conf, xBStorage ).assembleAdd( integralGauss );
    isometryInfo.setGaussCurvatureL1( integralGauss );
    GaussCurvVec = GaussCurvVector;
 
    RealType integralRelShapeOp = 0.;
    RelativeShapeOperatorL2<ConfiguratorType> ( matOptConf._conf, shellHandler.getChartToUndeformedShell_Cache(), xBStorage ).assembleAdd( integralRelShapeOp );
    RealType maxRelShapeOp = avRelShapeOpVector.maxCoeff();
    isometryInfo.setRelativeShapeOpL2( sqrt(integralRelShapeOp) );
    
    //compute average value of |D^2 u|^2 on each element and plot affine regions
    EvaluationHelper LinElastHelper ( matOptConf, shellHandler.getChartToUndeformedShell_Cache() );
    VectorType constZeroMaterial ( material.size() ); for( int i=0; i<constZeroMaterial.size(); ++i ) constZeroMaterial[i] = -1.;
    VectorType averagedMembraneStressVecWithoutMaterialFac ( totalStressVec.size() ), averagedBendingStressVecWithoutMaterialFac ( totalStressVec.size() ), averagedTotalStressVecWithoutMaterialFac ( totalStressVec.size() );
    LinElastHelper.evaluateAveragedStressOnElements( constZeroMaterial, solDisp, shellHandler.getChartToUndeformedShell(),
                                                     averagedMembraneStressVecWithoutMaterialFac, averagedBendingStressVecWithoutMaterialFac, averagedTotalStressVecWithoutMaterialFac );
    
    //Alternative Evaluation of Gauss Curvature with Lagrange Interpolation Operators
    typedef RefTriangMeshConfiguratorP1<DataTypeContainer, MeshType, typename ConfiguratorType::QuadRuleType> ConfiguratorTypeP1;
    ConfiguratorTypeP1 confP1 ( mesh );
    VectorType GaussCurvVector_Alternative ( totalStressVec.size() ); // MeanCurvVector ( totalStressVec.size() );
    VectorType avGaussCurvVector_Alternative ( totalStressVec.size() ); // avMeanCurvVector ( totalStressVec.size() );
    std::vector<Matrix32> DxBAtNodes ( numVertices );
    VectorType normalDeformed_LagrangeInterpolationP1 ( 3 * numVertices );
    VectorType D1xB_LagrangeInterpolationP1 ( 3 * numVertices ), D2xB_LagrangeInterpolationP1 ( 3 * numVertices );
    for( int nodeIdx = 0; nodeIdx < numVertices; ++nodeIdx ){
        TangentVecType tangentVec1, tangentVec2;
        for( int comp = 0; comp < 3; ++comp ) {
            tangentVec1[comp] = xB[nodeIdx +     numVertices + comp * numGlobalDofsDeform];
            tangentVec2[comp] = xB[nodeIdx + 2 * numVertices + comp * numGlobalDofsDeform];
        }
        TangentVecType normalVec = tangentVec1.cross( tangentVec2 );
        RealType norm = normalVec.norm();
        normalVec /= norm;
        for( int comp = 0; comp<3; ++comp ){
            normalDeformed_LagrangeInterpolationP1[nodeIdx + comp * numVertices] = normalVec[comp];
            D1xB_LagrangeInterpolationP1[nodeIdx + comp * numVertices] = tangentVec1[comp];
            D2xB_LagrangeInterpolationP1[nodeIdx + comp * numVertices] = tangentVec2[comp];
        }
        DxBAtNodes[nodeIdx].col(0) = tangentVec1; DxBAtNodes[nodeIdx].col(1) = tangentVec2;
    }
    DKTFEVectorFunctionEvaluator<ConfiguratorTypeP1> normalDeformed_LagrangeInterpolationP1_DFD( confP1, normalDeformed_LagrangeInterpolationP1, 3 ),
    D1xB_LagrangeInterpolationP1_DFD( confP1, D1xB_LagrangeInterpolationP1, 3 ),
    D2xB_LagrangeInterpolationP1_DFD( confP1, D2xB_LagrangeInterpolationP1, 3 );
    for( int elementIdx=0; elementIdx<mesh.getNumElements(); ++elementIdx ){
        const typename MeshType::ElementType& El ( mesh.getTriang( elementIdx ) );
        const typename ConfiguratorTypeP1::BaseFuncSetType &bfs = confP1.getBaseFunctionSet ( El );
        const int numQuadPoints = bfs.numQuadPoints( );
        RealType a = 0.;
        for ( int q = 0; q < numQuadPoints; ++q ){
            Matrix32 Dn;
            normalDeformed_LagrangeInterpolationP1_DFD.evaluateGradientAtQuadPoint ( El, q, Dn  );
            typename DataTypeContainer::Point3DType D1xB, D2xB;
            D1xB_LagrangeInterpolationP1_DFD.evaluateAtQuadPoint ( El, q, D1xB  );
            D2xB_LagrangeInterpolationP1_DFD.evaluateAtQuadPoint ( El, q, D2xB  );
            Matrix32 DxB; DxB.col(0) = D1xB; DxB.col(1) = D2xB;
            a += (DxB.transpose() * Dn ).determinant() * bfs.getWeight ( q ); 
        }
        GaussCurvVector_Alternative[elementIdx] = El.getAreaOfRefTriangle() * std::abs(a);
    }
    for( int elIdx=0; elIdx < mesh.getNumElements(); ++elIdx ) 
        avGaussCurvVector_Alternative[elIdx] = GaussCurvVector_Alternative[elIdx] / areaOnElementsVec[elIdx];
    
    //Evaluation of GaussMap: n: M_B \to S^2
    //elementwise \int_El sqrt( det( Dn^T Dn)
    //use Lagrangeinterpolation of n
    VectorType GaussMapVector( totalStressVec.size() );
    VectorType avGaussMapVector( totalStressVec.size() );
    for( int elementIdx=0; elementIdx<mesh.getNumElements(); ++elementIdx ){
        const typename MeshType::ElementType& El ( mesh.getTriang( elementIdx ) );
        const typename ConfiguratorTypeP1::BaseFuncSetType &bfs = confP1.getBaseFunctionSet ( El );
        const int numQuadPoints = bfs.numQuadPoints( );
        RealType a = 0.;
        for ( int q = 0; q < numQuadPoints; ++q ){
            Matrix32 Dn;
            normalDeformed_LagrangeInterpolationP1_DFD.evaluateGradientAtQuadPoint ( El, q, Dn  );
            a += (Dn.transpose() * Dn ).determinant() * bfs.getWeight ( q ); 
        }
        GaussMapVector[elementIdx] = El.getAreaOfRefTriangle() * std::abs(a);
    }
    for( int elIdx=0; elIdx < mesh.getNumElements(); ++elIdx ) 
        avGaussMapVector[elIdx] = GaussMapVector[elIdx] / areaOnElementsVec[elIdx];
    
    //Alternative Evaluation of Gauss Curvature by convolution wiht nodalwise hat function
    std::vector<RealType> hatKernelWeightVec ( numVertices );
    std::vector<Matrix22> convolutedShapeOpVec ( numVertices ); 
    for( int nodeIdx=0; nodeIdx<numVertices; ++nodeIdx ){
        hatKernelWeightVec[nodeIdx] = 0.;
        convolutedShapeOpVec[nodeIdx].setZero();
    }
    for( int elementIdx=0; elementIdx<mesh.getNumElements(); ++elementIdx ){
        
        const typename MeshType::ElementType& El ( mesh.getTriang( elementIdx ) );
        const typename ConfiguratorTypeP1::BaseFuncSetType &bfs = confP1.getBaseFunctionSet ( El );
        const int numQuadPoints = bfs.numQuadPoints( );
        
        for ( int q = 0; q < numQuadPoints; ++q ){
            
            const Matrix22& gInv = xBStorage.getFirstFFInv( El.getGlobalElementIdx(), q );
            const Matrix22& h = xBStorage.getSecondFF( El.getGlobalElementIdx(), q );
            const Matrix22 ShapeOpB = gInv * h; 
            
            for( int bfNum=0; bfNum<3; ++bfNum ){    
                const int globalNodeIdx = confP1.localToGlobal ( El, bfNum );
                const RealType bf = bfs.evaluateOnRefTriang ( bfNum, q );
                convolutedShapeOpVec[globalNodeIdx] += El.getAreaOfRefTriangle() * bf * bfs.getWeight ( q ) * ShapeOpB;
                hatKernelWeightVec[globalNodeIdx] += El.getAreaOfRefTriangle() * bf * bfs.getWeight ( q );
            }
        }
    }
    for( int nodeIdx=0; nodeIdx<numVertices; ++nodeIdx ){
        convolutedShapeOpVec[nodeIdx] /= hatKernelWeightVec[nodeIdx];
    }
    VectorType convolutedGaussCurvVec ( numVertices );
    for( int nodeIdx=0; nodeIdx<numVertices; ++nodeIdx ){
        convolutedGaussCurvVec[nodeIdx] = convolutedShapeOpVec[nodeIdx].determinant();
    }
    //EigenVectors and Eigenvalues of ShapeOperator
    VectorType EigenValues1ShapeOp (numVertices), EigenValues2ShapeOp (numVertices), EigenValuesSmallestAbsShapeOp (numVertices);
    std::vector<TangentVecType> EigenVecs1ShapeOp (numVertices), EigenVecs2ShapeOp(numVertices), EigenVecsSmallestAbsShapeOp (numVertices); 
    for( int nodeIdx=0; nodeIdx<numVertices; ++nodeIdx ){
        Eigen::SelfAdjointEigenSolver<Matrix22> eigenSolver;
        eigenSolver.compute( convolutedShapeOpVec[nodeIdx] );
        EigenValues1ShapeOp[nodeIdx] = eigenSolver.eigenvalues()[0];
        EigenValues2ShapeOp[nodeIdx] = eigenSolver.eigenvalues()[1];
        DomVecType eigenVec1 = eigenSolver.eigenvectors().col(0);
        DomVecType eigenVec2 = eigenSolver.eigenvectors().col(1);
        EigenVecs1ShapeOp[nodeIdx] = DxBAtNodes[nodeIdx] * eigenVec1;
        EigenVecs2ShapeOp[nodeIdx] = DxBAtNodes[nodeIdx] * eigenVec2;
        if( std::abs( EigenValues1ShapeOp[nodeIdx] ) < std::abs( EigenValues1ShapeOp[nodeIdx] ) ){
            EigenValuesSmallestAbsShapeOp[nodeIdx] = std::abs(EigenValues1ShapeOp[nodeIdx]);
            EigenVecsSmallestAbsShapeOp[nodeIdx] = EigenVecs1ShapeOp[nodeIdx];
        }else{
            EigenValuesSmallestAbsShapeOp[nodeIdx] = std::abs(EigenValues2ShapeOp[nodeIdx]);
            EigenVecsSmallestAbsShapeOp[nodeIdx] = EigenVecs2ShapeOp[nodeIdx];
        }
    }
    //convert
    VectorType EigenVecs1ShapeOp_convert (3 * numVertices), EigenVecs2ShapeOp_convert( 3 * numVertices), EigenVecsSmallestAbsShapeOp_convert ( 3 * numVertices); 
    for( int nodeIdx=0; nodeIdx<numVertices; ++nodeIdx ){
      for(int comp=0; comp<3;++comp){
          EigenVecs1ShapeOp_convert[nodeIdx + comp*numVertices] = EigenVecs1ShapeOp[nodeIdx][comp];
          EigenVecs2ShapeOp_convert[nodeIdx + comp*numVertices] = EigenVecs2ShapeOp[nodeIdx][comp];
          EigenVecsSmallestAbsShapeOp_convert[nodeIdx + comp*numVertices] = EigenVecsSmallestAbsShapeOp[nodeIdx][comp];
      }
    }
    //compute L1 norm of convoluted Gauss curvature 
    RealType convolutedGaussCurvL1 = 0.;
    ScalarFctLpNormToP<ConfiguratorType,ConfiguratorTypeP1> ( confP1, shellHandler.getChartToUndeformedShell_Cache(), convolutedGaussCurvVec, 1. ).assembleAdd( convolutedGaussCurvL1 );
    isometryInfo.setConvGaussCurvatureL1( convolutedGaussCurvL1 );
    
    
    //Second Derivative of xB on Elements
    VectorType D2xBSqrVector ( totalStressVec.size() ), avD2xBSqrVector ( totalStressVec.size() );
    RealType integralD2uSqr = 0.;
    for( int elementIdx=0; elementIdx<mesh.getNumElements(); ++elementIdx ){
        const typename MeshType::ElementType& El ( mesh.getTriang( elementIdx ) );
        const typename ConfiguratorType::BaseFuncSetType &bfs = conf.getBaseFunctionSet ( El );
        const int numQuadPoints = bfs.numQuadPoints( );
        RealType a = 0.;
        for ( int q = 0; q < numQuadPoints; ++q ){
            typename DataTypeContainer::Matrix22 D2xB_1, D2xB_2, D2xB_3; D2xB_1.setZero(); D2xB_2.setZero(); D2xB_3.setZero();
            for ( int bfNum = 0; bfNum < static_cast<int> ( conf.getNumLocalDofs ( El ) ); ++bfNum ) {
                typename DataTypeContainer::Matrix22 v;
                v = bfs.evaluateHessianOnRefTriang ( bfNum, q );
                D2xB_1 += xB[ conf.localToGlobal ( El, bfNum ) ] * v; 
                D2xB_2 += xB[ conf.localToGlobal ( El, bfNum ) + numGlobalDofsDeform] * v; 
                D2xB_3 += xB[ conf.localToGlobal ( El, bfNum ) + 2*numGlobalDofsDeform] * v; 
            }
            
            a += (pesopt::ddProd<RealType,typename DataTypeContainer::Matrix22>( D2xB_1, D2xB_1 ) + pesopt::ddProd<RealType,typename DataTypeContainer::Matrix22>( D2xB_2, D2xB_2 ) + pesopt::ddProd<RealType,typename DataTypeContainer::Matrix22>( D2xB_3, D2xB_3 ) ) * bfs.getWeight ( q );
        }
        D2xBSqrVector[elementIdx] = El.getAreaOfRefTriangle() * std::abs(a);
        integralD2uSqr += D2xBSqrVector[elementIdx];
    }
    for( int elIdx=0; elIdx < mesh.getNumElements(); ++elIdx ) avD2xBSqrVector[elIdx] = D2xBSqrVector[elIdx] / areaOnElementsVec[elIdx];
    isometryInfo.setD2uL2( sqrt(integralD2uSqr) );
    
    //Approx Second Derivative of xB on Elements
    VectorType ApproxD2xBSqrVector ( totalStressVec.size() ), avApproxD2xBSqrVector ( totalStressVec.size() );
    RealType integralApproxD2uSqr = 0.;
    for( int elementIdx=0; elementIdx<mesh.getNumElements(); ++elementIdx ){
        const typename MeshType::ElementType& El ( mesh.getTriang( elementIdx ) );
        const typename ConfiguratorType::ApproxGradientBaseFuncSetType &bfs = conf.getApproxGradientBaseFunctionSet ( El );
        const int numQuadPoints = bfs.numQuadPoints( );
        RealType a = 0.;
        for ( int q = 0; q < numQuadPoints; ++q ){
            typename DataTypeContainer::Matrix22 D2xB_1, D2xB_2, D2xB_3; D2xB_1.setZero(); D2xB_2.setZero(); D2xB_3.setZero();
            for ( int bfNum = 0; bfNum < static_cast<int> ( conf.getNumLocalDofs ( El ) ); ++bfNum ) {
                typename DataTypeContainer::Matrix22 v;
                v = bfs.evaluateApproxHessianOnRefTriang ( bfNum, q );
                D2xB_1 += xB[ conf.localToGlobal ( El, bfNum ) ] * v; 
                D2xB_2 += xB[ conf.localToGlobal ( El, bfNum ) + numGlobalDofsDeform] * v; 
                D2xB_3 += xB[ conf.localToGlobal ( El, bfNum ) + 2*numGlobalDofsDeform] * v; 
            }
            
            a += (pesopt::ddProd<RealType,typename DataTypeContainer::Matrix22>( D2xB_1, D2xB_1 ) + pesopt::ddProd<RealType,typename DataTypeContainer::Matrix22>( D2xB_2, D2xB_2 ) + pesopt::ddProd<RealType,typename DataTypeContainer::Matrix22>( D2xB_3, D2xB_3 ) ) * bfs.getWeight ( q );
        }
        ApproxD2xBSqrVector[elementIdx] = El.getAreaOfRefTriangle() * std::abs(a);
        integralApproxD2uSqr += ApproxD2xBSqrVector[elementIdx];
    }
    for( int elIdx=0; elIdx < mesh.getNumElements(); ++elIdx )  avApproxD2xBSqrVector[elIdx] = ApproxD2xBSqrVector[elIdx] / areaOnElementsVec[elIdx];
    isometryInfo.setApproxD2uL2( sqrt(integralApproxD2uSqr) );
    
    
    //for each Element: n_av = 1/3 sum_i=1,2,3 n_i, e = sum_i=1,2,3 |n_i - n_av|^2
    VectorType DiffNormalVector ( totalStressVec.size() ), avDiffNormalVector( totalStressVec.size() );
    for( int elementIdx=0; elementIdx<mesh.getNumElements(); ++elementIdx ){
        const typename MeshType::ElementType& El ( mesh.getTriang( elementIdx ) );
        
        int nodeIdx0 = El.getGlobalNodeIdx(0);
        typename DataTypeContainer::TangentVecType tangentVec1_nodeIdx0, tangentVec2_nodeIdx0;
        for( int comp = 0; comp < 3; ++comp ) {
                tangentVec1_nodeIdx0[comp] = xB[nodeIdx0 +     numVertices + comp * numGlobalDofsDeform];
                tangentVec2_nodeIdx0[comp] = xB[nodeIdx0 + 2 * numVertices + comp * numGlobalDofsDeform];
        }
        typename DataTypeContainer::TangentVecType normalVec_nodeIdx0 = tangentVec1_nodeIdx0.cross( tangentVec2_nodeIdx0 );
        RealType norm_nodeIdx0 = normalVec_nodeIdx0.norm();
        normalVec_nodeIdx0 /= norm_nodeIdx0;
        
        int nodeIdx1 = El.getGlobalNodeIdx(1);
        typename DataTypeContainer::TangentVecType tangentVec1_nodeIdx1, tangentVec2_nodeIdx1;
        for( int comp = 0; comp < 3; ++comp ) {
                tangentVec1_nodeIdx1[comp] = xB[nodeIdx1 +     numVertices + comp * numGlobalDofsDeform];
                tangentVec2_nodeIdx1[comp] = xB[nodeIdx1 + 2 * numVertices + comp * numGlobalDofsDeform];
        }
        typename DataTypeContainer::TangentVecType normalVec_nodeIdx1 = tangentVec1_nodeIdx1.cross( tangentVec2_nodeIdx1 );
        RealType norm_nodeIdx1 = normalVec_nodeIdx1.norm();
        normalVec_nodeIdx1 /= norm_nodeIdx1;
        
        int nodeIdx2 = El.getGlobalNodeIdx(2);
        typename DataTypeContainer::TangentVecType tangentVec1_nodeIdx2, tangentVec2_nodeIdx2;
        for( int comp = 0; comp < 3; ++comp ) {
                tangentVec1_nodeIdx2[comp] = xB[nodeIdx2 +     numVertices + comp * numGlobalDofsDeform];
                tangentVec2_nodeIdx2[comp] = xB[nodeIdx2 + 2 * numVertices + comp * numGlobalDofsDeform];
        }
        typename DataTypeContainer::TangentVecType normalVec_nodeIdx2 = tangentVec1_nodeIdx2.cross( tangentVec2_nodeIdx2 );
        RealType norm_nodeIdx2 = normalVec_nodeIdx2.norm();
        normalVec_nodeIdx2 /= norm_nodeIdx2;
        
        typename DataTypeContainer::TangentVecType normalVec_avg = 1./3. * ( normalVec_nodeIdx0 + normalVec_nodeIdx1 + normalVec_nodeIdx2 );
        
        RealType e = 0.;
        e += (normalVec_nodeIdx0 - normalVec_avg).squaredNorm();
        e += (normalVec_nodeIdx1 - normalVec_avg).squaredNorm();
        e += (normalVec_nodeIdx2 - normalVec_avg).squaredNorm();
        
        DiffNormalVector[elementIdx] = e; 
    }
    for( int elIdx=0; elIdx < mesh.getNumElements(); ++elIdx )  avDiffNormalVector[elIdx] = DiffNormalVector[elIdx] / areaOnElementsVec[elIdx];
    
    
    //save and plot
    std::vector<VectorType> CurvVector = { isoErrorL2SqrOnElements, avIsoErrorL2SqrOnElements,
                                           isoErrorL1OnElements, avIsoErrorL1OnElements,
                                           GaussCurvVector, MeanCurvVector, avGaussCurvVector, avMeanCurvVector,
                                           GaussCurvVector_Alternative, avGaussCurvVector_Alternative,
                                           GaussMapVector, avGaussMapVector,
                                           RelShapeOpVector, avRelShapeOpVector,
                                           averagedBendingStressVecWithoutMaterialFac,
                                           D2xBSqrVector, avD2xBSqrVector, ApproxD2xBSqrVector, avApproxD2xBSqrVector,
                                           DiffNormalVector, avDiffNormalVector
    };
    std::vector<string> nameCurvVector = { "IsometryErrorL2Sqr", "avIsometryErrorL2Sqr",
                                           "IsometryErrorL1", "avIsometryErrorL1",
                                           "GaussCurv", "MeanCurv", "avGaussCurv", "avMeanCurv",
                                           "GaussCurvAlternative", "avGaussCurvAlternative",
                                           "GaussMap", "avGaussMap",
                                           "RelShapeOp", "avRelShapeOp",
                                           "averagedBendingStressWithoutMaterialFacOnElements",
                                           "D2uSqr", "avD2uSqr", "ApproxD2uSqr", "avApproxD2uSqr",
                                           "DiffNormalVector", "avDiffNormalVector"
    };
    ShellWithMaterialPlotter<ConfiguratorType,ConfiguratorTypePf::_ShellFEType> shellPlotter( matOptConf._conf, shellHandler.getChartToUndeformedShell(), shellHandler.getDirichletMask(), saveDirectory, parser.template get<string>("saving.VTKFileType") );
//  TODO instead of update mesh use new ShellWithMaterialPlotter
//     this->saveStressOnUndeformedShell( "Isometry", shellPlotter, material, shellHandler.getChartToUndeformedShell(), solDisp, CurvVector, nameCurvVector );
//     this->plotStress( saveDirectory, "Isometry", nameCurvVector );
//     //save nodalwise values
//     shellPlotter.clearData();
//     shellPlotter.updateMesh( "undeform", shellHandler.getChartToUndeformedShell() );
//     shellPlotter.addScalarData ( convolutedGaussCurvVec, "convolutedGaussCurv", VERTEX_DATA );
//     shellPlotter.addScalarData ( EigenValues1ShapeOp, "EigenValues1ShapeOp", VERTEX_DATA );
//     shellPlotter.addScalarData ( EigenValues2ShapeOp, "EigenValues2ShapeOp", VERTEX_DATA );
//     shellPlotter.addScalarData ( EigenValuesSmallestAbsShapeOp, "EigenValuesSmallestAbsShapeOp", VERTEX_DATA );
//     shellPlotter.saveShellToFile( "IsometryNodalWise" );
//     //
//     shellPlotter.clearData();
//     shellPlotter.updateMesh( mesh );
//     shellPlotter.updateMesh( "disp", solDisp );
//     shellPlotter.addScalarData ( convolutedGaussCurvVec, "convolutedGaussCurv", VERTEX_DATA );
//     shellPlotter.addScalarData ( EigenValues1ShapeOp, "EigenValues1ShapeOp", VERTEX_DATA );
//     shellPlotter.addScalarData ( EigenValues2ShapeOp, "EigenValues2ShapeOp", VERTEX_DATA );
//     shellPlotter.addScalarData ( EigenValuesSmallestAbsShapeOp, "EigenValuesSmallestAbsShapeOp", VERTEX_DATA );
//     shellPlotter.addVectorData ( EigenVecs1ShapeOp_convert, 3, "EigenVecs1ShapeOp", VERTEX_DATA );
//     shellPlotter.addVectorData ( EigenVecs2ShapeOp_convert, 3, "EigenVecs2ShapeOp", VERTEX_DATA );
//     shellPlotter.addVectorData ( EigenVecsSmallestAbsShapeOp_convert, 3, "EigenVecsSmallestAbsShapeOp", VERTEX_DATA );
//     shellPlotter.saveShellToFile ( "IsometryNodalWiseDeformed" );    
    
    
    //image of gauss map 
    auto data = vtkSmartPointer<vtkPolyData>::New();
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto vertices = vtkSmartPointer<vtkCellArray>::New();
    for( int nodeIter=0; nodeIter<numVertices; ++nodeIter ){ 
        vtkIdType pid[1];
        pid[0] = points->InsertNextPoint ( normalDeformed_LagrangeInterpolationP1[nodeIter],
                                  normalDeformed_LagrangeInterpolationP1[nodeIter + 1 * numVertices], 
                                  normalDeformed_LagrangeInterpolationP1[nodeIter + 2 * numVertices] );
        vertices->InsertNextCell ( 1, pid );
    }
    data->SetPoints ( points );
    data->SetVerts(vertices);
    
    string vtkExtension;
    VTKSaver vtkWriter;
    vtkWriter.saveVTKDataSet( pesopt::strprintf( "%s/imageGaussMapAtVertices.vtk", parser.template get<string> ("saving.saveDirectory").c_str() ).c_str(), data, vtkExtension, VTKPOLYDATA );
        
    //END TEST CURVATURE

    
    //!plot results
    if( parser.template get<int> ( "saving.plotResults" ) == 1 ){
        cout << endl << "plot results" << endl;
        this->plotResults( saveDirectory, nameStressOnElementVec );
        runPdflatexWithBashFile( "AllResults.tex", saveDirectory );
    }
    
  }//end computeOptimalDeformation_OnCurrentMesh


  
  
  
  void computeOptimalDeformation_Adaptive ( MeshType &mesh,
                                            VectorType& material,
                                            VectorType &initDisp,
                                            const int numAdaptiveRefinementSteps ) const {
                                                
    pesopt::consoleOutput( "compute optimal design with adaptive refinement" );
    
    //! informations to store for all refinementsteps
    std::vector<VectorType> SolutionDisplamentVec;
    std::vector<RealType> gridSizeVec;
    std::vector<VectorType> membraneStressVector, bendingStressVector, totalStressVector;
    std::vector<VectorType> GaussCurvVector;
    std::vector<DeformationOptimizationShellEnergyInfo<RealType>> energyInfoVec;
    std::vector<IsometryInfo<RealType>> isometryInfoVec;
    std::vector<string> saveDirectoryVec;
    
    
    ConfiguratorType conf ( mesh );
    ShellHandlerInterface<ConfiguratorType> shellHandler( conf, this->_parser.template get<string>( "InputMesh.chartXAType" ), static_cast<ShellBoundaryType>( this->_parser.template get<int>( "InputMesh.ShellType" ) ), this->_parser.template get<bool> ( "InputMesh.ClampedBoundaryCondition" ) );
    MaskType DirichletMask ( shellHandler.getDirichletMask() );
    if( this->_parser.template get<int> ( "InputMesh.tangentSpaceType" ) == 2 ){
        mesh.generateApproximativeTangentSpaceAtNodes( DirichletMask ); //! \todo should be done directly in triangleMesh?
        mesh.updateAllProjectionCoefficients(); 
    }
    
    
    //! use adaptive refinement
    for( int refinementStep=0; refinementStep <= numAdaptiveRefinementSteps; ++refinementStep ) {
        
      const string saveDirectoryRefinementStep = this->_parser.createSubDirectory( "refinementStep" + std::to_string(refinementStep) );
      saveDirectoryVec.push_back( saveDirectoryRefinementStep );
        
      VectorType InitDisplacement, SolutionDisplacement;
      MeshType oldMesh( mesh ); VectorType markedElements;  
      
      //! MARK and REFINE
      if( refinementStep > 0 ){
           cout << endl << "start mark and refine" << endl;
           auto startTime = std::chrono::high_resolution_clock::now();         
           InitDisplacement.resize( SolutionDisplamentVec[refinementStep - 1].size() ); InitDisplacement = SolutionDisplamentVec[refinementStep - 1];
           SolutionDisplacement.resize( InitDisplacement.size() ); SolutionDisplacement = InitDisplacement;
          if( this->_parser.template get<bool> ( "ConstraintProblem.adaptiveMarkingUseTotalStressVec") ){
              this->template markAndRefineForGivenElementErrorVector<LinElastEnergyType::_DiscreteFunctionCacheType> ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask, totalStressVector[refinementStep-1], markedElements  ); 
          }else{
          if( this->_parser.template get<bool> ( "ConstraintProblem.adaptiveMarkingUseGaussCurvatureVec") ){
              this->template markAndRefineForGivenElementErrorVector<LinElastEnergyType::_DiscreteFunctionCacheType> ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask, GaussCurvVector[refinementStep-1], markedElements  ); 
          }else{
              this->template markAndRefine<LinElastEnergyType::_DiscreteFunctionCacheType> ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask, markedElements );
          }
          }
          std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
          std::cout << endl << "duration for markAndRefine = " << diff.count() << " sec" << endl;
      }else{
          InitDisplacement.resize( initDisp.size() ); InitDisplacement = initDisp;
          SolutionDisplacement.resize( InitDisplacement.size() ); SolutionDisplacement = InitDisplacement;
      }

      //! set initialmaterial on refined mesh
      VectorType initMaterial ( material.size() );
      if( refinementStep == 0 ) this->setInitialMaterial( this->_parser.template get<int> ( "ConstraintProblem.firstStepMaterialInitializationType" ), mesh, initMaterial, material  );
      else                      this->setInitialMaterial( this->_parser.template get<int> ( "ConstraintProblem.adaptiveMaterialInitializationType" ), mesh, initMaterial, material  );
      material = initMaterial;
      
      //! SOLVE
      VectorType membraneStressVec ( mesh.getNumElements() ), bendingStressVec ( mesh.getNumElements() ), totalStressVec( mesh.getNumElements() );
      VectorType GaussCurvVec ( mesh.getNumElements() );
      DeformationOptimizationShellEnergyInfo<RealType> energyInfo;
      IsometryInfo<RealType> isometryInfo;
      computeOptimalDeformation_OnCurrentMesh ( mesh, material, InitDisplacement, SolutionDisplacement, DirichletMask, 
                                                membraneStressVec, bendingStressVec, totalStressVec,
                                                GaussCurvVec,
                                                energyInfo, isometryInfo, saveDirectoryRefinementStep, refinementStep ); 
      
      //! SAVE
      SolutionDisplamentVec.push_back( SolutionDisplacement );
      gridSizeVec.push_back( DKTTriangMeshInfo<MeshType> ( mesh ).getMinAreaSqrt() ); 
      membraneStressVector.push_back( membraneStressVec ); bendingStressVector.push_back( bendingStressVec ); totalStressVector.push_back( totalStressVec );
      GaussCurvVector.push_back( GaussCurvVec );
      energyInfoVec.push_back( energyInfo );
      isometryInfo.setGridSize( gridSizeVec[refinementStep] );
      isometryInfoVec.push_back( isometryInfo );
      isometryInfo.template saveToFile<ParameterParserType>( "isometryInfo", saveDirectoryRefinementStep );
      
    }// end for refinement step
 
 
 
    initDisp.resize( SolutionDisplamentVec[numAdaptiveRefinementSteps].size() ); initDisp = SolutionDisplamentVec[numAdaptiveRefinementSteps];
   
    this->template plotConvergenceOfApdaptiveRefinement<IsometryConstraint>( energyInfoVec );
    
    VectorType fineSolution = SolutionDisplamentVec[numAdaptiveRefinementSteps];
    std::vector<VectorType> SolutionDisplamentVecProlongated;
    std::vector<RealType> ErrorD2L2_FineSolutionVec(numAdaptiveRefinementSteps+1), EOCD2L2_FineSolutionVec(numAdaptiveRefinementSteps);
    for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
            VectorType currentSolutionProlongated = SolutionDisplamentVec[refinementStep]; 
            mesh.prolongateVectorFct_DKTLinearly( currentSolutionProlongated );
            SolutionDisplamentVecProlongated.push_back( currentSolutionProlongated );
    }
    ConfiguratorType confFine ( mesh );
    for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
        DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> fineSolutionDFD ( conf, fineSolution, 3);
        DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> solDFD( conf, SolutionDisplamentVecProlongated[refinementStep], 3 );
        RealType tmp = 0.;
        SecondDerivativeEnergy<ConfiguratorType> ( conf, fineSolutionDFD, solDFD ).assembleAdd( tmp );
        ErrorD2L2_FineSolutionVec[refinementStep] = sqrt( tmp );
        isometryInfoVec[refinementStep].setErrorApproxD2uToFineSolutionL2( sqrt(tmp) );
    } 
    this->plotConvergenceIsometryOfApdaptiveRefinement( isometryInfoVec );
    
    
    if( this->_parser.template get<bool> ("saving.removeOldRefinementSteps") ){
        for( int refinementStep = 0; refinementStep < saveDirectoryVec.size() - 1; ++refinementStep ) pesopt::deleteDirectory( saveDirectoryVec[refinementStep] );
        const string saveDirectoryFinal = this->_parser.template get<string> ("saving.saveDirectory") + "/FinalResult";
        pesopt::moveDirectory( saveDirectoryVec[saveDirectoryVec.size()-1], saveDirectoryFinal );
    }else{
        const string saveDirectoryFinal = this->_parser.template get<string> ("saving.saveDirectory") + "/FinalResult";
        pesopt::copyDirectory( saveDirectoryVec[saveDirectoryVec.size()-1], saveDirectoryFinal );
    }
    
    
  }//end computeOptimalDeformation_Adaptive
    
};
    


#endif
