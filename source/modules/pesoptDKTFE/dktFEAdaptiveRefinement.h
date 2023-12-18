#ifndef __DKTFEADAPTIVEREFINEMENT_H
#define __DKTFEADAPTIVEREFINEMENT_H

#include <pesopt_IO.h>
#include <dktFEHandler.h>
#include <dktFEFunctionEvaluator.h>
#include <dktFEQuadrature.h>
#include <dktFEConfigurators.h>

// computes vector with \int_El 1
template <typename ConfiguratorType, DiscreteFunctionCacheType _DiscreteVectorFunctionCacheType >
void evaluateAreaOfElements( const ConfiguratorType &conf,
                             const ShellHandler<ConfiguratorType,_DiscreteVectorFunctionCacheType> &shellHandler, 
                             typename ConfiguratorType::VectorType &areaElementsVec ) {
    const int numElements = conf.getInitializer().getNumTriangs();
    areaElementsVec.resize( numElements ); areaElementsVec.setZero();
    for( int elementIdx=0; elementIdx<numElements; ++elementIdx ){
           const typename ConfiguratorType::ElementType& El ( conf.getInitializer().getTriang( elementIdx ) );
           const typename ConfiguratorType::BaseFuncSetType &bfs = conf.getBaseFunctionSet ( El );
           const int numQuadPoints = bfs.numQuadPoints( );
           typename ConfiguratorType::RealType areaElement = 0.;
           for ( int q = 0; q < numQuadPoints; ++q ){
               areaElement += shellHandler.getChartToUndeformedShell_Cache().getArea(elementIdx,q) * bfs.getWeight ( q );
           }
           areaElementsVec[elementIdx] = El.getAreaOfRefTriangle() * areaElement;
    }
}
    
    
template<typename RealType, typename PointType>
void computeBarycentricCoordinates( const PointType & pM, const PointType &p1, const PointType &p2, RealType &l1, RealType &l2 ) {
        PointType p2M = pM - p2;
        PointType p21 = p1 - p2;
        l1 = p2M.norm() / p21.norm();
        l2 = 1. - l1;
}
    
template<typename RealType, typename PointType>
void computeBarycentricCoordinates( const PointType & pM, const PointType &p1, const PointType &p2, const PointType &p3, RealType &l1, RealType &l2, RealType &l3 ) {
        PointType p12; p12 = p2 - p1; PointType p13; p13 = p3 - p1; PointType p1M; p1M = pM - p1;
        RealType d00 = p12.dot(p12);
        RealType d01 = p12.dot(p13);
        RealType d11 = p13.dot(p13);
        RealType d20 = p1M.dot(p12);
        RealType d21 = p1M.dot(p13);
        RealType invDenom = 1.0 / (d00 * d11 - d01 * d01);
        l2 = (d11 * d20 - d01 * d21) * invDenom;
        l3 = (d00 * d21 - d01 * d20) * invDenom;
        l1 = 1.0 - l2 - l3;
}


template <typename MeshType, typename VectorType = typename MeshType::VectorType>
void markElements_ElementErrorVec( MeshType &mesh, const VectorType &errorOnElements, const typename MeshType::RealType thresholdPercentil = 0.5 ) { 
  typedef typename MeshType::RealType RealType;
  std::vector<RealType> errorOnElementsSorted ( errorOnElements.size() );
  for( int i=0; i<errorOnElements.size(); ++i ) errorOnElementsSorted[i] = errorOnElements[i];
  std::sort( errorOnElementsSorted.begin(), errorOnElementsSorted.end() );
  RealType threshold = errorOnElementsSorted[ static_cast<int>( thresholdPercentil * errorOnElements.size() ) ];  
  for ( int elementIdx = 0; elementIdx < mesh.getNumTriangs(); ++elementIdx  ){
        if( errorOnElements[elementIdx] > threshold ) mesh.mark( elementIdx );
  }
}


template <typename ConfiguratorTypePf>
void markElements_L2gradient( ConfiguratorTypePf& confpf, typename ConfiguratorTypePf::InitType& mesh, typename ConfiguratorTypePf::VectorType& v, 
                              typename ConfiguratorTypePf::RealType threshold = 0.5 ){
  
  DKTFEScalarFunctionEvaluator<ConfiguratorTypePf> vDFD ( confpf, v );
  typename ConfiguratorTypePf::DomVecType pfgrad;
  for ( int elementIdx = 0; elementIdx < mesh.getNumTriangs(); ++elementIdx  ){
        vDFD.evaluateGradientAtQuadPoint ( mesh.getTriang(elementIdx), 0, pfgrad ); //v linear on each triangle, so can evaluate gradient at quadpoint 0
        if( pfgrad.squaredNorm() > threshold ) mesh.mark( elementIdx );
  }
}



template <typename ConfiguratorType>
void markElements_L2ApproxHessian( ConfiguratorType& conf, typename ConfiguratorType::InitType& mesh, 
                                   const DiscreteVectorFunctionStorage<ConfiguratorType, FirstAndSecondOrder> &xAStorage,
                                   typename ConfiguratorType::VectorType& disp,
                                   typename ConfiguratorType::RealType thresholdPercentil = 0.5 ){
    
  DKTFEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( conf, disp, 3 );
  typename ConfiguratorType::Matrix33 hessianAsVec;
  typename ConfiguratorType::VectorType errorOnElements ( mesh.getNumTriangs() );
  typename ConfiguratorType::VectorType areaOfElements ( mesh.getNumTriangs() );
  for ( int elementIdx = 0; elementIdx < mesh.getNumTriangs(); ++elementIdx  ){
      typename ConfiguratorType::RealType normOnElement = 0.;
      typename ConfiguratorType::RealType areaOnElement = 0.;
      const typename ConfiguratorType::ElementType& El ( mesh.getTriang( elementIdx ) );
      const typename ConfiguratorType::BaseFuncSetType &bfs = conf.getBaseFunctionSet ( El );
      for( int q = 0; q < conf.maxNumQuadPoints( ); ++q ){
        dispDFD.evaluateApproxHessianAsVecAtQuadPoint ( mesh.getTriang(elementIdx), q, hessianAsVec );
        normOnElement += xAStorage.getArea( elementIdx, q ) * pesopt::ddProd<typename ConfiguratorType::RealType, typename ConfiguratorType::Matrix33>(hessianAsVec,hessianAsVec) * bfs.getWeight(q);
        areaOnElement += xAStorage.getArea( elementIdx, q ) * pesopt::ddProd<typename ConfiguratorType::RealType, typename ConfiguratorType::Matrix33>(hessianAsVec,hessianAsVec) * bfs.getWeight(q);
      }
      errorOnElements[elementIdx] = El.getAreaOfRefTriangle() * normOnElement;
      areaOfElements[elementIdx] = El.getAreaOfRefTriangle() * areaOnElement;
  }
  markElements_ElementErrorVec<typename ConfiguratorType::InitType>( mesh, errorOnElements, thresholdPercentil );
}


template <typename MeshType>
void markClampedElements( MeshType& mesh, const typename MeshType::MaskType& mask ){
  
  for ( int elementIdx = 0; elementIdx < mesh.getNumTriangs(); ++elementIdx  ){
        const bool n1OnBdr = mask[mesh.getTriangNodeIdx(elementIdx, 0)];
        const bool n2OnBdr = mask[mesh.getTriangNodeIdx(elementIdx, 1)];
        const bool n3OnBdr = mask[mesh.getTriangNodeIdx(elementIdx, 2)];
        if( n1OnBdr && n2OnBdr && n3OnBdr ) mesh.mark( elementIdx );
  }
}


template< typename AdaptiveMeshType >
void findParentsInOldMesh( const AdaptiveMeshType &oldMesh,                     
                           const AdaptiveMeshType &newMesh,  
                           const int globalNodeIdx,
                           std::vector<int> &nodeIdxVecParents ) {
    typename std::map< int, typename AdaptiveMeshType::ParentInformation >::const_iterator iter;
    iter = newMesh.getInterpolationMap().find( globalNodeIdx );
    for(int p=0; p<2; ++p ){
        const int globalNodeIdxParent = iter->second.globalIndices[p];
        if( globalNodeIdxParent >= oldMesh.getNumVertices() ){
            findParentsInOldMesh<AdaptiveMeshType>( oldMesh, newMesh, globalNodeIdxParent, nodeIdxVecParents );
        }else{
            bool notContainedInParentList = true;
            for(int i=0;i<nodeIdxVecParents.size();i++){
                if( globalNodeIdxParent == nodeIdxVecParents[i]  ) notContainedInParentList = false;
            }
            if(notContainedInParentList) nodeIdxVecParents.push_back(globalNodeIdxParent);
        }
    }
}
    
    
    
    
    
    
    
template< typename AdaptiveMeshType>
void updateRefinedPositions_LinearInterpolationOfNormals( const AdaptiveMeshType &oldMesh, 
                                                          AdaptiveMeshType &newMesh ){
    typedef typename AdaptiveMeshType::RealType RealType;
    typedef typename AdaptiveMeshType::DomVecType DomVecType;
    typedef typename AdaptiveMeshType::Point3DType Point3DType;
    typedef typename AdaptiveMeshType::TangentVecType TangentVecType;
    
    //update of positions and tangent-,normal-space in refined mesh
    const int oldSize = oldMesh.getNumVertices();
    const int newSize = newMesh.getNumVertices();
    typename std::map< int, typename AdaptiveMeshType::ParentInformation >::const_iterator iter;
    
    //find ElementIndex and Refcoords in oldMesh
    std::vector<int> ElementIndexVec ( newSize - oldSize );
    std::vector<DomVecType> RefCoordVec ( newSize - oldSize );
    for( int nodeIdx = oldSize; nodeIdx < newSize; ++nodeIdx ){
        std::vector<int> nodeIdxVecParents;
        findParentsInOldMesh<AdaptiveMeshType>( oldMesh, newMesh, nodeIdx, nodeIdxVecParents );
        std::vector<int> commonElements;
        oldMesh.getCommonElements( nodeIdxVecParents, commonElements );
        //compute barycentric coords
        std::vector<RealType> lambda ( 3 );
        if(nodeIdxVecParents.size() == 2  ){
            computeBarycentricCoordinates<RealType,Point3DType>( newMesh.getVertex(nodeIdx), newMesh.getVertex(nodeIdxVecParents[0]), newMesh.getVertex(nodeIdxVecParents[1]), lambda[0], lambda[1] );
            lambda[2] = 0.;
        }
        if(nodeIdxVecParents.size() == 3 ){
            computeBarycentricCoordinates<RealType,Point3DType>( newMesh.getVertex(nodeIdx), newMesh.getVertex(nodeIdxVecParents[0]), newMesh.getVertex(nodeIdxVecParents[1]), newMesh.getVertex(nodeIdxVecParents[2]), lambda[0], lambda[1], lambda[2] );
        }
        //compute refcoord in commonElements[0]
        DomVecType RefCoordNewPoint; RefCoordNewPoint.setZero();
        for(int i=0; i<nodeIdxVecParents.size(); ++i ){
           //int localIdxParent = oldMesh.getLocalIdxWrtAdjacentElementIdx( nodeIdxVecParents[i], 0 );
           int localIdxParent = oldMesh.getLocalIdxWrtGlobalElementIdx( nodeIdxVecParents[i], commonElements[0] );
           DomVecType RefCoordParent; RefCoordParent.setZero();
           if( localIdxParent == 1 ) RefCoordParent[0] = 1.;
           if( localIdxParent == 2 ) RefCoordParent[1] = 1.;
           RefCoordNewPoint += lambda[i] * RefCoordParent;
        }
        ElementIndexVec[nodeIdx - oldSize] = commonElements[0];
        RefCoordVec[nodeIdx - oldSize] = RefCoordNewPoint;
    }
        
    for( int nodeIdx = oldSize; nodeIdx < newSize; ++nodeIdx ){
        // find new node index in map
        iter = newMesh.getInterpolationMap().find( nodeIdx );
        if( iter == newMesh.getInterpolationMap().end() ) 
            throw std::invalid_argument ( pesopt::strprintf ( "unknown vertex in file %s at line %d.", __FILE__, __LINE__ ).c_str() );
        // use linear interpolation of normals
        if( (iter->second.globalIndices[0] > nodeIdx ) || (iter->second.globalIndices[1] > nodeIdx ) )
            throw std::invalid_argument ( pesopt::strprintf ( "parent indices should be smaller then child??? in file %s at line %d.", __FILE__, __LINE__ ).c_str() );
        TangentVecType normalInterpolated = 0.5 * ( newMesh.getNormalVec(iter->second.globalIndices[0]) + newMesh.getNormalVec(iter->second.globalIndices[1]) );

        normalInterpolated.normalize();
        TangentVecType tv1Int, tv2Int;
        newMesh.getTangentSpaceFromNormal( tv1Int, tv2Int, normalInterpolated );
        // update
        newMesh.setTangentVec1( nodeIdx, tv1Int );  newMesh.setTangentVec2( nodeIdx, tv2Int ); newMesh.setNormalVec( nodeIdx, normalInterpolated );
    }//end for old size
    newMesh.updateAllTriangles();
}
    
    
    
    
    
    
    
    
    
template< typename ConfiguratorType>
void updateRefinedPositions_InterpolationOfNormalsAndFunctionOnNodes( const ConfiguratorType &oldConf, 
                                                                      const typename ConfiguratorType::VectorType &oldChartToUndeformedShell,
                                                                      typename ConfiguratorType::InitType &newMesh ){
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::DomVecType DomVecType;
    typedef typename ConfiguratorType::InitType AdaptiveMeshType;
    typedef typename AdaptiveMeshType::Point3DType Point3DType;
    typedef typename AdaptiveMeshType::TangentVecType TangentVecType;
    
    
#ifdef __DEBUGINTERPOLATIONADAPTIVEMESH
     cout << endl 
          << "================================" << endl
          << "start updateRefinedPositions_InterpolationOfNormalsAndFunctionOnNodes" << endl
          << endl;
#endif
    
    //update of positions and tangent-,normal-space in refined mesh
    const int oldSize = oldConf.getInitializer().getNumVertices();
    const int newSize = newMesh.getNumVertices();
    typename std::map< int, typename AdaptiveMeshType::ParentInformation >::const_iterator iter;
    
#ifdef   __DEBUGINTERPOLATIONADAPTIVEMESHFORPLATE
    for( int nodeIdx = 0; nodeIdx < newSize; ++nodeIdx ){
        const Point3DType coords = newMesh.getVertex(nodeIdx);
        if( std::abs( coords(2) ) > 1.e-16 )
            cout << endl
                 << "ERROR in adaptve refinement : " << endl 
                 << "nodeIdx = " << nodeIdx << endl 
                 << "coords = " << coords.transpose() << endl;
    }
#endif
    
    //find ElementIndex and Refcoords in oldMesh
    std::vector<int> ElementIndexVec ( newSize - oldSize );
    std::vector<DomVecType> RefCoordVec ( newSize - oldSize );
    for( int nodeIdx = oldSize; nodeIdx < newSize; ++nodeIdx ){
        std::vector<int> nodeIdxVecParents;
        findParentsInOldMesh<AdaptiveMeshType>( oldConf.getInitializer(), newMesh, nodeIdx, nodeIdxVecParents );
        std::vector<int> commonElements;
        oldConf.getInitializer().getCommonElements( nodeIdxVecParents, commonElements );
        //compute barycentric coords
        std::vector<RealType> lambda ( 3 );
        if(nodeIdxVecParents.size() == 2  ){
            computeBarycentricCoordinates<RealType,Point3DType>( newMesh.getVertex(nodeIdx), newMesh.getVertex(nodeIdxVecParents[0]), newMesh.getVertex(nodeIdxVecParents[1]), lambda[0], lambda[1] );
            lambda[2] = 0.;
        }
        if(nodeIdxVecParents.size() == 3 ){
            computeBarycentricCoordinates<RealType,Point3DType>( newMesh.getVertex(nodeIdx), newMesh.getVertex(nodeIdxVecParents[0]), newMesh.getVertex(nodeIdxVecParents[1]), newMesh.getVertex(nodeIdxVecParents[2]), lambda[0], lambda[1], lambda[2] );
        }
        //compute refcoord in commonElements[0]
        DomVecType RefCoordNewPoint; RefCoordNewPoint.setZero();
        for(int i=0; i<nodeIdxVecParents.size(); ++i ){
           //int localIdxParent = oldConf.getInitializer().getLocalIdxWrtAdjacentElementIdx( nodeIdxVecParents[i], 0 );
           int localIdxParent = oldConf.getInitializer().getLocalIdxWrtGlobalElementIdx( nodeIdxVecParents[i], commonElements[0] );
           DomVecType RefCoordParent; RefCoordParent.setZero();
           if( localIdxParent == 1 ) RefCoordParent[0] = 1.;
           if( localIdxParent == 2 ) RefCoordParent[1] = 1.;
           RefCoordNewPoint += lambda[i] * RefCoordParent;
        }
        ElementIndexVec[nodeIdx - oldSize] = commonElements[0];
        RefCoordVec[nodeIdx - oldSize] = RefCoordNewPoint;
    }
        
    for( int nodeIdx = oldSize; nodeIdx < newSize; ++nodeIdx ){
        // find new node index in map
        iter = newMesh.getInterpolationMap().find( nodeIdx );
        if( iter == newMesh.getInterpolationMap().end() ) 
            throw std::invalid_argument ( pesopt::strprintf ( "unknown vertex in file %s at line %d.", __FILE__, __LINE__ ).c_str() );
        // use linear interpolation of normals
        if( (iter->second.globalIndices[0] > nodeIdx ) || (iter->second.globalIndices[1] > nodeIdx ) )
            throw std::invalid_argument ( pesopt::strprintf ( "parent indices should be smaller then child??? in file %s at line %d.", __FILE__, __LINE__ ).c_str() );
        TangentVecType normalInterpolated = 0.5 * ( newMesh.getNormalVec(iter->second.globalIndices[0]) + newMesh.getNormalVec(iter->second.globalIndices[1]) );
#ifdef   __DEBUGINTERPOLATIONADAPTIVEMESHFORPLATE
        if( ( std::abs( normalInterpolated(0) ) > 1.e-16) ||  ( std::abs( normalInterpolated(1) ) > 1.e-16) || ( std::abs( normalInterpolated(2) - 1. ) > 1.e-16) )
            cout << endl  << "ERROR in normal " << endl;
#endif
        normalInterpolated.normalize();
        TangentVecType tv1Int, tv2Int;
        newMesh.getTangentSpaceFromNormal( tv1Int, tv2Int, normalInterpolated );
#ifdef __DEBUGINTERPOLATIONADAPTIVEMESH
        cout << endl
             << "==================" << endl;
        cout << "generate normal and tangentspace at new node with index " << nodeIdx << endl;
        cout << "parent0 = " << iter->second.globalIndices[0] << endl
             << "parent1 = " << iter->second.globalIndices[1] << endl;
        cout << "normal parent0        = " << newMesh.getNormalVec(iter->second.globalIndices[0]).transpose() << endl;
        cout << "normal parent1        = " << newMesh.getNormalVec(iter->second.globalIndices[1]).transpose() << endl;
        cout << "(normal0 + normal1)/2 = " << 0.5*( newMesh.getNormalVec(iter->second.globalIndices[0]).transpose() + newMesh.getNormalVec(iter->second.globalIndices[1]).transpose() ) << endl;
        cout << "normalInterpolated    = " << normalInterpolated.transpose() << endl;
        cout << "tangentVec1 = " << tv1Int.transpose() << endl;
        cout << "tangentVec2 = " << tv2Int.transpose() << endl;
#endif  
        //use exact Position in Old Mesh:
        DKTFEVectorFunctionEvaluator<ConfiguratorType> dfd ( oldConf, oldChartToUndeformedShell, 3 );
        Point3DType newVertex; 
        int ElementIndexOldMesh = ElementIndexVec[nodeIdx - oldSize];
        DomVecType RefCoordsOldMesh; RefCoordsOldMesh = RefCoordVec[nodeIdx - oldSize];
        dfd.evaluate( oldConf.getInitializer().getTriang(ElementIndexOldMesh), RefCoordsOldMesh, newVertex );
#ifdef   __DEBUGINTERPOLATIONADAPTIVEMESHFORPLATE
        if( std::abs( newVertex(2) ) > 1.e-16 ) 
            cout << endl 
                 << "ERROR in evaluation of new position: " << endl 
                 << "nodeIdx = " << nodeIdx << endl 
                 << "coords = " << newVertex.transpose() << endl 
                 << "normal = " << normalInterpolated.transpose() << endl
                 << "ElementInOldMesh = " << ElementIndexOldMesh << endl
                 << "RefCoordsInOldMesh = " << RefCoordsOldMesh.transpose() << endl
                 << endl;
#endif 
        newMesh.setVertex( nodeIdx, newVertex );
        // update
        newMesh.setTangentVec1( nodeIdx, tv1Int );  newMesh.setTangentVec2( nodeIdx, tv2Int ); newMesh.setNormalVec( nodeIdx, normalInterpolated );
    }//end for old size
    newMesh.updateAllTriangles();
    
#ifdef __DEBUGINTERPOLATIONADAPTIVEMESH
     cout << endl 
          << "finished updateRefinedPositions_InterpolationOfNormalsAndFunctionOnNodes" << endl
          << "================================" << endl
          << endl;
#endif
}




template< typename ConfiguratorType>
void updateRefinedPositions_ProjectOntoSphere( const ConfiguratorType &oldConf, 
                                               const typename ConfiguratorType::VectorType &oldChartToUndeformedShell,
                                               typename ConfiguratorType::InitType &newMesh ){
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::DomVecType DomVecType;
    typedef typename ConfiguratorType::InitType AdaptiveMeshType;
    typedef typename AdaptiveMeshType::Point3DType Point3DType;
    typedef typename AdaptiveMeshType::TangentVecType TangentVecType;
    
    //update of positions and tangent-,normal-space in refined mesh
    const AdaptiveMeshType &oldMesh = oldConf.getInitializer();
    const int oldSize = oldMesh.getNumVertices();
    const int newSize = newMesh.getNumVertices();
    typename std::map< int, typename AdaptiveMeshType::ParentInformation >::const_iterator iter;
    
    //find ElementIndex and Refcoords in oldMesh
    std::vector<int> ElementIndexVec ( newSize - oldSize );
    std::vector<DomVecType> RefCoordVec ( newSize - oldSize );
    for( int nodeIdx = oldSize; nodeIdx < newSize; ++nodeIdx ){
        std::vector<int> nodeIdxVecParents;
        findParentsInOldMesh<AdaptiveMeshType>( oldConf.getInitializer(), newMesh, nodeIdx, nodeIdxVecParents );
        std::vector<int> commonElements;
        oldConf.getInitializer().getCommonElements( nodeIdxVecParents, commonElements );
        //compute barycentric coords
        std::vector<RealType> lambda ( 3 );
        if(nodeIdxVecParents.size() == 2  ){
            computeBarycentricCoordinates<RealType,Point3DType>( newMesh.getVertex(nodeIdx), newMesh.getVertex(nodeIdxVecParents[0]), newMesh.getVertex(nodeIdxVecParents[1]), lambda[0], lambda[1] );
            lambda[2] = 0.;
        }
        if(nodeIdxVecParents.size() == 3 ){
            computeBarycentricCoordinates<RealType,Point3DType>( newMesh.getVertex(nodeIdx), newMesh.getVertex(nodeIdxVecParents[0]), newMesh.getVertex(nodeIdxVecParents[1]), newMesh.getVertex(nodeIdxVecParents[2]), lambda[0], lambda[1], lambda[2] );
        }
        //compute refcoord in commonElements[0]
        DomVecType RefCoordNewPoint; RefCoordNewPoint.setZero();
        for(int i=0; i<nodeIdxVecParents.size(); ++i ){
           //int localIdxParent = oldConf.getInitializer().getLocalIdxWrtAdjacentElementIdx( nodeIdxVecParents[i], 0 );
           int localIdxParent = oldConf.getInitializer().getLocalIdxWrtGlobalElementIdx( nodeIdxVecParents[i], commonElements[0] );
           DomVecType RefCoordParent; RefCoordParent.setZero();
           if( localIdxParent == 1 ) RefCoordParent[0] = 1.;
           if( localIdxParent == 2 ) RefCoordParent[1] = 1.;
           RefCoordNewPoint += lambda[i] * RefCoordParent;
        }
        ElementIndexVec[nodeIdx - oldSize] = commonElements[0];
        RefCoordVec[nodeIdx - oldSize] = RefCoordNewPoint;
    }
        
    for( int nodeIdx = oldSize; nodeIdx < newSize; ++nodeIdx ){
        // find new node index in map
        iter = newMesh.getInterpolationMap().find( nodeIdx );
        if( iter == newMesh.getInterpolationMap().end() ) throw std::invalid_argument ( pesopt::strprintf ( "unknown vertex in file %s at line %d.", __FILE__, __LINE__ ).c_str() );
        // use linear interpolation of normals
        TangentVecType normalInterpolated = 0.5 * ( newMesh.getNormalVec(iter->second.globalIndices[0]) + newMesh.getNormalVec(iter->second.globalIndices[1]) );
        normalInterpolated.normalize();
        TangentVecType tv1Int, tv2Int;
        newMesh.getTangentSpaceFromNormal( tv1Int, tv2Int, normalInterpolated );  
        //use exact Position of Old Mesh:
        DKTFEVectorFunctionEvaluator<ConfiguratorType> dfd ( oldConf, oldChartToUndeformedShell, 3 );
        Point3DType newVertex; 
        int ElementIndexOldMesh = ElementIndexVec[nodeIdx - oldSize];
        DomVecType RefCoordsOldMesh; RefCoordsOldMesh = RefCoordVec[nodeIdx - oldSize];
        dfd.evaluate( oldConf.getInitializer().getTriang(ElementIndexOldMesh), RefCoordsOldMesh, newVertex );
        newVertex.normalize();
        newMesh.setVertex( nodeIdx, newVertex );
        // update
        newMesh.setTangentVec1( nodeIdx, tv1Int );  newMesh.setTangentVec2( nodeIdx, tv2Int ); newMesh.setNormalVec( nodeIdx, normalInterpolated );
    }//end for old size
    newMesh.updateAllTriangles();
}
     
            
template< typename ConfiguratorType >
void switchProlongationOfMesh( const int ProlongationType,
                               const typename ConfiguratorType::InitType &oldMesh, 
                               typename ConfiguratorType::InitType & newMesh ){    
    
    cout << endl << "start prolongation of mesh" << endl;
    auto startTime = std::chrono::high_resolution_clock::now();      
    
     typedef typename ConfiguratorType::RealType RealType;
     typedef typename ConfiguratorType::InitType MeshType;
     
     switch( ProlongationType ){
          case 0:{
             cout << "do nothing" << endl;
             break;
          }
         //use chart map to interpolate
         case 1:{
             cout << "interpolate with FEM function" << endl;
             const ConfiguratorType oldConf ( oldMesh );
             const ShellHandlerInterface<ConfiguratorType> oldShellHandler( oldConf, "id" );
             updateRefinedPositions_InterpolationOfNormalsAndFunctionOnNodes<ConfiguratorType>( oldConf, oldShellHandler.getChartToUndeformedShell(), newMesh );
             break;
          }
          //use DKT chart to interpolate
//           case 2:{
//               cout << "interpolate with DKT function wrt UnitTriang" << endl;
//               typedef DuodecQuadrature<RealType, typename ConfiguratorType::DomVecType>                                    QuadType;
//               typedef UnitTriangMeshConfiguratorDKT<typename ConfiguratorType::DTContainer, MeshType, QuadType >           ConfiguratorTypeUnitTriangDKT;
//               const ConfiguratorTypeUnitTriangDKT oldConfDKT ( oldMesh );
//               const ShellHandlerInterface<ConfiguratorTypeUnitTriangDKT> oldShellHandler( oldConfDKT, "id" );
//               updateRefinedPositions_InterpolationOfNormalsAndFunctionOnNodes<ConfiguratorTypeUnitTriangDKT>( oldConfDKT, oldShellHandler.getChartToUndeformedShell(), newMesh );
//               break;
//           }
          //project onto sphere
          case 11:{
              //updateRefinedPositions_ProjectOntoSphere<MeshType>( oldMesh, newMesh );
              const ConfiguratorType oldConf ( oldMesh );
              const ShellHandlerInterface<ConfiguratorType> oldShellHandler( oldConf, "id" );
              updateRefinedPositions_ProjectOntoSphere<ConfiguratorType>( oldConf, oldShellHandler.getChartToUndeformedShell(), newMesh );
              break;   
          }
          //use chart map to interpolate and project boundary on unit circle
         case 12:{
             cout << "interpolate with FEM function" << endl;
             const ConfiguratorType oldConf ( oldMesh );
             const ShellHandlerInterface<ConfiguratorType> oldShellHandler( oldConf, "id" );
             updateRefinedPositions_InterpolationOfNormalsAndFunctionOnNodes<ConfiguratorType>( oldConf, oldShellHandler.getChartToUndeformedShell(), newMesh );
             typename MeshType::MaskType topologicalBoundaryMask ( newMesh.getNumVertices() );
             newMesh.findTopologicalBoundaryMask( topologicalBoundaryMask ); 
             for( int nodeIdx = 0; nodeIdx < newMesh.getNumVertices(); ++nodeIdx ){
                   if( topologicalBoundaryMask[nodeIdx] ){
                    typename MeshType::Point3DType newVertex = newMesh.getVertex( nodeIdx ); 
//                     cout << "newVertex.norm = " << newVertex.norm() << endl;
                    newVertex.normalize();
                    newMesh.setVertex( nodeIdx, newVertex );
                   }
             }
             newMesh.updateAllTriangles();
             break;
          }
          default : 
              throw std::invalid_argument ( pesopt::strprintf ( "wrong method %d in file %s at line %d.", ProlongationType,  __FILE__, __LINE__ ).c_str() );
              break;
     }
     
    std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << endl << "duration for prolongation of mesh = " << diff.count() << " sec" << endl;
}

        
template< typename ConfiguratorType >
void switchProlongationOfBoundaryMask( const bool ProlongationType,
                                       const typename ConfiguratorType::InitType & newMesh,
                                       typename ConfiguratorType::MaskType &mask ){  
    const ShellFEType feType = ConfiguratorType::_ShellFEType;
    if( feType == NodalValuedDofs )  newMesh.prolongateBoundaryMask( mask, ProlongationType  );
    if( feType == C1Dofs )           newMesh.prolongateBoundaryMaskDKT( mask, ProlongationType );
}
      
      
      
template< typename ConfiguratorTypePf >
void switchProlongationOfPhaseField( const int ProlongationType,
                                     const typename ConfiguratorTypePf::InitType & newMesh,
                                     const typename ConfiguratorTypePf::MaskType & oldMask,
                                     const typename ConfiguratorTypePf::MaskType & newMask,
                                     typename ConfiguratorTypePf::VectorType &material ){  
    
    const ShellFEType feType = ConfiguratorTypePf::_ShellFEType;
    
    switch( ProlongationType ){
        
        //set all (new???) values to 0
        case 0:{
            if( feType == ElementValuedDofs ){
                const int newSize = newMesh.getNumElements();
                material.resize(newSize);
                for(int i=0; i<newSize; ++i ) material[i] = 0.0;
            }
            if( feType == NodalValuedDofs ){
                const int newSize = newMesh.getNumVertices();
                material.resize(newSize);
                for(int i=0; i<newSize; ++i ) material[i] = 0.0;
            }
            if( feType == C1Dofs ){
                const int newSize = 3 * newMesh.getNumVertices();
                material.resize(newSize);
                for(int i=0; i<newSize; ++i ) material[i] = 0.0;
            }
        }break;
        
        //prolongate linearly
        case 1:{
            if( feType == ElementValuedDofs ) newMesh.prolongateScalarFct_P0( material);
            if( feType == NodalValuedDofs ) newMesh.prolongateScalarFct_P1Linearly( material ); 
            if( feType == C1Dofs ) newMesh.prolongateScalarFct_DKTLinearly( material ); 
        }break;
            
        //prolongate linearly + info at boundary
        case 2:{
                if( feType == NodalValuedDofs ){
                    newMesh.prolongateScalarFct_P1Linearly( material ); 
                    for( int i=oldMask.size(); i<newMask.size(); ++i ){
                        if( (newMask[i]) &&  (material[i] > -0.5) ) material[i] = 1.;    
                    }
                }
                if( feType == C1Dofs ) newMesh.prolongateScalarFct_DKTLinearly( material ); 
        }break;
            
        default : 
            throw std::invalid_argument ( pesopt::strprintf ( "wrong method in file %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
}
 
 
template< typename ConfiguratorType >
void switchProlongationOfDisplacement( const int ProlongationType,
                                       const typename ConfiguratorType::InitType & newMesh,
                                       typename ConfiguratorType::VectorType &disp ){  
    cout << endl << "start prolongation of displacement" << endl;
    auto startTime = std::chrono::high_resolution_clock::now();  
    
    const ShellFEType feType = ConfiguratorType::_ShellFEType;
    switch( ProlongationType ){
            case 0:
                //only to resize
                if( feType == NodalValuedDofs ) newMesh.prolongateVectorFct_P1Linearly( disp ); 
                if( feType == C1Dofs ) newMesh.prolongateVectorFct_DKTLinearly( disp ); 
                disp.setZero();
                break;
            case 1:
                if( feType == NodalValuedDofs ) newMesh.prolongateVectorFct_P1Linearly( disp ); 
                if( feType == C1Dofs ) newMesh.prolongateVectorFct_DKTLinearly( disp ); 
                break;
            default : 
                throw std::invalid_argument ( pesopt::strprintf ( "wrong method in file %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    std::chrono::duration<typename ConfiguratorType::RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << endl << "duration for prolongation of displacement = " << diff.count() << " sec" << endl;
}
            
            
            
template <typename ConfiguratorType>
void FEMProlongationOfDKTDisplacement(
    const ConfiguratorType &oldConf,
    const ConfiguratorType &newConf,
    typename ConfiguratorType::VectorType &disp
  ) {
            
            
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::DomVecType DomVecType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::InitType AdaptiveMeshType;
    typedef typename ConfiguratorType::Point3DType Point3DType;
    typedef typename ConfiguratorType::Matrix32 Matrix32;
            
    //update of positions and tangent-,normal-space in refined mesh
    const int oldSize = oldConf.getInitializer().getNumVertices();
    const int newSize = newConf.getInitializer().getNumVertices();
    typename std::map< int, typename AdaptiveMeshType::ParentInformation >::const_iterator iter;
    
    VectorType oldDisp = disp;
    disp.resize( 3 * 3 * newSize );
    for( int i = 0; i < oldSize; ++i ){
        for(int k=0; k < 3 * 3; ++k )
            disp[i + k * newSize] = oldDisp[i+k*oldSize];
    }
    
    
    //find ElementIndex and Refcoords in oldMesh
    std::vector<int> ElementIndexVec ( newSize - oldSize );
    std::vector<DomVecType> RefCoordVec ( newSize - oldSize );
    for( int nodeIdx = oldSize; nodeIdx < newSize; ++nodeIdx ){
        std::vector<int> nodeIdxVecParents;
        findParentsInOldMesh<AdaptiveMeshType>( oldConf.getInitializer(), newConf.getInitializer(), nodeIdx, nodeIdxVecParents );
        std::vector<int> commonElements;
        oldConf.getInitializer().getCommonElements( nodeIdxVecParents, commonElements );
        //compute barycentric coords
        std::vector<RealType> lambda ( 3 );
        if(nodeIdxVecParents.size() == 2  ){
            computeBarycentricCoordinates<RealType,Point3DType>( newConf.getInitializer().getVertex(nodeIdx), newConf.getInitializer().getVertex(nodeIdxVecParents[0]), newConf.getInitializer().getVertex(nodeIdxVecParents[1]), lambda[0], lambda[1] );
            lambda[2] = 0.;
        }
        if(nodeIdxVecParents.size() == 3 ){
            computeBarycentricCoordinates<RealType,Point3DType>( newConf.getInitializer().getVertex(nodeIdx), newConf.getInitializer().getVertex(nodeIdxVecParents[0]), newConf.getInitializer().getVertex(nodeIdxVecParents[1]), newConf.getInitializer().getVertex(nodeIdxVecParents[2]), lambda[0], lambda[1], lambda[2] );
        }
        //compute refcoord in commonElements[0]
        DomVecType RefCoordNewPoint; RefCoordNewPoint.setZero();
        for(int i=0; i<nodeIdxVecParents.size(); ++i ){
           int localIdxParent = oldConf.getInitializer().getLocalIdxWrtGlobalElementIdx( nodeIdxVecParents[i], commonElements[0] );
           DomVecType RefCoordParent; RefCoordParent.setZero();
           if( localIdxParent == 1 ) RefCoordParent[0] = 1.;
           if( localIdxParent == 2 ) RefCoordParent[1] = 1.;
           RefCoordNewPoint += lambda[i] * RefCoordParent;
        }
        ElementIndexVec[nodeIdx - oldSize] = commonElements[0];
        RefCoordVec[nodeIdx - oldSize] = RefCoordNewPoint;
    }
        
    DKTFEVectorFunctionEvaluator<ConfiguratorType> dfd ( oldConf, oldDisp, 3 );
    for( int nodeIdx = oldSize; nodeIdx < newSize; ++nodeIdx ){
        // find new node index in map
        iter = newConf.getInitializer().getInterpolationMap().find( nodeIdx );
        if( iter == newConf.getInitializer().getInterpolationMap().end() ) 
            throw std::invalid_argument ( pesopt::strprintf ( "unknown vertex in file %s at line %d.", __FILE__, __LINE__ ).c_str() );
        //use exact Position in Old Mesh:
        int ElementIndexOldMesh = ElementIndexVec[nodeIdx - oldSize];
        DomVecType RefCoordsOldMesh; RefCoordsOldMesh = RefCoordVec[nodeIdx - oldSize];
       
        Point3DType newVertex; 
        dfd.evaluate( oldConf.getInitializer().getTriang(ElementIndexOldMesh), RefCoordsOldMesh, newVertex );
        for ( int comp = 0; comp<3; ++comp )
            disp(nodeIdx + comp*newConf.getNumGlobalDofs() )  = newVertex(comp);
        
        Matrix32 Dx;
        dfd.evaluateApproxGradient( oldConf.getInitializer().getTriang(ElementIndexOldMesh), RefCoordsOldMesh, Dx );
        for ( int comp = 0; comp<3; ++comp ) {
            disp(nodeIdx + 1*newConf.getInitializer().getNumVertices() + comp*newConf.getNumGlobalDofs() ) = Dx(comp, 0);
            disp(nodeIdx + 2*newConf.getInitializer().getNumVertices() + comp*newConf.getNumGlobalDofs() ) = Dx(comp, 1);
        }
    }//end for old size
            
            
}
            
            
            
            
            
            

#endif
