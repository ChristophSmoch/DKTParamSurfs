#ifndef __QUOCFEADAPTIVEREFINEMENT_H
#define __QUOCFEADAPTIVEREFINEMENT_H


#include <feAdaptiveRefinement.h>
#include <feFunctionEvaluator.h>
#include "quocFEConfigurators.h"



template<typename QuocMeshType >
class FEAdaptiveMesh<QuocMeshType, QUOCMESH> 
: public FEAdaptiveMeshBase<QuocMeshType>{
    
  typedef typename QuocMeshType::ElementType              ElementType;
  typedef typename ElementType::ElementNodeIndicesType    ElementNodeIndicesType;
  typedef typename QuocMeshType::RealType                 RealType;
  typedef typename QuocMeshType::RealVecChart             RealVecChart;
  typedef typename QuocMeshType::IntVecChart              IntVecChart;
  typedef typename QuocMeshType::PointType                PointType;
  typedef typename QuocMeshType::VectorType               VectorType;
  typedef typename QuocMeshType::MaskType                 MaskType;
    
public:
    
    
  FEAdaptiveMesh( const QuocMeshType &mesh ) :
    FEAdaptiveMeshBase<QuocMeshType> ( mesh ) { }
    
  void refineAll( ) const override{
      refineByFactor( 2. );
  }
  
//   void refineMarkedElements( ) {
//      throw std::logic_error ( pesopt::strprintf ( "AdaptiveQuocMesh::refineMarked is not defined in %s at line %d.", __FILE__, __LINE__ ).c_str() );
//   }
  
  void refineByFactor( const RealType factor ) const override {
    this->_refinementLevel++;
     PointType lengthVec; this->_initMesh.getLenghtVec( lengthVec );
     IntVecChart numDofVecStartLevel; this->_initMesh.getNumDofVec( numDofVecStartLevel );
     IntVecChart numDofVecRefined; 
     for( int i=0; i<numDofVecRefined.size(); ++i ){
         numDofVecRefined[i] =  std::round( std::pow(factor, this->_refinementLevel ) * ( numDofVecStartLevel[i] - 1.) + 1. );
      }
      QuocMeshType refinedMesh ( numDofVecRefined, lengthVec );
      this->_meshVector.push_back( refinedMesh );
  }
    
    
  template <typename ConfiguratorType>
  void prolongateScalarFunctionAtNodes ( VectorType &vec ) const {
     VectorType vecOld(vec);
     this->template prolongateScalarFunctionAtNodes<ConfiguratorType> ( vecOld, vec );
 }
 
 
  template <typename ConfiguratorType>
  void prolongateScalarFunctionAtNodes ( const VectorType &vecOld, VectorType &vecNew  ) const {
     ConfiguratorType confOld ( this->_meshVector[this->_refinementLevel-1] );
     FEScalarFunctionEvaluator<ConfiguratorType> discreteFctOld ( confOld, vecOld ); 
     vecNew.resize ( this->_meshVector[this->_refinementLevel].getNumVertices() );
     for( int nodeIdxCurrentLevel=0; nodeIdxCurrentLevel <  this->_meshVector[this->_refinementLevel].getNumVertices(); nodeIdxCurrentLevel++ ){
         const PointType& GlobalCoords =  this->_meshVector[this->_refinementLevel].getVertex ( nodeIdxCurrentLevel );
         int elementNumberOld; PointType LocalCoordOld;
         confOld.getLocalCoords ( GlobalCoords, elementNumberOld, LocalCoordOld );
         vecNew[nodeIdxCurrentLevel] = discreteFctOld.evaluate(  this->_meshVector[this->_refinementLevel-1].getElement(elementNumberOld), LocalCoordOld );
     }
 }
 
 
 
 
  template <typename ConfiguratorType>
  void prolongateScalarFunctionOnElements ( VectorType &vec ) const {
     VectorType vecOld(vec);
     this->template prolongateScalarFunctionOnElements<ConfiguratorType> ( vecOld, vec );
 }
 
 
  template <typename ConfiguratorType>
  void prolongateScalarFunctionOnElements ( const VectorType &vecOld, VectorType &vecNew  ) const {
     ConfiguratorType confOld ( this->_meshVector[this->_refinementLevel-1] );
     vecNew.resize ( this->_meshVector[this->_refinementLevel].getNumElements() );
     for( int elementIdxCurrentLevel=0; elementIdxCurrentLevel <  this->_meshVector[this->_refinementLevel].getNumElements(); elementIdxCurrentLevel++ ){
         
         PointType LocalCoordMidPoint;
         for (int i = 0; i < LocalCoordMidPoint.size(); ++i ) LocalCoordMidPoint(i) = 0.5;
         PointType GlobalCoordsNew;
         this->_meshVector[this->_refinementLevel].getElement(elementIdxCurrentLevel).getGlobalCoord ( LocalCoordMidPoint, GlobalCoordsNew );
         
         int elementNumberOld; PointType LocalCoordOld;
         confOld.getLocalCoords ( GlobalCoordsNew, elementNumberOld, LocalCoordOld );

         vecNew[elementIdxCurrentLevel] = vecOld[elementNumberOld];
     }
 }
 
 
 
 
  template <typename ConfiguratorType>
  void prolongateVectorFunctionOnElements ( VectorType &vec ) const {
     VectorType vecOld(vec);
     this->template prolongateScalarFunctionOnElements<ConfiguratorType> ( vecOld, vec );
 }
 
 
  template <typename ConfiguratorType>
  void prolongateVectorFunctionOnElements ( const VectorType &vecOld, VectorType &vecNew, 
                                            const int numComponents = ConfiguratorType::dimChartDomain ) const {

    const int numElementsOld = this->_meshVector[this->_refinementLevel - 1].getNumElements();
    std::vector<Eigen::Ref<const VectorType> > refsOld;
    refsOld.reserve( numComponents );
    for( int c=0; c < numComponents; ++c )
        refsOld.push_back ( vecOld.segment( c * numElementsOld, numElementsOld) );
  
    const int numElementsNew = this->_meshVector[this->_refinementLevel].getNumElements();
    vecNew.resize ( numComponents * numElementsNew );
    // TODO see https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
    // Version with Eigen::Ref
//     std::vector<Eigen::Ref<VectorType> > refsNew;
//     refsNew.reserve( numComponents );
//     for( int c=0; c < numComponents; ++c )
//         refsNew.push_back ( vecNew.segment( c * numElementsNew, numElementsNew) );
//     for( int c=0; c < numComponents; ++c )
//         prolongateScalarFunctionOnElements<ConfiguratorType>(refsOld[c], refsNew[c] );
    // Version with hard copy
     std::vector<VectorType> refsNew;
     refsNew.reserve( numComponents );
     for( int c=0; c < numComponents; ++c )
        refsNew.push_back ( VectorType(numElementsNew) );
     for( int c=0; c < numComponents; ++c )
        prolongateScalarFunctionOnElements<ConfiguratorType>(refsOld[c], refsNew[c] );
     for (int i = 0; i<numElementsNew; ++i )
         for( int c=0; c < numComponents; ++c )
             vecNew(i + c*numElementsNew) = refsNew[c](i);
  }
    
};


#endif
