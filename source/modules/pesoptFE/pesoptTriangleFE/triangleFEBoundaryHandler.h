#ifndef __TRIANGLEFEBOUNDARYHANDLER_H
#define __TRIANGLEFEBOUNDARYHANDLER_H

#include <feBoundaryHandler.h>


// template< typename ConfiguratorType >
// class TriangleFEDirichletBoundaryConditionHandler
// : public FEBoundaryConditionHandler<ConfiguratorType>{
//   
// protected:
//   
//   typedef typename ConfiguratorType::RealType               RealType;
//   typedef typename ConfiguratorType::InitType               MeshType;
//   typedef typename ConfiguratorType::MaskType               MaskType;
//   typedef typename ConfiguratorType::PointType              PointType;
//   typedef typename ConfiguratorType::VectorType             VectorType;
//   typedef typename ConfiguratorType::DTContainer            DataTypeContainer;
//   typedef pesopt::BoostParser ParameterParserType;
//   
//   const ParameterParserType &_parser;
//   
// public:
//   
//   TriangleFEDirichletBoundaryConditionHandler( const ParameterParserType & parser, const ConfiguratorType &conf ) :
//     FEBoundaryConditionHandler<ConfiguratorType> ( conf ),
//     _parser ( parser )
//   {
//       generateDirichletBoundaryMask( this->_DirichletMask );
//   }
//   
//   bool hasDirichletBoundary() const override { return true; }
//   bool hasPeriodicBoundary() const override { return false; }
//   
//   void generateDirichletBoundaryMask ( MaskType &mask ) const{
//     mask.resize( this->_numGlobalDofs, false );
//     //TODO here for example left boundary
//     for ( int nodeIdx=0; nodeIdx < this->_numVertices; ++nodeIdx ) {
//         const PointType& coords ( this->_mesh.getVertex(nodeIdx) );
//         if ( coords [0] == 0. ){
//             mask[nodeIdx] = true;
//         }
//     }
//   }
//   
// };



#endif 
