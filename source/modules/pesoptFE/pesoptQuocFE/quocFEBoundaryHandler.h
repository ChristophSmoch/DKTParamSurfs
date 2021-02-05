#ifndef __QUOCFEBOUNDARYHANDLER_H
#define __QUOCFEBOUNDARYHANDLER_H

#include <feBoundaryHandler.h>


//TODO depending on dimChartDomain
// template< typename ConfiguratorType >
// void generateDirichletBoundaryMaskUponBoundaryType ( const QuocBoundaryType boundaryType, 
//                                                      const typename ConfiguratorType::InitType &mesh,
//                                                      typename ConfiguratorType::MaskType & mask ) {
//   switch ( boundaryType ){
//       
//       case NOBOUNDARY : {
//       } break;
//       
//       case LEFT : {
//             mask = mesh._boundaryLeft;
//       } break;
//       case RIGHT : {
//             mask = mesh._boundaryRight;
//       }break;
// //       case TOP : {
// //             mask = mesh._boundaryTop;
// //       }break;
// //       case BOTTOM : {
// //             mask = mesh._boundaryBottom;
// //       }break;
// //       case FRONT : {
// //             mask = mesh._boundaryFront;
// //       }break;
// //       case BACK : {
// //             mask = mesh._boundaryBack;
// //       }break;
//       case ALL : {
//            mask = mesh.getBoundaryMask();
//       }break;
//   
//       default :
//         throw std::invalid_argument( pesopt::strprintf ( "Wrong boundary condition. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
//         break;
//     }
// }


// template< typename ConfiguratorType >
// class QuocDirichletBoundaryConditionHandler
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
//   QuocDirichletBoundaryConditionHandler( const ParameterParserType & parser, const ConfiguratorType &conf ) : 
//     FEBoundaryConditionHandler<ConfiguratorType> ( conf ),
//     _parser ( parser )
//   {
//       generateDirichletBoundaryMask( this->_DirichletMask );
//   }
//   
//   bool hasDirichletBoundary() const override { return true; }
//   bool hasPeriodicBoundary() const override { return false; }
//   
//   void generateDirichletBoundaryMask ( MaskType & mask ) const{
//     mask.resize( this->_numGlobalDofs, false );
//     generateDirichletBoundaryMaskUponBoundaryType<ConfiguratorType> ( static_cast<QuocBoundaryType>( _parser.template get<int>( "InputMesh.DirichletBoundaryType" ) ), this->_mesh, mask );
//   }
// };



template< typename ConfiguratorType >
class QuocPeriodicBoundaryConditionHandler : 
public FEBoundaryConditionHandler<ConfiguratorType>{
  
protected:
  
  typedef typename ConfiguratorType::RealType               RealType;
  typedef typename ConfiguratorType::InitType               MeshType;
  typedef typename ConfiguratorType::MaskType               MaskType;
  typedef typename ConfiguratorType::PointType              PointType;
  typedef typename ConfiguratorType::VectorType             VectorType;
  typedef typename ConfiguratorType::DTContainer            DataTypeContainer;
  typedef pesopt::BoostParser ParameterParserType;
  
  const ParameterParserType &_parser;
  
public:
  
  QuocPeriodicBoundaryConditionHandler( const ParameterParserType & parser, const ConfiguratorType &conf ) : 
    FEBoundaryConditionHandler<ConfiguratorType> ( conf ),
    _parser ( parser )
   {
      this->_PeriodicMask = this->_mesh._boundaryPeriodic;
      this->_PeriodicIndices = this->_mesh._periodicIdentificationIndices;
   }
  
  bool hasDirichletBoundary() const override { return true; }
  bool hasPeriodicBoundary() const override { return false; }
  
};



template< typename ConfiguratorType >
class QuocDirichletAndPeriodicBoundaryConditionHandler
: public FEBoundaryConditionHandler<ConfiguratorType>{
  
protected:
  
  typedef typename ConfiguratorType::RealType               RealType;
  typedef typename ConfiguratorType::InitType               MeshType;
  typedef typename ConfiguratorType::MaskType               MaskType;
  typedef typename ConfiguratorType::PointType              PointType;
  typedef typename ConfiguratorType::VectorType             VectorType;
  typedef typename ConfiguratorType::DTContainer            DataTypeContainer;  
  typedef pesopt::BoostParser ParameterParserType;
  
  const ParameterParserType &_parser;
  
public:
  
  QuocDirichletAndPeriodicBoundaryConditionHandler( const ParameterParserType & parser, const ConfiguratorType &conf ) : 
    FEBoundaryConditionHandler<ConfiguratorType> ( conf ),
    _parser ( parser )
   { 
      this->_DirichletMask = this->_mesh.getBoundaryMask();
      this->_PeriodicMask = this->_mesh._boundaryPeriodic;
      this->_PeriodicIndices = this->_mesh._periodicIdentificationIndices;
   }
  
  bool hasDirichletBoundary() const override { return true; }
  bool hasPeriodicBoundary() const override { return true; }
  
};










// template< typename ConfiguratorType >
// class QuocDirichletAndMirrorPeriodicBoundaryConditionHandler
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
//   typedef typename DataTypeContainer::IntVecChart            IntVecChart;
//   
//   const ParameterParserType &_parser;
//   const ConfiguratorType &_conf;
//   const MeshType &_mesh;
//   const int _numVertices, _numGlobalDofs, _numElements;
//   mutable MaskType _DirichletMask;
//   mutable MaskType _PeriodicMask;
//   mutable std::vector<int> _PeriodicIndices;
//   
// public:
//   
//   QuocDirichletAndMirrorPeriodicBoundaryConditionHandler( const ParameterParserType & parser, const ConfiguratorType &conf ) : 
//   _parser ( parser ),
//   _conf ( conf ),
//   _mesh( conf.getMesh() ), _numElements ( _mesh.getNumElements() ),
//   _numVertices( conf.getNumGlobalDofs() ), _numGlobalDofs ( conf.getNumGlobalDofs() ) {
//       
//           generateDirichletBoundaryMask( _DirichletMask );
//           
//           _PeriodicMask.resize( _numGlobalDofs, false ); 
//           _PeriodicIndices.resize( _numGlobalDofs, -1 );
//           
//           _PeriodicMask = _mesh._boundaryPeriodic;
//           _PeriodicIndices = _mesh._periodicIdentificationIndices;
//           
//            //TEST for 7x7
//            IntVecChart nodeIndices; nodeIndices[0] = 4; nodeIndices[1] = 3;
//            IntVecChart mirrorNodeIndices; mirrorNodeIndices[0] = 2; mirrorNodeIndices[1] = 3;
//            _PeriodicMask[_mesh.getGlobalNodeIndex(nodeIndices)] = true;
//            _PeriodicIndices[_mesh.getGlobalNodeIndex(nodeIndices)] = _mesh.getGlobalNodeIndex(mirrorNodeIndices);
//           
// //          for( int nodeIdx=0; nodeIdx<_mesh.getNumVertices(); ++nodeIdx ){
// //              if( _mesh._boundaryPeriodic[nodeIdx] ){
// //                  //TODO
// //              }else{
// //                 IntVecChart nodeIndices = _mesh.getNodeIndices( nodeIdx );
// //                 IntVecChart mirrorNodeIndices = nodeIndices;
// //                 if( nodeIndices[0] > _mesh.getNumDofs(0) / 2 ){
// //                     _PeriodicMask[nodeIdx] = true;
// //                     mirrorNodeIndices[0] = _mesh.getNumDofs(0) - 1 - nodeIndices[0];
// //                     cout << "xIdx = " << nodeIndices[0] << ", mirrorXIdx = " << mirrorNodeIndices[0] << endl;
// //                 }
// //                 _PeriodicIndices[nodeIdx] = _mesh.getGlobalNodeIndex(mirrorNodeIndices);
// //              }
// //           }
// //           
// //           //TODO for 5x5
// //           _PeriodicIndices[23] = 1;
// //           
// //           for(int i=0; i<_PeriodicIndices.size(); ++i) 
// //               cout << "i=" << i << ", index = " << _PeriodicIndices[i] << ", mask = " << _PeriodicMask[i] << endl;
//       
//   }
//   
//   bool hasDirichletBoundary() const override { return true; }
//   bool hasPeriodicBoundary() const override { return true; }
//   
//   const MaskType &getDirichletMask ( ) const override{ return _DirichletMask;}
//   const MaskType &getPeriodicMask ( ) const override{ return _PeriodicMask;}
//   const std::vector<int> &getPeriodicIndices ( ) const override{ return _PeriodicIndices;}
//   
// };


#endif //__QUOCBOUNDARYHANDLER_H
