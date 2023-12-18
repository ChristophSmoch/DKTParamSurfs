#ifndef __QUOCFEBOUNDARYHANDLER_H
#define __QUOCFEBOUNDARYHANDLER_H

#include <feBoundaryHandler.h>





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













#endif //__QUOCBOUNDARYHANDLER_H
