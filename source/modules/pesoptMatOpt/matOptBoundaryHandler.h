#ifndef __MATOPTBOUNDARYHANDLER_H
#define __MATOPTBOUNDARYHANDLER_H

#include <pesopt_IO.h>
#include <feBoundaryHandler.h>


template <typename ConfiguratorType>
class MaterialOptimizationBoundaryConditionHandler {
 
 typedef typename ConfiguratorType::DTContainer DataTypeContainer;
 typedef pesopt::BoostParser ParameterParserType;

 FEBoundaryConditionHandler<ConfiguratorType> *bdryHandlerMaterial;
 FEBoundaryConditionHandler<ConfiguratorType> *bdryHandlerDisplacement;

public:

 MaterialOptimizationBoundaryConditionHandler ( const ParameterParserType & parser, 
                                                const ConfiguratorType &conf 
   ) {
    
    //boundary condition for material
    const string boundaryConditionMaterial = parser.template get<string>("Boundary.boundaryConditionMaterial");
    if( boundaryConditionMaterial == "Zero" ){
       bdryHandlerMaterial = new FEZeroBoundaryConditionHandler<ConfiguratorType> ( conf );
    }else{
        throw std::invalid_argument( "boundaryConditonMaterial not known" );
    }
    
    //boundary condition for displacement
    const string boundaryConditionDisplacement = parser.template get<string>("Boundary.boundaryConditionDisplacement");
    // ! GENERAL
    if( boundaryConditionDisplacement ==  "Zero" ){
       bdryHandlerDisplacement = new FEZeroBoundaryConditionHandler<ConfiguratorType> ( conf );
    }else if ( boundaryConditionDisplacement ==  "Dirichlet" ){
       bdryHandlerDisplacement = new FEDirichletBoundaryConditionHandler<ConfiguratorType> ( parser, conf );   
    // ! UNIT CUBE
    }else if ( boundaryConditionDisplacement ==  "LeftSide" ){
       bdryHandlerDisplacement = new FEDirichletBoundaryConditionHandlerLeftSide<ConfiguratorType> ( conf );
    }else if ( boundaryConditionDisplacement ==  "DownSide" ){
       bdryHandlerDisplacement = new FEDirichletBoundaryConditionHandlerDownSide<ConfiguratorType> ( parser, conf );
    }else if ( boundaryConditionDisplacement ==  "DirichletRectangleStruts" ){
       bdryHandlerDisplacement = new FEDirichletBoundaryConditionHandlerRectangleStruts<ConfiguratorType> ( parser, conf );  
    // ! RIGHT TRIANGLE 30-60
    }else if ( boundaryConditionDisplacement ==  "DirichletRightTriangle3060Struts" ){
       bdryHandlerDisplacement = new FEDirichletBoundaryConditionHandlerRightTriangle3060Struts<ConfiguratorType> ( parser, conf );   
    // ! EQUILATERAL TRIANGLE (1.0), (0, sqrt(3)),  (0, -1)
    }else if ( boundaryConditionDisplacement ==  "DirichletEquilateralTriangleStruts" ){
       bdryHandlerDisplacement = new FEDirichletBoundaryConditionHandlerEquilateralTriangleStruts<ConfiguratorType> ( parser, conf );    
    }else{
        throw std::invalid_argument( "boundaryConditonDisplacement not known" );
    }
 }
 
 ~MaterialOptimizationBoundaryConditionHandler () {
   delete bdryHandlerMaterial;
   delete bdryHandlerDisplacement;
 }

    
 const FEBoundaryConditionHandler<ConfiguratorType> &getBoundaryHandlerMaterial ( ) { return *bdryHandlerMaterial;}
 const FEBoundaryConditionHandler<ConfiguratorType> &getBoundaryHandlerDisplacement ( ) { return *bdryHandlerDisplacement;}

};


#endif
