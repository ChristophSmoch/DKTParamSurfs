#ifndef __FEMATERIALGENERATOR_H
#define __FEMATERIALGENERATOR_H


# include <pesopt_VTK.h>

template< typename ConfiguratorType >
class FEMaterialGenerator{
  
public:
  
  typedef typename ConfiguratorType::RealType               RealType;
  typedef typename ConfiguratorType::InitType               MeshType;
  typedef typename ConfiguratorType::MaskType               MaskType;
  typedef typename ConfiguratorType::PointType              PointType;
  typedef typename ConfiguratorType::VectorType             VectorType;
  typedef typename ConfiguratorType::DTContainer            DataTypeContainer;
  typedef pesopt::BoostParser ParameterParserType;
  
  const ParameterParserType &_parser;
  const ConfiguratorType &_conf;
  const MeshType &_mesh;
  const int _numVertices, _numGlobalDofs, _numElements;
  const RealType _pfLowerBound, _pfUpperBound;
  
public:
  
  FEMaterialGenerator( const ParameterParserType & Parser, const ConfiguratorType &conf ) : 
  _parser ( Parser),
  _conf ( conf ),
  _mesh( conf.getMesh() ), _numElements ( _mesh.getNumElements() ),
  _numVertices( conf.getNumGlobalDofs() ), _numGlobalDofs ( conf.getNumGlobalDofs() ),
  _pfLowerBound ( _parser.template get<RealType>("PhaseField.BoxConstraint_l") ),
  _pfUpperBound ( _parser.template get<RealType>("PhaseField.BoxConstraint_u" ) )
  {  }
  
 //==========================================================================================================================
 //==================================   Material     ========================================================================
 //==========================================================================================================================
  void constructConstantMaterial( VectorType &material, RealType materialConstant ) const{
    for( int i=0; i<_numVertices; ++i ) material[i] = materialConstant;
  }
  
  void constructLayeredMaterial( VectorType &material, 
                                 const RealType startLayer, const RealType endLayer, 
                                 const int direction = 0
                               ) const{
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const PointType& coords ( _mesh.getVertex(nodeIdx) );
        if( (coords[direction] <= endLayer) && (coords[direction] >= startLayer) ) material[nodeIdx] = _pfUpperBound;
        else material[nodeIdx] = _pfLowerBound;
    }      
  }
  
//   void constructHoles( VectorType &material, VTKDataSupp supp  ) const{ 
//       
//     const RealType pi = 4 * atan ( 1.0 );  
//       
//     RealType lx = 0.0, ly = 0.0;
//     for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
//                 const PointType& coords ( _mesh.getVertex(nodeIdx) );
//                 if( coords[0] > lx ) lx = coords[0];
//                 if( coords[1] > ly ) ly = coords[1];
//     }
//     
//     switch (supp){
//         
//         case VERTEX_DATA:            
//             for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
//                 const PointType& coords ( _mesh.getVertex(nodeIdx) );
//                 RealType tmp = cos( 6. * pi * coords[0] / lx ) * cos( 4. * pi * coords[1] ) + 0.6;
//                 tmp -= std::max<RealType>( 200.*(0.01 - coords[0] * coords[0] - (coords[1] - 0.5*ly) * (coords[1] - 0.5*ly) ) , 0. );
//                 tmp -= std::max<RealType>( 100. * ( coords[0] + coords[1] - lx - ly + 0.1), 0. );
//                 tmp -= std::max<RealType>( 100. * ( coords[0] - coords[1] - lx + 0.1 ), 0. );
//             
//                 if( tmp > 1.0 ) tmp = _pfUpperBound;
//                 if( tmp < -1.0 ) tmp = _pfLowerBound;
//                 material[nodeIdx] = tmp;
//             }
//         break;
//             
//         case FACE_DATA:
//             for ( int elementIdx=0; elementIdx < _numElements; ++elementIdx ) {
//                 PointType midpoint; midpoint.setZero();                
//                 for( int localNodeIdx=0; localNodeIdx < 4; ++localNodeIdx ){
//                     int globalNodeIdx = _mesh.getElementNodeIdx ( elementIdx, localNodeIdx );
//                     PointType coord = _mesh.getVertex( globalNodeIdx );
//                     midpoint += coord;
//                 }
//                 midpoint *= 0.25;
//                 
//                 RealType tmp = cos( 6. * pi * midpoint[0] / lx ) * cos( 4. * pi * midpoint[1] ) + 0.6;
//                 tmp -= std::max<RealType>( 200.*(0.01 - midpoint[0] * midpoint[0] - (midpoint[1] - 0.5*ly) * (midpoint[1] - 0.5*ly) ) , 0. );
//                 tmp -= std::max<RealType>( 100. * ( midpoint[0] + midpoint[1] - lx - ly + 0.1), 0. );
//                 tmp -= std::max<RealType>( 100. * ( midpoint[0] - midpoint[1] - lx + 0.1 ), 0. );
//             
//                 if( tmp > 1.0 ) tmp = 1.0;
//                 if( tmp < -1.0 ) tmp = -1.0;
//                 material[elementIdx] = tmp;
//             }            
//         break;
//     }
//   }
  
  
  void constructSchwarzPSurface( VectorType &material ) const{ 
    RealType lx = 0.0, ly = 0.0, lz = 0.0;
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const PointType& coords ( _mesh.getVertex(nodeIdx) );
        if( coords[0] > lx ) lx = coords[0];
        if( coords[1] > ly ) ly = coords[1];
        if( coords[2] > lz ) lz = coords[2];
    }
    const RealType pi = 4 * atan ( 1.0 );
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const PointType& coords ( _mesh.getVertex(nodeIdx) );
        RealType tmp = cos( 2. * pi * coords[0] / lx ) + cos( 2. * pi *  coords[1] / ly ) + cos( 2. * pi *  coords[2] / lz );
        if( tmp > 0.0 ) material[nodeIdx] = _pfUpperBound;
        if( tmp < 0.0 ) material[nodeIdx] = _pfLowerBound;
        if( tmp == 0.0 ) material[nodeIdx] = 0.5 * (_pfLowerBound + _pfUpperBound );
    }
  }
  
  void constructGyroid( VectorType &material ) const{ 
      RealType lx = 0.0, ly = 0.0, lz = 0.0;
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const PointType& coords ( _mesh.getVertex(nodeIdx) );
        if( coords[0] > lx ) lx = coords[0];
        if( coords[1] > ly ) ly = coords[1];
        if( coords[2] > lz ) lz = coords[2];
    }
    const RealType pi = 4 * atan ( 1.0 );
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const PointType& coords ( _mesh.getVertex(nodeIdx) );
        RealType tmp = sin( 2. * pi * coords[0] / lx ) * cos( 2. * pi * coords[1] / ly ) 
                     + sin( 2. * pi * coords[1] / ly ) * cos( 2. * pi * coords[2] / lz ) 
                     + sin( 2. * pi * coords[2] / lz ) * cos( 2. * pi * coords[0] / lx );
        if( tmp > 0.0 ) material[nodeIdx] = _pfUpperBound;
        if( tmp < 0.0 ) material[nodeIdx] = _pfLowerBound;
        if( tmp == 0.0 ) material[nodeIdx] = 0.5 * (_pfLowerBound + _pfUpperBound );
    }
  }
  
  void constructRandomMaterial( VectorType &material ) const{
      //this gives random values in [-1,1]
      VectorType tmp ( material.size() );
      tmp = VectorType::Random( material.size() ); 
      //scale to [_pfLowerBound,_pfUpperBound]
      pesopt::scaleVector<RealType,VectorType> ( tmp, -1, 1, material, _pfLowerBound, _pfUpperBound );
  }
  
//TODO 
//   void switchMaterialType( VectorType &material, VTKDataSupp dataSupp = VERTEX_DATA ) const{
 void switchMaterialType( VectorType &material ) const{
      switch( _parser.template get<int>( "Material.initMaterialType" ) ){
          
         case -2:{
             cout  << endl << "load vector from input mesh file " << endl;
             VTKSaver vtkReader;
             std::vector<RealType> dataVec;
//              vtkReader.template getScalarPointDataVec<RealType> ( pesopt::strprintf ( "%s", _parser.template get<string> ( "Material.materialFile" ).c_str(), "material", dataVec );
             vtkReader.template getScalarPointDataVec<RealType> ( 
                 _parser.template get<string> ( "InputMesh.fileName" ), "material", dataVec );
             for ( int i = 0; i<dataVec.size(); i++ ) material[i] = dataVec[i];
             cout  << endl << "finished load vector from input mesh file " << endl;      
         } break;
//         case -1:
//           pesopt::loadVectorFromFile<VectorType>( material, pesopt::strprintf ( "%s", _parser.template get<string> ( "Material.materialFile" ).c_str () ) ); 
//           material *= -1.;
//           break;
        case 0:
          cout  << endl << "load vector from file " << _parser.template get<string> ( "Material.materialFile" ).c_str() << endl;
          pesopt::loadVectorFromFile<VectorType>( material, pesopt::strprintf ( "%s", _parser.template get<string> ( "Material.materialFile" ).c_str() ) ); 
          break;
        case 1:
          constructConstantMaterial( material, _parser.template get<double>( "Material.materialConstant" ) );
          break;
        case 2:
          constructLayeredMaterial( material, _parser.template get<double>( "Material.materialStartLayer" ), _parser.template get<double>( "Material.materialEndLayer" ), _parser.template get<double>( "Material.materialDirectionLayer" ) );
          break;
//         case 10:
//             switch(dataSupp){
//                 case VERTEX_DATA:
//                     constructHoles( material, VERTEX_DATA );
//                     break;
//                 case FACE_DATA:
//                     constructHoles( material, FACE_DATA );
//                     break;
//             }
//             break;
        case 20:
            constructSchwarzPSurface( material );
            break;
//         case 21:
//             constructSchwarzPSurface( material );
//             material *= -1.;
//             break;
        case 22:
            constructGyroid( material );
            break;
//         case 23:
//             constructGyroid( material );
//             material *= -1.;
//             break;
        case 1000 :
            constructRandomMaterial( material );
            break;
        default:
          throw std::invalid_argument( pesopt::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      }
  }
  
 //==========================================================================================================================
 //==================================   Boundary     ========================================================================
 //==========================================================================================================================
 // TODO difference: periodic or not
 void switchMaterialBoundaryType( VectorType& material_l, VectorType& material_u, 
                                  VectorType& initMaterial, 
                                  const bool hasPeriodicBoundary ){ 
        
     const int dimChartDomain = _conf.dimChartDomain;
       
      switch( _parser.template get<int>("Boundary.boundaryConditionMaterialType") ){
          
       //=====================================================
       //        Unit cube [0,1]^d
       // ====================================================
          
          //simple strut (in 2d: segment, in 3d: square)
          case 1: {
              const RealType strutSize = _parser.template get<RealType>("Boundary.strutSizeBoundaryConditionMaterial");  
              for( int nodeIdx=0; nodeIdx<initMaterial.size(); ++nodeIdx ){
                  const PointType& coord = _mesh.getVertex( nodeIdx );
                  for( int coordDir=0; coordDir < dimChartDomain; ++coordDir ){
                     
                     bool setBoundaryConstraints = false;
                     if ( hasPeriodicBoundary ) {
                         if( coord[coordDir] == 0 ) setBoundaryConstraints = true;
                     }else{
                         if( (coord[coordDir] == 0 ) ||  ( coord[coordDir] == 1 ) ) setBoundaryConstraints = true; 
                     }
                     
                     if( setBoundaryConstraints ){
                         bool use_u = true;
                         for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                           int dir = (coordDir + otherDir) % dimChartDomain;
                           if( std::abs( coord[dir] - 0.5  ) >= strutSize ) use_u = false;
                         }
                         if( use_u ){  material_l[nodeIdx] = _pfUpperBound; initMaterial[nodeIdx] = _pfUpperBound;
                         }else{        material_u[nodeIdx] = _pfLowerBound; initMaterial[nodeIdx] = _pfLowerBound; }
                     }
                    
                  }
                }
            } break;
            
            
          //strut with distance (in 2d: segment x segment, in 3d: square x segment)
          case 2: {
              const RealType strutSize = _parser.template get<RealType>("Boundary.strutSizeBoundaryConditionMaterial");
              const RealType thicknessStrut = _parser.template get<RealType>("MaterialOptimization.thicknessStrut");
              for( int nodeIdx=0; nodeIdx<initMaterial.size(); ++nodeIdx ){
                  const PointType& coord = _mesh.getVertex( nodeIdx );
                  for( int coordDir=0; coordDir < dimChartDomain; ++coordDir ){
                    
                     bool setBoundaryConstraints = false;
                     if ( hasPeriodicBoundary ) {
                         if( coord[coordDir] <= thicknessStrut  ) setBoundaryConstraints = true;
                     }else{
                         if( ( coord[coordDir] <= thicknessStrut  ) ||  ( coord[coordDir] >= 1. - thicknessStrut  ) ) setBoundaryConstraints = true; 
                     }
                    
                    if( setBoundaryConstraints ){
                         bool use_u = true;
                         for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                           int dir = (coordDir + otherDir) % dimChartDomain;
                           if( std::abs( coord[dir] - 0.5  ) >= strutSize ) use_u = false;
                         }
                         if( use_u ){ material_l[nodeIdx] = _pfUpperBound;  initMaterial[nodeIdx] = _pfUpperBound;
                         }else{       material_u[nodeIdx] = _pfLowerBound;  initMaterial[nodeIdx] = _pfLowerBound; }
                    }
                  }
                }
            } break;
            
            
          //strut with distance depending on mesh size (in 2d: segment x segment, in 3d: square x segment)
          case 4: {
              const RealType strutSize = _parser.template get<RealType>("Boundary.strutSizeBoundaryConditionMaterial");
              const RealType thicknessStrutFactor = _parser.template get<RealType>("Boundary.thicknessStrutFactor");
              for( int nodeIdx=0; nodeIdx<initMaterial.size(); ++nodeIdx ){
                  const PointType& coord = _mesh.getVertex( nodeIdx );
                  for( int coordDir=0; coordDir < dimChartDomain; ++coordDir ){
                      
                    bool setBoundaryConstraints = false;
                     if ( hasPeriodicBoundary ) {
                         if( coord[coordDir] <= thicknessStrutFactor * _mesh.getMeshSize( coordDir )  ) setBoundaryConstraints = true;
                     }else{
                         if( ( coord[coordDir] <= thicknessStrutFactor * _mesh.getMeshSize( coordDir )  ) ||  ( coord[coordDir] >= 1. - thicknessStrutFactor * _mesh.getMeshSize( coordDir )  ) ) setBoundaryConstraints = true; 
                     }
                    
                    if( setBoundaryConstraints ){
                         bool use_u = true;
                         for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                           int dir = (coordDir + otherDir) % dimChartDomain;
                           if( std::abs( coord[dir] - 0.5  ) >= strutSize ) use_u = false;
                         }
                         if( use_u ){ material_l[nodeIdx] = _pfUpperBound;  initMaterial[nodeIdx] = _pfUpperBound;
                         }else{       material_u[nodeIdx] = _pfLowerBound;  initMaterial[nodeIdx] = _pfLowerBound; }
                    }
                  }
                }
            } break;
            
          // strut with distance depending on mesh size (in 2d: segment x segment, in 3d: square x segment) 
          // plus interface with of phasefield
          case 5: {
              const RealType strutSize = _parser.template get<RealType>("Boundary.strutSizeBoundaryConditionMaterial");
              const RealType thicknessStrutFactor = _parser.template get<RealType>("Boundary.thicknessStrutFactor");
              const RealType eps_factor = _parser.template get<RealType>("PhaseField.eps_factor");
              for( int nodeIdx=0; nodeIdx<initMaterial.size(); ++nodeIdx ){
                  const PointType& coord = _mesh.getVertex( nodeIdx );
                  for( int coordDir=0; coordDir < dimChartDomain; ++coordDir ){
                      
                     bool setBoundaryConstraints = false;
                     if ( hasPeriodicBoundary ) {
                         if( coord[coordDir] <= thicknessStrutFactor * _mesh.getMeshSize( coordDir )  ) setBoundaryConstraints = true;
                     }else{
                         if( ( coord[coordDir] <= thicknessStrutFactor * _mesh.getMeshSize( coordDir )  ) ||  ( coord[coordDir] >= 1. - thicknessStrutFactor * _mesh.getMeshSize( coordDir )  ) ) setBoundaryConstraints = true; 
                     }
                    
                    if( setBoundaryConstraints ){
                         int phasefieldType = 1; //0 - lowerBound, 1 - upperbound, 2 - between
                         for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                           const int dir = (coordDir + otherDir) % dimChartDomain;
                           if( std::abs( coord[dir] - 0.5  ) >= strutSize + eps_factor * _mesh.getMeshSize(dir) ){
                               phasefieldType = 0;
                               break;
                           }
                           if( std::abs( coord[dir] - 0.5  ) > strutSize ) phasefieldType = 2;
                         }
                         switch( phasefieldType ) {
                            case 0:{   material_u[nodeIdx] = _pfLowerBound;  initMaterial[nodeIdx] = _pfLowerBound;  }break;
                            case 1:{   material_l[nodeIdx] = _pfUpperBound;  initMaterial[nodeIdx] = _pfUpperBound;  }break;
                            case 2:{   initMaterial[nodeIdx] = 0.5 * (_pfLowerBound + _pfUpperBound);  }break;
                            default: break;
                         }
                    }
                  }
                }
            } break;
            
            
            
          //circle (in 2d: segment, in 3d: circle)
//           case 11: {
//               const RealType strutSize = _parser.template get<RealType>("Boundary.strutSizeBoundaryConditionMaterial");  
//               for( int nodeIdx=0; nodeIdx<initMaterial.size(); ++nodeIdx ){
//                   const PointType& coord = _mesh.getVertex( nodeIdx );
//                       
//                   for( int coordDir=0; coordDir < dimChartDomain; ++coordDir ){
//                   
//                      if( coord[coordDir] == 0 ){
//                          
//                          bool use_u = true;
//                          RealType distSqr = 0.;
//                          for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
//                            int dir = (coordDir + otherDir) % dimChartDomain;
//                            distSqr += pesopt::Sqr( coord[dir] - 0.5  );
//                          }
//                          if( distSqr >= strutSize * strutSize ) use_u = false;
//                          
//                          if( use_u ){
//                            material_l[nodeIdx] = _pfUpperBound;
//                            initMaterial[nodeIdx] = _pfUpperBound;
//                         }else{
//                             material_u[nodeIdx] = _pfLowerBound;
//                             initMaterial[nodeIdx] = _pfLowerBound;
//                         }
//                     }
//                     
//                   }
//                 }
//             } break;
            
            
          // TODO periodic boundary
          //Schwarz P
          case 20:{
              VectorType SchwarzPSurface ( initMaterial.size() );
              this->constructSchwarzPSurface( SchwarzPSurface );
              for( int nodeIdx=0; nodeIdx<initMaterial.size(); ++nodeIdx ){
                  if( _mesh.getBoundaryMask()[nodeIdx] ){
                    if( SchwarzPSurface[nodeIdx] == _pfUpperBound ){
                          material_l[nodeIdx] = _pfUpperBound;
                          initMaterial[nodeIdx] = _pfUpperBound;
                    }else{
                        material_u[nodeIdx] = _pfLowerBound;
                        initMaterial[nodeIdx] = _pfLowerBound;
                    }
                 }
            }
          } break;
          
          
       //=====================================================
       //        right triangle (0, 0), (1, 0),  (0, 1/sqrt(3)) 
       // ====================================================       
           
          // strut with distance depending on mesh size (in 2d: segment x segment, in 3d: TODO) 
          // plus interface with of phasefield
            case 105: {
              const RealType strutSize = _parser.template get<RealType>("Boundary.strutSizeBoundaryConditionMaterial");
              const RealType thicknessStrutFactor = _parser.template get<RealType>("Boundary.thicknessStrutFactor");
              const RealType eps_factor = _parser.template get<RealType>("PhaseField.eps_factor");
              for( int nodeIdx=0; nodeIdx<initMaterial.size(); ++nodeIdx ){
                  const PointType& coord = _mesh.getVertex( nodeIdx );
                  int coordDir=1;
                      
                    bool setBoundaryConstraints = false;
                    if( ( coord[coordDir] <= thicknessStrutFactor * _mesh.getInterfaceWith( )  ) ) setBoundaryConstraints = true; 
                    
                    if( setBoundaryConstraints ){
                         int phasefieldType = 1; //0 - lowerBound, 1 - upperbound, 2 - between
                         for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                           const int dir = (coordDir + otherDir) % dimChartDomain;
                           if( std::abs( coord[dir] - 0.5  ) >= strutSize + eps_factor * _mesh.getInterfaceWith() ){
                               phasefieldType = 0;
                               break;
                           }
                           if( std::abs( coord[dir] - 0.5  ) > strutSize ) phasefieldType = 2;
                         }
                         switch( phasefieldType ) {
                            case 0:{   material_u[nodeIdx] = _pfLowerBound;  initMaterial[nodeIdx] = _pfLowerBound;  }break;
                            case 1:{   material_l[nodeIdx] = _pfUpperBound;  initMaterial[nodeIdx] = _pfUpperBound;  }break;
                            case 2:{   initMaterial[nodeIdx] = 0.5 * (_pfLowerBound + _pfUpperBound);  }break;
                            default: break;
                         }
                    }
                }
            } break;
           
       //=====================================================
       //        equilateral triangle (1, 0), (0, sqrt(3)),  (-1, 0) 
       // ====================================================       
           
          // strut with distance depending on mesh size (in 2d: segment x segment, in 3d: TODO) 
          // plus interface with of phasefield
            case 205: {
              const RealType strutSize = _parser.template get<RealType>("Boundary.strutSizeBoundaryConditionMaterial");
              const RealType thicknessStrutFactor = _parser.template get<RealType>("Boundary.thicknessStrutFactor");
              const RealType eps_factor = _parser.template get<RealType>("PhaseField.eps_factor");
              for( int nodeIdx=0; nodeIdx<initMaterial.size(); ++nodeIdx ){
                  const PointType& coord = _mesh.getVertex( nodeIdx );
                  
                  int coordDir=1;
                  bool setBoundaryConstraints = false;
                  if( ( coord[coordDir] <= thicknessStrutFactor * _mesh.getInterfaceWith( )  ) ) setBoundaryConstraints = true; 
                  if( setBoundaryConstraints ){
                         int phasefieldType = 1; //0 - lowerBound, 1 - upperbound, 2 - between
                         for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                           const int dir = (coordDir + otherDir) % dimChartDomain;
                           if( ( std::abs( coord[dir] - 0.5  ) >= strutSize + eps_factor * _mesh.getInterfaceWith() )
                          &&   ( std::abs( coord[dir] + 0.5  ) >= strutSize + eps_factor * _mesh.getInterfaceWith() )
                             ){
                               phasefieldType = 0;
                               break;
                           }
                           if(   ( std::abs( coord[dir] - 0.5  ) > strutSize )
                              && ( std::abs( coord[dir] + 0.5  ) > strutSize )
                             )
                               phasefieldType = 2;
                         }
                         switch( phasefieldType ) {
                            case 0:{   material_u[nodeIdx] = _pfLowerBound;  initMaterial[nodeIdx] = _pfLowerBound;  }break;
                            case 1:{   material_l[nodeIdx] = _pfUpperBound;  initMaterial[nodeIdx] = _pfUpperBound;  }break;
                            case 2:{   initMaterial[nodeIdx] = 0.5 * (_pfLowerBound + _pfUpperBound);  }break;
                            default: break;
                         }
                   }
                   
                  coordDir=1; 
                  setBoundaryConstraints = false;
                  typename DataTypeContainer::DerivativeVectorValuedType RotationMatrix;
                  RotationMatrix(0, 0) = -0.5; RotationMatrix(0, 1) = 0.5 * sqrt(3);
                  RotationMatrix(1, 0) = -0.5 * sqrt(3); RotationMatrix(1, 1) = -0.5;
                  PointType d; d(0) = -0.5; d(1) = 0.5 * sqrt(3);
                  PointType rotatedCoords = RotationMatrix * coord + d;
                  if( ( rotatedCoords[coordDir] <= thicknessStrutFactor * _mesh.getInterfaceWith( )  ) ) setBoundaryConstraints = true; 
                  if( setBoundaryConstraints ){
                         int phasefieldType = 1; //0 - lowerBound, 1 - upperbound, 2 - between
                         for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                           const int dir = (coordDir + otherDir) % dimChartDomain;
                           if( ( std::abs( rotatedCoords[dir] - 0.5  ) >= strutSize + eps_factor * _mesh.getInterfaceWith() )
                          &&   ( std::abs( rotatedCoords[dir] + 0.5  ) >= strutSize + eps_factor * _mesh.getInterfaceWith() )
                             ){
                               phasefieldType = 0;
                               break;
                           }
                           if(   ( std::abs( rotatedCoords[dir] - 0.5  ) > strutSize )
                              && ( std::abs( rotatedCoords[dir] + 0.5  ) > strutSize )
                             )
                               phasefieldType = 2;
                         }
                         switch( phasefieldType ) {
                            case 0:{   material_u[nodeIdx] = _pfLowerBound;  initMaterial[nodeIdx] = _pfLowerBound;  }break;
                            case 1:{   material_l[nodeIdx] = _pfUpperBound;  initMaterial[nodeIdx] = _pfUpperBound;  }break;
                            case 2:{   initMaterial[nodeIdx] = 0.5 * (_pfLowerBound + _pfUpperBound);  }break;
                            default: break;
                         }
                   }
                   
                   
                   
                  coordDir=1; 
                  setBoundaryConstraints = false;
                  typename DataTypeContainer::DerivativeVectorValuedType RotationMatrix2;
                  RotationMatrix2(0, 0) = -0.5; RotationMatrix2(0, 1) = -0.5 * sqrt(3);
                  RotationMatrix2(1, 0) = 0.5 * sqrt(3); RotationMatrix2(1, 1) = -0.5;
                  PointType d2; d2(0) = 0.5; d2(1) = 0.5 * sqrt(3);
                  PointType rotatedCoords2 = RotationMatrix2 * coord + d2;
                  if( ( rotatedCoords2[coordDir] <= thicknessStrutFactor * _mesh.getInterfaceWith( )  ) ) setBoundaryConstraints = true; 
                  if( setBoundaryConstraints ){
                         int phasefieldType = 1; //0 - lowerBound, 1 - upperbound, 2 - between
                         for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                           const int dir = (coordDir + otherDir) % dimChartDomain;
                           if( ( std::abs( rotatedCoords2[dir] - 0.5  ) >= strutSize + eps_factor * _mesh.getInterfaceWith() )
                          &&   ( std::abs( rotatedCoords2[dir] + 0.5  ) >= strutSize + eps_factor * _mesh.getInterfaceWith() )
                             ){
                               phasefieldType = 0;
                               break;
                           }
                           if(   ( std::abs( rotatedCoords2[dir] - 0.5  ) > strutSize )
                              && ( std::abs( rotatedCoords2[dir] + 0.5  ) > strutSize )
                             )
                               phasefieldType = 2;
                         }
                         switch( phasefieldType ) {
                            case 0:{   material_u[nodeIdx] = _pfLowerBound;  initMaterial[nodeIdx] = _pfLowerBound;  }break;
                            case 1:{   material_l[nodeIdx] = _pfUpperBound;  initMaterial[nodeIdx] = _pfUpperBound;  }break;
                            case 2:{   initMaterial[nodeIdx] = 0.5 * (_pfLowerBound + _pfUpperBound);  }break;
                            default: break;
                         }
                   }    
                   
              }
            } break;
           
           
           
           
          default:
              throw std::invalid_argument( pesopt::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
              break;

        }
   }
  
};



#endif
