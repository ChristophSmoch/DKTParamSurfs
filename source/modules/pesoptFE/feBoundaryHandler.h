#ifndef __FEBOUNDARYHANDLER_H
#define __FEBOUNDARYHANDLER_H

 


//virtual base class
template< typename ConfiguratorType >
class FEBoundaryConditionHandler{
  
protected:
  
  typedef typename ConfiguratorType::RealType               RealType;
  typedef typename ConfiguratorType::InitType               MeshType;
  typedef typename ConfiguratorType::MaskType               MaskType;
  typedef typename ConfiguratorType::PointType              PointType;
  typedef typename ConfiguratorType::VectorType             VectorType;
  typedef typename ConfiguratorType::DTContainer            DataTypeContainer;
  typedef pesopt::BoostParser ParameterParserType;
  
  const ConfiguratorType &_conf;
  const MeshType &_mesh;
  const int _numVertices, _numGlobalDofs, _numElements;
  MaskType _DirichletMask;
  MaskType _PeriodicMask;
  std::vector<int> _PeriodicIndices;
  
public:
  
  FEBoundaryConditionHandler( const ConfiguratorType &conf ) :
    _conf ( conf ),
    _mesh( conf.getMesh() ), 
    _numElements ( _mesh.getNumElements() ),
    _numVertices( conf.getNumGlobalDofs() ),
    _numGlobalDofs ( conf.getNumGlobalDofs() ) {
        this->_DirichletMask.resize( this->_numGlobalDofs, false );
        this->_PeriodicMask.resize( this->_numGlobalDofs, false );
        this->_PeriodicIndices.resize( this->_numGlobalDofs ); //TODO if no periodic mask, then _PeriodicIndices(i) = i for all i;
    }
  
  
  virtual bool hasDirichletBoundary() const = 0;
  virtual bool hasPeriodicBoundary() const = 0;
  
  const MaskType &getDirichletMask ( ) const { return _DirichletMask; };
  const MaskType &getPeriodicMask ( ) const { return _PeriodicMask; };
  const std::vector<int> &getPeriodicIndices ( ) const { return _PeriodicIndices; };
  
  
   //! \note the following operations only delete the periodic value (i.e. set it to zero) 
  void collabseVector ( VectorType & vec ) const  {
    for ( int nodeIdx=0; nodeIdx<this->_mesh.getNumVertices(); nodeIdx++ ) {
        if( this->_PeriodicMask[nodeIdx] ) {
            vec[nodeIdx] = 0.0;
        }
    }
  }
  
  //! \note the following operations is additive: it adds the periodic values to the corresponding entry 
  void collabseVectorAdditive ( VectorType &vec ) const  {
    for ( int nodeIdx=0; nodeIdx<this->_mesh.getNumVertices(); nodeIdx++ ) {
        if( this->_PeriodicMask[nodeIdx] ){
            vec[this->_PeriodicIndices[nodeIdx]] += vec[nodeIdx];
            vec[nodeIdx] = 0.0;
        }
    }
  }
  
  void extendVector ( VectorType &vec ) const  {
    for ( int nodeIdx=0; nodeIdx<this->_mesh.getNumVertices(); nodeIdx++ ) {
      if( this->_PeriodicMask[nodeIdx] ) {
          vec[nodeIdx] = vec[this->_PeriodicIndices[nodeIdx]];
      }
    }
  }
  
  void extendVector ( const VectorType & vec, VectorType & vecPeriodicallyExtended  ) const  {
    vecPeriodicallyExtended = vec;
    this->extendVector( vecPeriodicallyExtended );
  }
  
  
  void collabseMultiVector ( VectorType & dispPeriodic ) const  {
    for ( int nodeIdx=0; nodeIdx<this->_mesh.getNumVertices(); nodeIdx++ ) {
      if( this->_PeriodicMask[nodeIdx] ){
          for( int comp=0; comp<this->_conf.dimChartDomain; ++comp) 
              dispPeriodic[nodeIdx + comp * this->_mesh.getNumVertices()] = 0.0;
      }
    }
  }
  
  void collabseMultiVectorAdditive ( VectorType & dispPeriodic ) const  {
    for ( int nodeIdx=0; nodeIdx<this->_mesh.getNumVertices(); nodeIdx++ ) {
      if( this->_PeriodicMask[nodeIdx] ){
          for( int comp=0; comp<this->_conf.dimChartDomain; ++comp){
              dispPeriodic[this->_PeriodicIndices[nodeIdx] + comp * this->_mesh.getNumVertices()] += dispPeriodic[nodeIdx + comp * this->_mesh.getNumVertices()];
              dispPeriodic[nodeIdx + comp * this->_mesh.getNumVertices()] = 0.0;
          }
      }
    }
  }
  
  void extendMultiVector ( VectorType & dispPeriodicallyExtended ) const  {
    for ( int nodeIdx=0; nodeIdx<this->_mesh.getNumVertices(); nodeIdx++ ) {
      if( this->_PeriodicMask[nodeIdx] ){
          for( int comp=0; comp<this->_conf.dimChartDomain; ++comp) 
              dispPeriodicallyExtended[nodeIdx + comp * this->_mesh.getNumVertices()] = dispPeriodicallyExtended[this->_PeriodicIndices[nodeIdx] + comp * this->_mesh.getNumVertices()];
      }
    }
  }
  
  void extendMultiVector ( const VectorType & dispPeriodic, VectorType & dispPeriodicallyExtended ) const  {
    dispPeriodicallyExtended = dispPeriodic;
    this->extendMultiVector( dispPeriodicallyExtended );
  }
  
  
  
};





template< typename ConfiguratorType >
class FEZeroBoundaryConditionHandler
: public FEBoundaryConditionHandler<ConfiguratorType>{
  
protected:
  
  typedef typename ConfiguratorType::RealType               RealType;
  typedef typename ConfiguratorType::InitType               MeshType;
  typedef typename ConfiguratorType::MaskType               MaskType;
  typedef typename ConfiguratorType::PointType              PointType;
  typedef typename ConfiguratorType::VectorType             VectorType;
  typedef typename ConfiguratorType::DTContainer            DataTypeContainer;
  typedef pesopt::BoostParser ParameterParserType;
  
public:
  
  FEZeroBoundaryConditionHandler( const ConfiguratorType &conf ) : 
    FEBoundaryConditionHandler<ConfiguratorType> ( conf ){
      this->_DirichletMask.resize( this->_numGlobalDofs, false );
      this->_PeriodicMask.resize( this->_numGlobalDofs, false );
      this->_PeriodicIndices.resize( this->_numGlobalDofs );
    }
  
  bool hasDirichletBoundary() const override { return false; }
  bool hasPeriodicBoundary() const override { return false; }
  
};





template< typename ConfiguratorType >
class FEDirichletBoundaryConditionHandler
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
  
  FEDirichletBoundaryConditionHandler( const ParameterParserType &parser, const ConfiguratorType &conf ) : 
    FEBoundaryConditionHandler<ConfiguratorType> ( conf ),
    _parser ( parser ) {
        this->_DirichletMask.resize( this->_numGlobalDofs, false );
        this->_DirichletMask = this->_mesh.getBoundaryMask();
    }
  
  bool hasDirichletBoundary() const override { return true; }
  bool hasPeriodicBoundary() const override { return false; }
  
};



template< typename ConfiguratorType >
class FEDirichletBoundaryConditionHandlerLeftSide
: public FEBoundaryConditionHandler<ConfiguratorType>{
  
protected:
  
  typedef typename ConfiguratorType::RealType               RealType;
  typedef typename ConfiguratorType::InitType               MeshType;
  typedef typename ConfiguratorType::MaskType               MaskType;
  typedef typename ConfiguratorType::PointType              PointType;
  typedef typename ConfiguratorType::VectorType             VectorType;
  typedef typename ConfiguratorType::DTContainer            DataTypeContainer;
  typedef pesopt::BoostParser ParameterParserType;
  
//   const ParameterParserType &_parser;
  
public:
  
//   FEDirichletBoundaryConditionHandlerLeftSide( const ParameterParserType &parser, const ConfiguratorType &conf ) : 
    FEDirichletBoundaryConditionHandlerLeftSide( const ConfiguratorType &conf ) : 
    FEBoundaryConditionHandler<ConfiguratorType> ( conf )
//     _parser ( parser )
    {
        this->_DirichletMask.resize( this->_numGlobalDofs, false );
        for ( int nodeIdx=0; nodeIdx < this->_numVertices; ++nodeIdx ) {
            const PointType& coords ( this->_mesh.getVertex(nodeIdx) );
            if ( coords [0] == 0. ) this->_DirichletMask[nodeIdx] = true;
        }
    }
  
  bool hasDirichletBoundary() const override { return true; }
  bool hasPeriodicBoundary() const override { return false; }
  
};


template< typename ConfiguratorType >
class FEDirichletBoundaryConditionHandlerDownSide
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
  
  FEDirichletBoundaryConditionHandlerDownSide( const ParameterParserType &parser, const ConfiguratorType &conf ) : 
    FEBoundaryConditionHandler<ConfiguratorType> ( conf ),
    _parser ( parser ) {
        this->_DirichletMask.resize( this->_numGlobalDofs, false );
        for ( int nodeIdx=0; nodeIdx < this->_numVertices; ++nodeIdx ) {
            const PointType& coords ( this->_mesh.getVertex(nodeIdx) );
            if ( coords [1] == 0. ) this->_DirichletMask[nodeIdx] = true;
        }
    }
  
  bool hasDirichletBoundary() const override { return true; }
  bool hasPeriodicBoundary() const override { return false; }
  
};



template< typename ConfiguratorType >
class FEDirichletBoundaryConditionHandlerRectangleStruts
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
  
  FEDirichletBoundaryConditionHandlerRectangleStruts( const ParameterParserType &parser, const ConfiguratorType &conf ) : 
    FEBoundaryConditionHandler<ConfiguratorType> ( conf ),
    _parser ( parser ) {
        
        this->_DirichletMask.resize( this->_numGlobalDofs, false );
        
         const RealType strutSize = _parser.template get<RealType>("Boundary.strutSizeBoundaryConditionDisplacement"); 
         const int dimChartDomain = static_cast<int> ( ConfiguratorType::dimChartDomain );
         
         for( int nodeIdx=0; nodeIdx<this->_numVertices; ++nodeIdx ){
            const PointType& coord = this->_mesh.getVertex( nodeIdx );
            for( int coordDir=0; coordDir < dimChartDomain; ++coordDir ){
                bool setBoundaryConstraints = false;
                if( (coord[coordDir] == 0 ) ||  ( coord[coordDir] == 1 ) ) setBoundaryConstraints = true; 
                if( setBoundaryConstraints ){
                    bool use_u = true;
                    for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                        int dir = (coordDir + otherDir) % dimChartDomain;
                        if( std::abs( coord[dir] - 0.5  ) >= strutSize ) use_u = false;
                    }
                    if( use_u ) this->_DirichletMask[nodeIdx] = true;
                }
                    
            }
        }
    }
  
  bool hasDirichletBoundary() const override { return true; }
  bool hasPeriodicBoundary() const override { return false; }
  
};



template< typename ConfiguratorType >
class FEDirichletBoundaryConditionHandlerRightTriangle3060Struts
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
  
  FEDirichletBoundaryConditionHandlerRightTriangle3060Struts( const ParameterParserType &parser, const ConfiguratorType &conf ) : 
    FEBoundaryConditionHandler<ConfiguratorType> ( conf ),
    _parser ( parser ) {
        
        this->_DirichletMask.resize( this->_numGlobalDofs, false );
        
         const RealType strutSize = _parser.template get<RealType>("Boundary.strutSizeBoundaryConditionDisplacement"); 
         const int dimChartDomain = static_cast<int> ( ConfiguratorType::dimChartDomain );
         
         for( int nodeIdx=0; nodeIdx<this->_numVertices; ++nodeIdx ){
            const PointType& coord = this->_mesh.getVertex( nodeIdx );
            int coordDir = 1;
            bool setBoundaryConstraints = false;
            if( (coord[coordDir] == 0 ) ) setBoundaryConstraints = true; 
            if( setBoundaryConstraints ){
                bool use_u = true;
                for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                    int dir = (coordDir + otherDir) % dimChartDomain;
                    if( std::abs( coord[dir] - 0.5  ) > strutSize ) use_u = false;
                }
                if( use_u ) this->_DirichletMask[nodeIdx] = true;
            }
                    
         }
    }
  
  bool hasDirichletBoundary() const override { return true; }
  bool hasPeriodicBoundary() const override { return false; }
  
};




template< typename ConfiguratorType >
class FEDirichletBoundaryConditionHandlerEquilateralTriangleStruts
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
  
  FEDirichletBoundaryConditionHandlerEquilateralTriangleStruts( const ParameterParserType &parser, const ConfiguratorType &conf ) : 
    FEBoundaryConditionHandler<ConfiguratorType> ( conf ),
    _parser ( parser ) {
        
        this->_DirichletMask.resize( this->_numGlobalDofs, false );
        
         const RealType strutSize = _parser.template get<RealType>("Boundary.strutSizeBoundaryConditionDisplacement"); 
         const int dimChartDomain = static_cast<int> ( ConfiguratorType::dimChartDomain );
         
         for( int nodeIdx=0; nodeIdx<this->_numVertices; ++nodeIdx ){
            const PointType& coord = this->_mesh.getVertex( nodeIdx );
            
            int coordDir = 1;
            bool setBoundaryConstraints = false;
            if( ( coord[coordDir] <= 0.01 * this->_mesh.getInterfaceWith( )  ) ) setBoundaryConstraints = true; 
            if( setBoundaryConstraints ){
                bool use_u = true;
                for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                    int dir = (coordDir + otherDir) % dimChartDomain;
                    if(  ( std::abs( coord[dir] - 0.5  ) > strutSize )
                      && ( std::abs( coord[dir] + 0.5  ) > strutSize )
                      ) use_u = false;
                }
                if( use_u ) this->_DirichletMask[nodeIdx] = true;
            }
            
            
            coordDir=1; 
            setBoundaryConstraints = false;
            typename DataTypeContainer::DerivativeVectorValuedType RotationMatrix;
            RotationMatrix(0, 0) = -0.5; RotationMatrix(0, 1) = 0.5 * sqrt(3);
            RotationMatrix(1, 0) = -0.5 * sqrt(3); RotationMatrix(1, 1) = -0.5;
            PointType d; d(0) = -0.5; d(1) = 0.5 * sqrt(3);
            PointType rotatedCoords = RotationMatrix * coord + d;
            if( ( rotatedCoords[coordDir] <= 0.01 * this->_mesh.getInterfaceWith( )  ) ) setBoundaryConstraints = true; 
            if( setBoundaryConstraints ){
                bool use_u = true;
                for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                    int dir = (coordDir + otherDir) % dimChartDomain;
                    if(   ( std::abs( rotatedCoords[dir] - 0.5  ) > strutSize )
                       && ( std::abs( rotatedCoords[dir] + 0.5  ) > strutSize )
                      ) use_u = false;
                }
                if( use_u ) this->_DirichletMask[nodeIdx] = true;
            }
            
            coordDir=1; 
            setBoundaryConstraints = false;
            typename DataTypeContainer::DerivativeVectorValuedType RotationMatrix2;
            RotationMatrix2(0, 0) = -0.5; RotationMatrix2(0, 1) = -0.5 * sqrt(3);
            RotationMatrix2(1, 0) = 0.5 * sqrt(3); RotationMatrix2(1, 1) = -0.5;
            PointType d2; d2(0) = 0.5; d2(1) = 0.5 * sqrt(3);
            PointType rotatedCoords2 = RotationMatrix2 * coord + d2;
            if( ( rotatedCoords2[coordDir] <= 0.01 * this->_mesh.getInterfaceWith( )  ) ) setBoundaryConstraints = true; 
            if( setBoundaryConstraints ){
                bool use_u = true;
                for( int otherDir=1; otherDir < dimChartDomain; ++otherDir ){
                    int dir = (coordDir + otherDir) % dimChartDomain;
                    if(   ( std::abs( rotatedCoords2[dir] - 0.5  ) > strutSize )
                       && ( std::abs( rotatedCoords2[dir] + 0.5  ) > strutSize )
                      ) use_u = false;
                }
                if( use_u ) this->_DirichletMask[nodeIdx] = true;
            }
                    
         }
    }
  
  bool hasDirichletBoundary() const override { return true; }
  bool hasPeriodicBoundary() const override { return false; }
  
};




// TODO clamped wrt one ring of boundary nodes
//      maybe distinguish if C1 dofs or tangent space derivatives are considered
template< typename ConfiguratorType >
class FEClampedBoundaryConditionHandler
: public FEBoundaryConditionHandler<ConfiguratorType>{
  
protected:
  
  typedef typename ConfiguratorType::RealType               RealType;
  typedef typename ConfiguratorType::InitType               MeshType;
  typedef typename ConfiguratorType::MaskType               MaskType;
  typedef typename ConfiguratorType::PointType              PointType;
  typedef typename ConfiguratorType::VectorType             VectorType;
  typedef typename ConfiguratorType::DTContainer            DataTypeContainer;
  typedef typename ConfiguratorType::ElementType            ElementType;
  typedef pesopt::BoostParser ParameterParserType;
  
  const ParameterParserType &_parser;
  
public:
  
  FEClampedBoundaryConditionHandler( const ParameterParserType &parser, const ConfiguratorType &conf ) : 
    FEBoundaryConditionHandler<ConfiguratorType> ( conf ),
    _parser ( parser ) {
        this->_DirichletMask.resize( this->_numGlobalDofs, false );
        this->_DirichletMask = this->_mesh.getBoundaryMask();
        MaskType DirichletMask( this->_DirichletMask );
        const int numNodesPerElement = static_cast<int>  ( ElementType::numNodes );
        for( int nodeIdx=0; nodeIdx< this->_numVertices; ++nodeIdx ){
            if( DirichletMask[nodeIdx] ){
                std::vector<int> commonElements;
                this->_mesh.getCommonElements( nodeIdx, commonElements );
                for( int i=0; i<commonElements.size(); ++i ){
                    int elementIdx = commonElements[i];
                    for( int localIndex = 0; localIndex < numNodesPerElement; ++localIndex ) 
                        this->_DirichletMask[this->_mesh.getElement( elementIdx ).getGlobalNodeIdx(localIndex)] = true;
                }
            }
          }
    }
  
  bool hasDirichletBoundary() const override { return true; }
  bool hasPeriodicBoundary() const override { return false; }
  
};







#endif //__FEBOUNDARYHANDLER_H
