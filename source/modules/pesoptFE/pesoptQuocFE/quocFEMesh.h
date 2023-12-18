#ifndef __QUOCFEMESH_H
#define __QUOCFEMESH_H

#include <pesopt_IO.h>
#include <pesopt_VTK.h>
#include "quocFEDefines.h"
#include "quocFEElement.h"
    


//Line (0) to (lx)
// DOFs in x: Nx
// meshsize: (hx) = (lx/(Nx-1))
template< typename DataTypeContainer = Quoc1DDataTypeContainer,
          typename QuocElementType = QuocElement1D<DataTypeContainer> >
class QuocMesh1D {
public:
  typedef QuocElementType ElementType;
  typedef typename ElementType::ElementNodeIndicesType      ElementNodeIndicesType;
//   typedef QuocBoundaryElementType                   BoundaryElementType;
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType              RealType;
  typedef typename DataTypeContainer::PointType             PointType;
  typedef typename DataTypeContainer::VectorType            VectorType;
  typedef typename DataTypeContainer::MaskType              MaskType;
  typedef typename DataTypeContainer::IntVecChart           IntVecChart;
  typedef typename DataTypeContainer::RealVecChart          RealVecChart; 
  typedef pesopt::BoostParser ParameterParserType;
  
  static const int dimChartDomain = 1;
  static const int _VTKCELLTYPE = 3; // VTK_LINE (=3)
  static const VTKDATATYPE _VTKDATATYPEUNDEFORMED = VTKSTRUCTUREDPOINTS;
  static const VTKDATATYPE _VTKDATATYPEDEFORMED =  VTKUNSTRUCTUREDGRID;
  static const PESOPT_FEMeshType _FEMeshType = QUOCMESH;
  
protected:
  int _Nx;
  IntVecChart _NumDofVec;
  RealType _lx;
  PointType _LengthVec;
  RealType _hx;
  std::vector< PointType > _vertexIterator;
  std::vector< ElementType > _elementIterator;
//   std::vector< BoundaryElementType > _boundaryElementIterator;
  mutable std::vector< std::vector<int>  > _vertexIteratorAdjacentElements;
  mutable std::vector< std::vector<int>  > _vertexIteratorAdjacentElementsLocalIdx;
  
public :
  MaskType _boundaryLeft, _boundaryRight;
  MaskType _boundaryMask;
  MaskType _boundaryPeriodic;
  std::vector<int> _periodicIdentificationIndices;
  
  
public:
  //! Create empty QuocMesh
//  QuocMesh1D ( ) : _vertexIterator(), _elementIterator()
// //   , _boundaryElementIterator() 
//   { }
  
  QuocMesh1D ( const int Nx, const RealType lx ) : 
    _Nx (Nx), _lx( lx ), _hx( _lx / static_cast<RealType>(_Nx - 1 ) )
  {
        _NumDofVec[0] = _Nx;
        _LengthVec[0] = _lx;
      
        const int numGlobalDofs = _Nx;
        _boundaryLeft.resize( numGlobalDofs, false ); _boundaryRight.resize( numGlobalDofs, false );
        _boundaryMask.resize( numGlobalDofs, false );
        _boundaryPeriodic.resize( numGlobalDofs, false );
                
        for( int xIdx=0; xIdx<_Nx; ++xIdx ){
              PointType vertex ( _lx * xIdx / static_cast<RealType>( _Nx - 1 ) );
              int nodeIdx = this->pushBackVertex( vertex );
              if( xIdx == 0 ) _boundaryLeft[nodeIdx] = true;
              if( xIdx == _Nx - 1 ) _boundaryRight[nodeIdx] = true;
              if( (xIdx == 0) || ( xIdx == _Nx - 1 )  ) _boundaryMask[nodeIdx] = true;
              if( xIdx == _Nx - 1 ) _boundaryPeriodic[nodeIdx] = true;
        }
          
      for( int nodeIdx = 0; nodeIdx < numGlobalDofs; ++nodeIdx ){
          if( (!_boundaryRight[nodeIdx]) ){
              //! \note ordering of elements has to coincide with ordering of base functions
              const int numNewElement = this->pushBackElement( ElementNodeIndicesType(nodeIdx, nodeIdx + 1 ) );  
//               if( this->getVertex(nodeIdx) == 0 ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::LEFT, _hx );
//               if( this->getVertex(nodeIdx + 1) == _lx ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::RIGHT, _hx );
          }
      }
      
      //Periodic indices to identify facing presenting node
      _periodicIdentificationIndices.resize( numGlobalDofs );
      for( int xIdx=0; xIdx<_Nx; ++xIdx ){
            int periodicXIdx = xIdx;
            if( xIdx == _Nx - 1 ) periodicXIdx = 0;
            _periodicIdentificationIndices[this->getGlobalNodeIndex(xIdx)] = this->getGlobalNodeIndex(periodicXIdx);
      }
      
      this->updateAllAdjacentElementsForNodes();
      
  }
  
  QuocMesh1D ( const IntVecChart &NumDofVec, const PointType &LenghtVec ) : QuocMesh1D ( NumDofVec[0], LenghtVec[0] ) {} 
    
    
   QuocMesh1D ( const string& fileName ){
      throw std::logic_error("constructor of QuocMesh1D is not defined for filename!");
   }
    
protected:
  void updateAllAdjacentElementsForNodes () const {
    _vertexIteratorAdjacentElements.clear();
    _vertexIteratorAdjacentElementsLocalIdx.clear();
    _vertexIteratorAdjacentElements.resize( this->getNumVertices() );
    _vertexIteratorAdjacentElementsLocalIdx.resize( this->getNumVertices() );
    for ( int elementIndex = 0; elementIndex < getNumElements(); ++elementIndex  )
        for(int localNodeIdx=0; localNodeIdx<2; ++localNodeIdx){
           const int globalNodeIdx = _elementIterator[elementIndex].getGlobalNodeIdx(localNodeIdx);
           _vertexIteratorAdjacentElements[globalNodeIdx].push_back(_elementIterator[elementIndex].getGlobalElementIdx());
           _vertexIteratorAdjacentElementsLocalIdx[globalNodeIdx].push_back(localNodeIdx);
        }
  }

public:
  void getCommonElements( const int nodeIdx1, std::vector<int> &commonElements ) const{
      for(int i=0;i<_vertexIteratorAdjacentElements[nodeIdx1].size(); ++i ) commonElements.push_back(_vertexIteratorAdjacentElements[nodeIdx1][i]); 
  }
    
    
public:
  int getNumVertices ( ) const { return ( _vertexIterator.size() );}
  int getNumElements ( ) const { return ( static_cast<int> ( _elementIterator.size() ) );}
//   int getNumBoundaryElements ( ) const { return ( static_cast<int> ( _boundaryElementIterator.size() ) );}
  
  //! insert new vertex and return global index
  int pushBackVertex ( const PointType newVertex ) {
    _vertexIterator.push_back ( newVertex );
    return getNumVertices() - 1;
  }
  //! insert new element and return global index
  int pushBackElement ( const ElementNodeIndicesType nodeIdx ) {
    int globalIdx = getNumElements();
    _elementIterator.push_back ( QuocElementType( globalIdx, nodeIdx, _vertexIterator ) );
    return globalIdx;
  }
  
//   int pushBackBoundaryElement ( const int elementIndex, QuocBoundaryType boundaryType, const RealType h ) {
//     int globalBoundaryIdx = getNumBoundaryElements();
// //     _boundaryElementIterator.push_back ( BoundaryElementType( globalBoundaryIdx, this->getElement(elementIndex), boundaryType, h ) );
//     return globalBoundaryIdx;
//   }

  const PointType& getVertex ( const int num ) const { return _vertexIterator[num];}
  void setVertex ( const int num, const PointType Arg ) { _vertexIterator[num] = Arg;}

  int getGlobalNodeIndex ( const int xIdx ) const {  return xIdx;}
  //int getGlobalIndex ( const CoordType &Coord ) const { return ( getGlobalNodeIndex ( Coord[0], Coord[1] ) );}
  int getGlobalElementIndex ( const int xFac ) const { return xFac;} 
  
  const QuocElementType& getElement ( const int num ) const {return _elementIterator[num];}
  void setElement ( const int num, const QuocElementType Arg ) { _elementIterator[num] = Arg;}

//   const BoundaryElementType& getBoundaryElement ( const int num ) const {return _boundaryElementIterator[num];}
//   BoundaryElementType & getBoundaryElement ( const int num ) {return _boundaryElementIterator[num];}
  
  int getElementNodeIdx  ( const int num, const int localNode ) const { return _elementIterator[num].getGlobalNodeIdx(localNode);}
  void setElementNodeIdx ( const int num, const int localNode, const int newIdx ) {_elementIterator[num].setGlobalNodeIdx(localNode, newIdx );}
  
  const RealType getMeshSize( const int direction ) const{ 
    if( direction == 0 ) return _hx;
    if( direction == 1 ) return 0.;
    if( direction == 2 ) return 0.;
    throw std::invalid_argument( pesopt::strprintf ( "Wrong direction in getMeshSize in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  const RealType getInterfaceWith( ) const{ return _hx; }
  const RealType getWidth( const int direction ) const {
   switch( direction ){
        case 0 :    return _lx;    break;
        default: throw std::invalid_argument( pesopt::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
    }
  }
  
  const int getNumDofs( const int direction ) const {
   switch( direction ){
        case 0 : return _Nx; break;
        case 1 : return 1; break; //for compatibility with vtk
        case 2 : return 1; break; //for compatibility with vtk
        default: throw std::invalid_argument( pesopt::strprintf ( "Wrong channel in QuocMesh1D get num dofs. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
    }
  }
  
  void getNumDofVec( IntVecChart &vec ) const { vec = _NumDofVec;}
  void getLenghtVec( PointType &vec ) const { vec = _LengthVec;}
  
  int getNumNodesOfElement(const int elIdx ) const { return 2; }
  
  const MaskType &getBoundaryMask( ) const { return _boundaryMask; }
  
  void saveToParser( ParameterParserType & parser ) const {
      parser.set ( "InputMesh.dimChartDomain", static_cast<int>( dimChartDomain ) );
      parser.set ( "InputMesh.FEType", "Line" );
      parser.set ( "InputMesh.MeshType", "Line" );
      parser.template setFixSizeVector<IntVecChart>("InputMesh.NumDofVec", _NumDofVec, true );
      parser.template setFixSizeVector<PointType>("InputMesh.LengthVec", _LengthVec, false );
  }
  
  void saveToParser( ParameterParserType & parser, const string fileNameMesh ) const {
//       parser.set ( "InputMesh.dimChartDomain", static_cast<int>( dimChartDomain ) );
//       parser.set ( "InputMesh.FEType", "Line" );
//       parser.set ( "InputMesh.MeshType", "File" );
//       parser.set ( "InputMesh.fileName", fileNameMesh ); 
      this->saveToParser( parser );
  }
  
  
};
    
    
    
    



//Rectangle (0,0) to (lx,ly)
// DOFs in x: Nx, y:Ny
// meshsize: (hx,hy) = (lx/(Nx-1), ly/(Ny-1))
template< typename DataTypeContainer = Quoc2DDataTypeContainer,
          typename QuocElementType = QuocElement2D<DataTypeContainer>, 
          typename QuocBoundaryElementType = QuocBoundaryElement2D<DataTypeContainer> >
class QuocMesh2D {
public:
  typedef QuocElementType ElementType;
  typedef typename ElementType::ElementNodeIndicesType          ElementNodeIndicesType;
  typedef QuocBoundaryElementType                               BoundaryElementType;
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType                  RealType;
  typedef typename DataTypeContainer::PointType                 PointType;
  typedef typename DataTypeContainer::VectorType                VectorType;
  typedef typename DataTypeContainer::MaskType                  MaskType;
  typedef typename DataTypeContainer::IntVecChart               IntVecChart;
  typedef typename DataTypeContainer::RealVecChart              RealVecChart;  
  typedef pesopt::BoostParser ParameterParserType;
  
  static const int dimChartDomain = 2;
  static const int _VTKCELLTYPE = 8; //VTK_PIXEL
  static const VTKDATATYPE _VTKDATATYPEUNDEFORMED = VTKSTRUCTUREDPOINTS;
  static const VTKDATATYPE _VTKDATATYPEDEFORMED =  VTKUNSTRUCTUREDGRID;
  static const PESOPT_FEMeshType _FEMeshType = QUOCMESH;
  
protected:
  int _Nx, _Ny;
  IntVecChart _NumDofVec;
  RealType _lx, _ly;
  PointType _LengthVec;
  RealType _hx, _hy;
  std::vector< PointType > _vertexIterator;
  std::vector< ElementType > _elementIterator;
  std::vector< BoundaryElementType > _boundaryElementIterator;
  mutable std::vector< std::vector<int>  > _vertexIteratorAdjacentElements;
  mutable std::vector< std::vector<int>  > _vertexIteratorAdjacentElementsLocalIdx;
public :
  MaskType _boundaryTop, _boundaryBottom, _boundaryLeft, _boundaryRight;
  MaskType _boundaryMask;
  MaskType _boundaryPeriodic;
  std::vector<int> _periodicIdentificationIndices;
  
public:
  //! Create empty QuocMesh
// !  QuocMesh2D ( ) : _vertexIterator(), _elementIterator(), _boundaryElementIterator()  { }
  
  QuocMesh2D ( const int Nx, const int Ny, const RealType lx, const RealType ly ) : 
    _Nx (Nx), _Ny(Ny), 
    _lx( lx ), _ly( ly ),
    _hx( _lx / static_cast<RealType>(_Nx - 1 ) ), _hy( _ly / static_cast<RealType>(_Ny - 1 ) )
  {
        _NumDofVec[0] = _Nx; _NumDofVec[1] = _Ny;
        _LengthVec[0] = _lx; _LengthVec[1] = _ly;
      
        const int numGlobalDofs = _Nx * _Ny;
        _boundaryTop.resize( numGlobalDofs, false ); _boundaryBottom.resize( numGlobalDofs, false );
        _boundaryLeft.resize( numGlobalDofs, false ); _boundaryRight.resize( numGlobalDofs, false );
        _boundaryMask.resize( numGlobalDofs, false );
        _boundaryPeriodic.resize( numGlobalDofs, false );
        
        for( int yIdx=0; yIdx<_Ny; ++yIdx )
            for( int xIdx=0; xIdx<_Nx; ++xIdx ){
              PointType vertex ( _lx * xIdx / static_cast<RealType>( _Nx - 1 ),  _ly * yIdx / static_cast<RealType>( _Ny - 1 ) );
              int nodeIdx = this->pushBackVertex( vertex );
              if( xIdx == 0 ) _boundaryLeft[nodeIdx] = true;
              if( xIdx == _Nx - 1 ) _boundaryRight[nodeIdx] = true;
              if( yIdx == 0 ) _boundaryBottom[nodeIdx] = true;
              if( yIdx == _Ny - 1 ) _boundaryTop[nodeIdx] = true;
              if( (xIdx == 0) || ( xIdx == _Nx - 1 ) || ( yIdx == 0 ) || ( yIdx == _Ny - 1 )  ) _boundaryMask[nodeIdx] = true;
              if( ( xIdx == _Nx - 1 ) || ( yIdx == _Ny - 1 ) ) _boundaryPeriodic[nodeIdx] = true;
          }
          
      for( int nodeIdx = 0; nodeIdx < numGlobalDofs; ++nodeIdx ){
          if( (!_boundaryRight[nodeIdx]) && (! _boundaryTop[nodeIdx]) ){
              //! \note ordering of elements has to coincide with ordering of base functions
              const int numNewElement = this->pushBackElement( ElementNodeIndicesType(nodeIdx, nodeIdx + 1, nodeIdx + _Nx, nodeIdx + _Nx + 1 ) );  
              if( this->getVertex(nodeIdx)[0] == 0 ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::LEFT, _hy );
              if( this->getVertex(nodeIdx + 1)[0] == _lx ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::RIGHT, _hy );
              if( this->getVertex(nodeIdx)[1] == 0 ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::BOTTOM, _hx );
              if( this->getVertex(nodeIdx + _Nx)[1] == _ly ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::TOP, _hx );
          }
      }
      
      //Periodic indices to identify facing presenting node
      _periodicIdentificationIndices.resize( numGlobalDofs );
      for( int yIdx=0; yIdx<_Ny; ++yIdx )
          for( int xIdx=0; xIdx<_Nx; ++xIdx ){
                int periodicXIdx = xIdx, periodicYIdx = yIdx;
                if( xIdx == _Nx - 1 ) periodicXIdx = 0;
                if( yIdx == _Ny - 1 ) periodicYIdx = 0;
                _periodicIdentificationIndices[this->getGlobalNodeIndex(xIdx,yIdx)] = this->getGlobalNodeIndex(periodicXIdx,periodicYIdx);
          }
      
      this->updateAllAdjacentElementsForNodes();
  }
  
  QuocMesh2D ( const IntVecChart &NumDofVec, const PointType &LenghtVec ) : QuocMesh2D ( NumDofVec[0], NumDofVec[1], LenghtVec[0], LenghtVec[1] ) {} 
  
  QuocMesh2D ( const string& fileName ){
      throw std::logic_error("constructor of QuocMesh2D is not defined for filename!");
   }
   
protected:
  void updateAllAdjacentElementsForNodes () const {
    _vertexIteratorAdjacentElements.clear();
    _vertexIteratorAdjacentElementsLocalIdx.clear();
    _vertexIteratorAdjacentElements.resize( this->getNumVertices() );
    _vertexIteratorAdjacentElementsLocalIdx.resize( this->getNumVertices() );
    for ( int elementIndex = 0; elementIndex < getNumElements(); ++elementIndex  )
        for(int localNodeIdx=0; localNodeIdx<2; ++localNodeIdx){
           const int globalNodeIdx = _elementIterator[elementIndex].getGlobalNodeIdx(localNodeIdx);
           _vertexIteratorAdjacentElements[globalNodeIdx].push_back(_elementIterator[elementIndex].getGlobalElementIdx());
           _vertexIteratorAdjacentElementsLocalIdx[globalNodeIdx].push_back(localNodeIdx);
        }
  }

public:
  void getCommonElements( const int nodeIdx1, std::vector<int> &commonElements ) const{
      for(int i=0;i<_vertexIteratorAdjacentElements[nodeIdx1].size(); ++i ) commonElements.push_back(_vertexIteratorAdjacentElements[nodeIdx1][i]); 
  }

public:
  int getNumVertices ( ) const { return ( _vertexIterator.size() );}
  int getNumElements ( ) const { return ( static_cast<int> ( _elementIterator.size() ) );}
  int getNumBoundaryElements ( ) const { return ( static_cast<int> ( _boundaryElementIterator.size() ) );}
  
  //! insert new vertex and return global index
  int pushBackVertex ( const PointType newVertex ) {
    _vertexIterator.push_back ( newVertex );
    return getNumVertices() - 1;
  }
  //! insert new element and return global index
  int pushBackElement ( const ElementNodeIndicesType nodeIdx ) {
    int globalIdx = getNumElements();
    _elementIterator.push_back ( QuocElementType( globalIdx, nodeIdx, _vertexIterator ) );
    return globalIdx;
  }
  
  int pushBackBoundaryElement ( const int elementIndex, QuocBoundaryType boundaryType, const RealType h ) {
    int globalBoundaryIdx = getNumBoundaryElements();
    _boundaryElementIterator.push_back ( BoundaryElementType( globalBoundaryIdx, this->getElement(elementIndex), boundaryType, h ) );
    return globalBoundaryIdx;
  }

  const PointType& getVertex ( const int num ) const { return _vertexIterator[num];}
  void setVertex ( const int num, const PointType Arg ) { _vertexIterator[num] = Arg;}

  int getGlobalNodeIndex ( const int xIdx, const int yIdx ) const {  return ( xIdx + yIdx *_Nx );}
  int getGlobalNodeIndex ( const IntVecChart &indices ) const { return this->getGlobalNodeIndex( indices[0], indices[1] ); }
  //int getGlobalIndex ( const CoordType &Coord ) const { return ( getGlobalNodeIndex ( Coord[0], Coord[1] ) );}
  int getGlobalElementIndex ( const int xFac, const int yFac ) const { return ( xFac + yFac * (_Nx - 1 ) );} 
  
  const QuocElementType& getElement ( const int num ) const {return _elementIterator[num];}
  void setElement ( const int num, const QuocElementType Arg ) { _elementIterator[num] = Arg;}

  const BoundaryElementType& getBoundaryElement ( const int num ) const {return _boundaryElementIterator[num];}
  
  int getElementNodeIdx  ( const int num, const int localNode ) const { return _elementIterator[num].getGlobalNodeIdx(localNode);}
  void setElementNodeIdx ( const int num, const int localNode, const int newIdx ) {_elementIterator[num].setGlobalNodeIdx(localNode, newIdx );}
  
  const RealType getMeshSize( const int direction ) const{ 
    if( direction == 0 ) return _hx;
    if( direction == 1 ) return _hy;
    if( direction == 2 ) return 0.;
    throw std::invalid_argument( pesopt::strprintf ( "Wrong direction in getMeshSize in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  const RealType getInterfaceWith( ) const{ return std::sqrt(_hx * _hy ); }
  const RealType getWidth( const int direction ) const {
   switch( direction ){
        case 0 :    return _lx;    break;
        case 1 : return _ly; break;
        default: throw std::invalid_argument( pesopt::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
    }
  }
  
  const int getNumDofs( const int direction ) const {
   switch( direction ){
        case 0 : return _Nx; break;
        case 1 : return _Ny; break;
        case 2 : return 1; break; //for compatibility with vtk
        default: throw std::invalid_argument( pesopt::strprintf ( "Wrong channel in QuocMesh2D get num dofs. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
    }
  }
  
  void getNumDofVec( IntVecChart &vec ) const { vec = _NumDofVec;}
  void getLenghtVec( PointType &vec ) const { vec = _LengthVec;}
  
  int getNumNodesOfElement(const int elIdx ) const { return 4; }
  
  
  IntVecChart getNodeIndices (int nodeIdx) const{
      int yIdx = nodeIdx / _Nx;
      int xIdx = nodeIdx % _Nx;
      IntVecChart indices (xIdx, yIdx); 
      return indices; 
  }
  
  const MaskType &getBoundaryMask( ) const { return _boundaryMask; }
  
  void saveToParser( ParameterParserType & parser ) const {
      parser.set ( "InputMesh.dimChartDomain", static_cast<int>( dimChartDomain ) );
      parser.set ( "InputMesh.FEType", "Rectangle" );
      parser.set ( "InputMesh.MeshType", "Rectangle" );
      parser.template setFixSizeVector<IntVecChart>("InputMesh.NumDofVec", _NumDofVec, true );
      parser.template setFixSizeVector<PointType>("InputMesh.LengthVec", _LengthVec, false );
  }
  
  void saveToParser( ParameterParserType & parser, const string /*fileNameMesh*/ ) const {
//     parser.set ( "InputMesh.dimChartDomain", static_cast<int>( dimChartDomain ) );
//       parser.set ( "InputMesh.FEType", "Rectangle" );
//       parser.set ( "InputMesh.MeshType", "File" );
//       parser.set ( "InputMesh.fileName", fileNameMesh ); 
     this->saveToParser( parser );
  }
  
};


 



//Rectangle (0,0,0) to (lx,ly,ly)
// DOFs in x: Nx, y:Ny, z:Nz
// meshsize: (hx,hy,hz) = (lx/(Nx-1), ly/(Ny-1), lz/(Nz-1))
template< typename DataTypeContainer = Quoc3DDataTypeContainer,
          typename QuocElementType = QuocElement3D<DataTypeContainer>, 
          typename QuocBoundaryElementType = QuocBoundaryElement3D<DataTypeContainer> >
class QuocMesh3D {
public:
  typedef QuocElementType ElementType;
  typedef typename ElementType::ElementNodeIndicesType      ElementNodeIndicesType;
  typedef QuocBoundaryElementType                           BoundaryElementType;
  typedef DataTypeContainer DTContainer;
  typedef typename DataTypeContainer::RealType              RealType;
  typedef typename DataTypeContainer::PointType             PointType;
  typedef typename DataTypeContainer::VectorType            VectorType;
  typedef typename DataTypeContainer::MaskType              MaskType;
  typedef typename DataTypeContainer::IntVecChart           IntVecChart;
  typedef typename DataTypeContainer::RealVecChart          RealVecChart;
  typedef pesopt::BoostParser ParameterParserType;
  
  static const int dimChartDomain = 3;
  static const int _VTKCELLTYPE = 11; //VTK_VOXEL
  static const VTKDATATYPE _VTKDATATYPEUNDEFORMED = VTKSTRUCTUREDPOINTS;
  static const VTKDATATYPE _VTKDATATYPEDEFORMED =  VTKUNSTRUCTUREDGRID;
  static const PESOPT_FEMeshType _FEMeshType = QUOCMESH;
  
protected:
  int _Nx, _Ny, _Nz;
  IntVecChart _NumDofVec;
  RealType _lx, _ly, _lz;
  PointType _LengthVec;
  RealType _hx, _hy, _hz;
  std::vector< PointType > _vertexIterator;
  std::vector< ElementType > _elementIterator;
  std::vector< BoundaryElementType > _boundaryElementIterator;
  mutable std::vector< std::vector<int>  > _vertexIteratorAdjacentElements;
  mutable std::vector< std::vector<int>  > _vertexIteratorAdjacentElementsLocalIdx;
public:
  MaskType _boundaryFront, _boundaryBack, _boundaryTop, _boundaryBottom, _boundaryLeft, _boundaryRight;
  MaskType _boundaryMask;
  MaskType _boundaryPeriodic;
  std::vector<int> _periodicIdentificationIndices;
  
public:
  //! Create empty QuocMesh
//   QuocMesh3D ( ) : _vertexIterator(), _elementIterator()  { }
  
  QuocMesh3D ( const int Nx, const int Ny, const int Nz, const RealType lx, const RealType ly, const RealType lz ) : 
    _Nx (Nx), _Ny(Ny), _Nz(Nz),
    _lx( lx ), _ly( ly ), _lz( lz ),
    _hx( _lx / static_cast<RealType>(_Nx - 1 ) ), _hy( _ly / static_cast<RealType>(_Ny - 1 ) ), _hz( _lz / static_cast<RealType>(_Nz - 1 ) )
    {
        
        _NumDofVec[0] = _Nx; _NumDofVec[1] = _Ny; _NumDofVec[2] = _Nz;
        _LengthVec[0] = _lx; _LengthVec[1] = _ly; _LengthVec[2] = _lz;
        
        const int numGlobalDofs = _Nx * _Ny * _Nz;
        _boundaryFront.resize( numGlobalDofs, false ); _boundaryBack.resize( numGlobalDofs, false );
        _boundaryTop.resize( numGlobalDofs, false ); _boundaryBottom.resize( numGlobalDofs, false );
        _boundaryLeft.resize( numGlobalDofs, false ); _boundaryRight.resize( numGlobalDofs, false );
        _boundaryMask.resize( numGlobalDofs, false );
        _boundaryPeriodic.resize( numGlobalDofs, false );
        
        for( int zIdx=0; zIdx<_Nz; ++zIdx )
            for( int yIdx=0; yIdx<_Ny; ++yIdx )
                for( int xIdx=0; xIdx<_Nx; ++xIdx ){
                    PointType vertex ( _lx * xIdx / static_cast<RealType>( _Nx - 1 ), 
                                       _ly * yIdx / static_cast<RealType>( _Ny - 1 ),
                                       _lz * zIdx / static_cast<RealType>( _Nz - 1 ) 
                                     );
                    int nodeIdx = this->pushBackVertex( vertex );
                    if( xIdx == 0 ) _boundaryLeft[nodeIdx] = true;
                    if( xIdx == _Nx - 1 ) _boundaryRight[nodeIdx] = true;
                    if( yIdx == 0 ) _boundaryFront[nodeIdx] = true;
                    if( yIdx == _Ny - 1 ) _boundaryBack[nodeIdx] = true;
                    if( zIdx == 0 ) _boundaryBottom[nodeIdx] = true;
                    if( zIdx == _Nz - 1 ) _boundaryTop[nodeIdx] = true;
                    
                    if( (xIdx == 0) || ( xIdx == _Nx - 1 ) || ( yIdx == 0 ) || ( yIdx == _Ny - 1 )  || ( zIdx == 0 ) ||  ( zIdx == _Nz - 1 ) ) _boundaryMask[nodeIdx] = true;
                    
                    if( ( xIdx == _Nx - 1 ) || ( yIdx == _Ny - 1 ) || ( zIdx == _Nz - 1 ) ) _boundaryPeriodic[nodeIdx] = true;
          }
          
      for( int nodeIdx = 0; nodeIdx < numGlobalDofs; ++nodeIdx ){
          if( (!_boundaryRight[nodeIdx]) && (!_boundaryBack[nodeIdx]) && (!_boundaryTop[nodeIdx]) ){
              ElementNodeIndicesType indices;
              indices << nodeIdx, nodeIdx + 1, nodeIdx + _Nx, nodeIdx + _Nx + 1,
                         nodeIdx + _Nx * _Ny, nodeIdx + _Nx * _Ny + 1, nodeIdx + _Nx * _Ny + _Nx, nodeIdx + _Nx * _Ny + _Nx + 1;
              const int numNewElement = this->pushBackElement( indices );   
              
              if( this->getVertex(nodeIdx)[0] == 0 ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::LEFT, _hy * _hz );
              if( this->getVertex(nodeIdx + 1)[0] == _lx ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::RIGHT, _hy * _hz );
              if( this->getVertex(nodeIdx)[1] == 0 ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::FRONT, _hx * _hz );
              if( this->getVertex(nodeIdx + _Nx)[1] == _ly ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::BACK, _hx * _hz );
              if( this->getVertex(nodeIdx)[2] == 0 ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::BOTTOM, _hx * _hy );
              if( this->getVertex(nodeIdx + _Nx * _Ny)[2] == _lz ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::TOP, _hx * _hy );
          }
      }

      //Periodic indices
      _periodicIdentificationIndices.resize( numGlobalDofs, -1 );
      for( int zIdx=0; zIdx<_Nz; ++zIdx )
          for( int yIdx=0; yIdx<_Ny; ++yIdx )
             for( int xIdx=0; xIdx<_Nx; ++xIdx ){
                int periodicXIdx = xIdx, periodicYIdx = yIdx, periodicZIdx = zIdx;
                if( xIdx == _Nx - 1 ) periodicXIdx = 0;
                if( yIdx == _Ny - 1 ) periodicYIdx = 0;
                if( zIdx == _Nz - 1) periodicZIdx = 0;
                _periodicIdentificationIndices[this->getGlobalNodeIndex(xIdx,yIdx,zIdx)] = this->getGlobalNodeIndex(periodicXIdx,periodicYIdx,periodicZIdx);
             }
             
       this->updateAllAdjacentElementsForNodes();
  }
  
  
  QuocMesh3D ( const IntVecChart &NumDofVec, const PointType &LenghtVec ) : 
  QuocMesh3D ( NumDofVec[0], NumDofVec[1], NumDofVec[2], LenghtVec[0], LenghtVec[1], LenghtVec[2] ) {}  
  
  QuocMesh3D ( const string& fileName ){
      throw std::logic_error("constructor of QuocMesh1D is not defined for filename!");
      this->loadFromVTK( fileName );
  }
  
  
protected:
  void updateAllAdjacentElementsForNodes () const {
    _vertexIteratorAdjacentElements.clear();
    _vertexIteratorAdjacentElementsLocalIdx.clear();
    _vertexIteratorAdjacentElements.resize( this->getNumVertices() );
    _vertexIteratorAdjacentElementsLocalIdx.resize( this->getNumVertices() );
    for ( int elementIndex = 0; elementIndex < getNumElements(); ++elementIndex  )
        for(int localNodeIdx=0; localNodeIdx<2; ++localNodeIdx){
           const int globalNodeIdx = _elementIterator[elementIndex].getGlobalNodeIdx(localNodeIdx);
           _vertexIteratorAdjacentElements[globalNodeIdx].push_back(_elementIterator[elementIndex].getGlobalElementIdx());
           _vertexIteratorAdjacentElementsLocalIdx[globalNodeIdx].push_back(localNodeIdx);
        }
  }

public:
  void getCommonElements( const int nodeIdx1, std::vector<int> &commonElements ) const{
      for(int i=0;i<_vertexIteratorAdjacentElements[nodeIdx1].size(); ++i ) commonElements.push_back(_vertexIteratorAdjacentElements[nodeIdx1][i]); 
  }
    
public:
  int getNumVertices ( ) const { return ( _vertexIterator.size() );}
  int getNumElements ( ) const { return ( static_cast<int> ( _elementIterator.size() ) );}
  int getNumBoundaryElements ( ) const { return ( static_cast<int> ( _boundaryElementIterator.size() ) );}
  
  //! insert new vertex and return global index
  int pushBackVertex ( const PointType newVertex ) {
    _vertexIterator.push_back ( newVertex );
    return getNumVertices() - 1;
  }
  //! insert new element and return global index
  int pushBackElement ( const ElementNodeIndicesType nodeIdx ) {
    int globalIdx = getNumElements();
    _elementIterator.push_back ( ElementType( globalIdx, nodeIdx, _vertexIterator ) );
    return globalIdx;
  }
  
  int pushBackBoundaryElement ( const int elementIndex, QuocBoundaryType boundaryType, const RealType h ) {
    int globalBoundaryIdx = getNumBoundaryElements();
    _boundaryElementIterator.push_back ( BoundaryElementType( globalBoundaryIdx, this->getElement(elementIndex), boundaryType, h) );
    return globalBoundaryIdx;
  }

  const PointType& getVertex ( const int num ) const { return _vertexIterator[num];}
  void setVertex ( const int num, const PointType Arg ) { _vertexIterator[num] = Arg;}

  int getGlobalNodeIndex ( const int xIdx, const int yIdx, const int zIdx ) const { return ( xIdx + yIdx * _Nx + zIdx * _Nx * _Ny );}
  int getGlobalNodeIndex ( const IntVecChart &indices ) const { return this->getGlobalNodeIndex( indices[0], indices[1], indices[2] ); }
  int getGlobalElementIndex ( const int xFac, const int yFac, const int zFac ) const { return ( xFac + yFac * (_Nx - 1 ) + zFac * (_Nx - 1 ) * ( _Ny - 1 ) );}  

  const ElementType& getElement ( const int num ) const {return _elementIterator[num];}
  void setElement ( const int num, const ElementType Arg ) { _elementIterator[num] = Arg;}

  const QuocBoundaryElementType& getBoundaryElement ( const int num ) const {return _boundaryElementIterator[num];}
  
  int getElementNodeIdx  ( const int num, const int localNode ) const { return _elementIterator[num].getGlobalNodeIdx(localNode);}
  void setElementNodeIdx ( const int num, const int localNode, const int newIdx ) {_elementIterator[num].setGlobalNodeIdx(localNode, newIdx );}
  
  const RealType getMeshSize( const int direction ) const{ 
    if( direction == 0 ) return _hx;
    if( direction == 1 ) return _hy;
    if( direction == 2 ) return _hz;
    throw std::invalid_argument( pesopt::strprintf ( "Wrong direction in getMeshSize in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  const RealType getInterfaceWith( ) const{ return std::cbrt(_hx * _hy * _hz);}
  const RealType getWidth( const int direction ) const {
     switch( direction ){
        case 0 : return _lx; break;
        case 1 : return _ly; break;
        case 2 : return _lz; break;
        default: throw std::invalid_argument( pesopt::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
    }
  }
  
  const int getNumDofs( const int direction ) const {
   switch( direction ){
        case 0 : return _Nx; break;
        case 1 : return _Ny; break;
        case 2 : return _Nz; break;
        default: throw std::invalid_argument( pesopt::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
    }
  }
  
  void getNumDofVec( IntVecChart &vec ) const { vec = _NumDofVec;}
  void getLenghtVec( PointType &vec ) const { vec = _LengthVec;}
  
  int getNumNodesOfElement(const int elIdx ) const { return 8; }
  
//   IntVecChart getElementDirectionIndices(int elementIdx) const{
//       int zIdx = elementIdx / ((_Nx - 1) * (_Ny - 1));
//       int rest = elementIdx % ((_Nx - 1) * (_Ny - 1));
//       int yIdx = rest / (_Nx - 1);
//       int xIdx = rest % (_Nx - 1); 
//       IntVecChart indices (xIdx,yIdx,zIdx);
//       return indices;
//   }
      
  IntVecChart getNodeIndices (int nodeIdx) const{
      int zIdx = nodeIdx / (_Nx * _Ny);
      int rest = nodeIdx % (_Nx * _Ny);
      int yIdx = rest / _Nx;;
      int xIdx = rest % _Nx;
      IntVecChart indices (xIdx, yIdx, zIdx); 
      return indices; 
  }

//   void findClosestVertex( const PointType &coords, int &vertexIndex ) const{
//       vertexIndex = 0;
//       RealType distSqr = ( coords - _vertexIterator[0] ).squaredNorm();
//       for( int i=0; i<this->getNumVertices(); ++i ){
//           RealType tmp = ( coords - _vertexIterator[i] ).squaredNorm();
//           if( tmp < distSqr ){
//               vertexIndex = i; distSqr = tmp;
//           }
//       }
//   }
  
  
  
protected: 
    
   void loadFromVTK( const string& filename ) {
    VTKSaver vtkReader;
    int numVertices;
     
//     std::vector<PointType> pointVec;
//     vtkReader.getPoints<PointType>( filename, pointVec );
//     numVertices = pointVec.size();    
//     _vertexIterator.resize( numVertices );
//     for( int i=0; i<numVertices; ++i ) this->setVertex ( i, pointVec[i] );
// 
//     std::vector<ElementNodeIndicesType> triangleVec;
//     vtkReader.getCells<ElementNodeIndicesType>( filename, triangleVec );
//     numTriangles = triangleVec.size();    
//     _elementIterator.resize( numTriangles );
//     for( int i = 0; i<numTriangles; ++i ) this->setElement ( i, triangleVec[i] ); 
  } 
  
public:
  
    
  const MaskType &getBoundaryMask( ) const { return _boundaryMask; }
    
  void saveToParser( ParameterParserType & parser ) const {
      parser.set ( "InputMesh.dimChartDomain", static_cast<int>( dimChartDomain ) );
      parser.set ( "InputMesh.FEType", "Cuboid" );
      parser.set ( "InputMesh.MeshType", "Cuboid" );
      parser.template setFixSizeVector<IntVecChart>("InputMesh.NumDofVec", _NumDofVec, true );
      parser.template setFixSizeVector<PointType>("InputMesh.LengthVec", _LengthVec, false );
  }
  
  void saveToParser( ParameterParserType & parser, const string fileNameMesh ) const {
/*      parser.set ( "InputMesh.dimChartDomain", static_cast<int>( dimChartDomain ) );
      parser.set ( "InputMesh.FEType", "Cuboid" );
      parser.set ( "InputMesh.MeshType", "File" );
      parser.set ( "InputMesh.fileName", fileNameMesh );*/ 
      this->saveToParser( parser );
  }
  
};


#endif
