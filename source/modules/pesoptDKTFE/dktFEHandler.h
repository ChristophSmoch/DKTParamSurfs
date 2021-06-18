#ifndef __DKTFEHANDLER_H
#define __DKTFEHANDLER_H

#include <dktFEBoundary.h>
#include <dktFEFunctionEvaluator.h>
#include <dktFEConfigurators.h>

#ifdef PESOPT_WITH_VTK
 #include <VTKMeshSaver.h>
#endif


template< typename ConfiguratorType >
class ShellHandlerInterface{

public:

  typedef typename ConfiguratorType::RealType       RealType;
  typedef typename ConfiguratorType::InitType       MeshType;
  typedef typename ConfiguratorType::MaskType       MaskType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Point3DType    Point3DType;
  typedef typename ConfiguratorType::VectorType     VectorType;
  typedef typename ConfiguratorType::DTContainer    DataTypeContainer;

  const ConfiguratorType &_conf;
  const MeshType &_mesh;
  const int _numVertices, _numGlobalDofs;
  ShellBoundaryType _shellBdrType;
  MaskType _DirichletMask;
  int _numBoundaryNodes;
  bool _clampedBoundaryCondition;
  mutable VectorType _xA;

public:


 ShellHandlerInterface( const ConfiguratorType &conf, const string chartXAType ) :
  _conf ( conf ), _mesh( conf.getInitializer() ),
  _numVertices ( _mesh.getNumVertices() ),  _numGlobalDofs ( conf.getNumGlobalDofs() ),
  _xA ( 3 * conf.getNumGlobalDofs() )
  {
    cout << "warning: ShellHandlerInterface without boundary mask!!!!!!!!!!!!!!!!!!!" << endl;
    generateChart_xA ( chartXAType );
  }


  ShellHandlerInterface( const ConfiguratorType &conf, const string chartXAType, const ShellBoundaryType &shellBdrType, const bool clampedBdr ) :
  _conf ( conf ),
  _mesh( conf.getInitializer() ),
  _numVertices ( _mesh.getNumVertices() ),
  _numGlobalDofs ( conf.getNumGlobalDofs() ),
  _shellBdrType ( shellBdrType ), _clampedBoundaryCondition ( clampedBdr ),
  _xA ( 3 * conf.getNumGlobalDofs() )
  {
    generateChart_xA ( chartXAType );
    generateDirichletBoundaryMask( _DirichletMask, _numBoundaryNodes );
  }


  ShellHandlerInterface( const ConfiguratorType &conf, const string chartXAType, const ShellBoundaryType &shellBdrType, const MaskType &DirichletMask, const bool clampedBdr ) :
  _conf ( conf ),
  _mesh( conf.getInitializer() ),
  _numVertices ( _mesh.getNumVertices() ),
  _numGlobalDofs ( conf.getNumGlobalDofs() ),
  _shellBdrType ( shellBdrType ), _DirichletMask( DirichletMask ), _clampedBoundaryCondition ( clampedBdr ),
  _xA ( 3 * conf.getNumGlobalDofs() )
  {
    generateChart_xA ( chartXAType );
    _numBoundaryNodes = 0; for( int i=0; i<_mesh.getNumVertices(); ++i) if(_DirichletMask[i]) _numBoundaryNodes++;
  }


private:
    //from BeamHandler
      void constructLeftRightClampedInit_PWLin( VectorType &K, const RealType alpha ) const{
          int sizeFirstPart = K.size() / 4;
          int sizeSecondPart =  3 * K.size() / 4;
          for(int i=0; i < sizeFirstPart; i++)            K[i] = 4. * alpha * static_cast<RealType> ( i ) / static_cast<RealType> ( K.size() - 1. );
          for(int i=sizeFirstPart; i<sizeSecondPart; i++) K[i] = -4. * alpha * static_cast<RealType> ( i ) / static_cast<RealType> ( K.size() - 1. ) + 2. * alpha;
          for(int i=sizeSecondPart; i<K.size(); i++)      K[i] = 4. * alpha * static_cast<RealType> ( i ) / static_cast<RealType> ( K.size() - 1. ) - 4. * alpha;
      }

     void computeCurveFromPhase( const VectorType &K, VectorType &CurveReal, VectorType &CurveImaginary ) const {
        CurveReal[0] = 0.0;
        CurveImaginary[0] = 0.0;
        RealType tau = 1.0 / (K.size() - 1.0 );
        for(int i=1; i< K.size(); i++){
            CurveReal[i] = CurveReal[i-1] + tau * cos( K[i-1] );
            CurveImaginary[i] = CurveImaginary[i-1] + tau * sin( K[i-1] );
        }
      }


public:

  void generateChart_xA( const string chartXAType) const{

      //identity
      if( chartXAType == "id" ){
        for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
            const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
            for( int comp=0; comp<3; ++comp )  _xA[ nodeIdx + _numGlobalDofs * comp ] = coords[comp];
            if( ConfiguratorType::_ShellFEType == C1Dofs ){
                TangentVecType firstTangentVecAtNode  ( _mesh.getTangentVec1 (nodeIdx) );
                TangentVecType secondTangentVecAtNode ( _mesh.getTangentVec2 (nodeIdx) );
                for( int comp=0; comp<3; ++comp ){
                  _xA[ nodeIdx +     _numVertices + _numGlobalDofs * comp ] = firstTangentVecAtNode  [comp];
                  _xA[ nodeIdx + 2 * _numVertices + _numGlobalDofs * comp ] = secondTangentVecAtNode [comp];
                }
            }
        }
      }

    //stereographic projection
    if( chartXAType == "CircToSphere" ){
        for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
            const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
            const RealType normSqr = coords[0] * coords[0] + coords[1] * coords[1];
            Point3DType coordsOnSphere;
            coordsOnSphere(0) = 2. * coords(0) / (normSqr + 1.);
            coordsOnSphere(1) = 2. * coords(1) / (normSqr + 1.);
            coordsOnSphere(2) = (1. - normSqr) / (normSqr + 1.);
            for( int comp=0; comp<3; ++comp ){
                _xA[ nodeIdx + _numGlobalDofs * comp ] = coordsOnSphere[comp];
            }
            if( ConfiguratorType::_ShellFEType == C1Dofs ){
                RealType fac = 2. / pesopt::Sqr( normSqr + 1. );
                TangentVecType firstTangentVecAtNode;
                firstTangentVecAtNode(0) = fac * ( - coords(0) * coords(0) + coords(1) * coords(1) + 1. );
                firstTangentVecAtNode(1) = fac * -2. * coords(0) * coords(1);
                firstTangentVecAtNode(2) = -2. * coords(0);
                TangentVecType secondTangentVecAtNode;
                secondTangentVecAtNode(0) = fac * -2. * coords(0) * coords(1);
                secondTangentVecAtNode(1) = fac * ( coords(0) * coords(0) - coords(1) * coords(1) + 1. );
                secondTangentVecAtNode(2) = -2. * coords(1);
                for( int comp=0; comp<3; ++comp ){
                  _xA[ nodeIdx +     _numVertices + _numGlobalDofs * comp ] = firstTangentVecAtNode  [comp];
                  _xA[ nodeIdx + 2 * _numVertices + _numGlobalDofs * comp ] = secondTangentVecAtNode [comp];
                }
            }
        }
      }


      if( chartXAType == "PlateToCylinder" ){
          const RealType pi = 4 * atan ( 1.0 );
          const RealType radius = 1. / pi;
          for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
            const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
            Point3DType coordsOnCylinder;
            coordsOnCylinder(0) = radius * sin( coords(0) * pi );
            coordsOnCylinder(1) = coords(1);
            coordsOnCylinder(2) = radius * cos( coords(0) * pi );
            for( int comp=0; comp<3; ++comp )_xA[ nodeIdx + _numGlobalDofs * comp ] = coordsOnCylinder[comp];
            if( ConfiguratorType::_ShellFEType == C1Dofs ){
                TangentVecType firstTangentVecAtNode;
                firstTangentVecAtNode(0) = radius * pi * cos( coords(0) * pi );
                firstTangentVecAtNode(1) = 0.;
                firstTangentVecAtNode(2) = - radius * pi * sin( coords(0) * pi );
                TangentVecType secondTangentVecAtNode; secondTangentVecAtNode(0) = 0.; secondTangentVecAtNode(1) = 1.; secondTangentVecAtNode(2) = 0.;
                for( int comp=0; comp<3; ++comp ){
                  _xA[ nodeIdx +     _numVertices + _numGlobalDofs * comp ] = firstTangentVecAtNode  [comp];
                  _xA[ nodeIdx + 2 * _numVertices + _numGlobalDofs * comp ] = secondTangentVecAtNode [comp];
                }
            }
        }
      }


      // TODO Christoph: Karte zur undeformierten Schale (nicht isometrisch)
      if( chartXAType == "PlateToNonIsometricTest" ){
          const RealType pi = 4 * atan ( 1.0 );
          const RealType radius = 1. / pi;
          for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
            const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
            Point3DType coordsOnNonIsometry;
            coordsOnNonIsometry(0) = coords(0);
            coordsOnNonIsometry(1) = coords(1);
            coordsOnNonIsometry(2) = 0.5 * ( (coords(0) - 0.5)*(coords(0) - 0.5) - (coords(1) - 0.5)*(coords(1) - 0.5));
            // coordsOnNonIsometry(2) = 0.5 * ( (coords(0) - 0.5)*(coords(0) - 0.5) + (coords(1) - 0.5)*(coords(1) - 0.5));
            for( int comp=0; comp<3; ++comp )_xA[ nodeIdx + _numGlobalDofs * comp ] = coordsOnNonIsometry[comp];
            if( ConfiguratorType::_ShellFEType == C1Dofs ){
                TangentVecType firstTangentVecAtNode;
                firstTangentVecAtNode(0) = 1.;
                firstTangentVecAtNode(1) = 0.;
                firstTangentVecAtNode(2) = coords(0) - 0.5;
                TangentVecType secondTangentVecAtNode; secondTangentVecAtNode(0) = 0.; secondTangentVecAtNode(1) = 1.; secondTangentVecAtNode(2) = - coords(1) + 0.5;
                // TangentVecType secondTangentVecAtNode; secondTangentVecAtNode(0) = 0.; secondTangentVecAtNode(1) = 1.; secondTangentVecAtNode(2) =  coords(1) - 0.5;
                for( int comp=0; comp<3; ++comp ){
                  _xA[ nodeIdx +     _numVertices + _numGlobalDofs * comp ] = firstTangentVecAtNode  [comp];
                  _xA[ nodeIdx + 2 * _numVertices + _numGlobalDofs * comp ] = secondTangentVecAtNode [comp];
                }
            }
        }
      }

      // Half Cone Christoph
      if( chartXAType == "PlateToHalfCone" ){
          const RealType pi = 4 * atan ( 1.0 );
          for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
            const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
            const RealType radius = 1. / pi * (1 - coords(1)/2);
            Point3DType coordsOnHalfCone;
            coordsOnHalfCone(0) = radius * cos( coords(0) * pi );
            coordsOnHalfCone(1) = coords(1);
            coordsOnHalfCone(2) = radius * sin( coords(0) * pi );
            for( int comp=0; comp<3; ++comp )_xA[ nodeIdx + _numGlobalDofs * comp ] = coordsOnHalfCone[comp];
            if( ConfiguratorType::_ShellFEType == C1Dofs ){
                TangentVecType firstTangentVecAtNode;
                firstTangentVecAtNode(0) = -1. * pi * radius * sin(coords(0) * pi);
                firstTangentVecAtNode(1) = 0.;
                firstTangentVecAtNode(2) = pi * radius * cos(coords(0) * pi);
                TangentVecType secondTangentVecAtNode;
                secondTangentVecAtNode(0) = -1. * 1. / 2. * 1. / pi * cos( coords(0) * pi );
                secondTangentVecAtNode(1) = 1.;
                secondTangentVecAtNode(2) = -1. * 1. / 2. * 1. / pi * sin( coords(0) * pi );
                // TangentVecType secondTangentVecAtNode; secondTangentVecAtNode(0) = 0.; secondTangentVecAtNode(1) = 1.; secondTangentVecAtNode(2) =  coords(1) - 0.5;
                for( int comp=0; comp<3; ++comp ){
                  _xA[ nodeIdx +     _numVertices + _numGlobalDofs * comp ] = firstTangentVecAtNode  [comp];
                  _xA[ nodeIdx + 2 * _numVertices + _numGlobalDofs * comp ] = secondTangentVecAtNode [comp];
                }
            }
        }
      }
      // Half Cone Christoph
      if( chartXAType == "PlateToWobbly" ){
          const RealType pi = 4 * atan ( 1.0 );
          for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
            const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
            const RealType radius = 1. / pi ;
            const RealType oszfac = 10. ;
            Point3DType coordsOnHalfCone;
            coordsOnHalfCone(0) = coords(0);
            coordsOnHalfCone(1) = coords(1);
            coordsOnHalfCone(2) =  1/oszfac -1. * radius * sin( coords(1) * pi * oszfac ) * cos( coords(0) * pi * oszfac);
            for( int comp=0; comp<3; ++comp )_xA[ nodeIdx + _numGlobalDofs * comp ] = coordsOnHalfCone[comp];
            if( ConfiguratorType::_ShellFEType == C1Dofs ){
                TangentVecType firstTangentVecAtNode;
                firstTangentVecAtNode(0) = 1.;
                firstTangentVecAtNode(1) = 0.;
                firstTangentVecAtNode(2) =   pi * radius * sin( coords(1) * pi* oszfac ) * sin( coords(0) * pi* oszfac);
                TangentVecType secondTangentVecAtNode;
                secondTangentVecAtNode(0) = 0.;
                secondTangentVecAtNode(1) = 1.;
                secondTangentVecAtNode(2) =  -1. * pi * radius * cos( coords(1) * pi* oszfac ) * cos( coords(0) * pi)* oszfac;
                // TangentVecType secondTangentVecAtNode; secondTangentVecAtNode(0) = 0.; secondTangentVecAtNode(1) = 1.; secondTangentVecAtNode(2) =  coords(1) - 0.5;
                for( int comp=0; comp<3; ++comp ){
                  _xA[ nodeIdx +     _numVertices + _numGlobalDofs * comp ] = firstTangentVecAtNode  [comp];
                  _xA[ nodeIdx + 2 * _numVertices + _numGlobalDofs * comp ] = secondTangentVecAtNode [comp];
                }
            }
        }
      }


      // chart(x,y) = (x,y, 1/4(2x-1)^2 - 1/4(2y-1)^2)
      // D chart(x,y) =
      // ( 1             0
      //   0             1
      //   (2x-1) -(2y-1)    )
      if( chartXAType == "PlateToQuadShape" ){
          for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
            const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
            Point3DType coordsOnShell;
            coordsOnShell(0) = coords(0);
            coordsOnShell(1) = coords(1);
            coordsOnShell(2) = 0.25 * ( pesopt::Sqr( 2. * coords(0) - 1. ) -  pesopt::Sqr( 2. * coords(1) - 1. ) );
            for( int comp=0; comp<3; ++comp )_xA[ nodeIdx + _numGlobalDofs * comp ] = coordsOnShell[comp];
            if( ConfiguratorType::_ShellFEType == C1Dofs ){
                TangentVecType firstTangentVecAtNode; firstTangentVecAtNode(0) = 1.; firstTangentVecAtNode(1) = 0.; firstTangentVecAtNode(2) = ( 2. * coords(0) - 1. );
                TangentVecType secondTangentVecAtNode; secondTangentVecAtNode(0) = 0.; secondTangentVecAtNode(1) = 1.; secondTangentVecAtNode(2) = - ( 2. * coords(1) - 1. );
                for( int comp=0; comp<3; ++comp ){
                  _xA[ nodeIdx +     _numVertices + _numGlobalDofs * comp ] = firstTangentVecAtNode  [comp];
                  _xA[ nodeIdx + 2 * _numVertices + _numGlobalDofs * comp ] = secondTangentVecAtNode [comp];
                }
            }
        }
      }




      //push plate together from left and right:
      // 1) construct phase K : [0,1] \to [-pi,pi]
      // 2) compute curve gamma(t) = \int_0^t e^(i K(t) ) dt
      // 3) xA(x,y,z=0) = ( gamma_1(x), y, gamma_2(x) )

      if( chartXAType == "PlateLeftRight" ){

        const RealType alpha = 2.0; //TODO optional
        VectorType phase ( 1000 ); //TODO require that mesh has N_x \times N_y nodes
        this->constructLeftRightClampedInit_PWLin( phase, alpha );
        VectorType gammaReal ( phase.size() ), gammaIm ( phase.size() );
        this->computeCurveFromPhase( phase, gammaReal, gammaIm );

        for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
            const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
            const RealType t = coords[0];
            const int index = static_cast<int>( t * (phase.size()-1) );
            Point3DType newCoords ( gammaReal[index], coords[1], gammaIm[index] );
            for( int comp=0; comp<3; ++comp )  _xA[ nodeIdx + _numGlobalDofs * comp ] = newCoords[comp];
            if( ConfiguratorType::_ShellFEType == C1Dofs ){
                TangentVecType firstTangentVecAtNode  ( cos( phase[index] ), 0., sin( phase[index] ) );
                TangentVecType secondTangentVecAtNode ( 0., 1., 0. );
                for( int comp=0; comp<3; ++comp ){
                  _xA[ nodeIdx +     _numVertices + _numGlobalDofs * comp ] = firstTangentVecAtNode  [comp];
                  _xA[ nodeIdx + 2 * _numVertices + _numGlobalDofs * comp ] = secondTangentVecAtNode [comp];
                }
            }
        }
      }

  }

  const VectorType &getChartToUndeformedShell (  ) const { return _xA;}

  void generateDirichletBoundaryMask ( MaskType & mask, int & numBoundaryNodes ) const{
    cout << "generateDirichletBoundaryMask with boundary type = " << _shellBdrType << endl;
    mask.resize( _numGlobalDofs, false );
    const bool ShellFETypeC1Dofs =  (ConfiguratorType::_ShellFEType == C1Dofs ) ? true : false;
    generateDirichletBoundaryMaskUponShellBoundaryType<MeshType, ShellFETypeC1Dofs> ( _shellBdrType, _mesh, mask, numBoundaryNodes, _clampedBoundaryCondition );
    cout << "generated with numBoundaryNodes = " << numBoundaryNodes << endl;
  }

  void setDirichletMask ( const MaskType &mask ) {
      _DirichletMask = mask;
      _numBoundaryNodes = 0; for( int i=0; i<_mesh.getNumVertices(); ++i) if(_DirichletMask[i]) _numBoundaryNodes++;
  }
  const MaskType & getDirichletMask ( ) const { return _DirichletMask;}
  const int & getNumBoundaryNodes() const { return _numBoundaryNodes; }

  bool hasClampedBoundary ( ) const { return _clampedBoundaryCondition; }

};



template< typename ConfiguratorType >
class ShellDeformationGenerator{

protected:

  typedef typename ConfiguratorType::RealType       RealType;
  typedef typename ConfiguratorType::InitType       MeshType;
  typedef typename ConfiguratorType::MaskType       MaskType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Point3DType    Point3DType;
  typedef typename ConfiguratorType::VectorType     VectorType;
  typedef typename ConfiguratorType::DTContainer    DataTypeContainer;

  const ConfiguratorType &_conf;
  const MeshType &_mesh;
  const int _numVertices, _numGlobalDofs;
  const MaskType &_DirichletMask;

public:
    ShellDeformationGenerator( const ConfiguratorType &conf, const MaskType &DirichletMask ) :
    _conf( conf ), _mesh( conf.getInitializer() ),
  _numVertices ( _mesh.getNumVertices() ), _numGlobalDofs( conf.getNumGlobalDofs() ),
  _DirichletMask( DirichletMask ) {}



//! construct displacements
  //===============================================================================================================================================================================
//   // generates identity
//   void generateIdentityDeformation ( pesopt::MultiVector<RealType> &deformation ) const{
//     int numVertices = _mesh.getNumVertices();
//     deformation.reallocate( 3, 3 * numVertices );
//     for( int i = 0; i < numVertices; i++ ){
//         pesopt::Vec3<RealType> coords = _mesh.getVertex( i );
//         (deformation[0]).set( i, coords[0] );
//         (deformation[1]).set( i, coords[1] );
//         (deformation[2]).set( i, coords[0] );
//         (deformation[0]).set( i + numVertices  , 1.0 );
//         (deformation[0]).set( i + 2*numVertices, 0.0 );
//         (deformation[1]).set( i + numVertices  , 0.0 );
//         (deformation[1]).set( i + 2*numVertices, 1.0 );
//         (deformation[2]).set( i + numVertices  , 0.0 );
//         (deformation[2]).set( i + 2*numVertices, 0.0 );
//     }
//   }

//===============================================================================================================================================================================
// constructDisplacementForm1DCase

  // linear interpolation between init and final slope
  RealType computeKForSimpleCurve( const RealType coords_x, const RealType initSlope, const RealType finalSlope ) const{
    RealType init = atan( initSlope );
    RealType fin = atan ( finalSlope );
    return ( 1.0 - coords_x ) * init + coords_x * fin;
  }

  // computes int_0^x exp( i (K(t) + K(0) ) d t
  void computeCurveFromPhaseForSimpleCurve( const RealType & coords_x, RealType & gamma_x, RealType & gamma_z, const RealType initSlope, const RealType finalSlope ) const{

    gamma_x = 0., gamma_z = 0.;

    int numSteps = 50; //TODO numQuadPoints
    RealType tau = coords_x / static_cast<RealType> (numSteps);

    for(int i=0; i< numSteps; ++i){
      RealType coord_i = (static_cast<RealType> (i) + 0.5 ) * tau;
      RealType K_i = computeKForSimpleCurve ( coord_i, initSlope, finalSlope );
      gamma_x += tau * cos( K_i );
      gamma_z += tau * sin( K_i );
    }
  }

  void constructSimpleDisplacement( VectorType &disp, const RealType initSlope, const RealType finalSlope ) const{
    disp.setZero();
    const int numVertices = _mesh.getNumVertices();
    const int numGlobalDofs = 3 * numVertices;

    // construct displacement w = \phi - (x,y,0) with
    // w(x,y) = ( \int_0^x cos(K(t) + K_0) dt - x, 0, \int_0^x sin(K(t) + K_0) dt ) )
    // D w(x,y) = ( cos(K(x) + K_0) - 1       0
    //                   0                    0
    //              sin(K(x) + K_0)           0
    for( int nodeIdx = 0; nodeIdx < numVertices; ++nodeIdx ){
        const Point3DType& coords = _mesh.getVertex( nodeIdx );
        RealType K_x = computeKForSimpleCurve( coords[0], initSlope, finalSlope );
        RealType gamma_x, gamma_z;
        computeCurveFromPhaseForSimpleCurve( coords[0], gamma_x, gamma_z, initSlope, finalSlope );
        disp[nodeIdx] = gamma_x - coords[0];
        disp[nodeIdx + numVertices] = cos( K_x ) - 1.0;
        disp[nodeIdx + 2 * 3 * numVertices] = gamma_z;
        disp[nodeIdx + 2 * 3 * numVertices + numVertices] =  sin( K_x );
    }

    for(int nodeIdx=0; nodeIdx< numVertices; ++nodeIdx){
      if( _DirichletMask[nodeIdx] ){
        for( int comp=0; comp<3; ++comp){
          disp[ nodeIdx + comp * numGlobalDofs] = 0.0;
          disp[nodeIdx + comp * numGlobalDofs + 2 * numVertices] = 0.0;
        }
        disp[nodeIdx + numVertices] = cos ( atan(initSlope) ) - 1.0;
        disp[nodeIdx + numVertices + numGlobalDofs] = 0.0;
        disp[nodeIdx + numVertices + 2 * numGlobalDofs] = sin ( atan(initSlope) );
      }
    }

  }


  RealType computeKForTwist( const RealType & coords_x, const RealType initSlope = 0.0 ) const {

    RealType init = atan( initSlope );
    //RealType final = atan ( finalSlope );
    RealType fin = 6 * atan ( 1.0 ); //=3pi/2

    RealType sizeFirstPart = 1. / 8.;
    RealType endTwist = 1. / 2.;

    if( coords_x < sizeFirstPart ) return init;
    if( (coords_x >= sizeFirstPart) && (coords_x <= endTwist) ){
      return  (coords_x - sizeFirstPart) * (fin - init) / (endTwist - sizeFirstPart ) + init;
    }else return fin;

  }

  // computes int_0^x exp( i (K(t) + K(0) ) d t
  void computeCurveFromPhaseForTwist( const RealType coords_x, RealType & gamma_x, RealType & gamma_z, const RealType initSlope ) const{

    gamma_x = 0.0; gamma_z = 0.0;

    int numSteps = 50; //TODO numQuadPoints
    RealType tau = coords_x / static_cast<RealType> (numSteps);

    for(int i=0; i< numSteps; ++i){
      RealType coord_i = ( static_cast<RealType> (i) + 0.5 ) * tau;
      RealType K_i = computeKForTwist ( coord_i, initSlope );
      gamma_x += tau * cos( K_i );
      gamma_z += tau * sin( K_i );
    }
  }

  void construcTwistDisplacement( VectorType &disp, const RealType initSlope ) const{
    disp.setZero();
    const int numVertices = _mesh.getNumVertices();
    const int numGlobalDofs = 3 * numVertices;
    // construct displacement w = \phi - (x,y,0) with
    // w(x,y) = ( \int_0^x cos(K(t) + K_0) dt - x, 0, \int_0^x sin(K(t) + K_0) dt ) )
    // D w(x,y) = ( cos(K(x) + K_0) - 1       0
    //                   0                    0
    //              sin(K(x) + K_0)           0
    for( int nodeIdx = 0; nodeIdx < numVertices; ++nodeIdx ){
        const Point3DType& coords = _mesh.getVertex( nodeIdx );
        RealType K_index_x = computeKForTwist( coords[0], initSlope );
        RealType gamma_x, gamma_z;
        computeCurveFromPhaseForTwist( coords[0], gamma_x, gamma_z, initSlope );
        disp[nodeIdx] = gamma_x - coords[0];
        disp[nodeIdx + numVertices ] = cos( K_index_x ) - 1.0;
        disp[nodeIdx + 2 * numGlobalDofs ] = gamma_z;
        disp[nodeIdx + 2 * numGlobalDofs + numVertices] = sin( K_index_x );
    }

    for(int nodeIdx=0; nodeIdx< numVertices; ++nodeIdx){
      if( _DirichletMask[nodeIdx] ){
        for( int comp=0; comp<3; ++comp){
          disp[ nodeIdx + comp * numGlobalDofs] = 0.0;
          disp[nodeIdx + comp * numGlobalDofs + 2 * numVertices] = 0.0;
        }
        disp[nodeIdx + numVertices] = cos ( atan(initSlope) ) - 1.0;
        disp[nodeIdx + numVertices + numGlobalDofs] = 0.0;
        disp[nodeIdx + numVertices + 2 * numGlobalDofs] = sin ( atan(initSlope) );
      }
    }

  }

  template<typename ParameterParserType>
  void switchDisplacementType( VectorType &disp, const int type, const ParameterParserType &parser ) const{
      switch( type ){
        case 0: disp.setZero(); break;
        case 1 : constructSimpleDisplacement( disp, parser.template get<RealType>( "ConstraintProblem.initSlope") , parser.template get<RealType>( "ConstraintProblem.finalSlope" ) ); break;
        case 2 : construcTwistDisplacement( disp, parser.template get<RealType>( "ConstraintProblem.initSlope" ) ); break;
        case 101:{
            for( int nodeIdx=0; nodeIdx<_DirichletMask.size(); ++nodeIdx ){
              if( _DirichletMask[nodeIdx] ){
                  for( int c=0; c<3; ++c ) disp[nodeIdx + c * _DirichletMask.size() ] = 0.;
              }
              else{
//                   for( int c=0; c<3; ++c ) disp[nodeIdx + c * _DirichletMask.size() ] = -1.e-5;
                  TangentVecType tmp = TangentVecType::Random( 3 );
                  disp[nodeIdx + 0 * _DirichletMask.size() ] = 1.e-5 * tmp(0);
                  disp[nodeIdx + 1 * _DirichletMask.size() ] = 1.e-5 * tmp(1);
                  disp[nodeIdx + 2 * _DirichletMask.size() ] = 1.e-5 * tmp(2);
              }
        }
        }break;
//         case 102:{
//             for( int nodeIdx=0; nodeIdx<_DirichletMask.size(); ++nodeIdx ){
//               if( _DirichletMask[nodeIdx] ){
//                   for( int c=0; c<3; ++c ) disp[nodeIdx + c * _DirichletMask.size() ] = 0.;
//               }
//               else{
//                   for( int c=0; c<3; ++c ) disp[nodeIdx + c * _DirichletMask.size() ] = 1.e-5;
//               }
//             }
//         }break;
        default: throw std::invalid_argument( pesopt::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
  }

//   void switchDisplacementType( VectorType &disp, const int type, const ParameterParserType &parser ) const{
//       switch( type ){
//         case 0: disp.setZero(); break;
//         case 1 : constructSimpleDisplacement( disp, parser.template get<RealType>( "ConstraintProblem.initSlope") , parser.template get<RealType>( "ConstraintProblem.finalSlope" ) ); break;
//         case 2 : construcTwistDisplacement( disp, parser.template get<RealType>( "ConstraintProblem.initSlope" ) ); break;
//         case 101:{
//
//         }break;
//         default: throw std::invalid_argument( pesopt::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
//       }
//   }

};




template< typename ConfiguratorType, DiscreteFunctionCacheType _DiscreteVectorFunctionCacheType>
class ShellHandler :
 public ShellHandlerInterface<ConfiguratorType>,
 public ShellDeformationGenerator<ConfiguratorType>
 {

public:

  typedef typename ConfiguratorType::RealType             RealType;
  typedef typename ConfiguratorType::InitType             MeshType;
  typedef typename ConfiguratorType::MaskType             MaskType;
  typedef typename ConfiguratorType::TangentVecType       TangentVecType;
  typedef typename ConfiguratorType::DomVecType           DomVecType;
  typedef typename ConfiguratorType::Point3DType          Point3DType;
  typedef typename ConfiguratorType::Matrix22             Matrix22;
  typedef typename ConfiguratorType::Matrix32             Matrix32;
  typedef typename ConfiguratorType::Matrix33             Matrix33;
  typedef typename ConfiguratorType::Tensor322Type        Tensor322Type;
  typedef typename ConfiguratorType::VectorType           VectorType;
  typedef typename ConfiguratorType::DTContainer          DataTypeContainer;
  typedef pesopt::BoostParser ParameterParserType;

  const ConfiguratorType &_conf;
  const MeshType &_mesh;
  const int _numVertices, _numGlobalDofs;
  const DiscreteVectorFunctionStorage <ConfiguratorType, _DiscreteVectorFunctionCacheType> *_xACachePtr;

public:

   ShellHandler( const ConfiguratorType &conf, const string chartXAType, const string saveDirectory ) :
    ShellHandlerInterface<ConfiguratorType> ( conf, chartXAType,  static_cast<ShellBoundaryType>( 100 ), 0 ),
    ShellDeformationGenerator<ConfiguratorType> ( conf, this->getDirichletMask() ),
    _conf ( conf ), _mesh( conf.getInitializer() ), _numVertices ( _mesh.getNumVertices() ), _numGlobalDofs ( conf.getNumGlobalDofs() )
  {
    _xACachePtr = new DiscreteVectorFunctionStorage<ConfiguratorType, _DiscreteVectorFunctionCacheType> ( _conf, this->getChartToUndeformedShell(), 3 );
  }


  ShellHandler( const ParameterParserType &Parser, const ConfiguratorType &conf ) :
    ShellHandlerInterface<ConfiguratorType> ( conf, Parser.template get<string>( "InputMesh.chartXAType" ),  static_cast<ShellBoundaryType>( Parser.template get<int>( "InputMesh.ShellType" ) ), Parser.template get<bool> ( "InputMesh.ClampedBoundaryCondition" ) ),
    ShellDeformationGenerator<ConfiguratorType> ( conf, this->getDirichletMask() ),
    _conf ( conf ), _mesh( conf.getInitializer() ), _numVertices ( _mesh.getNumVertices() ), _numGlobalDofs ( conf.getNumGlobalDofs() )
  {
    _xACachePtr = new DiscreteVectorFunctionStorage<ConfiguratorType, _DiscreteVectorFunctionCacheType> ( _conf, this->getChartToUndeformedShell(), 3 );
  }

  ShellHandler( const ParameterParserType &Parser, const ConfiguratorType &conf, const MaskType &DirichletMask ) :
    ShellHandlerInterface<ConfiguratorType> ( conf, Parser.template get<string>( "InputMesh.chartXAType" ), static_cast<ShellBoundaryType>( Parser.template get<int>( "InputMesh.ShellType" ) ), DirichletMask, Parser.template get<bool> ( "InputMesh.ClampedBoundaryCondition" ) ),
    ShellDeformationGenerator<ConfiguratorType> ( conf, DirichletMask ),
    _conf ( conf ),  _mesh( conf.getInitializer() ), _numVertices ( _mesh.getNumVertices() ), _numGlobalDofs ( conf.getNumGlobalDofs() )
  {
    _xACachePtr = new DiscreteVectorFunctionStorage<ConfiguratorType, _DiscreteVectorFunctionCacheType> ( _conf, this->getChartToUndeformedShell(), 3 );
  }

  ~ShellHandler() {  delete _xACachePtr;};

  const DiscreteVectorFunctionStorage <ConfiguratorType, _DiscreteVectorFunctionCacheType> &getChartToUndeformedShell_Cache () const { return *_xACachePtr;}

};





















template<typename ConfiguratorType, ShellFEType MaterialFEType>
class ShellMaterialGenerator{ };


template<typename ConfiguratorType>
class ShellMaterialGenerator<ConfiguratorType,NodalValuedDofs>{

 protected:

  typedef typename ConfiguratorType::RealType       RealType;
  typedef typename ConfiguratorType::InitType       MeshType;
  typedef typename ConfiguratorType::MaskType       MaskType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Point3DType    Point3DType;
  typedef typename ConfiguratorType::VectorType     VectorType;
  typedef typename ConfiguratorType::DTContainer    DataTypeContainer;

  const MeshType &_mesh;
  const int _numVertices;

 public:
   ShellMaterialGenerator( const ConfiguratorType &conf ) :  _mesh( conf.getInitializer() ), _numVertices ( _mesh.getNumVertices() ) {}

  ShellMaterialGenerator( const MeshType &mesh ) : _mesh( mesh ), _numVertices ( _mesh.getNumVertices() ) {}

  void constructConstantMaterial( VectorType &material, RealType materialConstant ) const{
      material.setZero();
      for( int i=0; i<_numVertices; ++i ) material[i] = materialConstant;
   }

  //! constructs material in a box spanned by start and end point
  void constructLayerMaterial( VectorType &material,  const TangentVecType startPoint, const TangentVecType endPoint,
                               const RealType mhard = 1.0, const RealType msoft = -1.0 ) const{
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
        if( (coords[0] >= startPoint[0]) && (coords[1] >= startPoint[1]) && (coords[2] >= startPoint[2])
         && (coords[0] <= endPoint[0]) && (coords[1] <= endPoint[1]) && (coords[2] <= endPoint[2]) ){
            material[nodeIdx] = mhard;
        }
    }
  }

 //! construct ray spanned by start + t rayVector with thickness
 void constructRayMaterial( VectorType &material,  const TangentVecType startPoint, const TangentVecType rayVector, const RealType thickness,
                               const RealType mhard = 1.0, const RealType msoft = -1.0 ) const{
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
        TangentVecType projection;
        RealType projStep = rayVector.dot( coords - startPoint ); projStep /= rayVector.squaredNorm();
        projection = startPoint + projStep * rayVector;
        RealType distance = (projection - coords ).norm();
        if( distance < thickness )  material[nodeIdx] = mhard;
    }
  }

 //! construct ray spanned by start + t rayVector with thickness
 void constructTriangleMaterial( VectorType &material,  const RealType length, const RealType area,
                               const RealType mhard = 1.0, const RealType msoft = -1.0 ) const{
    const RealType height = 2. * area / length; const RealType startY = 0.5 * (1. - length );
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
        bool rightBound = false; if( coords[0] <= height ) rightBound = true;
        bool upperBound = false; RealType up = -0.5 * length / height * coords[0] + ( 1. - startY ); if( coords[1] <= up ) upperBound = true;
        bool lowerBound = false; RealType low = 0.5 * length / height * coords[0] + startY; if( coords[1] >= low ) lowerBound = true;
        if( rightBound && upperBound && lowerBound )  material[nodeIdx] = mhard;
    }
  }

  void constructHoles( VectorType &material ) const{
    RealType lx = 0.0, ly = 0.0;
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
        if( coords[0] > lx ) lx = coords[0];
        if( coords[1] > ly ) ly = coords[1];
    }
    const RealType pi = 4 * atan ( 1.0 );
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
        RealType tmp = cos( 6. * pi * coords[0] / lx ) * cos( 4. * pi * coords[1] ) + 0.6;
        tmp -= std::max<RealType>( 200.*(1.e-2 - coords[0] * coords[0] - (coords[1] - 0.5*ly) * (coords[1] - 0.5*ly) ) , 0. );
        tmp -= std::max<RealType>( 100. * ( coords[0] + coords[1] - lx - ly + 0.1), 0. );
        tmp -= std::max<RealType>( 100. * ( coords[0] - coords[1] - lx + 0.1 ), 0. );

        if( tmp > 1.0 ) tmp = 1.0;
        if( tmp < -1.0 ) tmp = -1.0;
        material[nodeIdx] = tmp;
    }
  }

  void constructDisc( VectorType &material, const RealType radius, const Point3DType& center ) const{
    const RealType radiusSqr = pesopt::Sqr( radius );
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const Point3DType& coords ( _mesh.getVertex(nodeIdx) );
        const RealType distToCenter = pesopt::Sqr( coords(0) - center(0) ) + pesopt::Sqr( coords(1) - center(1) ) + pesopt::Sqr( coords(2) - center(2) );
        if( distToCenter < radiusSqr ) material[nodeIdx] = 1.0;
    }
  }

  void constructRandomMaterial( VectorType &material ) const{  material = VectorType::Random( material.size() );}

//   template<typename ParameterParserType>
//   void switchMaterialType( VectorType &material, const ParameterParserType &parser ) const{
//       switch( parser.template get<int>( "InitMaterial.initMaterialType" ) ){
//         case 1:
//           constructConstantMaterial( material, parser.template get<double>( "InitMaterial.materialConstant" ) );
//           break;
//         case 2:{
//            constructConstantMaterial( material, -1. );
//            const int numLayers = parser.template get<int>( "InitMaterial.numLayers" );
//            for( int layer=1; layer<=numLayers; ++layer ){
//              TangentVecType startPoint; parser.template getFixSizeVector<RealType,TangentVecType> ( pesopt::strprintf( "InitMaterial.StartRangeLayer%d", layer ).c_str(), startPoint );
//              TangentVecType endPoint; parser.template getFixSizeVector<RealType,TangentVecType> ( pesopt::strprintf( "InitMaterial.EndRangeLayer%d", layer ).c_str(), endPoint );
//              constructLayerMaterial( material, startPoint, endPoint );
//            }
//         }break;
//         case 3:{
//            constructConstantMaterial( material, -1. );
//            const int numRays = parser.template get<int>( "InitMaterial.numRays" );
//            for( int ray=1; ray<=numRays; ++ray ){
//              TangentVecType startPoint; parser.template getFixSizeVector<RealType,TangentVecType> ( pesopt::strprintf( "InitMaterial.StartRangeRay%d", ray ).c_str(), startPoint );
//              TangentVecType rayVector; parser.template getFixSizeVector<RealType,TangentVecType> ( pesopt::strprintf( "InitMaterial.VectorRay%d", ray ).c_str(), rayVector );
//              RealType thickness = parser.template get<RealType> ( pesopt::strprintf( "InitMaterial.ThicknessRay%d", ray ).c_str() );
//              constructRayMaterial( material, startPoint, rayVector, thickness );
//            }
//         }break;
//         case 4:{
//            constructConstantMaterial( material, -1. );
//            constructTriangleMaterial( material, parser.template get<RealType> ("InitMaterial.TriangleLenght" ), parser.template get<RealType> ("InitMaterial.TriangleArea" ) );
//         }break;
//         case 10:
//             constructHoles( material );
//             break;
//         case 1000 :
//             constructRandomMaterial( material );
//             break;
//         default:
//           throw std::invalid_argument( pesopt::strprintf ( "Wrong channel for initMaterialType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
//           break;
//       }
//   }

  void switchMaterialTypeForFixedArea( const int designType, const RealType areaHardMaterial, VectorType &material, string &designTypeName ) const{
        switch( designType ){
            case 1:{
                designTypeName = "BoundaryLayerMaterial";
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint(0.,0.,0.);
                TangentVecType endPoint (areaHardMaterial, 1.,0.);
                constructLayerMaterial( material, startPoint, endPoint );
            }break;
            case 2:{
                designTypeName = "OrthogonalLayerMiddle";
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint(0.,0.5 - 0.5 * areaHardMaterial,0.);
                TangentVecType endPoint (1.,0.5 + 0.5 * areaHardMaterial ,0.);
                constructLayerMaterial( material, startPoint, endPoint );
            }break;
            case 3:{
                designTypeName = "2OrthogonalLayersBoundary";
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1(0.,0.,0.);
                TangentVecType endPoint1 (1.,0.5 * areaHardMaterial ,0.);
                constructLayerMaterial( material, startPoint1, endPoint1 );
                TangentVecType startPoint2(0.,1 - 0.5 * areaHardMaterial, 0.);
                TangentVecType endPoint2 (1.,1. ,0.);
                constructLayerMaterial( material, startPoint2, endPoint2 );
            }break;
            case 4:{
                designTypeName = "2OrthogonalLayersMiddle";
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1(0., 0.25 - 0.25 * areaHardMaterial ,0.);
                TangentVecType endPoint1 (1., 0.25 + 0.25 * areaHardMaterial ,0.);
                constructLayerMaterial( material, startPoint1, endPoint1 );
                TangentVecType startPoint2(0., 0.75 - 0.25 * areaHardMaterial, 0.);
                TangentVecType endPoint2 (1., 0.75 + 0.25 * areaHardMaterial , 0.);
                constructLayerMaterial( material, startPoint2, endPoint2 );
            }break;
            case 5:{
                designTypeName = "BoundaryAndOrthogonalLayer";
                RealType delta = 1. - std::sqrt( 1. - areaHardMaterial );
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1 (0.,0.,0.);
                TangentVecType endPoint1 (delta, 1., 0. );
                constructLayerMaterial( material, startPoint1, endPoint1 );
                TangentVecType startPoint2 (delta ,0.5 - 0.5 * delta , 0.);
                TangentVecType endPoint2 (1., 0.5 + 0.5 * delta , 0. );
                constructLayerMaterial( material, startPoint2, endPoint2 );
            }break;
            case 6:{
                designTypeName = "BoundaryAnd2OrthogonalLayer";
                RealType delta = 3./4. - std::sqrt( 9./16. - 0.5 * areaHardMaterial );
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1 (0.,0.,0.);
                TangentVecType endPoint1 (delta, 1., 0. );
                constructLayerMaterial( material, startPoint1, endPoint1 );
                TangentVecType startPoint2 (delta , 0. , 0.);
                TangentVecType endPoint2 (1., delta , 0. );
                constructLayerMaterial( material, startPoint2, endPoint2 );
                TangentVecType startPoint3 (delta , 1. - delta , 0.);
                TangentVecType endPoint3 (1., 1. , 0. );
                constructLayerMaterial( material, startPoint3, endPoint3 );
            }break;
            case 7:{
                designTypeName = "BoundaryAnd4OrthogonalLayer";
                RealType delta = 3./4. - std::sqrt( 9./16. - 0.5 * areaHardMaterial );
                RealType alpha = (1. - 2. * delta) / 3.;
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1 (0.,0.,0.);
                TangentVecType endPoint1 (delta, 1., 0. );
                constructLayerMaterial( material, startPoint1, endPoint1 );
                TangentVecType startPoint2 (delta , 0. , 0.);
                TangentVecType endPoint2 (1., 0.5 * delta , 0. );
                constructLayerMaterial( material, startPoint2, endPoint2 );
                TangentVecType startPoint4 (delta, 0.5 * delta + alpha , 0.);
                TangentVecType endPoint4 (1., delta + alpha , 0. );
                constructLayerMaterial( material, startPoint4, endPoint4 );
                TangentVecType startPoint5 (delta , delta + 2 * alpha, 0.);
                TangentVecType endPoint5 (1., 1.5*delta + 2 * alpha , 0. );
                constructLayerMaterial( material, startPoint5, endPoint5 );
                TangentVecType startPoint3 (delta , 1. - 0.5 * delta , 0.);
                TangentVecType endPoint3 (1., 1. , 0. );
                constructLayerMaterial( material, startPoint3, endPoint3 );
            }break;
            case 8:{
                designTypeName = "SquareAtClampedBoundary";
                constructConstantMaterial( material, -1. );
                RealType side = sqrt( areaHardMaterial );
                TangentVecType startPoint(0.,0.5 - 0.5 * side,0.);
                TangentVecType endPoint (side,0.5 + 0.5 * side ,0.);
                constructLayerMaterial( material, startPoint, endPoint );
            }break;
            case 11:{
                designTypeName = "Cross";
                RealType delta = 1. - std::sqrt( 1. - areaHardMaterial );
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1 (0.5 - 0.5 * delta,0.,0.);
                TangentVecType endPoint1 ( 0.5 + 0.5 * delta, 1., 0. );
                constructLayerMaterial( material, startPoint1, endPoint1 );
                TangentVecType startPoint2 (0, 0.5 - 0.5 * delta , 0.);
                TangentVecType endPoint2 (1., 0.5 + 0.5 * delta , 0. );
                constructLayerMaterial( material, startPoint2, endPoint2 );
            }break;
            case 12:{
                designTypeName = "DiagonalCross";
                // area \approx 2 sqrt(2) thickness
                RealType thickness = areaHardMaterial / std::sqrt(8);
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1 ( 0.,0.,0. );
                TangentVecType rayVector1 ( 1.,1.,0.);
                constructRayMaterial( material, startPoint1, rayVector1, thickness );
                TangentVecType startPoint2 ( 0.,1.,0. );
                TangentVecType rayVector2(1.,-1.,0.);
                constructRayMaterial( material, startPoint2, rayVector2, thickness );
           }break;

            case 31:{
                designTypeName = "HolesDiffuse";
                //TODO area \approx 2 sqrt(2) thickness
                constructHoles( material );
           }break;
           case 32:{
                designTypeName = "Holes01";
                //TODO area \approx 2 sqrt(2) thickness
                constructHoles( material );
                for( int i=0; i<material.size(); ++i){
                  if( material[i] > 0. ) material[i] = 1.;
                  else material[i] = -1.;
                }
           }break;

            case 1000:{
                designTypeName = "Random";
                constructRandomMaterial( material );
                material *= 2. * areaHardMaterial;
                for( int i=0; i<material.size(); ++i ) material[i] += 2. * areaHardMaterial - 1.;
           }break;
            case 1001:{
                designTypeName = "Constant";
                RealType materialConstant = 2. * areaHardMaterial - 1.;
                constructConstantMaterial( material, materialConstant );
           }break;

            default:
                throw std::invalid_argument( pesopt::strprintf ( "Wrong channel: designType = %d. In File %s at line %d.", designType, __FILE__, __LINE__ ).c_str() );
                break;
         }
    }

};




template<typename ConfiguratorType>
class ShellMaterialGenerator<ConfiguratorType,ElementValuedDofs>{

 protected:

  typedef typename ConfiguratorType::RealType       RealType;
  typedef typename ConfiguratorType::InitType       MeshType;
  typedef typename ConfiguratorType::MaskType       MaskType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Point3DType    Point3DType;
  typedef typename ConfiguratorType::VectorType     VectorType;
  typedef typename ConfiguratorType::DTContainer    DataTypeContainer;

  const MeshType &_mesh;
  const int _numElements;

 public:

  ShellMaterialGenerator( const ConfiguratorType &conf ) :  _mesh( conf.getInitializer() ), _numElements ( _mesh.getNumElements() ) {}

  ShellMaterialGenerator( const MeshType &mesh ) : _mesh( mesh ), _numElements ( _mesh.getNumElements() ) {}

  void constructConstantMaterial( VectorType &material, RealType materialConstant ) const{
      material.setZero();
      for( int elementIdx=0; elementIdx<_numElements; ++elementIdx ) material[elementIdx] = materialConstant;
   }

  //! constructs material in a box spanned by start and end point
  void constructLayerMaterial( VectorType &material,  const TangentVecType startPoint, const TangentVecType endPoint,
                               const RealType mhard = 1.0, const RealType msoft = -1.0 ) const{
    for( int elementIdx=0; elementIdx<_numElements; ++elementIdx ) {
        const Point3DType& coords ( _mesh.getTriang(elementIdx).getMidPoint() );
        if( (coords[0] >= startPoint[0]) && (coords[1] >= startPoint[1]) && (coords[2] >= startPoint[2])
         && (coords[0] <= endPoint[0]) && (coords[1] <= endPoint[1]) && (coords[2] <= endPoint[2]) ){
            material[elementIdx] = mhard;
        }
    }
  }

 //! construct ray spanned by start + t rayVector with thickness
 void constructRayMaterial( VectorType &material,  const TangentVecType startPoint, const TangentVecType rayVector, const RealType thickness,
                               const RealType mhard = 1.0, const RealType msoft = -1.0 ) const{
    for( int elementIdx=0; elementIdx<_numElements; ++elementIdx ) {
        const Point3DType& coords ( _mesh.getTriang(elementIdx).getMidPoint() );
        TangentVecType projection;
        RealType projStep = rayVector.dot( coords - startPoint ); projStep /= rayVector.squaredNorm();
        projection = startPoint + projStep * rayVector;
        RealType distance = (projection - coords ).norm();
        if( distance < thickness )  material[elementIdx] = mhard;
    }
  }

 //! construct ray spanned by start + t rayVector with thickness
 void constructTriangleMaterial( VectorType &material,  const RealType length, const RealType area,
                               const RealType mhard = 1.0, const RealType msoft = -1.0 ) const{
    const RealType height = 2. * area / length; const RealType startY = 0.5 * (1. - length );
    for( int elementIdx=0; elementIdx<_numElements; ++elementIdx ) {
        const Point3DType& coords ( _mesh.getTriang(elementIdx).getMidPoint() );
        bool rightBound = false; if( coords[0] <= height ) rightBound = true;
        bool upperBound = false; RealType up = -0.5 * length / height * coords[0] + ( 1. - startY ); if( coords[1] <= up ) upperBound = true;
        bool lowerBound = false; RealType low = 0.5 * length / height * coords[0] + startY; if( coords[1] >= low ) lowerBound = true;
        if( rightBound && upperBound && lowerBound )  material[elementIdx] = mhard;
    }
  }

  void constructHoles( VectorType &material ) const{
    RealType lx = 0.0, ly = 0.0;
    for( int elementIdx=0; elementIdx<_numElements; ++elementIdx ) {
        const Point3DType& coords ( _mesh.getTriang(elementIdx).getMidPoint() );
        if( coords[0] > lx ) lx = coords[0];
        if( coords[1] > ly ) ly = coords[1];
    }
    const RealType pi = 4 * atan ( 1.0 );
    for( int elementIdx=0; elementIdx<_numElements; ++elementIdx ) {
        const Point3DType& coords ( _mesh.getTriang(elementIdx).getMidPoint() );
        RealType tmp = cos( 6. * pi * coords[0] / lx ) * cos( 4. * pi * coords[1] ) + 0.6;
        tmp -= std::max<RealType>( 200.*(1.e-2 - coords[0] * coords[0] - (coords[1] - 0.5*ly) * (coords[1] - 0.5*ly) ) , 0. );
        tmp -= std::max<RealType>( 100. * ( coords[0] + coords[1] - lx - ly + 0.1), 0. );
        tmp -= std::max<RealType>( 100. * ( coords[0] - coords[1] - lx + 0.1 ), 0. );

        if( tmp > 1.0 ) tmp = 1.0;
        if( tmp < -1.0 ) tmp = -1.0;
        material[elementIdx] = tmp;
    }
  }

  void constructDisc( VectorType &material, const RealType radius, const Point3DType& center ) const{
    const RealType radiusSqr = pesopt::Sqr( radius );
    for( int elementIdx=0; elementIdx<_numElements; ++elementIdx ) {
        const Point3DType& coords ( _mesh.getTriang(elementIdx).getMidPoint() );
        const RealType distToCenter = pesopt::Sqr( coords(0) - center(0) ) + pesopt::Sqr( coords(1) - center(1) ) + pesopt::Sqr( coords(2) - center(2) );
        if( distToCenter < radiusSqr ) material[elementIdx] = 1.0;
    }
  }

  void constructRandomMaterial( VectorType &material ) const{  material = VectorType::Random( material.size() );}

  void switchMaterialTypeForFixedArea( const int designType, const RealType areaHardMaterial, VectorType &material, string &designTypeName ) const{
        switch( designType ){
            case 1:{
                designTypeName = "BoundaryLayerMaterial";
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint(0.,0.,0.);
                TangentVecType endPoint (areaHardMaterial, 1.,0.);
                constructLayerMaterial( material, startPoint, endPoint );
            }break;
            case 2:{
                designTypeName = "OrthogonalLayerMiddle";
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint(0.,0.5 - 0.5 * areaHardMaterial,0.);
                TangentVecType endPoint (1.,0.5 + 0.5 * areaHardMaterial ,0.);
                constructLayerMaterial( material, startPoint, endPoint );
            }break;
            case 3:{
                designTypeName = "2OrthogonalLayersBoundary";
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1(0.,0.,0.);
                TangentVecType endPoint1 (1.,0.5 * areaHardMaterial ,0.);
                constructLayerMaterial( material, startPoint1, endPoint1 );
                TangentVecType startPoint2(0.,1 - 0.5 * areaHardMaterial, 0.);
                TangentVecType endPoint2 (1.,1. ,0.);
                constructLayerMaterial( material, startPoint2, endPoint2 );
            }break;
            case 4:{
                designTypeName = "2OrthogonalLayersMiddle";
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1(0., 0.25 - 0.25 * areaHardMaterial ,0.);
                TangentVecType endPoint1 (1., 0.25 + 0.25 * areaHardMaterial ,0.);
                constructLayerMaterial( material, startPoint1, endPoint1 );
                TangentVecType startPoint2(0., 0.75 - 0.25 * areaHardMaterial, 0.);
                TangentVecType endPoint2 (1., 0.75 + 0.25 * areaHardMaterial , 0.);
                constructLayerMaterial( material, startPoint2, endPoint2 );
            }break;
            case 5:{
                designTypeName = "BoundaryAndOrthogonalLayer";
                RealType delta = 1. - std::sqrt( 1. - areaHardMaterial );
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1 (0.,0.,0.);
                TangentVecType endPoint1 (delta, 1., 0. );
                constructLayerMaterial( material, startPoint1, endPoint1 );
                TangentVecType startPoint2 (delta ,0.5 - 0.5 * delta , 0.);
                TangentVecType endPoint2 (1., 0.5 + 0.5 * delta , 0. );
                constructLayerMaterial( material, startPoint2, endPoint2 );
            }break;
            case 6:{
                designTypeName = "BoundaryAnd2OrthogonalLayer";
                RealType delta = 3./4. - std::sqrt( 9./16. - 0.5 * areaHardMaterial );
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1 (0.,0.,0.);
                TangentVecType endPoint1 (delta, 1., 0. );
                constructLayerMaterial( material, startPoint1, endPoint1 );
                TangentVecType startPoint2 (delta , 0. , 0.);
                TangentVecType endPoint2 (1., delta , 0. );
                constructLayerMaterial( material, startPoint2, endPoint2 );
                TangentVecType startPoint3 (delta , 1. - delta , 0.);
                TangentVecType endPoint3 (1., 1. , 0. );
                constructLayerMaterial( material, startPoint3, endPoint3 );
            }break;
            case 7:{
                designTypeName = "BoundaryAnd4OrthogonalLayer";
                RealType delta = 3./4. - std::sqrt( 9./16. - 0.5 * areaHardMaterial );
                RealType alpha = (1. - 2. * delta) / 3.;
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1 (0.,0.,0.);
                TangentVecType endPoint1 (delta, 1., 0. );
                constructLayerMaterial( material, startPoint1, endPoint1 );
                TangentVecType startPoint2 (delta , 0. , 0.);
                TangentVecType endPoint2 (1., 0.5 * delta , 0. );
                constructLayerMaterial( material, startPoint2, endPoint2 );
                TangentVecType startPoint4 (delta, 0.5 * delta + alpha , 0.);
                TangentVecType endPoint4 (1., delta + alpha , 0. );
                constructLayerMaterial( material, startPoint4, endPoint4 );
                TangentVecType startPoint5 (delta , delta + 2 * alpha, 0.);
                TangentVecType endPoint5 (1., 1.5*delta + 2 * alpha , 0. );
                constructLayerMaterial( material, startPoint5, endPoint5 );
                TangentVecType startPoint3 (delta , 1. - 0.5 * delta , 0.);
                TangentVecType endPoint3 (1., 1. , 0. );
                constructLayerMaterial( material, startPoint3, endPoint3 );
            }break;
            case 8:{
                designTypeName = "SquareAtClampedBoundary";
                constructConstantMaterial( material, -1. );
                RealType side = sqrt( areaHardMaterial );
                TangentVecType startPoint(0.,0.5 - 0.5 * side,0.);
                TangentVecType endPoint (side,0.5 + 0.5 * side ,0.);
                constructLayerMaterial( material, startPoint, endPoint );
            }break;
            case 11:{
                designTypeName = "Cross";
                RealType delta = 1. - std::sqrt( 1. - areaHardMaterial );
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1 (0.5 - 0.5 * delta,0.,0.);
                TangentVecType endPoint1 ( 0.5 + 0.5 * delta, 1., 0. );
                constructLayerMaterial( material, startPoint1, endPoint1 );
                TangentVecType startPoint2 (0, 0.5 - 0.5 * delta , 0.);
                TangentVecType endPoint2 (1., 0.5 + 0.5 * delta , 0. );
                constructLayerMaterial( material, startPoint2, endPoint2 );
            }break;
            case 12:{
                designTypeName = "DiagonalCross";
                // area \approx 2 sqrt(2) thickness
                RealType thickness = areaHardMaterial / std::sqrt(8);
                constructConstantMaterial( material, -1. );
                TangentVecType startPoint1 ( 0.,0.,0. );
                TangentVecType rayVector1 ( 1.,1.,0.);
                constructRayMaterial( material, startPoint1, rayVector1, thickness );
                TangentVecType startPoint2 ( 0.,1.,0. );
                TangentVecType rayVector2(1.,-1.,0.);
                constructRayMaterial( material, startPoint2, rayVector2, thickness );
           }break;

            case 31:{
                designTypeName = "HolesDiffuse";
                //TODO area \approx 2 sqrt(2) thickness
                constructHoles( material );
           }break;
           case 32:{
                designTypeName = "Holes01";
                //TODO area \approx 2 sqrt(2) thickness
                constructHoles( material );
                for( int i=0; i<material.size(); ++i){
                  if( material[i] > 0. ) material[i] = 1.;
                  else material[i] = -1.;
                }
           }break;

            case 1000:{
                designTypeName = "Random";
                constructRandomMaterial( material );
                material *= 2. * areaHardMaterial;
                for( int i=0; i<material.size(); ++i ) material[i] += 2. * areaHardMaterial - 1.;
           }break;
            case 1001:{
                designTypeName = "Constant";
                RealType materialConstant = 2. * areaHardMaterial - 1.;
                constructConstantMaterial( material, materialConstant );
           }break;

            default:
                throw std::invalid_argument( pesopt::strprintf ( "Wrong channel: designType = %d. In File %s at line %d.", designType, __FILE__, __LINE__ ).c_str() );
                break;
         }
    }

};








template< typename ConfiguratorType >
class ShellPlotter{

protected:

  typedef typename ConfiguratorType::RealType       RealType;
  typedef typename ConfiguratorType::InitType       MeshType;
  typedef typename ConfiguratorType::MaskType       MaskType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Point3DType    Point3DType;
  typedef typename ConfiguratorType::VectorType     VectorType;
  typedef typename ConfiguratorType::DTContainer    DataTypeContainer;

  const ConfiguratorType &_conf;
  const MeshType &_mesh;
  const int _numVertices, _numElements, _numGlobalDofsDeform;
  const VectorType &_xA;
  const MaskType &_DirichletMask;
  const string _saveDirectory;
  const string _VTKFileType;

  mutable VTKMeshSaver<MeshType> _meshSaver;

public:
    ShellPlotter( const ConfiguratorType &conf, const VectorType &xA, const MaskType &DirichletMask, const string saveDirectory, const string VTKFileType ) :
  _conf ( conf ),  _mesh( conf.getInitializer() ),
  _numVertices ( _mesh.getNumVertices() ), _numElements( _mesh.getNumElements() ), _numGlobalDofsDeform ( conf.getNumGlobalDofs() ),
  _xA( xA ), _DirichletMask( DirichletMask ),
  _saveDirectory ( saveDirectory ), _VTKFileType ( VTKFileType ),
  _meshSaver( _mesh ) {}


protected:
 void getDeformedMesh( MeshType &mesh, VectorType &deformVec, const string dispOrDeform, const VectorType & dispOrDeformVec ) const{

    bool dispOrDeformType = false;
    if( dispOrDeform == "deform" ) { deformVec = dispOrDeformVec; dispOrDeformType = true; }
    if( dispOrDeform == "disp" ) { deformVec = dispOrDeformVec + _xA; dispOrDeformType = true; }
    if( dispOrDeform == "undeform" ) { deformVec = _xA; dispOrDeformType = true; }

    if( !dispOrDeformType ) throw std::invalid_argument( pesopt::strprintf ( "Wrong channel for dispOrDeformType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );

    for( int i = 0; i < mesh.getNumVertices(); ++i ){
        Point3DType coords;
        for( int comp = 0; comp < 3; ++comp ) coords[comp] = deformVec[i + comp * _numGlobalDofsDeform];
        mesh.setVertex( i, coords );
    }
 }

public:

  void clearData ( ) const { _meshSaver.clearData(); }
//   void updateMesh ( const MeshType & mesh ) const { _meshSaver.updateMesh( mesh ); }
//   void updateMesh ( const string dispOrDeform, const VectorType & dispOrDeformVec  ) const {
//     MeshType meshDeformed( _mesh ); VectorType deformVec ( dispOrDeformVec.size() );
//     this->getDeformedMesh( meshDeformed, deformVec, dispOrDeform, dispOrDeformVec );
//     this->updateMesh( meshDeformed );
//   }
  void addScalarData ( const VectorType & data, string dataDescr, VTKDataSupp supp ) const { _meshSaver.addScalarData( data, dataDescr, supp ); }
  void addVectorData ( const VectorType & data, const int numComponents, string dataDescr, VTKDataSupp supp ) const { _meshSaver.addVectorData( data, numComponents, dataDescr, supp ); }
  void addNormalData ( const VectorType & data, const int numComponents, string dataDescr, VTKDataSupp supp ) const { _meshSaver.addNormalData( data, numComponents, dataDescr, supp ); }

  void saveShellToFile ( const string outfile_base_name ) const{
    _meshSaver.save( pesopt::strprintf (  "%s/%s.%s", _saveDirectory.c_str(), outfile_base_name.c_str(), _VTKFileType.c_str() ), VTKPOLYDATA );
  }


  void saveShellToFile ( const string dispOrDeform,
                         const VectorType & dispOrDeformVec,
                         const string outfile_base_name,
                         const bool plotNormalField = true  ) const{

    // set deformed mesh
    MeshType meshDeformed ( this->_mesh ); VectorType deformVec ( dispOrDeformVec.size() );
    this->getDeformedMesh( meshDeformed, deformVec, dispOrDeform, dispOrDeformVec );

    VTKMeshSaver<MeshType> meshSaver ( meshDeformed );

    //boundary
    VectorType boundary ( this->_numVertices );
    for( int i=0; i<this->_numVertices; ++i ){
      if( this->_DirichletMask[i] ) boundary[i] = 1.;
      else                    boundary[i] = 0.;
    }
    meshSaver.addScalarData ( boundary, "boundary", VERTEX_DATA );

    //topological boundary
    MaskType topologicalBoundaryMask ( this->_numVertices );
    meshDeformed.findTopologicalBoundaryMask( topologicalBoundaryMask  );
    VectorType topologicalBoundary ( this->_numVertices );
    for( int i=0; i<this->_numVertices; ++i ){
      if( topologicalBoundaryMask[i] ) topologicalBoundary[i] = 1.;
      else                             topologicalBoundary[i] = 0.;
    }
    meshSaver.addScalarData ( topologicalBoundary, "topologicalBoundary", VERTEX_DATA );

    //! normals on deformed surface at nodes
    VectorType normalAtNodes ( 3 * this->_numVertices ), tangent1AtNodes ( 3 * this->_numVertices ), tangent2AtNodes( 3 * this->_numVertices );
    if( plotNormalField ){
        if( ConfiguratorType::_ShellFEType == C1Dofs ){
            for( int i = 0; i < meshDeformed.getNumVertices(); ++i ){
                    TangentVecType tangentVec1, tangentVec2;
                    for( int comp = 0; comp < 3; ++comp ) {
                        tangentVec1[comp] = deformVec[i + this->_numVertices + comp * this->_numGlobalDofsDeform];
                        tangentVec2[comp] = deformVec[i + 2 * this->_numVertices + comp * this->_numGlobalDofsDeform];
                    }
                    TangentVecType normalVec = tangentVec1.cross( tangentVec2 );
                    RealType norm = normalVec.norm();
                    normalVec /= norm;
                    for( int comp = 0; comp<3; ++comp ){
                        normalAtNodes[i + comp * this->_numVertices] = normalVec[comp];
                        tangent1AtNodes[i + comp * this->_numVertices] = tangentVec1[comp];
                        tangent2AtNodes[i + comp * this->_numVertices] = tangentVec2[comp];
                    }
            }
        }else{
            meshDeformed.updateAllTriangles();
            meshDeformed.generateApproximativeTangentSpaceAtNodes();
            //TODO optional with boundary mask
            for( int i = 0; i < meshDeformed.getNumVertices(); ++i ){
                    const TangentVecType & tangentVec1 = meshDeformed.getTangentVec1( i ),
                                           tangentVec2 = meshDeformed.getTangentVec2( i ),
                                           normalVec   = meshDeformed.getNormalVec( i );
                    for( int comp = 0; comp<3; ++comp ){
                        normalAtNodes[i + comp * this->_numVertices] = normalVec[comp];
                        tangent1AtNodes[i + comp * this->_numVertices] = tangentVec1[comp];
                        tangent2AtNodes[i + comp * this->_numVertices] = tangentVec2[comp];
                    }
            }
        }
        meshSaver.addNormalData ( normalAtNodes, 3, "normalField", VERTEX_DATA );
        meshSaver.addVectorData ( tangent1AtNodes, 3, "tangentVec1", VERTEX_DATA );
        meshSaver.addVectorData ( tangent2AtNodes, 3, "tangentVec2", VERTEX_DATA );
    }

    meshSaver.save( pesopt::strprintf (  "%s/%s.%s", this->_saveDirectory.c_str(), outfile_base_name.c_str(), this->_VTKFileType.c_str() ), VTKPOLYDATA );
  }


// TODO
//   void savePointCloudDeformedToFile (  const VectorType &disp ) const{
//
//     DKTFEVectorFunctionEvaluator<ConfiguratorType> xADFD ( _conf, _xA , 3 );
//     const VectorType _xB( _xA + disp );
//     DKTFEVectorFunctionEvaluator<ConfiguratorType> xBDFD ( _conf, _xB, 3 );
//     std::vector<Point3DType> pointVecUndeformed; pointVecUndeformed.reserve( _mesh.getNumTriangs() * _conf.maxNumQuadPoints () );
//     std::vector<Point3DType> pointVecDeformed; pointVecDeformed.reserve( _mesh.getNumTriangs() * _conf.maxNumQuadPoints () );
//
//     std::vector<RealType> pointVec_Material;
//     for ( int elementIdx = 0; elementIdx < _mesh.getNumTriangs(); ++elementIdx){
//         const typename ConfiguratorType::ElementType& El ( _mesh.getTriang( elementIdx ) );
//         for ( int localQuadPointIndex = 0; localQuadPointIndex < _conf.maxNumQuadPoints (); ++localQuadPointIndex){
//             Point3DType xA; xADFD.evaluateAtQuadPoint( El, localQuadPointIndex, xA );  pointVecUndeformed.push_back( xA );
//             Point3DType xB; xBDFD.evaluateAtQuadPoint( El, localQuadPointIndex, xB );  pointVecDeformed.push_back( xB );
//             const RealType pf = 0.; pointVec_Material.push_back( pf );
//         }
//     }
//     std::string saveDirectoryPointCloud = _saveDirectory + "/" + "PointCloud";
//     boost::filesystem::create_directory ( saveDirectoryPointCloud );
//
//     PointCloud<DataTypeContainer> pointCloud;
//     pointCloud.saveAsLegacyVTK( pesopt::strprintf ( "%s/PointCloud_Undeformed.%s", saveDirectoryPointCloud.c_str(), _VTKFileType.c_str() ), pointVecUndeformed, pointVec_Material );
//     pointCloud.saveAsLegacyVTK( pesopt::strprintf ( "%s/PointCloud_Deformed.%s", saveDirectoryPointCloud.c_str(), _VTKFileType.c_str() ), pointVecDeformed, pointVec_Material );
//   }

  void saveShellWithFaceDataToFile ( const string dispOrDeform, const VectorType & dispOrDeformVec,
                               const VectorType & faceDataVec, const string outfile_base_name ) const{
    MeshType mesh ( _mesh ); VectorType deformVec ( dispOrDeformVec.size() );
    this->getDeformedMesh( mesh, deformVec, dispOrDeform, dispOrDeformVec );
    VTKMeshSaver<MeshType> meshSaver ( mesh );
    meshSaver.addScalarData ( faceDataVec, "faceData", FACE_DATA );
    meshSaver.save( pesopt::strprintf ( "%s/%s.%s", _saveDirectory.c_str(), outfile_base_name.c_str(), _VTKFileType.c_str() ), VTKPOLYDATA );
  }


//   void plotShellFctAndDerivativesOnElements ( const string dispOrDeform, const VectorType & dispOrDeformVec, const string outfile_base_name  ) const{
//
//       DataOnElementsPlotter<MeshType> plotter;
//       //TODO compare with DFD( disp + xA )
//       DKTFEVectorFunctionEvaluator<ConfiguratorType> dispDFD ( _conf, dispOrDeformVec, 3 );
//       DKTFEVectorFunctionEvaluator<ConfiguratorType> xADFD ( _conf, _xA, 3 );
//
//       std::vector<Point3DType> pointVecUndeformed; pointVecUndeformed.reserve ( 3 * _mesh.getNumTriangs() );
//       for ( int elementIdx = 0; elementIdx < _mesh.getNumTriangs(); ++elementIdx){
//         const typename ConfiguratorType::ElementType& El ( _mesh.getTriang( elementIdx ) );
//         for( int localNodeIndex = 0; localNodeIndex < 3; ++localNodeIndex )
//             pointVecUndeformed.push_back( El.getNode(localNodeIndex) );
//       }
//       plotter.saveAsVTKPOLYDATA ( pointVecUndeformed, pesopt::strprintf ( "%s/%s_xA.%s", _saveDirectory.c_str(), outfile_base_name.c_str(), _VTKFileType.c_str() ) );
//
//       std::vector<Point3DType> pointVecDeformed; pointVecDeformed.reserve ( 3 * _mesh.getNumTriangs() );
// //       std::vector<Point3DType> pointVecD1Deformed; pointVecD1Deformed.reserve ( 3 * _mesh.getNumTriangs() );
// //       std::vector<Point3DType> pointVecD2Deformed; pointVecD2Deformed.reserve ( 3 * _mesh.getNumTriangs() );
// //       std::vector<Point3DType> pointVecD3Deformed; pointVecD3Deformed.reserve ( 3 * _mesh.getNumTriangs() );
//       std::vector<Point3DType> pointVecLaplaceUndeformed; pointVecLaplaceUndeformed.reserve ( 3 * _mesh.getNumTriangs() );
//       std::vector<Point3DType> pointVecLaplaceDeformed; pointVecLaplaceDeformed.reserve ( 3 * _mesh.getNumTriangs() );
//       for ( int elementIdx = 0; elementIdx < _mesh.getNumTriangs(); ++elementIdx){
//         const typename ConfiguratorType::ElementType& El ( _mesh.getTriang( elementIdx ) );
//         for( int localNodeIndex = 0; localNodeIndex < 3; ++localNodeIndex ){
//             DomVecType RefCoords;  El.getRefCoordsFromLocalIndex( localNodeIndex, RefCoords );
//             Point3DType coords, disp;
//             xADFD.evaluate( El, RefCoords, coords );
//             dispDFD.evaluate( El, RefCoords, disp );
//             pointVecDeformed.push_back( coords + disp );
//
// //             Matrix33 GradDisp;
// //             dispDFD.evaluateGradient( El, RefCoords, GradDisp );
//
//
//             //Laplacian
// //             Matrix22 gA, gAinv; xADFD.evaluateFirstFundamentalForm( El, RefCoords, gA ); gAinv = gA.inverse();
// //             Matrix32 dXA; xADFD.evaluateGradient( El, RefCoords, dXA ); //TODO dispDFD.evaluateApproxGradient( El, RefCoords, dX );
// //             Tensor322Type ddXA; xADFD.evaluateApproxHessianSym( El, RefCoords, ddXA );
// //             DomVecType vecForLaplaceA; xADFD.evaluateVectorForLaplacian ( gAinv, dXA, ddXA, vecForLaplaceA );
// //
// //             Matrix32 dX; dispDFD.evaluateGradient( El, RefCoords, dX ); //TODO dispDFD.evaluateApproxGradient( El, RefCoords, dX );
// //             Tensor322Type ddX; dispDFD.evaluateApproxHessianSym( El, RefCoords, ddX );
// //
// //             Point3DType laplaceXA, laplaceX;
// //             dispDFD.evaluateLaplaceBeltrami ( gAinv, dXA, ddXA, vecForLaplaceA, laplaceXA );
// //             dispDFD.evaluateLaplaceBeltrami ( gAinv, dX, ddX, vecForLaplaceA, laplaceX );
// //
// //
// //             pointVecLaplaceUndeformed.push_back( laplaceXA + coords );
// //             pointVecLaplaceDeformed.push_back( laplaceXA + laplaceX + coords ); //TODO????
//
//         }
//       }
//
//       plotter.saveAsVTKPOLYDATA ( pointVecDeformed, pesopt::strprintf ( "%s/%s_xB.%s", _saveDirectory.c_str(), outfile_base_name.c_str(), _VTKFileType.c_str() ) );
// //       plotter.saveAsVTKPOLYDATA ( pointVecLaplaceUndeformed, pesopt::strprintf ( "%s/%s_LapXA.vtk", _saveDirectory.c_str(), outfile_base_name.c_str() ) );
// //       plotter.saveAsVTKPOLYDATA ( pointVecLaplaceDeformed, pesopt::strprintf ( "%s/%s_LapXB.vtk", _saveDirectory.c_str(), outfile_base_name.c_str() ) );
//
//
// //     //! normals on deformed surface at nodes
// //     VectorType normalAtNodes ( 3 * _numVertices );
// //     for( int i = 0; i < mesh.getNumVertices(); ++i ){
// //         TangentVecType tangentVec1, tangentVec2;
// //         for( int comp = 0; comp < 3; ++comp ) {
// //             tangentVec1[comp] = deformVec[i + _numVertices + comp * _numGlobalDofsDeform];
// //             tangentVec2[comp] = deformVec[i + 2 * _numVertices + comp * _numGlobalDofsDeform];
// //         }
// //         TangentVecType normalVec = tangentVec1.cross( tangentVec2 );
// //         RealType norm = normalVec.norm();
// //         normalVec /= norm;
// //         for( int comp = 0; comp<3; ++comp ) normalAtNodes[i + comp * _numVertices] = normalVec[comp];
// //     }
//   }


//   void plotUndeformedShellWithForce ( const string outfile_base_name, const VectorType &force   ) const{
//     DiscreteVectorFunctionStorage<ConfiguratorType> forceCached ( _conf, force, 3 );
//     std::vector<Point3DType> pointVecAtQuadpoints; pointVecAtQuadpoints.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
//     std::vector<TangentVecType> forceVecAtQuadpoints; forceVecAtQuadpoints.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
//     for ( int elementIdx = 0; elementIdx < _conf.getInitializer().getNumTriangs(); ++elementIdx){
//         const typename ConfiguratorType::ElementType& El ( _conf.getInitializer().getTriang( elementIdx ) );
//         for ( int localQuadPointIndex = 0; localQuadPointIndex < _conf.maxNumQuadPoints (); ++localQuadPointIndex){
//             pointVecAtQuadpoints.push_back( _xACachePtr->_coords[elementIdx][localQuadPointIndex] );
//             forceVecAtQuadpoints.push_back ( forceCached._coords[elementIdx][localQuadPointIndex] );
//         }
//     }
//     PointCloud<DataTypeContainer> pointCloud;
//     pointCloud.saveAsLegacyVTK( pesopt::strprintf ( "%s/ForcePointCloud.%s", _saveDirectory.c_str(), _VTKFileType.c_str() ), pointVecAtQuadpoints, forceVecAtQuadpoints );
//   }



};







template< typename ConfiguratorType, ShellFEType MaterialFEType >
class ShellWithMaterialPlotter : public ShellPlotter<ConfiguratorType> {

protected:

  typedef typename ConfiguratorType::RealType       RealType;
  typedef typename ConfiguratorType::InitType       MeshType;
  typedef typename ConfiguratorType::MaskType       MaskType;
  typedef typename ConfiguratorType::TangentVecType TangentVecType;
  typedef typename ConfiguratorType::Point3DType    Point3DType;
  typedef typename ConfiguratorType::VectorType     VectorType;
  typedef typename ConfiguratorType::DTContainer    DataTypeContainer;

public:
    ShellWithMaterialPlotter( const ConfiguratorType &conf, const VectorType &xA, const MaskType &DirichletMask, const string saveDirectory, const string VTKFileType ) :
    ShellPlotter<ConfiguratorType>( conf, xA, DirichletMask, saveDirectory, VTKFileType ) {}


  void saveMeshWithMaterialToFile ( const VectorType & material, const string outfile_base_name ) const{
    VTKMeshSaver<MeshType> meshSaver ( this->_mesh );

    //material
    VectorType materialData;
    if( (MaterialFEType == NodalValuedDofs) || (MaterialFEType == C1Dofs) ){
        materialData.resize( this->_numVertices );
        for( int nodeIdx=0; nodeIdx<this->_numVertices; ++nodeIdx ) materialData[nodeIdx] = material[nodeIdx];
        meshSaver.addScalarData ( materialData, "material", VERTEX_DATA );
    }
    if( MaterialFEType == ElementValuedDofs ){
        materialData.resize( this->_numElements );
        for( int elementIdx=0; elementIdx<this->_numElements; ++elementIdx ) materialData[elementIdx] = material[elementIdx];
        meshSaver.addScalarData ( materialData, "material", FACE_DATA );
    }

    //boundary
    VectorType boundary ( this->_numVertices );
    for( int i=0; i< this->_numVertices; ++i ){
      if( this->_DirichletMask[i] ) boundary[i] = 1.;
      else                    boundary[i] = 0.;
    }
    meshSaver.addScalarData ( boundary, "boundary", VERTEX_DATA );

    //topological boundary
    MaskType topologicalBoundaryMask ( this->_numVertices );
    this->_mesh.findTopologicalBoundaryMask( topologicalBoundaryMask  );
    VectorType topologicalBoundary ( this->_numVertices );
    for( int i=0; i<this->_numVertices; ++i ){
      if( topologicalBoundaryMask[i] ) topologicalBoundary[i] = 1.;
      else                             topologicalBoundary[i] = 0.;
    }
    meshSaver.addScalarData ( topologicalBoundary, "topologicalBoundary", VERTEX_DATA );

    meshSaver.save( pesopt::strprintf (  "%s/%s.%s", this->_saveDirectory.c_str(), outfile_base_name.c_str(), this->_VTKFileType.c_str() ), VTKPOLYDATA );
  }


 void saveShellWithMaterialToFile ( const string dispOrDeform, const VectorType & dispOrDeformVec,
                                    const VectorType & material,
                                    const string outfile_base_name,
                                    const bool plotNormalField = true  ) const{

    // set deformed mesh
    MeshType meshDeformed ( this->_mesh ); VectorType deformVec ( dispOrDeformVec.size() );
    this->getDeformedMesh( meshDeformed, deformVec, dispOrDeform, dispOrDeformVec );

    VTKMeshSaver<MeshType> meshSaver ( meshDeformed );

    //material
    VectorType materialData;
    if( (MaterialFEType == NodalValuedDofs) || (MaterialFEType == C1Dofs) ){
        materialData.resize( this->_numVertices );
        for( int nodeIdx=0; nodeIdx<this->_numVertices; ++nodeIdx ) materialData[nodeIdx] = material[nodeIdx];
        meshSaver.addScalarData ( materialData, "material", VERTEX_DATA );
    }
    if( MaterialFEType == ElementValuedDofs ){
        materialData.resize( this->_numElements );
        for( int elementIdx=0; elementIdx<this->_numElements; ++elementIdx ) materialData[elementIdx] = material[elementIdx];
        meshSaver.addScalarData ( materialData, "material", FACE_DATA );
    }

    //boundary
    VectorType boundary ( this->_numVertices );
    for( int i=0; i<this->_numVertices; ++i ){
      if( this->_DirichletMask[i] ) boundary[i] = 1.;
      else                    boundary[i] = 0.;
    }
    meshSaver.addScalarData ( boundary, "boundary", VERTEX_DATA );

    //topological boundary
    MaskType topologicalBoundaryMask ( this->_numVertices );
    meshDeformed.findTopologicalBoundaryMask( topologicalBoundaryMask  );
    VectorType topologicalBoundary ( this->_numVertices );
    for( int i=0; i<this->_numVertices; ++i ){
      if( topologicalBoundaryMask[i] ) topologicalBoundary[i] = 1.;
      else                             topologicalBoundary[i] = 0.;
    }
    meshSaver.addScalarData ( topologicalBoundary, "topologicalBoundary", VERTEX_DATA );

    //! normals on deformed surface at nodes
    VectorType normalAtNodes ( 3 * this->_numVertices ), tangent1AtNodes ( 3 * this->_numVertices ), tangent2AtNodes( 3 * this->_numVertices );
    if( plotNormalField ){
        if( ConfiguratorType::_ShellFEType == C1Dofs ){
            for( int i = 0; i < meshDeformed.getNumVertices(); ++i ){
                    TangentVecType tangentVec1, tangentVec2;
                    for( int comp = 0; comp < 3; ++comp ) {
                        tangentVec1[comp] = deformVec[i + this->_numVertices + comp * this->_numGlobalDofsDeform];
                        tangentVec2[comp] = deformVec[i + 2 * this->_numVertices + comp * this->_numGlobalDofsDeform];
                    }
                    TangentVecType normalVec = tangentVec1.cross( tangentVec2 );
                    RealType norm = normalVec.norm();
                    normalVec /= norm;
                    for( int comp = 0; comp<3; ++comp ){
                        normalAtNodes[i + comp * this->_numVertices] = normalVec[comp];
                        tangent1AtNodes[i + comp * this->_numVertices] = tangentVec1[comp];
                        tangent2AtNodes[i + comp * this->_numVertices] = tangentVec2[comp];
                    }
            }
        }else{
            meshDeformed.updateAllTriangles();
            meshDeformed.generateApproximativeTangentSpaceAtNodes();
            //TODO optional with boundary mask
            for( int i = 0; i < meshDeformed.getNumVertices(); ++i ){
                    const TangentVecType & tangentVec1 = meshDeformed.getTangentVec1( i ),
                                           tangentVec2 = meshDeformed.getTangentVec2( i ),
                                           normalVec   = meshDeformed.getNormalVec( i );
                    for( int comp = 0; comp<3; ++comp ){
                        normalAtNodes[i + comp * this->_numVertices] = normalVec[comp];
                        tangent1AtNodes[i + comp * this->_numVertices] = tangentVec1[comp];
                        tangent2AtNodes[i + comp * this->_numVertices] = tangentVec2[comp];
                    }
            }
        }
        meshSaver.addNormalData ( normalAtNodes, 3, "normalField", VERTEX_DATA );
        meshSaver.addVectorData ( tangent1AtNodes, 3, "tangentVec1", VERTEX_DATA );
        meshSaver.addVectorData ( tangent2AtNodes, 3, "tangentVec2", VERTEX_DATA );
    }

    meshSaver.save( pesopt::strprintf (  "%s/%s.%s", this->_saveDirectory.c_str(), outfile_base_name.c_str(), this->_VTKFileType.c_str() ), VTKPOLYDATA );
  }

};



#endif //__DKTFEHANDLER_H
