#ifndef __QUOCFEBOUNDARYINTEGRATOR_H
#define __QUOCFEBOUNDARYINTEGRATOR_H

#include <feBoundaryIntegrator.h>
#include "quocFEDefines.h"

//!===========================================================================================================================
//! Boundary Integration Intefaces
//!===========================================================================================================================


//! Base class for computing \f$ \int_{\partial \Omega} f \cdot \psi_i \, da \f$. 
//The function \f$ f \f$ has to be implemented
// \psi_i are the (vector-valued!) basis functions
// template <typename ConfiguratorType, typename Imp>
// class QuocBoundaryIntegrationBase {
// public:
//   typedef typename ConfiguratorType::RealType RealType;
//   typedef typename ConfiguratorType::ElementType ElementType;
//   typedef typename ConfiguratorType::BoundaryElementType BoundaryElementType;
//   typedef typename ConfiguratorType::VectorType VectorType;
//   typedef typename ConfiguratorType::RealVecChart RealVecChart;
//   typedef typename ConfiguratorType::RealVecChartBoundary RealVecChartBoundary;
//   typedef typename ConfiguratorType::BoundaryQuadType BoundaryQuadType;
// 
//   const ConfiguratorType &_config;
//   
// public:
// 
//   QuocBoundaryIntegrationBase ( const ConfiguratorType &config ) : _config ( config ) {}
// 
//   void assembleAdd ( VectorType &dest ) const {
//     const int numGlobalDofs = _config.getNumGlobalDofs();
//     const int numLocalDofs = _config.getNumLocalBoundaryDofs();
//    
//     for ( int boundaryElementIdx = 0; boundaryElementIdx < _config.getMesh().getNumBoundaryElements(); ++boundaryElementIdx){
//       const BoundaryElementType& bdrEl ( _config.getMesh().getBoundaryElement( boundaryElementIdx ) );
//       BoundaryQuadType quad;
//       RealVecChart refCoord;
//       for ( int q = 0; q < _config.maxNumBoundaryQuadPoints(); q++ ) {
//         RealVecChartBoundary refCoordBoundary = quad.getRefCoord ( q );
//         bdrEl.getRefCoord( refCoordBoundary, refCoord );
//         RealVecChart aux; asImp().getNonlinearity ( bdrEl, refCoord, aux );
// 
//         for ( int locNodeIndex = 0; locNodeIndex < numLocalDofs; locNodeIndex++ ) {
//           int b = bdrEl.getNodeIndexOfElement( locNodeIndex );
//           for ( int comp = 0; comp < ConfiguratorType::dimChartDomain; ++comp ){
//             dest[ _config.localToGlobal ( bdrEl.getElement(), b ) + comp * numGlobalDofs ] += quad.getWeight ( q ) * bdrEl.getVolumeOfBoundaryElement() * _config.getBaseFunctionSet( ).evaluate ( b, refCoord ) * aux[comp];
//           }
//         }
//       }
//     } 
//     
//   }
//   
//   //! this method computes \f$ f \f$ and has to be implemented in the derived class
//  void getNonlinearity ( const BoundaryElementType &bdrEl, const RealVecChart refCoord, RealVecChart &aux ) const {
//     throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
//     return asImp().getNonlinearity ( bdrEl, refCoord, aux );
//   }
// 
// protected:
//   /** barton-nackman **/
//   Imp& asImp( ) { return static_cast<Imp&> ( *this ); }
//   const Imp& asImp( ) const { return static_cast<const Imp&> ( *this ); }
// 
// };


template <typename ConfiguratorType>
class IntegrateDispOverBoundary
  : public FEBoundaryIntegrationBase < ConfiguratorType, IntegrateDispOverBoundary<ConfiguratorType> > {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::BoundaryElementType BoundaryElementType;
  typedef typename ConfiguratorType::RealVecChart RealVecChart;

  const std::vector<RealVecChart> &_forceVec; // LEFT, _BOTTOM, (in 3d: _FRONT)
  
public:
  IntegrateDispOverBoundary( const ConfiguratorType &conf, const std::vector<RealVecChart> &f  ) 
  : FEBoundaryIntegrationBase< ConfiguratorType, IntegrateDispOverBoundary< ConfiguratorType> > ( conf ), _forceVec ( f ) {}

  void getNonlinearity( const BoundaryElementType &bdrEl, const RealVecChart & refCoord, RealVecChart &aux ) const {
      aux.setZero();
      switch( bdrEl.getBoundaryType() ){
             case LEFT:  { aux = _forceVec[0];}break;
             case BOTTOM:{ aux = _forceVec[1];}break;
             case RIGHT: { aux = _forceVec[2];}break;
             case TOP:   { aux = _forceVec[3];}break;
             case FRONT: { aux = _forceVec[4];}break;
             case BACK:  { aux = _forceVec[5];}break;
             default:  break;
         }
  }
};


#endif
