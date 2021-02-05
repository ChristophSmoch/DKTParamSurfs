#ifndef __EXPLICITEFUNCTION_H
#define __EXPLICITEFUNCTION_H


#include <feIntegrator.h>

template< typename DataTypeContainer>
class ExpliciteScalarFunctionBase {

protected:
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
  typedef typename DataTypeContainer::DerivativeVectorValuedType DerivativeVectorValuedType;
  
public:
   ExpliciteScalarFunctionBase() {}
   
  //!*************************************************************************************
  //!    functions that have to be provided in derived classes.
  //!*************************************************************************************
  virtual void evaluate ( const RealVecChart &coords, RealType &fct ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  virtual void evaluateDerivative ( const RealVecChart &coords, RealVecChart &Derivative ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  virtual void evaluateLaplace ( const RealVecChart &coords, RealType &lap ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  virtual void evaluateBiLaplace ( const RealVecChart &coords, RealType &lap ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }

  //!*************************************************************************************
  //!    functions that have to be provided in derived classes.
  //!*************************************************************************************
  template<typename MeshType>
  void evaluateAtNodes ( const MeshType &mesh, typename MeshType::VectorType &fctAtNodesVec ) const {
      for( int i=0; i<mesh.getNumVertices(); ++i ){
          const RealVecChart& coords = mesh.getVertex(i);
          RealType fct;
          this->evaluate( coords, fct );
          fctAtNodesVec[i] = fct;
      }
  }

};




template< typename DataTypeContainer>
class ExpliciteVectorFunctionBase {

protected:
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
  typedef typename DataTypeContainer::DerivativeVectorValuedType DerivativeVectorValuedType;
  
public:
   ExpliciteVectorFunctionBase() {}
   
  //!*************************************************************************************
  //!    functions that have to be provided in derived classes.
  //!*************************************************************************************
  virtual void evaluate ( const RealVecChart &coords, RealVecChart &fct ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  virtual void evaluateDerivative ( const RealVecChart &coords, DerivativeVectorValuedType &Derivative ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  virtual void evaluateInverse ( const RealVecChart &coords, RealVecChart &fct ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  virtual void evaluateDerivativeInverse ( const RealVecChart &coords, DerivativeVectorValuedType &Derivative ) const {
    throw std::invalid_argument( pesopt::strprintf ( "Called the base function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }

};




// //RealVecChart = RealVecChart, RealType = RealType, DerivativeType = RealVecChart
// //f(x,y) = x^3 + y^2 x + y^3
// //Df(x,y) = (3x^2 + y^2,  2xy + 3y^2 )
// template< typename RealVecChart, typename RealType, typename DerivativeType>
// class ExampleFunction_Simple 
// : public ExpliciteScalarFunctionBase<RealVecChart,RealType,DerivativeType,ExampleFunction_Simple<RealVecChart,RealType,DerivativeType> > {
//     
// public:
//     
//   ExampleFunction_Simple() : ExpliciteScalarFunctionBase<RealVecChart,RealType,DerivativeType,ExampleFunction_Simple<RealVecChart,RealType,DerivativeType> > () {}
//     
//    void evaluate ( const RealVecChart &coords, RealType &fct ) const {
//        fct = coords[0] * coords[0] * coords[0] + coords[0] * coords[1] * coords[1] + coords[1] * coords[1] * coords[1];
//    }
//    
//    void evaluateDerivative ( const RealVecChart &coords, DerivativeType &Derivative ) const {
//        Derivative.setZero();
//        Derivative(0) = 3. * coords[0] * coords[0] + coords[1] * coords[1];
//        Derivative(1) = 2. * coords[0] * coords[1] + 3. * coords[1] * coords[1];
//    }
//    
// //    void evaluateLaplace ( const RealVecChart &coords, RealType &fct ) const {
// //        fct.setZero();
// //        fct[2] = -2. * M_PI * M_PI * sin( M_PI * coords[0] ) * sin( M_PI * coords[1] ) ;
// //    }
//    
// //    void evaluateBiLaplace ( const RealVecChart &coords, RealType &fct ) const {
// //        fct.setZero();
// //        fct[2] = 4. * M_PI * M_PI * M_PI * M_PI * sin( M_PI * coords[0] ) * sin( M_PI * coords[1] );
// //    }
//    
//    //TODO
// //     exactLaplaceDKT[i+2*conf.getNumGlobalDofs()] = -2. * M_PI * M_PI * sin( M_PI * coords[0] ) * sin( M_PI * coords[1] );
// //     exactLaplaceDKT[i+2*conf.getNumGlobalDofs()+numVertices] = -2. * M_PI * M_PI * M_PI * cos( M_PI * coords[0] ) * sin( M_PI * coords[1] );
// //     exactLaplaceDKT[i+2*conf.getNumGlobalDofs()+2*numVertices] = -2. * M_PI * M_PI * M_PI * sin( M_PI * coords[0] ) * cos( M_PI * coords[1] );
// 
// /*
//     double L2Norm_UnitSquare ( ) const  { return 0.25; }
//     double L2NormDerivative_UnitSquare ( ) const  {return 4.9348; }
//     double L2NormLaplace_UnitSquare ( ) const  { return 97.409; }*/
//     
// };









// 1) function with simple Dirichlet boundary condition on (0,1)^n
//  u(x_1, .., x_n)   = \prod_{i=1}^n sin( pi x_i )  
// Du(x_1, .., x_n)   = \pi \left( cos(pi x_i) \prod_{j \neq i} sin(pi x_j)  \right)_{i=1, ..., n} 
// D^2u(x_1, .., x_n) = \pi^2 | - sin(pi x) sin(pi y)          cos(pi x) cos(pi y)   |
//             |                   cos(pi x) cos(pi y)      -  sin(pi x) sin(pi y)   |
// laplace(u)(x_1, .., x_n) = -n \pi^2 u(x_1, .., x_n) 
// bilaplace(u)(x_1, .., x_n) = n^2 pi^4 u(x_1, .., x_n) 
// TODO so far in 2D:
// int |u|^2 = 1/4
// int |D u|^2 = pi^2 / 2 \approx 4.9348
// int |D^2 u|^2 = pi^4 \approx 97.409
// int |laplace(u)|^2 = pi^4 \approx 97.409 
//RealVecChart = Vec3, RealType = Vec3, DerivativeType = Mat33
template< typename DataTypeContainer>
class ExampleFunction_SimplySupportedOnUnitCube 
: public ExpliciteScalarFunctionBase<DataTypeContainer > {
    
protected:
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
  typedef typename DataTypeContainer::DerivativeVectorValuedType DerivativeVectorValuedType;
    
public:
    
  ExampleFunction_SimplySupportedOnUnitCube () {}
    
   void evaluate ( const RealVecChart &coords, RealType &fct ) const override {
       fct = 1.;
       for ( int i = 0; i < coords.size(); ++i )
           fct *= sin( M_PI * coords[i] );
   }
   
   void evaluateDerivative ( const RealVecChart &coords, RealVecChart &Derivative ) const override {
       for ( int i = 0; i < coords.size(); ++i ) {
           RealType fac = 1.;
           for ( int j = 0; j < coords.size(); ++j ) {
             if ( j != i ) fac *= sin( M_PI * coords[j] );
           }
           Derivative(i) = M_PI * cos( M_PI * coords[i] ) * fac;
       }
   }
   
   void evaluateLaplace ( const RealVecChart &coords, RealType &lap ) const override {
       const int dim = coords.size();
       RealType fct; this->evaluate( coords,  fct );
       lap = - dim * M_PI * M_PI * fct;
   }
   
   void evaluateBiLaplace ( const RealVecChart &coords, RealType &bilap ) const override {
       RealType lap; this->evaluateLaplace( coords,  lap );
       bilap = lap * lap;
   }

   // TODO so far in 2D
   double L2Norm_UnitSquare ( ) const  { return 0.25; }
   double L2NormDerivative_UnitSquare ( ) const  {return 4.9348; }
   double L2NormLaplace_UnitSquare ( ) const  { return 97.409; }
    
};









// 2) u(x) = 256 \prod_i x_i^2 (x_i-1)^2
// Du(x) = 512 ( x_i (x_i-1) (2x_i-1) \prod_{j \neq i} x_j^2 (x_j-1)^2 )_i
// \trace D^2u(x) = \left( 512  (6x_i^2 - 6x_i +1) \prod_{j \neq i} x_j^2 (x_j-1)^2 \right)_i             
// laplace(u)(x,y) =  512 * \sum_i \left( (6x_i^2 - 6x_i +1) \prod_{j \neq i} x_j^2 (x_j-1)^2 \right)
// TODO in 2D
// bilaplace(u)(x,y) = 2048 ( 6x^2 - 6x + 1) * ( 6y^2 - 6y + 1) + 6144 ( ( x (x - 1) )^2 + (y (y-1))^2 ) ;
// int |u|^2 = 0.16512
// int |D u|^2 = 3.96288
// int |D^2 u|^2 = 166.441
// int |laplace(u)|^2 = 166.441
//RealVecChart = Vec3, RealType = Vec3, DerivativeType = Mat33
template< typename DataTypeContainer>
class ExampleFunction_ClampedBdryOnUnitCube 
: public ExpliciteScalarFunctionBase<DataTypeContainer > {
    
protected:
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::RealVecChart RealVecChart;
  typedef typename DataTypeContainer::DerivativeVectorValuedType DerivativeVectorValuedType;
    
public:
    
  ExampleFunction_ClampedBdryOnUnitCube () {}
    
   void evaluate ( const RealVecChart &coords, RealType &fct ) const {
       fct = 256.;
       for ( int i = 0; i < coords.size(); ++i )
         fct *= pesopt::Sqr( coords[i] * (coords[i] - 1.) );
   }
   
   void evaluateDerivative ( const RealVecChart &coords, RealVecChart &Derivative ) const override {
       for ( int i = 0; i < coords.size(); ++i ) {
           RealType fac = 512.;
           for ( int j = 0; j < coords.size(); ++j ) {
             if ( j != i ) fac *= pesopt::Sqr( coords[j] * (coords[j] - 1.) );
           }
           Derivative(i) = coords[i] * (coords[i] - 1.) * (2. * coords[i] - 1.) * fac;
       }
   }
   
   void evaluateLaplace ( const RealVecChart &coords, RealType &lap ) const {
       lap = 0.;
        for ( int i = 0; i < coords.size(); ++i ) {
           RealType fac = 512.;
           for ( int j = 0; j < coords.size(); ++j ) {
             if ( j != i ) fac *= pesopt::Sqr( coords[j] * (coords[j] - 1.) );
           }
           lap += ( 6. * pesopt::Sqr( coords[i] ) - 6. * coords[i] + 1.) * fac;
       }
   }
   
   // TODO so far in 2D
   void evaluateBiLaplace ( const RealVecChart &coords, RealType &bilap ) const {
       bilap = 0.;
       
       for ( int k = 0; k < coords.size(); ++k ) {
           RealType fac = 6144.;
           for ( int j = 0; j < coords.size(); ++j ) {
             if ( j != k ) fac *= pesopt::Sqr( coords[j] * (coords[j] - 1.) );
           }
           bilap += fac;
       }
       
       for ( int k = 0; k < coords.size(); ++k ) {
           RealType innersum = 0.;
           for ( int i = 0; i < coords.size(); ++i ) {
               RealType innerfac = 1024.;
               for ( int j = 0; j < coords.size(); ++j ) {
                  if ( (j != i) && (j !=  k) ) innerfac *= pesopt::Sqr( coords[j] * (coords[j] - 1.) );
               }
               innersum +=  ( 6. * pesopt::Sqr( coords[i] ) - 6. * coords[i] + 1.) * innerfac;
           }
           bilap += ( 6. * pesopt::Sqr( coords[k] ) - 6. * coords[k] + 1.) * innersum;
       }
      
//       bilap = 2048. * ( 6. * pesopt::Sqr( coords[0] ) - 6. * coords[0] + 1.) * ( 6. * pesopt::Sqr( coords[1] ) - 6. * coords[1] + 1.) 
//              + 6144. * ( pesopt::Sqr( coords[0] * (coords[0] - 1.) ) + pesopt::Sqr( coords[1] * (coords[1] - 1.) ) );
   }
   
   
   //TODO derivative of Laplace 
//     exactLaplaceDKT[i+2*conf.getNumGlobalDofs()]  = 512. * ( 6. * pesopt::Sqr( coords[0] ) - 6. * coords[0] + 1.) * pesopt::Sqr( coords[1] * (coords[1] - 1.) )
//         + 512. * ( 6. * pesopt::Sqr( coords[1] ) - 6. * coords[1] + 1.) * pesopt::Sqr( coords[0] * (coords[0] - 1.) );    
//     exactLaplaceDKT[i+2*conf.getNumGlobalDofs()+numVertices] = 
//                           3072. * ( 2. * coords[0] - 1. ) * pesopt::Sqr( coords[1] * (coords[1] - 1. ) )
//                         + 1024. * ( 6. * pesopt::Sqr( coords[1] ) - 6. * coords[1] + 1.) * coords[0] * ( coords[0] - 1.) * (2. * coords[0] - 1.);
//     exactLaplaceDKT[i+2*conf.getNumGlobalDofs()+2*numVertices] = 
//                           3072. * ( 2. * coords[1] - 1. ) * pesopt::Sqr( coords[0] * (coords[0] - 1. ) )
//                         + 1024. * ( 6. * pesopt::Sqr( coords[0] ) - 6. * coords[0] + 1.) * coords[1] * ( coords[1] - 1.) * (2. * coords[1] - 1.);
    
    // TODO so far in 2D
    double L2Norm_UnitSquare ( ) const  { return  0.16512;}
    double L2NormDerivative_UnitSquare ( ) const  {return 3.96288;}
    double L2NormLaplace_UnitSquare ( ) const  { return 213.995;}
};










template<typename ConfiguratorType>
class FEIntegratorVecConst
: public FEIntegratorVec< ConfiguratorType, FEIntegratorVecConst<ConfiguratorType> > {
protected: 

    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
    typedef typename ConfiguratorType::DTContainer DataTypeContainer;

    const RealType _c;
    
public:

    FEIntegratorVecConst ( const ConfiguratorType &conf, const RealType c  ) 
    : FEIntegratorVec<ConfiguratorType, FEIntegratorVecConst<ConfiguratorType> > ( conf ),
      _c ( c ) {  }

    RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
        return _c; 
    }

};



template<typename ConfiguratorType>
class FEIntegratorVecOfExpliciteFunction
: public FEIntegratorVec< ConfiguratorType, FEIntegratorVecOfExpliciteFunction<ConfiguratorType> > {
protected: 

    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::RealVecChart RealVecChart;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
    typedef typename ConfiguratorType::DTContainer DataTypeContainer;
    typedef typename ConfiguratorType::QuadRuleType QuadRuleType;
    QuadRuleType _quadRule;
    const ExpliciteScalarFunctionBase<DataTypeContainer> &_expliciteFct;
    
public:

    FEIntegratorVecOfExpliciteFunction ( const ConfiguratorType &conf, 
                                         const ExpliciteScalarFunctionBase<DataTypeContainer> &expliciteFct  ) 
    : FEIntegratorVec<ConfiguratorType, FEIntegratorVecOfExpliciteFunction<ConfiguratorType> > ( conf ),
      _expliciteFct( expliciteFct ) {  }

    RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int quadPoint ) const {
        RealType f;
        RealVecChart GlobalCoord;
        El.getGlobalCoord ( _quadRule.getRefCoord( quadPoint ), GlobalCoord );
        _expliciteFct.evaluate( GlobalCoord, f );
        return f;
    }

};



template<typename ConfiguratorType>
class FEIntegratorVecOfExpliciteFunctionLaplace
: public FEIntegratorVec< ConfiguratorType, FEIntegratorVecOfExpliciteFunctionLaplace<ConfiguratorType> > {
protected: 

    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::RealVecChart RealVecChart;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
    typedef typename ConfiguratorType::DTContainer DataTypeContainer;
    typedef typename ConfiguratorType::QuadRuleType QuadRuleType;
    QuadRuleType _quadRule;
    const ExpliciteScalarFunctionBase<DataTypeContainer> &_expliciteFct;
    
public:

    FEIntegratorVecOfExpliciteFunctionLaplace ( const ConfiguratorType &conf, 
                                         const ExpliciteScalarFunctionBase<DataTypeContainer> &expliciteFct  ) 
    : FEIntegratorVec<ConfiguratorType, FEIntegratorVecOfExpliciteFunctionLaplace<ConfiguratorType> > ( conf ),
      _expliciteFct( expliciteFct ) {  }

    RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int quadPoint ) const {
        RealType f;
        RealVecChart GlobalCoord;
        El.getGlobalCoord ( _quadRule.getRefCoord( quadPoint ), GlobalCoord );
        _expliciteFct.evaluateLaplace( GlobalCoord, f );
        return f;
    }

};


template<typename ConfiguratorType>
class FEIntegratorVecOfExpliciteFunctionBiLaplace
: public FEIntegratorVec< ConfiguratorType, FEIntegratorVecOfExpliciteFunctionBiLaplace<ConfiguratorType> > {
protected: 

    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::RealVecChart RealVecChart;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
    typedef typename ConfiguratorType::DTContainer DataTypeContainer;
    typedef typename ConfiguratorType::QuadRuleType QuadRuleType;
    QuadRuleType _quadRule;
    const ExpliciteScalarFunctionBase<DataTypeContainer> &_expliciteFct;
    
public:

    FEIntegratorVecOfExpliciteFunctionBiLaplace ( const ConfiguratorType &conf, 
                                         const ExpliciteScalarFunctionBase<DataTypeContainer> &expliciteFct  ) 
    : FEIntegratorVec<ConfiguratorType, FEIntegratorVecOfExpliciteFunctionBiLaplace<ConfiguratorType> > ( conf ),
      _expliciteFct( expliciteFct ) {  }

    RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int quadPoint ) const {
        RealType f;
        RealVecChart GlobalCoord;
        El.getGlobalCoord ( _quadRule.getRefCoord( quadPoint ), GlobalCoord );
        _expliciteFct.evaluateBiLaplace( GlobalCoord, f );
        return f;
    }

};




#endif
