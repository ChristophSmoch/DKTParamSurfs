#ifndef __MATOPTDEFINES__H
#define __MATOPTDEFINES__H

#include <unordered_map>
#include "matOptMultipleLoad.h"


template<typename RealType, const int dimChartDomain>
class MaterialPropertiesEvaluator { };


template<typename RealType>
class MaterialPropertiesEvaluator<RealType,1> { 
    
public:
    
 MaterialPropertiesEvaluator ( ) {}
    
 const RealType computeMuFromNuE( const RealType PoissonRatio,  const RealType ElastModulus  ) const { 
    return ElastModulus / ( 2. * ( 1. + PoissonRatio ) ) ;
 }
 
 const RealType computeMuFromNuE_DerivativeInNu( const RealType PoissonRatio,  const RealType ElastModulus  ) const { 
    return -1. * ElastModulus / ( 2. * ( 1. + PoissonRatio ) * ( 1. + PoissonRatio ) ) ;
 }
 
 const RealType computeMuFromNuE_DerivativeInE( const RealType PoissonRatio,  const RealType ElastModulus  ) const { 
    return 1. / ( 2. * ( 1. + PoissonRatio ) ) ;
 }
    

 const RealType computeLambdaFromNuE( const RealType PoissonRatio,  const RealType ElastModulus  ) const {
     return ElastModulus * PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - PoissonRatio ) );
 }
};


template<typename RealType>
class MaterialPropertiesEvaluator<RealType,2> { 
    
public:
    
 MaterialPropertiesEvaluator ( ) {}
    
 const RealType computeMuFromNuE( const RealType PoissonRatio,  const RealType ElastModulus  ) const { 
    return ElastModulus / ( 2. * ( 1. + PoissonRatio ) ) ;
 }
    
 const RealType computeMuFromNuE_DerivativeInNu( const RealType PoissonRatio,  const RealType ElastModulus  ) const { 
    return -1. * ElastModulus / ( 2. * ( 1. + PoissonRatio ) * ( 1. + PoissonRatio ) ) ;
 }
 
 const RealType computeMuFromNuE_DerivativeInE( const RealType PoissonRatio,  const RealType ElastModulus  ) const { 
    return 1. / ( 2. * ( 1. + PoissonRatio ) ) ;
 }
 
 const RealType computeLambdaFromNuE ( const RealType PoissonRatio,  const RealType ElastModulus  ) const {
  return ElastModulus * PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - PoissonRatio ) );
 }
 
 const RealType computeLambdaFromNuE_DerivativeInNu( const RealType PoissonRatio,  const RealType ElastModulus  ) const {
    return ElastModulus * ( 1. + PoissonRatio * PoissonRatio ) / std::pow( 1. - PoissonRatio * PoissonRatio,  2);
 }
 
 const RealType computeLambdaFromNuE_DerivativeInE( const RealType PoissonRatio,  const RealType ElastModulus  ) const {
    return PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - PoissonRatio ) );
 }
};


template<typename RealType>
class MaterialPropertiesEvaluator<RealType,3> { 
    
public:
    
 MaterialPropertiesEvaluator ( ) {}
    
 const RealType computeMuFromNuE( const RealType PoissonRatio,  const RealType ElastModulus  ) const{ 
    return ElastModulus / ( 2. * ( 1. + PoissonRatio ) ) ;
 }
  
 const RealType computeMuFromNuE_DerivativeInNu( const RealType PoissonRatio,  const RealType ElastModulus  ) const { 
    return -1. * ElastModulus / ( 2. * ( 1. + PoissonRatio ) * ( 1. + PoissonRatio ) ) ;
 }
 
 const RealType computeMuFromNuE_DerivativeInE( const RealType PoissonRatio,  const RealType ElastModulus  ) const { 
    return 1. / ( 2. * ( 1. + PoissonRatio ) ) ;
 }
    
 const RealType computeLambdaFromNuE( const RealType PoissonRatio,  const RealType ElastModulus  ) const {
    return ElastModulus * PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - 2. * PoissonRatio ) );
 }
 
  const RealType computeLambdaFromNuE_DerivativeInNu( const RealType PoissonRatio,  const RealType ElastModulus  ) const {
    return ElastModulus * ( 1. + 2. * PoissonRatio * PoissonRatio ) / std::pow( 1. - PoissonRatio - 2. * PoissonRatio * PoissonRatio, 2 );
 }
 
  const RealType computeLambdaFromNuE_DerivativeInE( const RealType PoissonRatio,  const RealType ElastModulus  ) const {
    return PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - 2. * PoissonRatio ) );
 }
};
  
  



template<typename RealType, const int dimChartDomain>
class MaterialProperties { };



template<typename RealType>
class MaterialProperties<RealType,1> {
protected:
  const int _dimChartDomain = 1;
  std::string _name;
  RealType _density;       // in kg/m^3
  RealType _diffusionCoefficient;
  RealType _elastModulus;  // in GPa
  RealType _poissonRatio;  // ( nu )
  RealType _lambda;        // 1st lame constant in GPa = N/m^2
  RealType _mu;            // 2nd lame constant (shear modulus) in GPa = N/m^2
  RealType _bulkModulus;   // in GPa 

public:
  MaterialProperties() :
    _name ( "DefaultMaterial" ),
    _density ( 0. ), 
    _diffusionCoefficient ( 0. ), 
    _elastModulus ( 0. ), _poissonRatio ( 0. ), 
    _lambda ( 0. ), _mu ( 0. ),
    _bulkModulus ( 0. ) 
    {}

  MaterialProperties ( std::string Name,
             RealType Density,
             RealType DiffusionCoefficient,
             RealType ElastModulus, RealType PoissonRatio ) :
    _name ( Name ),
    _density ( Density ),
    _diffusionCoefficient ( DiffusionCoefficient ),
    _elastModulus ( ElastModulus ),
    _poissonRatio ( PoissonRatio )
     { 
        this->computeLambda( _elastModulus, _poissonRatio );
        this->computeMu( _elastModulus, _poissonRatio );
        this->computeBulkModulus( _elastModulus, _poissonRatio );
     }
    
   MaterialProperties ( const MaterialProperties<RealType,2>& material ) :
    _name ( material.getName() ),
    _density ( material.getDensity() ),
    _diffusionCoefficient ( material.getDiffusionCoefficient() ),
    _elastModulus ( material.getElastModulus() ),
    _poissonRatio ( material.getPoissonRatio() ),
    _lambda ( material.getLambda() ),
    _mu ( material.getMu() ),
    _bulkModulus ( material.getBulkModulus( ) )
    { }

  // ************** inspectors ************** //
  const std::string getName() const { return _name;}
  RealType getDensity() const {return _density;}
  RealType getDiffusionCoefficient() const {return _diffusionCoefficient;}
  RealType getElastModulus() const {return _elastModulus;}
  RealType getPoissonRatio() const {return _poissonRatio;}
  RealType getLambda() const {return _lambda;}
  RealType getMu() const {return _mu;}
  RealType getBulkModulus() const {return _bulkModulus;}

  void print() const {
    std::cout << "-------- " << getName() << " --------" << std::endl
              << "dimChartDomain            : " << _dimChartDomain << std::endl
              << "density              : " << getDensity() << " kg/m^3" << std::endl
              << "diffusionCoefficient : " << getDiffusionCoefficient() << std::endl
              << "elast modulus   (E)  : " << getElastModulus() << " GPa" << std::endl
              << "poisson ratio   (nu) : " << getPoissonRatio()  << std::endl
              << "1st Lame    (lambda) : " << getLambda() << " GPa" << std::endl
              << "2nd Lame        (mu) : " << getMu() << " GPa" << std::endl
              << "bulk modulus    (K)  : " << getBulkModulus() << " GPa" << std::endl
              << endl;
  }
  
  template <typename ParameterParserType>
  void getInfo( ParameterParserType & parser ) const {
      parser.set( "Material.Name", getName() );
      parser.set( "Material.Density", getDensity() );
      parser.set( "Material.ElastModulus", getElastModulus() );
      parser.set( "Material.PoissonRatio", getPoissonRatio() );
      parser.set( "Material.lambda", getLambda() );
      parser.set( "Material.mu", getMu() );
      parser.set( "Material.bulkModulus", getBulkModulus() );
  }

  void setName ( const std::string name ) { _name = name; }
  void setDensity ( RealType value_KgPerM3 ){ _density = value_KgPerM3;}
  void setDiffusionCoefficient ( RealType value ) { _diffusionCoefficient = value; }
  void setElastModulus ( RealType value_GPa ) {_elastModulus = value_GPa;}
  void setPoissonRatio ( RealType value ){ _poissonRatio = value;}
  
  void computeLambda( const RealType ElastModulus, const RealType PoissonRatio) { 
      _lambda = ElastModulus * PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - PoissonRatio ) );
  }
  void computeMu( const RealType ElastModulus, const RealType PoissonRatio ) { _mu = ElastModulus / ( 2. * ( 1. + PoissonRatio ) ) ;}
  void computeBulkModulus( const RealType ElastModulus, const RealType PoissonRatio ) { 
      _bulkModulus = ElastModulus / ( 2. * ( 1. - PoissonRatio ) ) ;
  }
  
  void set( const MaterialProperties<RealType,2>& material ) {
    setName( material._name );
    setDensity( material._density ); setDiffusionCoefficient( material._diffusionCoefficient );
    setElastModulus( material._elastModulus ); setPoissonRatio( material._poissonRatio );
    computeLambda( _elastModulus, _poissonRatio ); 
    computeMu( _elastModulus, _poissonRatio );
    computeBulkModulus( _elastModulus, _poissonRatio );
  }
  
  template <typename VoigtVecType>
  void getVoigtTensorVec( VoigtVecType &vec ) {
      vec.setZero();
      //(C_1111)
      vec[0] = this->getLambda() + 2. * this->getMu();
  }
  
};




template<typename RealType>
class MaterialProperties<RealType,2> {
protected:
  const int _dimChartDomain = 2;
  std::string _name;
  RealType _density;       // in kg/m^3
  RealType _diffusionCoefficient;
  RealType _elastModulus;  // in GPa
  RealType _poissonRatio;  // ( nu )
  RealType _lambda;        // 1st lame constant in GPa = N/m^2
  RealType _mu;            // 2nd lame constant (shear modulus) in GPa = N/m^2
  RealType _bulkModulus;   // in GPa 

public:
  MaterialProperties() :
    _name ( "DefaultMaterial" ),
    _density ( 0. ), 
    _diffusionCoefficient ( 0. ), 
    _elastModulus ( 0. ), _poissonRatio ( 0. ), 
    _lambda ( 0. ), _mu ( 0. ),
    _bulkModulus ( 0. ) 
    {}

  MaterialProperties ( std::string Name,
             RealType Density,
             RealType DiffusionCoefficient,
             RealType ElastModulus, RealType PoissonRatio ) :
    _name ( Name ),
    _density ( Density ),
    _diffusionCoefficient ( DiffusionCoefficient ),
    _elastModulus ( ElastModulus ),
    _poissonRatio ( PoissonRatio )
     { 
        this->computeLambda( _elastModulus, _poissonRatio );
        this->computeMu( _elastModulus, _poissonRatio );
        this->computeBulkModulus( _elastModulus, _poissonRatio );
     }
    
   MaterialProperties ( const MaterialProperties<RealType,2>& material ) :
    _name ( material.getName() ),
    _density ( material.getDensity() ),
    _diffusionCoefficient ( material.getDiffusionCoefficient() ),
    _elastModulus ( material.getElastModulus() ),
    _poissonRatio ( material.getPoissonRatio() ),
    _lambda ( material.getLambda() ),
    _mu ( material.getMu() ),
    _bulkModulus ( material.getBulkModulus( ) )
    { }

  // ************** inspectors ************** //
  const std::string getName() const { return _name;}
  RealType getDensity() const {return _density;}
  RealType getDiffusionCoefficient() const {return _diffusionCoefficient;}
  RealType getElastModulus() const {return _elastModulus;}
  RealType getPoissonRatio() const {return _poissonRatio;}
  RealType getLambda() const {return _lambda;}
  RealType getMu() const {return _mu;}
  RealType getBulkModulus() const {return _bulkModulus;}

  void print() const {
    std::cout << "-------- " << getName() << " --------" << std::endl
              << "dimChartDomain            : " << _dimChartDomain << std::endl
              << "density              : " << getDensity() << " kg/m^3" << std::endl
              << "diffusionCoefficient : " << getDiffusionCoefficient() << std::endl
              << "elast modulus   (E)  : " << getElastModulus() << " GPa" << std::endl
              << "poisson ratio   (nu) : " << getPoissonRatio()  << std::endl
              << "1st Lame    (lambda) : " << getLambda() << " GPa" << std::endl
              << "2nd Lame        (mu) : " << getMu() << " GPa" << std::endl
              << "bulk modulus    (K)  : " << getBulkModulus() << " GPa" << std::endl
              << endl;
  }
  
  template <typename ParameterParserType>
  void getInfo( ParameterParserType & parser ) const {
      parser.set( "Material.Name", getName() );
      parser.set( "Material.Density", getDensity() );
      parser.set( "Material.ElastModulus", getElastModulus() );
      parser.set( "Material.PoissonRatio", getPoissonRatio() );
      parser.set( "Material.lambda", getLambda() );
      parser.set( "Material.mu", getMu() );
      parser.set( "Material.bulkModulus", getBulkModulus() );
  }

  void setName ( const std::string name ) { _name = name; }
  void setDensity ( RealType value_KgPerM3 ){ _density = value_KgPerM3;}
  void setDiffusionCoefficient ( RealType value ) { _diffusionCoefficient = value; }
  void setElastModulus ( RealType value_GPa ) {_elastModulus = value_GPa;}
  void setPoissonRatio ( RealType value ){ _poissonRatio = value;}
  
  void computeLambda( const RealType ElastModulus, const RealType PoissonRatio) { 
      _lambda = ElastModulus * PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - PoissonRatio ) );
  }
  void computeMu( const RealType ElastModulus, const RealType PoissonRatio ) { _mu = ElastModulus / ( 2. * ( 1. + PoissonRatio ) ) ;}
  void computeBulkModulus( const RealType ElastModulus, const RealType PoissonRatio ) { 
      _bulkModulus = ElastModulus / ( 2. * ( 1. - PoissonRatio ) ) ;
  }
  
  void set( const MaterialProperties<RealType,2>& material ) {
    setName( material._name );
    setDensity( material._density ); setDiffusionCoefficient( material._diffusionCoefficient );
    setElastModulus( material._elastModulus ); setPoissonRatio( material._poissonRatio );
    computeLambda( _elastModulus, _poissonRatio ); 
    computeMu( _elastModulus, _poissonRatio );
    computeBulkModulus( _elastModulus, _poissonRatio );
  }
  
  template <typename VoigtVecType>
  void getVoigtTensorVec( VoigtVecType &vec ) {
      vec.setZero();
      //(C_1111,C_2222,C_1212,C_1122,C_2212,C1112)
      vec[0] = this->getLambda() + 2. * this->getMu();
      vec[1] = this->getLambda() + 2. * this->getMu();
      vec[2] = this->getMu();
      vec[3] = this->getLambda();
  }
  
};



template<typename RealType>
class MaterialProperties<RealType,3> {
protected:
  const int _dimChartDomain = 3;
  std::string _name;
  RealType _density;       // in kg/m^3
  RealType _diffusionCoefficient;
  RealType _elastModulus;  // in GPa
  RealType _poissonRatio;  // ( nu )
  RealType _lambda;        // 1st lame constant in GPa = N/m^2
  RealType _mu;            // 2nd lame constant (shear modulus) in GPa = N/m^2
  RealType _bulkModulus;   // in GPa

public:
  MaterialProperties() :
    _name ( "DefaultMaterial" ),
    _density ( 0. ), 
    _diffusionCoefficient ( 0. ), 
    _elastModulus ( 0. ), _poissonRatio ( 0. ), 
    _lambda ( 0. ), _mu ( 0. ),
    _bulkModulus ( 0. ) 
    {}

  MaterialProperties ( std::string Name,
             RealType Density,
             RealType DiffusionCoefficient,
             RealType ElastModulus, RealType PoissonRatio ) :
    _name ( Name ),
    _density ( Density ),
    _diffusionCoefficient ( DiffusionCoefficient ),
    _elastModulus ( ElastModulus ),
    _poissonRatio ( PoissonRatio )
     { 
        this->computeLambda( _elastModulus, _poissonRatio );
        this->computeMu( _elastModulus, _poissonRatio );
        this->computeBulkModulus( _elastModulus, _poissonRatio );
     }
    
   MaterialProperties ( const MaterialProperties<RealType,3>& material ) :
    _name ( material.getName() ),
    _density ( material.getDensity() ),
    _diffusionCoefficient ( material.getDiffusionCoefficient() ),
    _elastModulus ( material.getElastModulus() ),
    _poissonRatio ( material.getPoissonRatio() ),
    _lambda ( material.getLambda() ),
    _mu ( material.getMu() ),
    _bulkModulus ( material.getBulkModulus( ) )
    { }

  // ************** inspectors ************** //
  const std::string getName() const { return _name;}
  RealType getDensity() const {return _density;}
  RealType getDiffusionCoefficient() const {return _diffusionCoefficient;}
  RealType getElastModulus() const {return _elastModulus;}
  RealType getPoissonRatio() const {return _poissonRatio;}
  RealType getLambda() const {return _lambda;}
  RealType getMu() const {return _mu;}
  RealType getBulkModulus() const {return _bulkModulus;}

  void print() const {
    std::cout << "-------- " << getName() << " --------" << std::endl
              << "dimChartDomain            : " << _dimChartDomain << std::endl
              << "density              : " << getDensity() << " kg/m^3" << std::endl
              << "diffusionCoefficient : " << getDiffusionCoefficient() << std::endl
              << "elast modulus   (E)  : " << getElastModulus() << " GPa" << std::endl
              << "poisson ratio   (nu) : " << getPoissonRatio()  << std::endl
              << "1st Lame    (lambda) : " << getLambda() << " GPa" << std::endl
              << "2nd Lame        (mu) : " << getMu() << " GPa" << std::endl
              << "bulk modulus    (K)  : " << getBulkModulus() << " GPa" << std::endl
              << endl;
  }
  
  template <typename ParameterParserType>
  void getInfo( ParameterParserType & parser ) const {
      parser.set( "Material.Name", getName() );
      parser.set( "Material.Density", getDensity() );
      parser.set( "Material.ElastModulus", getElastModulus() );
      parser.set( "Material.PoissonRatio", getPoissonRatio() );
      parser.set( "Material.lambda", getLambda() );
      parser.set( "Material.mu", getMu() );
      parser.set( "Material.bulkModulus", getBulkModulus() );
  }

  void setName ( const std::string name ) { _name = name; }
  void setDensity ( RealType value_KgPerM3 ){ _density = value_KgPerM3;}
  void setDiffusionCoefficient ( RealType value ) { _diffusionCoefficient = value; }
  void setElastModulus ( RealType value_GPa ) {_elastModulus = value_GPa;}
  void setPoissonRatio ( RealType value ){ _poissonRatio = value;}
  void computeLambda( const RealType ElastModulus, const RealType PoissonRatio) { 
      _lambda = ElastModulus * PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - 2. * PoissonRatio ) );
  }
  void computeMu( const RealType ElastModulus, const RealType PoissonRatio ) { _mu = ElastModulus / ( 2. * ( 1. + PoissonRatio ) ) ;}
  void computeBulkModulus( const RealType ElastModulus, const RealType PoissonRatio ) { 
      const int dimChartDomain = 3;
      _bulkModulus = ElastModulus / ( 3. * ( 1. - 2. * PoissonRatio ) ) ;
  }
  
  void set( const MaterialProperties<RealType,3>& material ) {
    setName( material._name );
    setDensity( material._density ); setDiffusionCoefficient( material._diffusionCoefficient );
    setElastModulus( material._elastModulus ); setPoissonRatio( material._poissonRatio );
    computeLambda( _elastModulus, _poissonRatio ); 
    computeMu( _elastModulus, _poissonRatio );
    computeBulkModulus( _elastModulus, _poissonRatio );
  }
  
  template <typename VoigtVecType>
  void getVoigtTensorVec( VoigtVecType &vec ) {
      vec.setZero();
      //(C_1111,C_2222,C_3333,C_1212,C_1122,C_2212,C1112)
      vec[0] = this->getLambda() + 2. * this->getMu();
      vec[1] = this->getLambda() + 2. * this->getMu();
      vec[2] = this->getLambda() + 2. * this->getMu();
      vec[3] = this->getMu();
      vec[4] = this->getMu();
      vec[5] = this->getMu();
      vec[6] = this->getLambda();
      vec[7] = this->getLambda();
      vec[11] = this->getLambda();
  }
  
};





#define QUOC_APPROXCHARFCT_FOURTHORDER

// here (a,b) = (0,1)
// also typical (a,b) = (-1,1)
template < typename RealType >
class PhaseFieldFunctions {
  public :  
      
//    typedef _ConfiguratorType ConfiguratorType;
//    typedef typename ConfiguratorType::RealType RealType;
//    typedef typename ConfiguratorType::DTContainer DataTypeContainer;
    
//    const ConfiguratorType &_conf;
   
//    private:
//   const RealType _factorDoubleWell;
  
   public:
   
   PhaseFieldFunctions ( ) {}
//     _factorDoubleWell ( 9.0 ) 
//     {}

#ifdef QUOC_APPROXCHARFCT_FOURTHORDER
    RealType approxCharFct_vol ( const RealType v ) const { return v; }
    RealType approxCharFct_vol_Derivative ( const RealType /*v*/ ) const { return 1.0;}
    RealType approxCharFct_vol_SecondDerivative ( const RealType /*v*/ ) const { return 0.0;}
    
    RealType approxCharFct_material ( const RealType v ) const { return std::pow( v,  4); }
    RealType approxCharFct_material_Derivative ( const RealType v ) const { return 4.0 * std::pow( v,  3 );}
    RealType approxCharFct_material_SecondDerivative ( const RealType v ) const { return 12. * std::pow( v,  2);}  
#endif


    RealType doubleWell ( const RealType v ) const { return 9.0 * std::pow( v * (v - 1.0),  2);}
    RealType doubleWellDerivative ( const RealType v ) const { return 9.0 * 2.0 * v * ( v - 1.0 ) * ( 2. * v - 1 );}
    RealType doubleWellSecondDerivative ( const RealType v ) const { return 9.0 * 2.0 * ( 6. * v * v - 6 * v + 1.0 );}
    
};






template < typename DataTypeContainer >
class MaterialOptimizationInfo {
    
public:
  typedef typename DataTypeContainer::RealType RealType;
  typedef pesopt::BoostParser ParameterParserType;
  
  const RealType _factorComplianceCost, _factorVolumeCost, _factorInterfaceCost;
  const RealType _epsInterfaceLength;

  //used e.g. if one is only interested in optimal deformation
  MaterialOptimizationInfo ( const ParameterParserType & parser ) :
    _factorComplianceCost ( parser.template get<double>( "MaterialOptimization.factorComplianceCost" ) ),
    _factorVolumeCost ( parser.template get<double>( "MaterialOptimization.factorVolumeCost" ) ),
    _factorInterfaceCost ( parser.template get<double> ( "MaterialOptimization.factorInterfaceCost" ) ),
    _epsInterfaceLength ( 1. ){
        cout << endl << endl 
             << "WARNING: this constructor of MaterialOptimizationInfo should not be used for material optimization with interface energy, since eps=1" 
             << endl << endl;
    }
    

  MaterialOptimizationInfo ( const ParameterParserType & parser, const RealType epsInterfaceLength ) :
    _factorComplianceCost ( parser.template get<double>( "MaterialOptimization.factorComplianceCost" ) ),
    _factorVolumeCost ( parser.template get<double>( "MaterialOptimization.factorVolumeCost" ) ),
    _factorInterfaceCost ( parser.template get<double> ( "MaterialOptimization.factorInterfaceCost" ) ),
    _epsInterfaceLength ( epsInterfaceLength ) { }

};




template < typename DataTypeContainer>
class MaterialOptimizationInfoMultipleLoad 
: public MaterialOptimizationInfo<DataTypeContainer> {
  
public :
  typedef typename DataTypeContainer::RealType RealType;
  typedef pesopt::BoostParser ParameterParserType;
  
protected:
    WeightFunctionMultipleLoad<RealType> _weightFctLoad;
 
public:
    
  MaterialOptimizationInfoMultipleLoad ( const ParameterParserType & parser, 
                                         const RealType epsInterfaceLength, 
                                         //const RealType epsFactor, 
                                         const std::vector<RealType> &weightFct_weightVec) :
    MaterialOptimizationInfo<DataTypeContainer> ( parser, epsInterfaceLength ),
    _weightFctLoad( parser.template get<RealType>("MaterialOptimization.weightFct_weight"), parser.template get<RealType>("MaterialOptimization.weightFct_p"), parser.template get<RealType>("MaterialOptimization.weightFct_q"), weightFct_weightVec ) {}

  const WeightFunctionMultipleLoad<RealType> &getWeightFunctionLoad( ) const{ return _weightFctLoad;}
    
};





#endif
