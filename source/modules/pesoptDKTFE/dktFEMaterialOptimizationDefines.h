#ifndef __DKTFEMATERIALOPTIMIZATIONDEFINES__H
#define __DKTFEMATERIALOPTIMIZATIONDEFINES__H


template<typename RealType>
class Material {
protected:
  const std::string _name;
  RealType _density;       // in kg/m^3
  RealType _diffusionCoefficient;
  RealType _elastModulus;  // in GPa
  RealType _poissonRatio;  // ( nu )
  RealType _lambda;        // 1st lame constant in GPa = N/m^2
  RealType _mu;            // 2nd lame constant (shear modulus) in GPa = N/m^2

public:
  Material() :
    _name ( "DefaultMaterial" ),
    _density ( 0 ), _diffusionCoefficient ( 0 ), _elastModulus ( 0 ), _poissonRatio ( 0 ), _lambda ( 0 ), _mu ( 0 ) {}

  Material ( std::string Name,
             RealType Density,
             RealType DiffusionCoefficient,
             RealType ElastModulus, RealType PoissonRatio ) :
    _name ( Name ),
    _density ( Density ),
    _diffusionCoefficient ( DiffusionCoefficient ),
    _elastModulus ( ElastModulus ),
    _poissonRatio ( PoissonRatio )
     {
        this->setLambda( _elastModulus, _poissonRatio );
        this->setMu( _elastModulus, _poissonRatio );
     }

   Material ( const Material<RealType>& material ) :
    _name ( material.getName() ),
    _density ( material.getDensity() ),
    _diffusionCoefficient ( material.getDiffusionCoefficient() ),
    _elastModulus ( material.getElastModulus() ),
    _poissonRatio ( material.getPoissonRatio() ),
    _lambda ( material.getLambda() ),
    _mu ( material.getMu() ){ }

  ~Material() {}

  // ************** inspectors ************** //
  const std::string getName() const { return _name;}
  RealType getDensity() const {return _density;}
  RealType getDiffusionCoefficient() const {return _diffusionCoefficient;}
  RealType getElastModulus() const {return _elastModulus;}
  RealType getPoissonRatio() const {return _poissonRatio;}
  RealType getLambda() const {return _lambda;}
  RealType getMu() const {return _mu;}

  void print() const {
    std::cout << "-------- " << getName() << " --------" << std::endl
              << "density              : " << getDensity() << " kg/m^3" << std::endl
              << "diffusionCoefficient : " << getDiffusionCoefficient() << " TODO" << std::endl
              << "elast modulus   (E)  : " << getElastModulus() << " GPa" << std::endl
              << "poisson ratio   (nu) : " << getPoissonRatio()  << std::endl
              << "1st Lame    (lambda) : " << getLambda() << " GPa" << std::endl
              << "2nd Lame        (mu) : " << getMu() << " GPa" << std::endl;
  }

  void setDensity ( RealType value_KgPerM3 ){ _density = value_KgPerM3;}
  void setDiffusionCoefficient ( RealType value ) { _diffusionCoefficient = value; }
  void setElastModulus ( RealType value_GPa ) {_elastModulus = value_GPa;}
  void setPoissonRatio ( RealType value ){ _poissonRatio = value;}
  void setLambda( const RealType ElastModulus, const RealType PoissonRatio) { _lambda = ElastModulus * PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - 2. * PoissonRatio ) );}
  void setMu( const RealType ElastModulus, const RealType PoissonRatio ) { _mu = ElastModulus / ( 2. * ( 1. + PoissonRatio ) ) ;}
  void setLambdaMu( const RealType lambda, const RealType mu ) {
    _lambda = lambda; _mu = mu;
    //TODO: E, nu
  }

  void set( const Material<RealType>& material ) {
    setDensity( material._density ); setDiffusionCoefficient( material._diffusionCoefficient );
    setElastModulus( material._elastModulus ); setPoissonRatio( material._poissonRatio );
    setLambda( _elastModulus, _poissonRatio ); setMu( _elastModulus, _poissonRatio );
  }
};





template<typename RealType>
class materialOptInfo{

public:
    Material<RealType> _HardMaterial, _SoftMaterial;
    const RealType _factorMembraneEnergy, _factorBendingEnergy;
    const RealType _thicknessHard, _thicknessSoft;

    template<typename ParameterParserType>
    materialOptInfo ( const ParameterParserType &parser ) :
    _HardMaterial ( ), _SoftMaterial (  ),
    _factorMembraneEnergy ( parser.template get<double>( "Material.factor_membraneEnergy" ) ),
    _factorBendingEnergy  ( parser.template get<double> ( "Material.factor_bendingEnergy" ) ),
    _thicknessHard ( parser.template get<double>( "Material.thickness_hard" ) ),
    _thicknessSoft ( parser.template get<double> ( "Material.thickness_soft" ) )
    {
       _HardMaterial.set(  Material<RealType> ( "Hard", 1.0, 1.0, parser.template get<RealType>( "Material.ElastModulus_hard" ), parser.template get<RealType>( "Material.PoissonRatio_hard" ) ) );
       _SoftMaterial.set(  Material<RealType> ( "Soft", 1.0, 1.0, parser.template get<RealType>( "Material.ElastModulus_soft" ), parser.template get<RealType>( "Material.PoissonRatio_soft" ) ) );
    }
};









template < typename _ConfiguratorType, typename _ConfiguratorTypePf >
class PhaseFieldFunctions {
  public :

   typedef _ConfiguratorType ConfiguratorType;
   typedef _ConfiguratorTypePf ConfiguratorTypePf;
   typedef typename ConfiguratorType::RealType RealType;
   typedef typename ConfiguratorType::DTContainer DataTypeContainer;

   const ConfiguratorType &_conf;
   const ConfiguratorTypePf &_confpf;
   const int _pfFunctionMaterialType; //TODO as template argument: 1 - linear, 2 - second order, 4- fourth order, -1 - hom
   const int _pfFunctionDoubleWellType; //TODO as template argument:  2 - second order, 4- fourth order
   mutable RealType _factorDoubleWell;

   PhaseFieldFunctions ( const ConfiguratorType &conf, const ConfiguratorTypePf &confpf,
                         const int pfFunctionMaterialType,
                         const int pfFunctionDoubleWellType,
                         const RealType factorDoubleWell = 1 //= 9. / 16.
                       ) :
    _conf ( conf ), _confpf ( confpf ), _pfFunctionMaterialType( pfFunctionMaterialType ), _pfFunctionDoubleWellType( pfFunctionDoubleWellType ),
    _factorDoubleWell ( factorDoubleWell )
    {
        if( _pfFunctionDoubleWellType == 2 ) _factorDoubleWell *= 0.40528473456935108577551;
        if( _pfFunctionDoubleWellType == 4 ) _factorDoubleWell *= 9. / 16.;
    }


    RealType approxCharFct_vol ( const RealType v ) const { return 0.5 * ( v + 1.0 ); }
    RealType approxCharFct_vol_Derivative ( const RealType /*v*/ ) const { return 0.5;}
    RealType approxCharFct_vol_SecondDerivative ( const RealType /*v*/ ) const { return 0.0;}


    RealType approxCharFct_material ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
        const RealType diff = c_hard - c_soft;
        RealType aux = 0.;
        if( std::abs( diff ) < 1.e-16 ){
            aux = c_hard;
        }else{
            switch( _pfFunctionMaterialType ){
                case -1:{
                    RealType tmp = (1. + v) / c_hard + (1. - v) / c_soft;
                    aux = 2. / tmp;
                }break;
                case 1:{
                    RealType chi = 0.5 * ( v + 1.0 );
                    aux = c_hard * chi + c_soft * (1.-chi);
                }break;
                case 4:{
                    RealType chi = pesopt::Sqr( pesopt::Sqr( v + 1.0 ) ) / 16.;
                    aux = c_hard * chi + c_soft * (1.-chi);
                }break;
                default :
                    throw std::invalid_argument( pesopt::strprintf ( "Wrong Type. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
                    break;
            }
        }
        return aux;
    }
    RealType approxCharFct_material_Derivative ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
        const RealType diff = c_hard - c_soft;
        RealType aux = 0.;
        if( std::abs( diff ) < 1.e-16 ){
            aux = 0.;
        }else{
            switch( _pfFunctionMaterialType ){
                case -1:{
                    RealType tmp = (1. + v) / c_hard + (1. - v) / c_soft;
                    aux = -2. * (1. / c_hard - 1. / c_soft ) / pesopt::Sqr( tmp );
                }break;
                case 1:{
                    RealType Dchi = 0.5;
                    aux = Dchi * diff;
                }break;
                case 4:{
                    RealType Dchi = 0.25 * pesopt::Cub( v + 1.0 );
                    aux = Dchi * diff;
                }break;
                default :
                    throw std::invalid_argument( pesopt::strprintf ( "Wrong Type. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
                    break;
            }
        }
        return aux;
    }


    RealType approxCharFct_thicknessSqr ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
        RealType chi = pesopt::Sqr( pesopt::Sqr( v + 1.0 ) ) / 16.;
        return pesopt::Sqr(c_hard) * chi + pesopt::Sqr(c_soft) * (1.-chi);
    }
    RealType approxCharFct_thicknessSqr_Derivative ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
        RealType Dchi = 0.25 * pesopt::Cub( v + 1.0 );
        return Dchi * ( pesopt::Sqr(c_hard) - pesopt::Sqr(c_soft) );
    }


// #ifdef SHELLFE_APPROXCHARFCT_FIRSTORDER
//
//     RealType approxCharFct_material ( const RealType v ) const { return 0.5 * ( v + 1.0 ); }
//     RealType approxCharFct_material_Derivative ( const RealType /*v*/ ) const { return 0.5;}
//     RealType approxCharFct_material_SecondDerivative ( const RealType /*v*/ ) const { return 0.0;}
//
//     RealType approxCharFct_material ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
//         RealType chi = approxCharFct_material( v );
//         return c_hard * chi + c_soft * (1.-chi);
//     }
//     RealType approxCharFct_material_Derivative ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
//         RealType tmp = approxCharFct_material_Derivative(v);
//         return (c_hard - c_soft) * tmp;
//     }
//
// #endif //SHELLFE_APPROXCHARFCT_SECONDORDER

// #ifdef SHELLFE_APPROXCHARFCT_SECONDORDER
//
//     RealType approxCharFct_material ( const RealType v ) const { return 0.25 * pesopt::Sqr( v + 1.0 ); }
//     RealType approxCharFct_material_Derivative ( const RealType v ) const { return 0.5 * ( v + 1.0 );}
//     RealType approxCharFct_material_SecondDerivative ( const RealType /*v*/ ) const { return 0.5;}
//
//     RealType approxCharFct_material ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
//         RealType chi = approxCharFct_material( v );
//         return c_hard * chi + c_soft * (1.-chi);
//     }
//     RealType approxCharFct_material_Derivative ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
//         RealType tmp = approxCharFct_material_Derivative(v);
//         return (c_hard - c_soft) * tmp;
//     }
//
// #endif //SHELLFE_APPROXCHARFCT_SECONDORDER

// #ifdef SHELLFE_APPROXCHARFCT_FOURTHORDER
//
//     RealType approxCharFct_material ( const RealType v ) const { return pesopt::Sqr( pesopt::Sqr( v + 1.0 ) ) / 16.; }
//     RealType approxCharFct_material_Derivative ( const RealType v ) const { return 0.25 * pesopt::Cub( v + 1.0 );}
//     RealType approxCharFct_material_SecondDerivative ( const RealType v ) const { return 0.75 * pesopt::Sqr( v + 1.0 );}
//
//     RealType approxCharFct_material ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
//         RealType chi = approxCharFct_material( v );
//         return c_hard * chi + c_soft * (1.-chi);
//     }
//     RealType approxCharFct_material_Derivative ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
//         RealType tmp = approxCharFct_material_Derivative(v);
//         return (c_hard - c_soft) * tmp;
//     }
// #endif //SHELLFE_APPROXCHARFCT_FOURTHORDER


// #ifdef SHELLFE_APPROXCHARFCT_1DHOMOGENIZATION
//
//     RealType approxCharFct_material ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
//         RealType tmp = (1. + v) / c_hard + (1. - v) / c_soft;
//         return 2. / tmp;
//     }
//     RealType approxCharFct_material_Derivative ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
//         RealType tmp = (1. + v) / c_hard + (1. - v) / c_soft;
//         return -2. * (1. / c_hard - 1. / c_soft ) / pesopt::Sqr( tmp );
//     }
//     RealType approxCharFct_material_SecondDerivative ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
//         RealType tmp = (1. + v) / c_hard + (1. - v) / c_soft;
//         return 4. * pesopt::Sqr(1. / c_hard - 1. / c_soft ) / pesopt::Cub( tmp );
//     }
//
//     //TODO: for more than one paramter:
// //     RealType approxCharFct_material ( const RealType v ) const { return 0.25 * pesopt::Sqr( v + 1.0 ); }
// //     RealType approxCharFct_material_Derivative ( const RealType v ) const { return 0.5 * ( v + 1.0 );}
// //     RealType approxCharFct_material_SecondDerivative ( const RealType /*v*/ ) const { return 0.5;}
// #endif //SHELLFE_APPROXCHARFCT_1DHOMOGENIZATION


// #ifdef SHELLFE_APPROXCHARFCT_WIRTH
//     RealType approxCharFct_vol ( const RealType v ) const { return 0.25 * pesopt::Sqr( v + 1.0 ); }
//     RealType approxCharFct_vol_Derivative ( const RealType v ) const { return 0.5 * ( v + 1.0 );}
//     RealType approxCharFct_vol_SecondDerivative ( const RealType /*v*/ ) const { return 0.5;}
//
//     RealType approxCharFct_material ( const RealType v ) const { return 0.25 * pesopt::Sqr( v + 1.0 ); }
//     RealType approxCharFct_material_Derivative ( const RealType v ) const { return 0.5 * ( v + 1.0 );}
//     RealType approxCharFct_material_SecondDerivative ( const RealType /*v*/ ) const { return 0.5;}
//
//     RealType approxCharFct_material ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
//         RealType chi = approxCharFct_material( v );
//         return c_hard * chi + c_soft * (1.-chi);
//     }
//     RealType approxCharFct_material_Derivative ( const RealType v, const RealType c_hard, const RealType c_soft ) const {
//         RealType tmp = approxCharFct_material_Derivative(v);
//         return (c_hard - c_soft) * tmp;
//     }
// #endif //SHELLFE_APPROXCHARFCT_WIRTH




    RealType doubleWell ( const RealType v ) const {
        RealType aux = 0.;
        switch( _pfFunctionDoubleWellType ){
            case 2:{
                aux = _factorDoubleWell * (1. - v * v);
            }break;
            case 4:{
                aux = _factorDoubleWell * pesopt::Sqr( v * v - 1.0 );
            }break;
            default :
                throw std::invalid_argument( pesopt::strprintf ( "Wrong Type. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
                break;
        }
        return aux;
    }
    RealType doubleWellDerivative  ( const RealType v ) const {
        RealType aux = 0.;
        switch( _pfFunctionDoubleWellType ){
            case 2:{
                aux = _factorDoubleWell * -2. * v;
            }break;
            case 4:{
                aux = _factorDoubleWell * 4.0 * ( v * v - 1.0 ) * v;
            }break;
            default :
                throw std::invalid_argument( pesopt::strprintf ( "Wrong Type. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
                break;
        }
        return aux;
    }
    RealType doubleWellSecondDerivative  ( const RealType v ) const {
        RealType aux = 0.;
        switch( _pfFunctionDoubleWellType ){
            case 2:{
                aux = _factorDoubleWell * -2.;
            }break;
            case 4:{
                aux = _factorDoubleWell * 4.0 * ( 3. * v * v - 1.0 );
            }break;
            default :
                throw std::invalid_argument( pesopt::strprintf ( "Wrong Type. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
                break;
        }
        return aux;
    }



//     RealType doubleWell ( const RealType v ) const { return _factorDoubleWell * pesopt::Sqr( v * v - 1.0 );}
//     RealType doubleWellDerivative ( const RealType v ) const { return _factorDoubleWell * 4.0 * ( v * v - 1.0 ) * v;}
//     RealType doubleWellSecondDerivative ( const RealType v ) const { return _factorDoubleWell * 4.0 * ( 3. * v * v - 1.0 );}

};



template < typename _ConfiguratorType, typename _ConfiguratorTypePf >
class ShellWithMaterialConfigurator : public PhaseFieldFunctions<_ConfiguratorType, _ConfiguratorTypePf> {

public :

  typedef _ConfiguratorType ConfiguratorType;
  typedef _ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::DTContainer DataTypeContainer;
  typedef typename ConfiguratorType::RealType RealType;
  typedef pesopt::BoostParser ParameterParserType;

  const materialOptInfo<RealType> _materialInfo;

  //used e.g. if one is only interested in optimal deformation
  ShellWithMaterialConfigurator ( const ParameterParserType &parser, const ConfiguratorType &conf, const ConfiguratorTypePf &confpf ) :
    PhaseFieldFunctions<_ConfiguratorType, _ConfiguratorTypePf> ( conf, confpf, parser.template get<int>("Material.pfFunctionMaterialType"), parser.template get<int>("Material.pfFunctionDoubleWellType"), parser.template get<double> ("MaterialOptimization.factorDoubleWell" ) ), _materialInfo ( parser ){ }
};



template < typename _ConfiguratorType, typename _ConfiguratorTypePf >
class MaterialOptimizationShellConfigurator : public PhaseFieldFunctions<_ConfiguratorType, _ConfiguratorTypePf> {

public :

  typedef _ConfiguratorType ConfiguratorType;
  typedef _ConfiguratorTypePf ConfiguratorTypePf;
  typedef typename ConfiguratorType::DTContainer DataTypeContainer;
  typedef typename ConfiguratorType::RealType RealType;
  typedef pesopt::BoostParser ParameterParserType;

  const materialOptInfo<RealType> _materialInfo;
  RealType _factorComplianceCost, _factorAreaCost, _factorInterfaceCost;
  RealType _areaConstraintLowerBound, _areaConstraintUpperBound;
  RealType _interfaceConstraintLowerBound, _interfaceConstraintUpperBound;
  const RealType _epsInterfaceLength;
  const RealType _epsFactor; //epsInterfaceLength = epsFactor * gridsize

  //used e.g. if one is only interested in optimal deformation
  MaterialOptimizationShellConfigurator ( const ParameterParserType &parser, const ConfiguratorType &conf, const ConfiguratorTypePf &confpf ) :
    PhaseFieldFunctions<_ConfiguratorType, _ConfiguratorTypePf> ( conf, confpf, parser.template get<int>("Material.pfFunctionMaterialType"), parser.template get<int>("Material.pfFunctionDoubleWellType"), parser.template get<double> ("MaterialOptimization.factorDoubleWell" ) ),
    _materialInfo ( parser ),
    _factorComplianceCost ( parser.template get<double>( "MaterialOptimization.factorComplianceCost" ) ),
    _factorAreaCost ( parser.template get<double>( "MaterialOptimization.factorAreaCost" ) ),
    _factorInterfaceCost ( parser.template get<double> ( "MaterialOptimization.factorInterfaceCost" ) ),
    _epsInterfaceLength ( 1. ), _epsFactor( 1.){
        cout << endl << endl
             << "WARNING: this constructor of MaterialOptimizationShellConfigurator should not be used for material optimization, since eps=1 and factorDoubleWell = default"
             << endl << endl;
    }

  MaterialOptimizationShellConfigurator ( const ParameterParserType &parser, const ConfiguratorType &conf, const ConfiguratorTypePf &confpf, const RealType epsInterfaceLength, const RealType epsFactor ) :
    PhaseFieldFunctions<_ConfiguratorType, _ConfiguratorTypePf> ( conf, confpf, parser.template get<int>("Material.pfFunctionMaterialType"), parser.template get<int>("Material.pfFunctionDoubleWellType"), parser.template get<double> ("MaterialOptimization.factorDoubleWell" ) ),
    _materialInfo ( parser ),
    _factorComplianceCost ( parser.template get<double>( "MaterialOptimization.factorComplianceCost" ) ),
    _factorAreaCost ( parser.template get<double>( "MaterialOptimization.factorAreaCost" ) ),
    _factorInterfaceCost ( parser.template get<double> ( "MaterialOptimization.factorInterfaceCost" ) ),
    _epsInterfaceLength ( epsInterfaceLength ), _epsFactor ( epsFactor ) {}
};




template <typename RealType>
class IsometryInfo {

public :

    RealType _gridSize;

    RealType _isometryErrorL1, _isometryErrorL2;
    RealType _GaussCurvatureL1, _ConvGaussCurvatureL1, _GaussCurvatureL1Diff, _GaussCurvatureInt;

    RealType _approxD2uL2, _D2uL2, _relativeShapeOpL2;

    RealType _errorApproxD2uToFineSolutionL2;
//     RealType _errorApproxD2uToFineSolutionL2_EOC;

    IsometryInfo ( ) {}

    void setGridSize( const RealType e ) { _gridSize = e; };

    void setIsometryErrorL1( const RealType e ) { _isometryErrorL1 = e; };
    void setIsometryErrorL2( const RealType e ) { _isometryErrorL2 = e; };
    void setGaussCurvatureL1( const RealType e ) { _GaussCurvatureL1 = e; };
    void setGaussCurvatureInt( const RealType e ) { _GaussCurvatureInt = e; };
    void setGaussCurvatureL1Diff( const RealType e ) { _GaussCurvatureL1Diff = e; };
    void setConvGaussCurvatureL1( const RealType e ) { _ConvGaussCurvatureL1 = e; };

    void setApproxD2uL2( const RealType e ) { _approxD2uL2 = e; }
    void setD2uL2( const RealType e ) { _D2uL2 = e; }
    void setRelativeShapeOpL2( const RealType e ) { _relativeShapeOpL2 = e; }

    void setErrorApproxD2uToFineSolutionL2( const RealType e ) { _errorApproxD2uToFineSolutionL2 = e; };

    template<typename ParameterParserType>
    void saveToFile( const string fileName, const string saveDirectory ) const {
        ParameterParserType resultsParser;

        resultsParser.set( "saving.saveDirectory", saveDirectory );

        resultsParser.set( "Isometry.gridSize", _gridSize );
        resultsParser.set( "Isometry.isometryErrorL1", _isometryErrorL1 );
        resultsParser.set( "Isometry.isometryErrorL2", _isometryErrorL2 );
        resultsParser.set( "Isometry.GaussCurvatureL1", _GaussCurvatureL1 );
        resultsParser.set( "Isometry.GaussCurvatureInt", _GaussCurvatureInt );
        resultsParser.set( "Isometry.GaussCurvatureL1Diff", _GaussCurvatureL1Diff );
        resultsParser.set( "Isometry.ConvGaussCurvatureL1", _ConvGaussCurvatureL1 );

        resultsParser.set( "Isometry.ApproxD2uL2", _approxD2uL2 );
        resultsParser.set( "Isometry.D2uL2", _D2uL2 );
        resultsParser.set( "Isometry.RelativeShapeOpL2", _relativeShapeOpL2 );

        resultsParser.set( "Isometry.ErrorApproxD2uToFineSolutionL2", _errorApproxD2uToFineSolutionL2 );

        resultsParser.saveToFile( fileName + ".ini" );
    }

};




template <typename RealType>
class DeformationOptimizationShellEnergyInfo {

public :

    RealType _potentialEnergy, _storedElasticEnergy, _dissipationEnergy;
    RealType _membraneEnergy, _bendingEnergy;
    RealType _residualConstraint;
    RealType _L2Norm, _LInfNormAtQuadPoints, _LInfNormAtNodes;

    DeformationOptimizationShellEnergyInfo ( ) {}

    void setComplianceEnergies( const RealType e1, const RealType e2, const RealType e3 ) { _potentialEnergy = e1; _storedElasticEnergy = e2; _dissipationEnergy = e3; };
    void setMembraneAndBendingEnergy( const RealType mem, const RealType ben ) { _membraneEnergy = mem; _bendingEnergy = ben; };

    void setResidualConstraint( const RealType tmp ) { _residualConstraint = tmp; };

    void setL2Norm( const RealType e ) { _L2Norm = e; };
    void setLInfNormAtQuadPoints( const RealType e ) { _LInfNormAtQuadPoints = e; };
    void setLInfNormAtNodes( const RealType e ) { _LInfNormAtNodes = e; };


    template<typename ParameterParserType>
    void saveToFile( const string fileName, const string saveDirectory ) const {
        ParameterParserType resultsParser;

        resultsParser.set( "saving.saveDirectory", saveDirectory );
        resultsParser.set( "Energy.potentialEnergy", _potentialEnergy );
        resultsParser.set( "Energy.storedElasticEnergy", _storedElasticEnergy );
        resultsParser.set( "Energy.dissipationEnergy", _dissipationEnergy );
        resultsParser.set( "Energy.membraneEnergy", _membraneEnergy );
        resultsParser.set( "Energy.bendingEnergy", _bendingEnergy );
        resultsParser.set( "Energy.residual", _residualConstraint );
        resultsParser.set( "Energy.L2Norm", _L2Norm );
        resultsParser.set( "Energy.LInfNormAtNodes", _LInfNormAtNodes );
        resultsParser.set( "Energy.LInfNormAtQuadPoints", _LInfNormAtQuadPoints );

        resultsParser.saveToFile( fileName + ".ini" );
    }


    template<typename ParameterParserType>
    void loadFromFile( const string fileName, const string saveDirectory )  {
        ParameterParserType parser ( pesopt::strprintf( "%s/%s.ini", saveDirectory.c_str(), fileName.c_str() ) );
        _potentialEnergy  = parser.template get<RealType>( "Energy.potentialEnergy" );
        _storedElasticEnergy = parser.template get<RealType>( "Energy.storedElasticEnergy" );
        _dissipationEnergy = parser.template get<RealType>( "Energy.dissipationEnergy" );
        _membraneEnergy = parser.template get<RealType>( "Energy.membraneEnergy" );
        _bendingEnergy = parser.template get<RealType>( "Energy.bendingEnergy" );
        _residualConstraint = parser.template get<RealType>( "Energy.residual" );
        _L2Norm = parser.template get<RealType>( "Energy.L2Norm" );
        _LInfNormAtNodes = parser.template get<RealType>( "Energy.LInfNormAtNodes" );
        _LInfNormAtQuadPoints = parser.template get<RealType>( "Energy.LInfNormAtQuadPoints" );
    }

};


template <typename MatOptConfigurator, bool IsometryConstraint>
class MaterialOptimizationShellEnergyInfo {

    typedef typename MatOptConfigurator::RealType RealType;

public :

    RealType _complianceEnergy, _area, _interfaceEnergy;
    RealType _residualCostFunctional;
    RealType _interfaceEnergyDirichletPart, _interfaceEnergyDoubleWellPart;
    RealType _potentialEnergy, _storedElasticEnergy, _dissipationEnergy;
    RealType _membraneEnergy, _bendingEnergy;
    RealType _residualConstraint;

    RealType _isometryConstraintIntegratedDiff, _isometryConstraintIntegratedL2;

    MaterialOptimizationShellEnergyInfo ( ) {}

    void setComplianceEnergy( const RealType complianceEnergy ) { _complianceEnergy = complianceEnergy; };
    void setInterfaceEnergy( const RealType interfaceEnergy ) { _interfaceEnergy = interfaceEnergy; };
    void setInterfaceEnergyDirichletPart( const RealType interfaceEnergy ) { _interfaceEnergyDirichletPart = interfaceEnergy; };
    void setInterfaceEnergyDoubleWellPart( const RealType interfaceEnergy ) { _interfaceEnergyDoubleWellPart = interfaceEnergy; };
    void setAreaEnergy( const RealType area ) { _area = area; };
    void setResidualCostFunctional( const RealType tmp ) { _residualCostFunctional = tmp; };

    void setComplianceEnergies( const RealType e1, const RealType e2, const RealType e3 ) { _potentialEnergy = e1; _storedElasticEnergy = e2; _dissipationEnergy = e3; };
    void setMembraneAndBendingEnergy( const RealType mem, const RealType ben ) { _membraneEnergy = mem; _bendingEnergy = ben; };

    void setIsometryConstaintIntegratedDiff( const RealType tmp ) { _isometryConstraintIntegratedDiff = tmp; };
    void setIsometryConstaintIntegratedL2( const RealType tmp ) { _isometryConstraintIntegratedL2 = tmp; };

    void setResidualConstraint( const RealType tmp ) { _residualConstraint = tmp; };
};


#endif
