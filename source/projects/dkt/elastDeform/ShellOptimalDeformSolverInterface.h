#ifndef __OPTIMALDEFORMSOLVERINTERFACE_H
#define __OPTIMALDEFORMSOLVERINTERFACE_H

#include <pesopt_IO.h>
#include <pesopt_DKT.h>
#include <energyDefines.h>


template <typename MatOptConfigurator>
class optimalDeformSolverInterface {
protected :
  typedef typename MatOptConfigurator::ConfiguratorType        ConfiguratorType;
  typedef typename MatOptConfigurator::ConfiguratorTypePf      ConfiguratorTypePf;
  typedef typename ConfiguratorType::DTContainer               DataTypeContainer;
  typedef typename ConfiguratorType::RealType                  RealType;
  typedef typename ConfiguratorType::MaskType                  MaskType;
  typedef typename ConfiguratorType::InitType                  MeshType;
  typedef typename ConfiguratorType::VectorType                VectorType;
  typedef typename ConfiguratorType::SparseMatrixType          SparseMatrixType;
  typedef pesopt::BoostParser ParameterParserType;
  
  const ParameterParserType &_parser;
  const MatOptConfigurator &_matOpConf;
  const ConfiguratorType & _conf;
  const int _numVertices, _numGlobalDofsDisp;

  mutable VectorType _pf;
  const MaskType & _mask;
  
  mutable VectorType _rhs;
  mutable VectorType _solDisp;
  
  mutable RealType _energy_Membrane, _energy_Bending, _energy_Pot;
  mutable int _outputLevel;

  public:

      optimalDeformSolverInterface( const ParameterParserType &parser,
                                    const MatOptConfigurator & matOpConf,
                                    const VectorType & Phasefield,
                                    const MaskType &mask,
                                    const VectorType & rhs_force ) :
         _parser ( parser ),
         _matOpConf ( matOpConf ),
         _conf ( matOpConf._conf ),
         _numVertices( matOpConf._conf.getInitializer().getNumVertices() ), _numGlobalDofsDisp( matOpConf._conf.getNumGlobalDofs() ), 
         _pf ( Phasefield ),
         _mask( mask ),
         _rhs ( rhs_force ),
         _solDisp ( 3*_conf.getNumGlobalDofs() ),
         _outputLevel( parser.template get<int> ( "ConstraintProblem.outputLevel" ) )
        { 
             this->_energy_Membrane = 0.0; this->_energy_Bending = 0.0; this->_energy_Pot = 0.; //this->_energy_Total = 0.;
        }

    const ParameterParserType& getParser( ) const { return _parser; }
    const MatOptConfigurator & getMatOptConfigurator() const { return _matOpConf;}

    const VectorType & getCurrentPf ( ) const { return _pf; };
    const VectorType& getRHS ( ) const{  return _rhs;  }
    const VectorType& getSolutionDisplacement() const { return _solDisp;}
    void setOutputLevel( const int outputLevel ) const {_outputLevel = outputLevel;};
    
    void getComplianceEnergies ( RealType & energy_P, RealType & energy_Membrane, RealType & energy_Bending, RealType & energy_Elastic, RealType &energy_dissipation ) const{
        energy_P = _energy_Pot;
        energy_Membrane = _energy_Membrane;
        energy_Bending = _energy_Bending;
        energy_Elastic = _energy_Membrane + _energy_Bending;
        energy_dissipation = 2. * ( _energy_Pot - energy_Elastic );
    }
    
    template<bool IsometryConstraint>
    void getPartialEnergyInfo(  DeformationOptimizationShellEnergyInfo<RealType> & energyInfo ) const {
        const RealType energy_Elastic = _energy_Membrane + _energy_Bending;
        const RealType energy_dissipation = 2. * ( _energy_Pot - energy_Elastic );
        energyInfo.setComplianceEnergies( _energy_Pot, energy_Elastic, energy_dissipation  );
        energyInfo.setMembraneAndBendingEnergy( _energy_Membrane, _energy_Bending );
   }
};




template <typename MatOptConfigurator, bool LinElastEnergy, typename ElastEnergyType, bool IsometryConstraint>
class optimalDeformSolver { };



// template <typename MatOptConfigurator, bool LinElastEnergy, typename ElastEnergyType, bool IsometryConstraint>
// class optimalDeformSolverAdaptive { };

#endif
