#ifndef __SHELLBUCKLINGOPTIMALDEFORMSOLVERINTERFACE_H
#define __SHELLBUCKLINGOPTIMALDEFORMSOLVERINTERFACE_H

#include <pesopt_IO.h>
#include <pesopt_DKT.h>
#include <pesopt_VTK.h>

#include "ShellIntegrals.h"
#include "ShellForces.h"
#include "ShellIsometryConstraint.h"

#include "ShellBucklingTikzPlotter.h"


//TODO so far for DKT-Deformations
template<typename MatOptConfigurator>
class EOCInterface {

 protected:
  typedef typename MatOptConfigurator::ConfiguratorType         ConfiguratorType;
  typedef typename MatOptConfigurator::ConfiguratorTypePf       ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType                   RealType;
  typedef typename ConfiguratorType::MaskType                   MaskType;
  typedef typename ConfiguratorType::InitType                   MeshType;
  typedef typename ConfiguratorTypePf::DTContainer              DataTypeContainer;
  typedef typename ConfiguratorType::VectorType                 VectorType;
  typedef pesopt::BoostParser ParameterParserType;

  const ParameterParserType &_parser;

public:
  EOCInterface ( const ParameterParserType &parser ) : _parser( parser ) { }

 template<DiscreteFunctionCacheType _DiscreteFunctionCacheType>
 void computeEOCInL2( const MeshType &fineMesh, const std::vector<RealType> &gridSizeVec, std::vector<VectorType> &solutionVector ) const{

        pesopt::consoleOutput( "EOCs in L2 compared to fine solution" );
        const int numAdaptiveRefinementSteps = solutionVector.size() - 1;
        VectorType fineSolution = solutionVector[numAdaptiveRefinementSteps];
        std::vector<VectorType> solutionVectorProlongated;
        std::vector<RealType> ErrorL2_FineSolutionVec(numAdaptiveRefinementSteps+1), EOCL2_FineSolutionVec(numAdaptiveRefinementSteps);
        for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
            VectorType currentSolutionProlongated = solutionVector[refinementStep];
            fineMesh.prolongateVectorFct_DKTLinearly( currentSolutionProlongated );
            solutionVectorProlongated.push_back( currentSolutionProlongated );
        }
        ConfiguratorType conf ( fineMesh );
        ConfiguratorTypePf confpf ( fineMesh );
        MatOptConfigurator matOptConf ( _parser, conf, confpf );
        ShellHandler<ConfiguratorType,_DiscreteFunctionCacheType> shellHandler ( _parser, conf );
        VectorType material( confpf.getNumGlobalDofs() ); for( int i=0; i<material.size(); ++i ) material[i] = 1.;
        for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
            VectorType diff ( fineSolution.size() ); diff = fineSolution - solutionVectorProlongated[refinementStep];
            RealType tmp = 0.;
            PfWeightedLpNormToP<MatOptConfigurator> ( matOptConf, shellHandler.getChartToUndeformedShell_Cache(), diff, material, 2. ).assembleAdd( tmp );
            ErrorL2_FineSolutionVec[refinementStep] = sqrt( tmp );
        }
        for( int refinementStep=1; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
            RealType facGridSize = gridSizeVec[refinementStep]/gridSizeVec[refinementStep-1];
            EOCL2_FineSolutionVec[refinementStep-1] = std::log( ErrorL2_FineSolutionVec[refinementStep] / ErrorL2_FineSolutionVec[refinementStep-1] ) / std::log(facGridSize);
        }
        pesopt::consoleOutput( "EOC in L^2" );
        for( int refinementStep=1; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ) cout << EOCL2_FineSolutionVec[refinementStep-1] << endl;
        pesopt::consoleOutput( "Diff in L^2" );
        for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ) cout << ErrorL2_FineSolutionVec[refinementStep] << endl;
 }

 template<DiscreteFunctionCacheType _DiscreteFunctionCacheType>
 void computeEOCInLInfty( const MeshType &fineMesh, const std::vector<RealType> &gridSizeVec, std::vector<VectorType> &solutionVector ) const{
        pesopt::consoleOutput( "EOCs in Linf compared to fine solution" );
        const int numAdaptiveRefinementSteps = solutionVector.size() - 1;
        VectorType fineSolution = solutionVector[numAdaptiveRefinementSteps];
        std::vector<VectorType> solutionVectorProlongated;
        std::vector<RealType> ErrorLInfty_FineSolutionVec(numAdaptiveRefinementSteps+1), EOCLinfty_FineSolutionVec(numAdaptiveRefinementSteps);
        for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
            VectorType currentSolutionProlongated = solutionVector[refinementStep];
            fineMesh.prolongateVectorFct_DKTLinearly( currentSolutionProlongated );
            solutionVectorProlongated.push_back( currentSolutionProlongated );
        }
        ConfiguratorType conf ( fineMesh ); ConfiguratorTypePf confpf ( fineMesh ); MatOptConfigurator matOptConf ( _parser, conf, confpf );
        ShellHandler<ConfiguratorType,_DiscreteFunctionCacheType> shellHandler ( _parser, conf );
        VectorType material( confpf.getNumGlobalDofs() ); for( int i=0; i<material.size(); ++i ) material[i] = 1.;
        for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
            VectorType diff ( fineSolution.size() ); diff = fineSolution - solutionVectorProlongated[refinementStep];
            RealType tmp = 0.;
            LInftyNorm<MatOptConfigurator> ( matOptConf, shellHandler.getChartToUndeformedShell_Cache(), diff  ).assembleAdd( tmp );
            ErrorLInfty_FineSolutionVec[refinementStep] = tmp;
        }
        for( int refinementStep=1; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
            RealType facGridSize = gridSizeVec[refinementStep]/gridSizeVec[refinementStep-1];
            EOCLinfty_FineSolutionVec[refinementStep-1] = std::log( ErrorLInfty_FineSolutionVec[refinementStep] / ErrorLInfty_FineSolutionVec[refinementStep-1] ) / std::log(facGridSize);
        }
        pesopt::consoleOutput( "EOC in L^Inf" );
        for( int refinementStep=1; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ) cout << EOCLinfty_FineSolutionVec[refinementStep-1] << endl;
        pesopt::consoleOutput( "Diff in L^Inf" );
        for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ) cout << ErrorLInfty_FineSolutionVec[refinementStep] << endl;
 }

 template<DiscreteFunctionCacheType _DiscreteFunctionCacheType>
 void computeEOCInLInftyNodalwise( const MeshType &fineMesh, const std::vector<RealType> &gridSizeVec, std::vector<VectorType> &solutionVector ) const{

        pesopt::consoleOutput( "EOCs in Linf nodalwise compared to fine solution" );
        const int numAdaptiveRefinementSteps = solutionVector.size() - 1;
        VectorType fineSolution = solutionVector[numAdaptiveRefinementSteps];
        std::vector<VectorType> solutionVectorProlongated;
        std::vector<RealType> ErrorLInftyNodalWise_FineSolutionVec(numAdaptiveRefinementSteps+1), EOCLinftyNodalWise_FineSolutionVec(numAdaptiveRefinementSteps);
        for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
            VectorType currentSolutionProlongated = solutionVector[refinementStep];
            fineMesh.prolongateVectorFct_DKTLinearly( currentSolutionProlongated );
            solutionVectorProlongated.push_back( currentSolutionProlongated );
        }
        ConfiguratorType conf ( fineMesh ); ConfiguratorTypePf confpf ( fineMesh ); MatOptConfigurator matOptConf ( _parser, conf, confpf );
        ShellHandler<ConfiguratorType,_DiscreteFunctionCacheType> shellHandler ( _parser, conf );
        VectorType material( confpf.getNumGlobalDofs() ); for( int i=0; i<material.size(); ++i ) material[i] = 1.;
        for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
            VectorType diff ( fineSolution.size() ); diff = fineSolution - solutionVectorProlongated[refinementStep];
            RealType dispLInftyAtNodes = 0.;
            for( int nodeIdx=0; nodeIdx<fineMesh.getNumVertices(); ++nodeIdx ){
                for( int comp=0; comp < 3; ++comp ){
                    RealType abs = std::abs( diff[nodeIdx + comp * conf.getNumGlobalDofs()] );
                    if( abs > dispLInftyAtNodes ) dispLInftyAtNodes = abs;
                }
            }
            ErrorLInftyNodalWise_FineSolutionVec[refinementStep] = dispLInftyAtNodes;
        }
        for( int refinementStep=1; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
            RealType facGridSize = gridSizeVec[refinementStep]/gridSizeVec[refinementStep-1];
            EOCLinftyNodalWise_FineSolutionVec[refinementStep-1] = std::log( ErrorLInftyNodalWise_FineSolutionVec[refinementStep] / ErrorLInftyNodalWise_FineSolutionVec[refinementStep-1] ) / std::log(facGridSize);
        }
        pesopt::consoleOutput( "EOC in L^Inf" );
        for( int refinementStep=1; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ) cout << EOCLinftyNodalWise_FineSolutionVec[refinementStep-1] << endl;
        pesopt::consoleOutput( "Diff in L^Inf" );
        for( int refinementStep=0; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ) cout << ErrorLInftyNodalWise_FineSolutionVec[refinementStep] << endl;
 }

};


template<typename DataTypeContainer,ShellFEType MaterialFEType>
class optimalDeformPlotter {

 protected:
  typedef typename DataTypeContainer::VectorType                VectorType;
  typedef pesopt::BoostParser ParameterParserType;

  const ParameterParserType &_parser;
  mutable VTKDataSupp _supp;
  const string _VTKFileType;
  const string _latexType;

public:
  optimalDeformPlotter ( const ParameterParserType &Parser ) :
  _parser( Parser ),
  _VTKFileType ( Parser.template get<string>("saving.VTKFileType") ),
  _latexType ( Parser.template get<string>("saving.LATEXType") ){
      if( (MaterialFEType == NodalValuedDofs) || (MaterialFEType == C1Dofs) ) _supp = VERTEX_DATA;
      if( MaterialFEType == ElementValuedDofs ) _supp = FACE_DATA;
 }


  void plotStress( const string saveDirectory, const string fileNameData, const std::vector<string> &nameStressOnElementVec ) const {
    ParameterParserType parserVTKStress     ( _parser.template get<string> ( "InputMesh.parserFileVTKPlotStress" ) );
    ShellVTKPlotter<DataTypeContainer> shellVtkPlotter( saveDirectory, this->_VTKFileType );
    shellVtkPlotter.plotStress( parserVTKStress, fileNameData, nameStressOnElementVec );
  }

  void plotResults( const string saveDirectory, const std::vector<string> &nameStressOnElementVec ) const {
    ParameterParserType parserVTKChart      ( _parser.template get<string> ( "InputMesh.parserFileVTKPlotChart" ) );
    ParameterParserType parserVTKUndeformed ( _parser.template get<string> ( "InputMesh.parserFileVTKPlotUndeformed" ) );
    ParameterParserType parserVTKDeformed   ( _parser.template get<string> ( "InputMesh.parserFileVTKPlotDeformed" ) );
    ParameterParserType parserVTKStress     ( _parser.template get<string> ( "InputMesh.parserFileVTKPlotStress" ) );
    ShellVTKPlotter<DataTypeContainer> shellVtkPlotter( saveDirectory, this->_VTKFileType );
    shellVtkPlotter.plotShellWithData( "chart", parserVTKChart, "material", _supp );
    shellVtkPlotter.plotShellWithData( "undeformedShell", parserVTKUndeformed, "material", _supp );
    shellVtkPlotter.plotShellWithData( "deformedShell", parserVTKDeformed, "material", _supp );
// TODO
//     const string fileNameData = "Stress";
//     shellVtkPlotter.plotStress( parserVTKStress, fileNameData, nameStressOnElementVec );
  }

};


// template<typename MatOptConfigurator, typename LinElastEnergyType, bool IsometryConstraint>
template<typename MatOptConfigurator>
class optimalDeformSolverInterfaceAdaptive
: public optimalDeformPlotter<typename MatOptConfigurator::DataTypeContainer,MatOptConfigurator::ConfiguratorTypePf::_ShellFEType> {

 protected:
  typedef typename MatOptConfigurator::ConfiguratorType         ConfiguratorType;
  typedef typename MatOptConfigurator::ConfiguratorTypePf       ConfiguratorTypePf;
  typedef typename ConfiguratorType::RealType                   RealType;
  typedef typename ConfiguratorType::MaskType                   MaskType;
  typedef typename ConfiguratorType::InitType                   MeshType;
  typedef typename ConfiguratorTypePf::DTContainer              DataTypeContainer;
  typedef typename ConfiguratorType::VectorType                 VectorType;
  typedef pesopt::BoostParser  ParameterParserType;

  const ParameterParserType &_parser;
  const string _VTKFileType;
  const string _latexType;

public:
  optimalDeformSolverInterfaceAdaptive ( const ParameterParserType &Parser ) :
  optimalDeformPlotter<DataTypeContainer,ConfiguratorTypePf::_ShellFEType> ( Parser ),
  _parser( Parser ),
  _VTKFileType ( Parser.template get<string>("saving.VTKFileType") ),
  _latexType ( Parser.template get<string>("saving.LATEXType") ){ }

protected:

void constructForceRHS( const ParameterParserType &parser,
                        const MatOptConfigurator &matOptConf,
                        const DiscreteVectorFunctionStorage<ConfiguratorType,FirstOrder> &xAStorage,
                        const MaskType &DirichletMask,
                        VectorType &rhs_force,
                        const RealType factorForce = 1.,
                        const bool setDirichletNodes = true ) const{
    createForceShellFE<ConfiguratorType> ( parser, matOptConf._conf, xAStorage, rhs_force, DirichletMask, setDirichletNodes );
    RealType facForce = parser.template get<RealType> ( "Force.factor_force" );
    RealType scaleForceByThickness = 1.;
    if( parser.template get<int> ( "Force.scaleForceWithThickness") == 1  ){
        const RealType thickness = std::sqrt( matOptConf._materialInfo._factorBendingEnergy );
        const RealType exponent =  parser.template get<RealType> ( "Force.scaleForceWithThicknessExponent" );
        scaleForceByThickness = std::pow( thickness, exponent );
    }
    rhs_force *= factorForce * facForce * scaleForceByThickness;

}


template<DiscreteFunctionCacheType _DiscreteFunctionCacheType>
void markAndRefine( MeshType &mesh, VectorType & solutionMaterial,
                    VectorType &initDisp, VectorType &solDisp,
                    MaskType &DirichletMask,
                    VectorType &markedElements ) const{
        const MeshType oldMesh ( mesh );
        const MaskType oldDirichletMask( DirichletMask );
        switch( _parser.template get<int> ( "ConstraintProblem.adaptiveMarkingType" ) ){
            case 1:{ mesh.markAll();} break;
            //2xmarkAll
            case 2:{ mesh.markAll();
                     mesh.refineMarkedTriangles( );
                     //! Prolongate i) mesh (with tangent space) ii) material iii) displacement
                     switchProlongationOfMesh<ConfiguratorType > ( _parser.template get<int> ( "ConstraintProblemAdaptiveProlongationTypes.ProlongationMeshType" ), oldMesh, mesh );
                     switchProlongationOfBoundaryMask<ConfiguratorType>( _parser.template get<bool> ( "ConstraintProblemAdaptiveProlongationTypes.ProlongationMeshBoundaryType" ), mesh, DirichletMask );
                     switchProlongationOfPhaseField<ConfiguratorTypePf> ( _parser.template get<int> ( "ConstraintProblemAdaptiveProlongationTypes.ProlongationPfType" ), mesh, oldDirichletMask, DirichletMask, solutionMaterial );
                     initDisp = solDisp;
                     switchProlongationOfDisplacement<ConfiguratorType> ( _parser.template get<int> ( "ConstraintProblemAdaptiveProlongationTypes.ProlongationDispType" ), mesh, initDisp );
                     mesh.markAll();
            } break;
//             case 11:{ ConfiguratorTypePf confpf ( mesh );
//                       const RealType thresholdMaterialGrad = 0.5;
//                       markElements_L2gradient<ConfiguratorTypePf> ( confpf, mesh, solutionMaterial, thresholdMaterialGrad );
//             }break;
//             case 12:{ ConfiguratorTypePf confpf ( mesh );
//                       const RealType thresholdMaterialGrad = 0.5;
//                       markElements_L2gradient<ConfiguratorTypePf> ( confpf, mesh, solutionMaterial, thresholdMaterialGrad );
//                       markClampedElements<MeshType> ( mesh, oldDirichletMask );
//             }break;
//             case 21:{ ConfiguratorTypePf confpf ( mesh );
//                       const RealType thresholdMaterialGrad = 0.5;
//                       markElements_L2gradient<ConfiguratorTypePf> ( confpf, mesh, solutionMaterial, thresholdMaterialGrad );
//                       const RealType thresholdPercentilIsometryError = 0.75;
//                       ConfiguratorType conf ( mesh );
//                       ShellHandler<ConfiguratorType,FirstOrder> shellHandler ( _parser, conf );
//                       IsometryConstraintIntegrated<MatOptConfigurator> isoConstraintIntegratedOp ( conf, shellHandler.getChartToUndeformedShell_Cache() );
//                       VectorType isoConstraintErrorOnElements( mesh.getNumElements() );
//                       isoConstraintIntegratedOp.applyOnElements ( solDisp, isoConstraintErrorOnElements );
//                       markElements_ElementErrorVec<MeshType> ( mesh, isoConstraintErrorOnElements, thresholdPercentilIsometryError );
//             }break;
//             case 22:{ ConfiguratorTypePf confpf ( mesh );
//                       const RealType thresholdMaterialGrad = 0.5;
//                       markElements_L2gradient<ConfiguratorTypePf> ( confpf, mesh, solutionMaterial, thresholdMaterialGrad );
//                       const RealType thresholdPercentilIsometryError = 0.75;
//                       ConfiguratorType conf ( mesh );
//                       ShellHandler<ConfiguratorType,FirstOrder> shellHandler ( _parser, conf );
//                       IsometryConstraintIntegrated<MatOptConfigurator> isoConstraintIntegratedOp ( conf, shellHandler.getChartToUndeformedShell_Cache() );
//                       VectorType isoConstraintErrorOnElements( mesh.getNumElements() );
//                       isoConstraintIntegratedOp.applyOnElements ( solDisp, isoConstraintErrorOnElements );
//                       markElements_ElementErrorVec<MeshType> ( mesh, isoConstraintErrorOnElements, thresholdPercentilIsometryError );
//                       markClampedElements<MeshType> ( mesh, oldDirichletMask );
//             }break;
            //IsometryConstraint
            case 23:{ ConfiguratorTypePf confpf ( mesh );
                      const RealType thresholdPercentilIsometryError = 0.75;
                      ConfiguratorType conf ( mesh );
                      ShellHandler<ConfiguratorType,FirstOrder> shellHandler ( _parser, conf );
                      IsometryConstraintIntegrated<MatOptConfigurator> isoConstraintIntegratedOp ( conf, shellHandler.getChartToUndeformedShell_Cache() );
                      VectorType isoConstraintErrorOnElements( mesh.getNumElements() );
                      isoConstraintIntegratedOp.applyOnElements ( solDisp, isoConstraintErrorOnElements );
                      markElements_ElementErrorVec<MeshType> ( mesh, isoConstraintErrorOnElements, thresholdPercentilIsometryError );
                      markClampedElements<MeshType> ( mesh, oldDirichletMask );
            }break;
            case 31:{ ConfiguratorType conf ( mesh );
                      ShellHandler<ConfiguratorType,_DiscreteFunctionCacheType> shellHandler ( _parser, conf );
                      const RealType thresholdPercentil = _parser.template get<RealType> ( "ConstraintProblem.adaptiveMarkingPercentil" );
                      markElements_L2ApproxHessian<ConfiguratorType> ( conf, mesh, shellHandler.getChartToUndeformedShell_Cache(), solDisp, thresholdPercentil );
            }break;
            case 32:{ ConfiguratorType conf ( mesh );
                      ShellHandler<ConfiguratorType,_DiscreteFunctionCacheType> shellHandler ( _parser, conf );
                      const RealType thresholdPercentil = _parser.template get<RealType> ( "ConstraintProblem.adaptiveMarkingPercentil" );
                      markElements_L2ApproxHessian<ConfiguratorType> ( conf, mesh, shellHandler.getChartToUndeformedShell_Cache(), solDisp, thresholdPercentil );
                      markClampedElements<MeshType> ( mesh, oldDirichletMask );
            }break;


            case 41:{ ConfiguratorType conf ( mesh );
                      ShellHandler<ConfiguratorType,_DiscreteFunctionCacheType> shellHandler ( _parser, conf );
                      const RealType thresholdPercentil = _parser.template get<RealType> ( "ConstraintProblem.adaptiveMarkingPercentil" );
                      markElements_L2ApproxHessian<ConfiguratorType> ( conf, mesh, shellHandler.getChartToUndeformedShell_Cache(), solDisp, thresholdPercentil );
                      markClampedElements<MeshType> ( mesh, oldDirichletMask );
            }break;
            default : throw std::invalid_argument ( pesopt::strprintf ( "wrong method in file %s at line %d.", __FILE__, __LINE__ ).c_str() );
        };

        mesh.getMarkedTrianglesMask( markedElements );
        mesh.refineMarkedTriangles( );
        //! Prolongate i) mesh (with tangent space) ii) material iii) displacement
        switchProlongationOfMesh<ConfiguratorType> ( _parser.template get<int> ( "ConstraintProblemAdaptiveProlongationTypes.ProlongationMeshType" ), oldMesh, mesh );
        switchProlongationOfBoundaryMask<ConfiguratorType>( _parser.template get<bool> ( "ConstraintProblemAdaptiveProlongationTypes.ProlongationMeshBoundaryType" ), mesh, DirichletMask );
        switchProlongationOfPhaseField<ConfiguratorTypePf> ( _parser.template get<int> ( "ConstraintProblemAdaptiveProlongationTypes.ProlongationPfType" ), mesh, oldDirichletMask, DirichletMask, solutionMaterial );
        initDisp = solDisp;
        switchProlongationOfDisplacement<ConfiguratorType> ( _parser.template get<int> ( "ConstraintProblemAdaptiveProlongationTypes.ProlongationDispType" ), mesh, initDisp );
  }



template<DiscreteFunctionCacheType _DiscreteFunctionCacheType>
void markAndRefineForGivenElementErrorVector( MeshType &mesh, VectorType & solutionMaterial, VectorType &initDisp, VectorType &solDisp, MaskType &DirichletMask,
                                              const VectorType &elementErrorVec, VectorType &markedElements
    ) const{
            const MeshType oldMesh ( mesh );
            const MaskType oldDirichletMask( DirichletMask );

            const RealType thresholdPercentil = _parser.template get<RealType> ( "ConstraintProblem.adaptiveMarkingPercentil" );
            markElements_ElementErrorVec<MeshType>( mesh, elementErrorVec, thresholdPercentil );

            mesh.getMarkedTrianglesMask( markedElements );
            mesh.refineMarkedTriangles( );
            //! Prolongate i) mesh (with tangent space) ii) material iii) displacement
            switchProlongationOfMesh<ConfiguratorType> ( _parser.template get<int> ( "ConstraintProblemAdaptiveProlongationTypes.ProlongationMeshType" ), oldMesh, mesh );
            switchProlongationOfBoundaryMask<ConfiguratorType>( _parser.template get<bool> ( "ConstraintProblemAdaptiveProlongationTypes.ProlongationMeshBoundaryType" ), mesh, DirichletMask );
            switchProlongationOfPhaseField<ConfiguratorTypePf> ( _parser.template get<int> ( "ConstraintProblemAdaptiveProlongationTypes.ProlongationPfType" ), mesh, oldDirichletMask, DirichletMask, solutionMaterial );
            initDisp = solDisp;
            switchProlongationOfDisplacement<ConfiguratorType> ( _parser.template get<int> ( "ConstraintProblemAdaptiveProlongationTypes.ProlongationDispType" ), mesh, initDisp );
}


  void setInitialMaterial( const int type, const MeshType &mesh, VectorType &initMaterial, const VectorType &solutionMaterial ) const {
     switch( type ){

            case 0:{ initMaterial = solutionMaterial; }break;

            case 1:{ //! load from file
                pesopt::loadVectorFromFile<VectorType>(  initMaterial, _parser.template get<string> ("CompareDesigns.materialFile") );
            }break;

            case 2:{ //! construct with shellHandler
                ConfiguratorType conf ( mesh );
                ShellMaterialGenerator<ConfiguratorType,ConfiguratorTypePf::_ShellFEType> shellHandler( conf );
                const RealType areaHardMaterial = _parser.template get<RealType> ("InitMaterial.areaHardMaterial");
                const int designType = _parser.template get<int> ("InitMaterial.designType");
                string designTypeName;
                shellHandler.switchMaterialTypeForFixedArea( designType, areaHardMaterial, initMaterial, designTypeName );
            }break;

            case 3:{ //! threshold to {-1,1}
                for( int i=0; i<initMaterial.size(); ++i ){
                    if( solutionMaterial[i] < 0 ) initMaterial[i] = -1.;
                    else                          initMaterial[i] = 1.;
                }
            }break;

            //BEGIN TODO only for P1
            case 11:{  //symmetrize in y direction: init(x,y,z) = 1/2 ( sol(x,y,z) + sol(x,1-y,z) )
                for( int nodeIdx=0; nodeIdx<initMaterial.size(); ++nodeIdx ){
                    int vertexIndex;
                    typename MeshType::Point3DType reflectedVertex ( mesh.getVertex(nodeIdx) );
                    reflectedVertex[1] = 1. - mesh.getVertex(nodeIdx)[1];
                    mesh.findClosestVertex( reflectedVertex, vertexIndex );
                    initMaterial[nodeIdx] = 0.5 * ( solutionMaterial[nodeIdx] + solutionMaterial[vertexIndex] );
                }
            }break;
            case 12:{  //symmetrize in x and y direction: init(x,y,z) = 1/4 ( sol(x,y,z) + sol(x,1-y,z) + sol(1-x,y,z) + sol(1-x,1-y,z) )
                for( int nodeIdx=0; nodeIdx<initMaterial.size(); ++nodeIdx ){
                    int reflectedVertexIndex_x, reflectedVertexIndex_y, reflectedVertexIndex_xy;
                    typename MeshType::Point3DType reflectedVertex_y ( mesh.getVertex(nodeIdx) ), reflectedVertex_x ( mesh.getVertex(nodeIdx) ), reflectedVertex_xy ( mesh.getVertex(nodeIdx) );
                    reflectedVertex_x[0] = 1. - mesh.getVertex(nodeIdx)[0]; mesh.findClosestVertex( reflectedVertex_x, reflectedVertexIndex_x );
                    reflectedVertex_y[1] = 1. - mesh.getVertex(nodeIdx)[1]; mesh.findClosestVertex( reflectedVertex_y, reflectedVertexIndex_y );
                    reflectedVertex_xy[0] = 1. - mesh.getVertex(nodeIdx)[0]; reflectedVertex_xy[1] = 1. - mesh.getVertex(nodeIdx)[1]; mesh.findClosestVertex( reflectedVertex_xy, reflectedVertexIndex_xy );
                    initMaterial[nodeIdx] = 0.25 * ( solutionMaterial[nodeIdx] + solutionMaterial[reflectedVertexIndex_x] + solutionMaterial[reflectedVertexIndex_y] + solutionMaterial[reflectedVertexIndex_xy]);
                }
            }break;

            case 15:{  //rotate in xy plane around (0,5,0,5) by 90 degree: init(x,y,z) = 1/2 ( sol(x,y,z) + sol(Q(x,y),z) )
                for( int nodeIdx=0; nodeIdx<initMaterial.size(); ++nodeIdx ){
                    int rotatedVertexIndex_xy;
                    typename MeshType::Point3DType rotatedVertex_xy ( mesh.getVertex(nodeIdx) );
                    rotatedVertex_xy[0] = mesh.getVertex(nodeIdx)[1]; rotatedVertex_xy[1] = 1. - mesh.getVertex(nodeIdx)[0]; mesh.findClosestVertex( rotatedVertex_xy, rotatedVertexIndex_xy );
                    initMaterial[nodeIdx] = 0.5 * ( solutionMaterial[nodeIdx] + solutionMaterial[rotatedVertexIndex_xy]);
                }
            }break;
            //BEGIN TODO only for P1

            case 1000:{ // initNew = ( solutionOld + Random) / 2
                VectorType rand = VectorType::Random( solutionMaterial.size() );
                initMaterial = 0.5 * ( solutionMaterial + rand );
            }break;
     };
 }

  void saveStressOnUndeformedShell( const string fileName,
                                    const ShellPlotter<ConfiguratorType> &shellPlotter,
                                    const VectorType &material, const VectorType &xA, const VectorType &displacement,
                                    const std::vector<VectorType> &StressOnElementVec, const std::vector<string> &nameStressOnElementVec ) const{
      //! averaged stress on elements;
// TODO instead of updateMesh use new ShellPlotter
//         shellPlotter.clearData();
//         shellPlotter.updateMesh( "undeform", xA );
//         for( int i=0; i<StressOnElementVec.size(); ++i ){
//             shellPlotter.addScalarData ( StressOnElementVec[i], nameStressOnElementVec[i], FACE_DATA );
//         }
//         shellPlotter.saveShellToFile( fileName );
  }

   template<DiscreteFunctionCacheType _DiscreteFunctionCacheType,bool IsometryConstraint>
   void saveNorms ( const MatOptConfigurator &matOptConf, const ShellHandler<ConfiguratorType,_DiscreteFunctionCacheType> &shellHandler, const VectorType &material,  const VectorType &solDisp, DeformationOptimizationShellEnergyInfo<RealType> &energyInfo ) const {
        RealType dispAtCenter = 0., dispLInftyAtNodes = 0., dispL2NormSqr = 0., dispLinftyNorm = 0.;
        for( int nodeIdx=0; nodeIdx< matOptConf._conf.getInitializer().getNumVertices(); ++nodeIdx ){
            for( int comp=0; comp < 3; ++comp ){
                RealType abs = std::abs( solDisp[nodeIdx + comp * matOptConf._conf.getNumGlobalDofs()] );
                if( abs > dispLInftyAtNodes ) dispLInftyAtNodes = abs;
            }
        }
        PfWeightedLpNormToP<MatOptConfigurator>  ( matOptConf, shellHandler.getChartToUndeformedShell_Cache(), solDisp, material, 2. ).assembleAdd( dispL2NormSqr );
        LInftyNorm<MatOptConfigurator> ( matOptConf, shellHandler.getChartToUndeformedShell_Cache(), solDisp               ).assembleAdd( dispLinftyNorm );
        energyInfo.setL2Norm( std::sqrt(dispL2NormSqr) );
        energyInfo.setLInfNormAtQuadPoints( dispLinftyNorm );
        energyInfo.setLInfNormAtNodes( dispLInftyAtNodes );
   }

   //plot energy diagramm: x - percent of energy, y - area
   void saveEnergyVsArea ( const string saveDirectory, const VectorType &totalStressVec, const VectorType &areaOnElementsVec ) const{
        std::map<int, RealType> map_ElementNumber_averagedEnergy;
        for( int elIdx=0; elIdx < totalStressVec.size(); ++elIdx ) map_ElementNumber_averagedEnergy.insert ( std::pair<int,RealType>(elIdx,totalStressVec[elIdx] / areaOnElementsVec[elIdx]) );
        // Defining a lambda function to compare two pairs. It will compare two pairs using second field
        typedef std::function<bool(std::pair<int,RealType>, std::pair<int,RealType>)> Comparator;
        Comparator compFunctor = [](std::pair<int,RealType> pair1 ,std::pair<int,RealType> pair2){ return pair1.second > pair2.second;  };
        // Declaring a set that will store the pairs using above comparision logic
        //WARNING set stores only one entry with same value, so have to allow multiset
        //         std::set<std::pair<int,RealType>, Comparator> sortedSet( map_ElementNumber_averagedEnergy.begin(), map_ElementNumber_averagedEnergy.end(), compFunctor);
        std::multiset<std::pair<int,RealType>, Comparator> sortedSet( map_ElementNumber_averagedEnergy.begin(), map_ElementNumber_averagedEnergy.end(), compFunctor);
        // Iterate over a set using range base for loop, it will display the items in sorted order of values
        std::vector<RealType> averagedTotalStressVecOnElementsSorted;
        std::vector<int> elementIdxSorted;
        for (std::pair<int,RealType> element : sortedSet){
            elementIdxSorted.push_back( element.first );
            averagedTotalStressVecOnElementsSorted.push_back( element.second );
        }
        std::vector<RealType> partOfEnergyVec, partOfAreaVec;
        RealType accumulatedEnergy = 0., accumulatedArea = 0.;
        for( int i=0; i<totalStressVec.size(); ++i ){
            accumulatedEnergy += totalStressVec[ elementIdxSorted[i] ];
            accumulatedArea += areaOnElementsVec[ elementIdxSorted[i] ];
            partOfEnergyVec.push_back( accumulatedEnergy );
            partOfAreaVec.push_back( accumulatedArea );
        }
        TikzPlotterCurve<RealType,std::vector<RealType>> ( saveDirectory, _latexType.c_str() ).plotCurve ( partOfEnergyVec, partOfAreaVec, "EnergyVsArea", "energy vs area" );
   }

   void saveResultsDeformationAndStress ( const string saveDirectory,
                      const ShellWithMaterialPlotter<ConfiguratorType,ConfiguratorTypePf::_ShellFEType> &shellPlotter,
                      const VectorType &material, const VectorType &xA,  const VectorType &solDisp,
                      const std::vector<VectorType> &StressOnElementVec, const std::vector<string> &nameStressOnElementVec ) const{
       pesopt::printVectorToFile<VectorType> ( solDisp, pesopt::strprintf( "%s/SolDisp.txt", saveDirectory.c_str() ).c_str(), 10 );
       shellPlotter.saveMeshWithMaterialToFile ( material,  pesopt::strprintf( "chart" ).c_str() );
       shellPlotter.saveShellWithMaterialToFile ( "undeform", xA, material, pesopt::strprintf( "undeformedShell" ).c_str()  );
       shellPlotter.saveShellWithMaterialToFile ( "disp", solDisp, material, pesopt::strprintf( "deformedShell" ).c_str() );

// TODO
//        const string fileNameData = "Stress";
//        this->saveStressOnUndeformedShell(  fileNameData, shellPlotter, material, xA, solDisp, StressOnElementVec, nameStressOnElementVec );
    }

    template<bool IsometryConstraint>
    void saveSummaryInTexFile( const string saveDirectory, const DeformationOptimizationShellEnergyInfo<RealType> &energyInfo ) const {
      TikzPlotterShellBucklingFileInfo fileInfo;
      fileInfo.setDeformationFiles( "chart_material.png",
                                    "undeformedShell_material.png",
                                    "deformedShell_material.png"   );
      fileInfo.setStressFiles( "Stress_averagedMembraneStressOnElements.png",
                               "Stress_averagedBendingStressOnElements.png",
                               "Stress_averagedTotalStressOnElements.png" );
      fileInfo.setEnergyVsAreaFile( "EnergyVsArea.pdf" );
      if( IsometryConstraint ){
       /*fileInfo.setStressForIsometryFiles( "Stress_averagedBendingStressWithoutMaterialFacOnElements.png", "Stress_affinePartOfDeformationOnElements.png" );*/
       fileInfo.setStressForIsometryFiles( "Stress_averagedBendingStressWithoutMaterialFacOnElements.png" );
      }
      TikzPlotterShellBuckling<DataTypeContainer> tikzPlotter( saveDirectory );
      tikzPlotter.template plotAll<IsometryConstraint>( "AllResults", energyInfo, fileInfo );
    }


   template<DiscreteFunctionCacheType _DiscreteFunctionCacheType, bool IsometryConstraint>
   void saveAllResults ( const string saveDirectory,
                      const MatOptConfigurator &matOptConf,
                      const ShellHandler<ConfiguratorType,_DiscreteFunctionCacheType> &shellHandler,
                      const VectorType &material,  const VectorType &solDisp,
                      const std::vector<VectorType> &StressOnElementVec, const std::vector<string> &nameStressOnElementVec,
                      const VectorType &totalStressVec, const VectorType &areaOnElementsVec,
                      DeformationOptimizationShellEnergyInfo<RealType> &energyInfo  ) const{

        ShellWithMaterialPlotter<ConfiguratorType,ConfiguratorTypePf::_ShellFEType> shellPlotter( matOptConf._conf, shellHandler.getChartToUndeformedShell(), shellHandler.getDirichletMask(), saveDirectory,  _VTKFileType );

        this->template saveNorms<_DiscreteFunctionCacheType,IsometryConstraint> ( matOptConf, shellHandler, material, solDisp, energyInfo );

        energyInfo.template saveToFile<ParameterParserType>( "energyInfo", saveDirectory );

        this->saveResultsDeformationAndStress(  saveDirectory, shellPlotter, material, shellHandler.getChartToUndeformedShell(), solDisp, StressOnElementVec, nameStressOnElementVec );

        this->saveEnergyVsArea( saveDirectory, totalStressVec, areaOnElementsVec );

        this->template saveSummaryInTexFile<IsometryConstraint>( saveDirectory, energyInfo );
    }



    template<bool IsometryConstraint>
    void plotConvergenceOfApdaptiveRefinement( const std::vector<DeformationOptimizationShellEnergyInfo<RealType>> &energyInfoVec ) const {
        const string saveDirectoryConvergence = this->_parser.createSubDirectory( "ConvergenceAdaptiveRef" );
        const int numAdaptiveRefinementSteps = energyInfoVec.size() - 1;
        VectorType membraneEnergyVec(numAdaptiveRefinementSteps+1), bendingEnergyVec(numAdaptiveRefinementSteps+1),
                   potEnergyVec(numAdaptiveRefinementSteps+1), storedElasticEnergyVec(numAdaptiveRefinementSteps+1), freeEnergyVec(numAdaptiveRefinementSteps+1);
        VectorType L2NormDispVec(numAdaptiveRefinementSteps+1), LInfNormDispAtQuadVec(numAdaptiveRefinementSteps+1), LInfNormDispAtNodesVec(numAdaptiveRefinementSteps+1);

        for( int refinementStep=0; refinementStep <= numAdaptiveRefinementSteps; ++refinementStep ){
            membraneEnergyVec[refinementStep] = energyInfoVec[refinementStep]._membraneEnergy;
            bendingEnergyVec[refinementStep] = energyInfoVec[refinementStep]._bendingEnergy;
            potEnergyVec[refinementStep] = energyInfoVec[refinementStep]._potentialEnergy;
            storedElasticEnergyVec[refinementStep] = energyInfoVec[refinementStep]._storedElasticEnergy;
            freeEnergyVec[refinementStep] = energyInfoVec[refinementStep]._dissipationEnergy;
            //           residualConstraintVec[refinementStep] = energyInfoVec[refinementStep]._residualConstraint;
            L2NormDispVec[refinementStep] = energyInfoVec[refinementStep]._L2Norm;
            LInfNormDispAtQuadVec[refinementStep] = energyInfoVec[refinementStep]._LInfNormAtQuadPoints;
            LInfNormDispAtNodesVec[refinementStep] = energyInfoVec[refinementStep]._LInfNormAtNodes;;
        }


        TikzPlotterCurve<RealType, VectorType> tikzPlotterCurve ( saveDirectoryConvergence );
        tikzPlotterCurve.plotCurve( membraneEnergyVec, "convergenceMembraneEnergy", "membrane energy" );
        tikzPlotterCurve.plotCurve( bendingEnergyVec, "convergenceBendingEnergy", "bending energy" );
        tikzPlotterCurve.plotCurve( potEnergyVec, "convergencePotentialEnergy", "potential energy" );
        tikzPlotterCurve.plotCurve( storedElasticEnergyVec, "convergenceStoredElasticEnergy", "stored elastic energy" );
        tikzPlotterCurve.plotCurve( freeEnergyVec, "convergenceFreeEnergy", "free energy" );
        tikzPlotterCurve.plotCurve( L2NormDispVec, "convergenceL2Norm", "L2 Norm" );
        tikzPlotterCurve.plotCurve( LInfNormDispAtQuadVec, "convergenceLInfNormQuad", "LInf at Quadpoints" );
        tikzPlotterCurve.plotCurve( LInfNormDispAtNodesVec, "convergenceLInfNormNodes", "LInf at Nodes" );

        TikzPlotterShellBuckling<DataTypeContainer> tikzPlotterShellBuckling( saveDirectoryConvergence );
        std::vector<string> convergenceFiles = { "convergenceMembraneEnergy.pdf", "convergenceBendingEnergy.pdf", "convergenceFreeEnergy.pdf", "convergenceLInfNormNodes.pdf" };
        tikzPlotterShellBuckling.plotConvegenceOfAdaptiveRefinement( "AllResultsConvAdRef", convergenceFiles);
        runPdflatexWithBashFile( "AllResultsConvAdRef.tex", saveDirectoryConvergence );
    }


    void plotConvergenceIsometryOfApdaptiveRefinement( const std::vector<IsometryInfo<RealType>> &isometryInfoVec ) const {
        const string saveDirectoryConvergenceIsometry = this->_parser.createSubDirectory( "ConvergenceIsometryAdaptiveRef" );
        const int numAdaptiveRefinementSteps = isometryInfoVec.size() - 1;
        VectorType gridSizeVec(numAdaptiveRefinementSteps+1);
        VectorType isometryErrorL1Vec(numAdaptiveRefinementSteps+1);
        VectorType isometryErrorL2Vec(numAdaptiveRefinementSteps+1);
        VectorType GaussCurvatureL1Vec(numAdaptiveRefinementSteps+1);
        VectorType GaussCurvatureL1DiffVec(numAdaptiveRefinementSteps+1);
        VectorType ConvGaussCurvatureL1Vec(numAdaptiveRefinementSteps+1);
        VectorType errorApproxD2uToFineSolutionL2Vec(numAdaptiveRefinementSteps+1);

        for( int refinementStep=0; refinementStep <= numAdaptiveRefinementSteps; ++refinementStep ){
            gridSizeVec[refinementStep] = isometryInfoVec[refinementStep]._gridSize;
            isometryErrorL1Vec[refinementStep] = isometryInfoVec[refinementStep]._isometryErrorL1;
            isometryErrorL2Vec[refinementStep] = isometryInfoVec[refinementStep]._isometryErrorL2;
            GaussCurvatureL1Vec[refinementStep] = isometryInfoVec[refinementStep]._GaussCurvatureL1;
            GaussCurvatureL1DiffVec[refinementStep] = isometryInfoVec[refinementStep]._GaussCurvatureL1Diff;
            ConvGaussCurvatureL1Vec[refinementStep] = isometryInfoVec[refinementStep]._ConvGaussCurvatureL1;
            errorApproxD2uToFineSolutionL2Vec[refinementStep] = isometryInfoVec[refinementStep]._errorApproxD2uToFineSolutionL2;
        }

        VectorType EOC_isometryErrorL1Vec(numAdaptiveRefinementSteps+1); EOC_isometryErrorL1Vec[0] = 0.; //no previous iteration step to compare
        VectorType EOC_isometryErrorL2Vec(numAdaptiveRefinementSteps+1); EOC_isometryErrorL2Vec[0] = 0.; //no previous iteration step to compare
        VectorType EOC_GaussCurvatureL1Vec(numAdaptiveRefinementSteps+1); EOC_GaussCurvatureL1Vec[0] = 0.; //no previous iteration step to compare
        VectorType EOC_GaussCurvatureL1DiffVec(numAdaptiveRefinementSteps+1); EOC_GaussCurvatureL1DiffVec[0] = 0.; //no previous iteration step to compare
        VectorType EOC_ConvGaussCurvatureL1Vec(numAdaptiveRefinementSteps+1); EOC_ConvGaussCurvatureL1Vec[0] = 0.; //no previous iteration step to compare
        VectorType EOC_errorApproxD2uToFineSolutionL2Vec(numAdaptiveRefinementSteps+1); EOC_errorApproxD2uToFineSolutionL2Vec[0] = 0.; //no previous iteration step to compare

        for( int refinementStep=1; refinementStep<=numAdaptiveRefinementSteps; ++refinementStep ){
          RealType facGridSize = gridSizeVec[refinementStep]/gridSizeVec[refinementStep-1];
          EOC_isometryErrorL1Vec[refinementStep] = std::log( isometryErrorL1Vec[refinementStep] / isometryErrorL1Vec[refinementStep-1] ) / std::log(facGridSize);
          EOC_isometryErrorL2Vec[refinementStep] = std::log( isometryErrorL2Vec[refinementStep] / isometryErrorL2Vec[refinementStep-1] ) / std::log(facGridSize);
          EOC_GaussCurvatureL1Vec[refinementStep] = std::log( GaussCurvatureL1Vec[refinementStep] / GaussCurvatureL1Vec[refinementStep-1] ) / std::log(facGridSize);
          EOC_GaussCurvatureL1DiffVec[refinementStep] = std::log( GaussCurvatureL1DiffVec[refinementStep] / GaussCurvatureL1DiffVec[refinementStep-1] ) / std::log(facGridSize);
          EOC_ConvGaussCurvatureL1Vec[refinementStep] = std::log( ConvGaussCurvatureL1Vec[refinementStep] / ConvGaussCurvatureL1Vec[refinementStep-1] ) / std::log(facGridSize);
          EOC_errorApproxD2uToFineSolutionL2Vec[refinementStep] = std::log( errorApproxD2uToFineSolutionL2Vec[refinementStep] / errorApproxD2uToFineSolutionL2Vec[refinementStep-1] ) / std::log(facGridSize);
        }

        TexGeneratorTable<RealType, VectorType> tablePlotter ( saveDirectoryConvergenceIsometry );
        std::vector<VectorType> vecs; std::vector<string> describtions;
        vecs.push_back(gridSizeVec), describtions.push_back( "$h$" );
        vecs.push_back(isometryErrorL1Vec), describtions.push_back( "$|g_B - g_A|_{L^1}$" );
        vecs.push_back(EOC_isometryErrorL1Vec), describtions.push_back( "EOC $|g_B - g_A|_{L^1}$" );
        vecs.push_back(isometryErrorL2Vec), describtions.push_back( "$|g_B - g_A|_{L^2}$" );
        vecs.push_back(EOC_isometryErrorL2Vec), describtions.push_back( "EOC $|g_B - g_A|_{L^2}$" );
        vecs.push_back(GaussCurvatureL1Vec), describtions.push_back( "$|det( \\nabla \\nabla_h u_h \\cdot n_h) |_{L^1}$" );
        vecs.push_back(EOC_GaussCurvatureL1Vec), describtions.push_back( "EOC $|det( \\nabla \\nabla_h u_h \\cdot n_h) |_{L^1}$" );
        vecs.push_back(GaussCurvatureL1DiffVec), describtions.push_back( "$|det( \\nabla \\nabla_h u_h \\cdot n_h) - det( D^2 u_{\\text{fine}} \\cdot n_{\\text{fine}})|_{L^1}$" );
        vecs.push_back(EOC_GaussCurvatureL1DiffVec), describtions.push_back( "EOC $|det( \\nabla \\nabla_h u_h \\cdot n_h) - det( D^2 u_{\\text{fine}} \\cdot n_{\\text{fine}})|_{L^1}$" );
        vecs.push_back(ConvGaussCurvatureL1Vec), describtions.push_back( "$|det( \\hat \\phi_i * \\nabla \\nabla_h u_h \\cdot n_h) |_{L^1}$" );
        vecs.push_back(EOC_ConvGaussCurvatureL1Vec), describtions.push_back( "EOC $|det( \\phi_i *  \\nabla \\nabla_h u_h \\cdot n_h) |_{L^1}$" );
        vecs.push_back(errorApproxD2uToFineSolutionL2Vec), describtions.push_back( "$| \\nabla \\nabla_h u_h - D^2 u_{\\text{fine}}|_{L^2}$" );
        vecs.push_back(EOC_errorApproxD2uToFineSolutionL2Vec), describtions.push_back( "EOC $| \\nabla \\nabla_h u_h - D^2 u_{\\text{fine}}|_{L^2}$" );
        tablePlotter.plotTable ( vecs, "ConvergenceIsometry", "Convergence Isometry", describtions);
    }



};


#endif
