#include <pesopt_IO.h>
#include <pesopt_DKT.h>
#include <pesopt_VTK.h>

#include "ShellOptimalDeformIsometrySemiNonLinAdaptive.h"


typedef DataTypeContainerDKTFE     DataTypeContainer;
typedef typename DataTypeContainer::RealType                                                                RealType;
typedef typename DataTypeContainer::VectorType                                                              VectorType;
typedef DuodecQuadrature<RealType, typename DataTypeContainer::DomVecType>                                  QuadType;
typedef DKTPlateElement<DataTypeContainer>                                                                  TriangleTypeRefTriang;
typedef AdaptiveTriangMeshWithTangentSpace<DataTypeContainer, TriangleTypeRefTriang >                       MeshTypeRefTriang;
typedef RefTriangMeshConfiguratorDKT<DataTypeContainer, MeshTypeRefTriang, QuadType >                       ConfiguratorTypePlateTriang;
typedef RefTriangMeshConfiguratorP1<DataTypeContainer, MeshTypeRefTriang, QuadType >                        ConfiguratorTypePlateTriangPf;
typedef MaterialOptimizationShellConfigurator<ConfiguratorTypePlateTriang,ConfiguratorTypePlateTriangPf>    MatOptConfTypePlateTriang;
typedef typename ConfiguratorTypePlateTriang::MaskType                                                      MaskType;
typedef pesopt::BoostParser ParameterParserType;

int main(int argc, char ** argv) {

    pesopt::consoleOutputStartProgramm( "Find Optimal Deformation (nonlinear elasticity)" );

    cout << " Usage of programm ShelOptimalDeformNonLin : " << endl
         << " 1 - num adaptive refinement steps [default 0]" << endl
         << " 2 - [ParameterFile]" << endl
         << " 3 - [saveDirectory]" << endl;


    int numAdaptiveRefinementSteps = 0;
    if( argc > 1 ) numAdaptiveRefinementSteps = stoi(argv[1]);

    ParameterParserType parser;
    if (argc <= 2){
        parser = ParameterParserType ( "../../../../parameters/dkt/OptDeformIsometryNonLin.ini" );
        parser.addCounterToSaveDirectory( "../../../../parameters/counter.txt", "/OptDeformIsometrySemiNonLin" );
        parser.saveToFile( "ParameterParser.ini" );
    }
    if( argc > 2 ) parser = ParameterParserType ( argv[2] );
    parser.set ( "ConstraintProblem.outputLevel", 3 );
    parser.set ( "ConstraintProblem.outputLevelIpopt", 5 );

    if( argc > 3 ) parser.set ( "saving.saveDirectory", argv[3] );





    //! start watch
    auto startTime = std::chrono::high_resolution_clock::now();

    const string triangleType = parser.get<string> ( "InputMesh.TriangleType" );

    if( ( triangleType == "PlateElement" )  || (triangleType == "PlateElementMappedToShell" )  ){

        const string fileName = parser.get<string> ( "InputMesh.file" );
        const int tangentSpaceType = parser.get<int> ( "InputMesh.tangentSpaceType" );
        MeshTypeRefTriang mesh ( fileName, tangentSpaceType );
        ConfiguratorTypePlateTriang conf ( mesh ); ConfiguratorTypePlateTriangPf confpf ( mesh );

        VectorType initDisp( 3 * conf.getNumGlobalDofs() );
        initDisp.setZero();

        VectorType material( confpf.getNumGlobalDofs() );
        ShellMaterialGenerator<ConfiguratorTypePlateTriang,ConfiguratorTypePlateTriangPf::_ShellFEType> shellHandler ( conf );
        string designTypeName;
        shellHandler.switchMaterialTypeForFixedArea( parser.template get<int> ("InitMaterial.designType" ),
                                                     parser.template get<RealType> ("InitMaterial.areaHardMaterial"),
                                                     material, designTypeName );

        optimalDeformSolverIsometrySemiNonLinAdaptive<MatOptConfTypePlateTriang, SemiNonlinearMembraneBendingEnergyOp<MatOptConfTypePlateTriang> > OptDeformFinderShell ( parser );
        OptDeformFinderShell.computeOptimalDeformation_Adaptive( mesh, material, initDisp, numAdaptiveRefinementSteps );


  }

    //! print elapsed time
    std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << "duration = " << diff.count() << " sec" << endl;

  return ( EXIT_SUCCESS );
}
