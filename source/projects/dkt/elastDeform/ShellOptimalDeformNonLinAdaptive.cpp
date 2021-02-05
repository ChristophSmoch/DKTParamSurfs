#include <pesopt_DKT.h>
#include <pesopt_IO.h>
#include <pesopt_VTK.h>

#include "ShellOptimalDeformNonLinAdaptive.h"

typedef DataTypeContainerDKTFE     DataTypeContainer;
typedef typename DataTypeContainer::RealType RealType;
typedef typename DataTypeContainer::VectorType  VectorType;
typedef pesopt::BoostParser ParameterParserType;

// typedef triangShellFE::CenterQuadrature<RealType, typename DataTypeContainer::DomVecType>  QuadType;
typedef DuodecQuadrature<RealType, typename DataTypeContainer::DomVecType> QuadType;


typedef DKTPlateElement<DataTypeContainer>                                                        TriangleTypeRefTriang;
typedef AdaptiveTriangMeshWithTangentSpace<DataTypeContainer, TriangleTypeRefTriang >             MeshTypeRefTriang;
typedef RefTriangMeshConfiguratorDKT<DataTypeContainer, MeshTypeRefTriang, QuadType >             ConfiguratorTypePlateTriang;
typedef RefTriangMeshConfiguratorP1<DataTypeContainer, MeshTypeRefTriang, QuadType >              ConfiguratorTypePlateTriangPf;
typedef MaterialOptimizationShellConfigurator<ConfiguratorTypePlateTriang,ConfiguratorTypePlateTriangPf> MatOptConfTypePlateTriang;

typedef typename ConfiguratorTypePlateTriang::MaskType                                            MaskType;




int main(int argc, char ** argv) {

    pesopt::consoleOutputStartProgramm( "Find Optimal Deformation (nonlinear elasticity)" );
      
    cout << " Usage of programm ShelOptimalDeformNonLin : " << endl
         << " 1 - [ParameterFile]" << endl
         << " 2 - [saveDirectory]" << endl
         << " 3 - [nonlinelastenergyType: 11 - Discrete Differential Geometry, 21 - Nonlinear Membrane + Linear Bending, 31 - Nonlinear Membrane and Bending]" << endl
         << " 4 - num adaptive refinement steps [default 0]" << endl;
         
    ParameterParserType parser;   
    if (argc == 1){
        parser = ParameterParserType ( "../../../../parameters/dkt/OptDeformNonLin.ini" );
        parser.addCounterToSaveDirectory( "../../../../parameters/counter.txt", "/OptDeform" );
        parser.saveToFile( "ParameterParser.ini" );
    }
    if( argc > 1 ) parser = ParameterParserType ( argv[1] );
    parser.set ( "ConstraintProblem.outputLevel", 3 );
    parser.set ( "ConstraintProblem.outputLevelIpopt", 5 );
    
    if( argc > 2 ) parser.set ( "saving.saveDirectory", argv[2] );
    
    
    int nonlinelastenergytype; 
    if( argc > 3 ) nonlinelastenergytype = stoi(argv[3]); 
    else nonlinelastenergytype = parser.get<int> ( "ConstraintProblem.ElastEnergy" );
    
    int numAdaptiveRefinementSteps = 0;
    if( argc > 4 ) numAdaptiveRefinementSteps = stoi(argv[4]);
         
    
    
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
        
        //TODO
        VectorType material( confpf.getNumGlobalDofs() );
        for( int i=0; i<material.size(); ++i ) material[i] = 1.;

        
        
        switch( nonlinelastenergytype ){
            
            //NonLinear Membrane + NonLinear Bending
            case 31: {
                optimalDeformSolverNonLinAdaptive<MatOptConfTypePlateTriang, NonlinearMembraneBendingEnergyOp<MatOptConfTypePlateTriang> > OptDeformFinderShell ( parser );
                OptDeformFinderShell.computeOptimalDeformation_Adaptive( mesh, material, initDisp, numAdaptiveRefinementSteps );
            }break;
            
            default: 
                    throw std::invalid_argument ( pesopt::strprintf ( "wrong ElastEnergy %d in File %s at line %d.", nonlinelastenergytype, __FILE__, __LINE__ ).c_str() );
                    break;
        }
        
        
  }//end if triangleType == 2
    
    //! print elapsed time
    std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << "duration = " << diff.count() << " sec" << endl;

  return ( EXIT_SUCCESS );
}
