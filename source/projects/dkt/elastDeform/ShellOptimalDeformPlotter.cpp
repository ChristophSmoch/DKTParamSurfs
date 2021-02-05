#include "ShellOptimalDeformSolverInterfaceAdaptive.h"

typedef DataTypeContainerDKTFE      DataTypeContainer;
typedef typename DataTypeContainer::RealType   RealType;
typedef typename DataTypeContainer::VectorType VectorType;
typedef pesopt::BoostParser ParameterParserType;




int main(int argc, char ** argv) {

    pesopt::consoleOutputStartProgramm( "Plot Deformation" );
         
    ParameterParserType parser;   
    if (argc < 3 ){
    cout << " Usage of programm: " << endl
         << " 1 - [ParameterFile]" << endl
         << " 2 - [saveDirectory]" << endl
         << " 3 - [Material Nodal (=1) or Element (=2)]" << endl;
         return 23;
    }
    
    parser = ParameterParserType ( argv[1] );
    parser.set ( "saving.saveDirectory", argv[2] );
    const string saveDirectory = parser.template get<string> ("saving.saveDirectory");
    
    const std::vector<string> nameStressOnElementVec = {"membraneStressOnElements", 
                                                       "bendingStressOnElements", 
                                                        "totalStressOnElements",
                                                        "averagedMembraneStressOnElements", "averagedBendingStressOnElements", 
                                                        "averagedTotalStressOnElements",
                                                        "areaOnElements"};
    
    const int VertexOrFaceData = stoi(argv[3]);
    if( VertexOrFaceData == 1 ){
        optimalDeformPlotter<DataTypeContainer,NodalValuedDofs> ( parser ).plotResults( saveDirectory, nameStressOnElementVec );
    }
    if( VertexOrFaceData == 2 ){
        optimalDeformPlotter<DataTypeContainer,ElementValuedDofs> ( parser ).plotResults( saveDirectory, nameStressOnElementVec );
    }
    

  return ( EXIT_SUCCESS );
}
