#include <pesopt_IO.h>
#include <pesopt_DKT.h>
#include <pesopt_VTK.h>

#include "ShellOptimalDeformIsometrySemiNonLinAdaptive.h"
#include "ShellDeformationEnergiesSemiNonLin.h"

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
typedef typename ConfiguratorTypePlateTriang::Point3DType Point3DType;
typedef typename ConfiguratorTypePlateTriang::TangentVecType TangentVecType;
typedef typename ConfiguratorTypePlateTriang::SparseMatrixType SparseMatrixType;

int main(int argc, char ** argv) {

    pesopt::consoleOutputStartProgramm( "Compare semilinear and nonlinear bending energy" );

    // cout << " Usage of programm ShelOptimalDeformNonLin : " << endl
    //      << " 1 - num adaptive refinement steps [default 0]" << endl
    //      << " 2 - [ParameterFile]" << endl
    //      << " 3 - [saveDirectory]" << endl;


    int numAdaptiveRefinementSteps = 0;
    if( argc > 1 ) numAdaptiveRefinementSteps = stoi(argv[1]);

    ParameterParserType parser;
    if (argc <= 2){
        parser = ParameterParserType ( "../../../../parameters/dkt/OptDeformIsometryNonLin.ini" );
        parser.addCounterToSaveDirectory( "../../../../parameters/counter.txt", "/OptDeformIsometryNonLin" );
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
        MatOptConfTypePlateTriang matOptConf ( parser, conf, confpf );

        VectorType initDisp( 3 * conf.getNumGlobalDofs() );
        initDisp.setZero();

        VectorType material( confpf.getNumGlobalDofs() );
        for ( int i =0; i<material.size(); ++i ) material[i] = 1.;
//         ShellMaterialGenerator<ConfiguratorTypePlateTriang,ConfiguratorTypePlateTriangPf::_ShellFEType> shellHandler ( conf );
//         string designTypeName;
//         shellHandler.switchMaterialTypeForFixedArea( parser.template get<int> ("InitMaterial.designType" ),
//                                                      parser.template get<RealType> ("InitMaterial.areaHardMaterial"),
//                                                      material, designTypeName );

        // TODO Christoph
        ShellHandler<ConfiguratorTypePlateTriang,FirstAndSecondOrder> shellHandlerXA ( parser, conf );
        VectorType xB ( shellHandlerXA.getChartToUndeformedShell().size() );

        const RealType pi = 4 * atan ( 1.0 );
        const RealType radius = 2. / pi;
        for ( int nodeIdx=0; nodeIdx < mesh.getNumVertices(); ++nodeIdx ) {
            const Point3DType& coords ( mesh.getVertex(nodeIdx) );
            Point3DType coordsOnCylinder;
            coordsOnCylinder(0) = coords(0);
            coordsOnCylinder(1) = radius * sin( coords(1) * pi/2 );
            coordsOnCylinder(2) = - radius * cos( coords(1) * pi/2 );
            // coordsOnCylinder(0) = coords(0) * cos(pi/2 );
            // coordsOnCylinder(1) = coords(1);
            // coordsOnCylinder(2) = coords(0) * sin( pi/2 );
            for( int comp=0; comp<3; ++comp )xB[ nodeIdx + conf.getNumGlobalDofs() * comp ] = coordsOnCylinder[comp];
            if( ConfiguratorTypePlateTriang::_ShellFEType == C1Dofs ){
                TangentVecType firstTangentVecAtNode;
                firstTangentVecAtNode(0) =  1.;
                firstTangentVecAtNode(1) = 0.;
                firstTangentVecAtNode(2) = 0.;
                // firstTangentVecAtNode(0) =  cos( pi/2 );
                // firstTangentVecAtNode(1) = 0.;
                // firstTangentVecAtNode(2) = sin( pi/2 );
                TangentVecType secondTangentVecAtNode;
                secondTangentVecAtNode(0) = 0.;
                secondTangentVecAtNode(1) =  pi/2 * radius * cos( coords(1) * pi/2);
                secondTangentVecAtNode(2) =  pi/2 * radius * sin( coords(1) * pi/2 );
                for( int comp=0; comp<3; ++comp ){
                  xB[ nodeIdx +     mesh.getNumVertices() + conf.getNumGlobalDofs() * comp ] = firstTangentVecAtNode  [comp];
                  xB[ nodeIdx + 2 * mesh.getNumVertices() + conf.getNumGlobalDofs() * comp ] = secondTangentVecAtNode [comp];
                }
            }
        }



        DiscreteVectorFunctionStorage<ConfiguratorTypePlateTriang,FirstAndSecondOrder> xBStorage ( matOptConf._conf, xB, 3 );

        SemiNonlinearBendingEnergyOp<MatOptConfTypePlateTriang> testSemiNonlinearEnergyOp ( matOptConf, shellHandlerXA.getChartToUndeformedShell_Cache(),
                             material,
                             matOptConf._materialInfo._factorBendingEnergy );
//       SemiNonlinearBendingEnergy<MatOptConfTypePlateTriang> testSemiNonlinearEnergyOp( matOptConf,
//                              shellHandlerXA.getChartToUndeformedShell_Cache(),
//                              xBStorage,
//                              material,
//                              matOptConf._materialInfo._factorBendingEnergy  );
        RealType testSemiNonlinearEnergy = 0.;
//         testSemiNonlinearEnergyOp.assembleAdd ( testSemiNonlinearEnergy );
        VectorType disp = xB - shellHandlerXA.getChartToUndeformedShell();
        testSemiNonlinearEnergyOp.evaluateEnergy( disp,  testSemiNonlinearEnergy );
        // full nonlin
        NonlinearBendingEnergy<MatOptConfTypePlateTriang> testNonlinearEnergyOp( matOptConf,
                             shellHandlerXA.getChartToUndeformedShell_Cache(),
                             xBStorage,
                             material,
                             matOptConf._materialInfo._factorBendingEnergy  );
        RealType testNonlinearEnergy = 0.;
        testNonlinearEnergyOp.assembleAdd ( testNonlinearEnergy );
        cout <<"============================ " <<  endl
            << "seminonlin energy  = " <<  testSemiNonlinearEnergy  <<  endl
            << "nonlin energy      = " <<  testNonlinearEnergy  <<  endl
            << "============================" <<  endl;


        // TEST Gradient
        cout << "test gradient" <<  endl;
        VectorType grad ( disp.size() );
        testSemiNonlinearEnergyOp.evaluateGradient( disp,  grad );
        cout <<  "grad.norm = " <<  grad.norm() <<  endl;

        // TEST Hessian
        // cout << "test hessian" << endl;
        // SparseMatrixType hess ( disp.size());
        // testSemiNonlinearEnergyOp.evaluateHessian(disp, hess);
        // cout << "hess.norm = " << hess.norm() << endl;
  }

    //! print elapsed time
    std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << "duration = " << diff.count() << " sec" << endl;

  return ( EXIT_SUCCESS );
}
