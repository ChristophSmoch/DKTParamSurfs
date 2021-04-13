#include <pesopt_IO.h>
#include <pesopt_DKT.h>
#include "ShellBucklingTikzPlotter.h"
#include "../dktFEMaterialOptimizationDefines.h"

#include "../elastDeform/ShellDeformationEnergiesWithMaterialInterfaces.h"

typedef double RealType;
typedef pesopt::BoostParser ParameterParserType;
typedef DataTypeContainerDKTFE  DataTypeContainer;
typedef typename DataTypeContainer::VectorType VectorType;
typedef typename DataTypeContainer::TangentVecType TangentVecType;


template<bool linearEnergyType, bool IsometryConstraint>
void executeProgrammOptDeformForSingleForce( const ParameterParserType &parserForceIter, const int ElastEnergy, const int numAdaptiveRefinementStepsSingleOpt ) {
    
    const string saveDirectoryForceIter = parserForceIter.template get<string>( "saving.saveDirectory" );
    
    string systemCommand;
    #ifdef PESOPT_WITH_OPENMP
    const string numThreadsSingleString = pesopt::strprintf( "%d", parserForceIter.get<int>( "BASH.numThreads" ) );
    systemCommand += "OMP_NUM_THREADS=" + numThreadsSingleString;
    #endif
                
    systemCommand += " nohup "; 
    if( parserForceIter.template get<bool> ("BASH.renewKerberosTicket") ) systemCommand += "krenew -- ";
    //     if( useNice )    systemCommand += "nice ";
    systemCommand += "xvfb-run -a --server-args \"-screen 0 1920x1080x24\" ";
                
    if( linearEnergyType ){
        if( IsometryConstraint ) systemCommand += "./ShellOptimalDeformationIsometryLinAdaptive ";
        //TODO else systemCommand += "./ShellOptimalDeformationLinAdaptive "                
    }else{
        if( IsometryConstraint ) systemCommand += "./ShellOptimalDeformationIsometryNonLinAdaptive ";
        else                     systemCommand += "./ShellOptimalDeformationNonLinAdaptive ";
    }
    
    systemCommand += pesopt::strprintf("%s/ParameterParser.ini", saveDirectoryForceIter.c_str() ) + " " + saveDirectoryForceIter;

    systemCommand += pesopt::strprintf( " %d", ElastEnergy );
    systemCommand += pesopt::strprintf( " %d", numAdaptiveRefinementStepsSingleOpt );
                
    systemCommand += " >" + saveDirectoryForceIter + "/ConsoleOutput.out 2>&1 ";
            
    cout << "systemCommand = " << systemCommand << endl;
    bool failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
    if ( failed ) cerr << endl 
                       << "programm" << endl 
                       << systemCommand << endl 
                       << "returned an error." << endl;
}


void executeProgrammPlotDeformForSingleForce(
    const ParameterParserType &parserForceIter, 
    const string saveDirectoryForceIter, 
    const string saveDirectoryForceIterFinal ) {
        
    string systemCommand;
    systemCommand += " nohup "; 
    if( parserForceIter.template get<bool> ("BASH.renewKerberosTicket") ) systemCommand += "krenew -- ";
    systemCommand += "xvfb-run -a --server-args \"-screen 0 1920x1080x24\" ";
    systemCommand += "./ShellOptimalDeformationPlotter ";
    systemCommand += pesopt::strprintf("%s/ParameterParser.ini", saveDirectoryForceIter.c_str() ) + " " + saveDirectoryForceIterFinal;
    systemCommand += " >" + saveDirectoryForceIter + "/ConsoleOutput.out 2>&1 ";     
    cout << "systemCommand = " << systemCommand << endl;
    bool failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
    if ( failed ) cerr << endl << "programm" << endl << systemCommand << endl  << "returned an error." << endl;
    runPdflatexWithBashFile( "AllResults.tex", saveDirectoryForceIterFinal );
}



void getSaveDirectoryForceIterVec( const ParameterParserType &parser,
                                   const std::vector<RealType> &compareFactorForceVec,
                                   std::vector<string> &saveDirectoryForceIterVec, std::vector<string> &subDirectoryForceIterVec ){
    for( int factorForceIter=0; factorForceIter<compareFactorForceVec.size(); ++factorForceIter ){
        ParameterParserType parserForceIter ( parser );
        string factorForceString = to_string( compareFactorForceVec[factorForceIter] );
        factorForceString.erase ( factorForceString.find_last_not_of('0') + 1, std::string::npos );
        string subDirectoryForceIter = parserForceIter.createSubDirectory( "FactorForce" + factorForceString );
        parserForceIter.set ( "saving.saveDirectory", subDirectoryForceIter );
        saveDirectoryForceIterVec[factorForceIter] = parserForceIter.get<string>( "saving.saveDirectory" );
        subDirectoryForceIterVec[factorForceIter] = "FactorForce" + factorForceString;
    }
}

void setSaveDirectoryThickness( ParameterParserType &parser, const RealType thickness ){
    string thicknessString = to_string( thickness );
    thicknessString.erase ( thicknessString.find_last_not_of('0') + 1, std::string::npos );
    const string saveDirectoryIter = parser.createSubDirectory( "Thickness" + thicknessString );
    parser.set ( "saving.saveDirectory", saveDirectoryIter );
}

void setSaveDirectoryForce( ParameterParserType &parser, const  RealType force ){
    string factorForceString = to_string( force );
    factorForceString.erase ( factorForceString.find_last_not_of('0') + 1, std::string::npos );
    const string saveDirectoryForceIter = parser.createSubDirectory( "FactorForce" + factorForceString );
    parser.set ( "saving.saveDirectory", saveDirectoryForceIter );
}

template<const bool IsometryConstraint>
void plotResultsToCompareThicknessAndForce( const string saveDirectory ) {  

    ParameterParserType parser ( saveDirectory + "/ParameterParser.ini" );
    std::vector<RealType> compareThicknessVec; parser.template getFreeSizeVector<RealType,std::vector<RealType>> ( "Compare.thicknessVec", compareThicknessVec );
    
    for( int thicknessIter=0; thicknessIter<compareThicknessVec.size(); ++thicknessIter ){
        ParameterParserType parserIterTmp ( parser );
        setSaveDirectoryThickness( parserIterTmp, compareThicknessVec[thicknessIter] );
        const string saveDirectoryIter = parserIterTmp.template get<string> ( "saving.saveDirectory" );
            
        ParameterParserType parserIter ( saveDirectoryIter + "/ParameterParser.ini" );
        std::vector<RealType> compareFactorForceVec; parserIter.template getFreeSizeVector<RealType,std::vector<RealType>> ( "Compare.factorForceVec", compareFactorForceVec );
            
        for( int factorForceIter=0; factorForceIter< compareFactorForceVec.size(); ++factorForceIter ){
            ParameterParserType parserForceIter ( parserIter );
            setSaveDirectoryForce( parserForceIter, compareFactorForceVec[factorForceIter] );
            const string saveDirectoryForceIter = parserForceIter.template get<string> ("saving.saveDirectory" );
            const string saveDirectoryForceIterFinal = saveDirectoryForceIter + "/FinalResult"; //TODO or numAdaptiveRefinementSteps?
            executeProgrammPlotDeformForSingleForce( parserForceIter, saveDirectoryForceIter, saveDirectoryForceIterFinal );
        }
        runPdflatexWithBashFile( "summary.tex", saveDirectoryIter );
    }
}


// template<const bool IsometryConstraint>
void summarizeResultsToCompareForces( const string saveDirectory, 
                                      const std::vector<string> &saveDirectoryForceIterVec, 
                                      const std::vector<string> &subDirectoryForceIterVec, 
                                      const std::vector<RealType> &compareFactorForceVec,
                                      const bool plotResults
                                    ) {  
        cout << endl << "start to summarize results to compare forces" << endl;
       
        ParameterParserType parser; parser.set ("saving.saveDirectory", saveDirectory );
        const string saveDirectoryCompareForce = parser.createSubDirectory ( "CompareScaledForce" );
         
        TikzPlotterCurve<RealType,std::vector<RealType>> tikzPlotterScaledForce( saveDirectoryCompareForce );
        std::vector<RealType> TikzEpotVec( compareFactorForceVec.size() ), TikzEstoredVec( compareFactorForceVec.size() ), TikzEfreeVec( compareFactorForceVec.size() ),
                              TikzEMemVec( compareFactorForceVec.size() ), TikzEBenVec( compareFactorForceVec.size() ), 
                              TikzL2Vec( compareFactorForceVec.size() ), TikzLInfNodesVec( compareFactorForceVec.size() ), TikzLInfQuadVec( compareFactorForceVec.size() );
        cout << endl << "get energy infos" << endl; 
        for( int factorForceIter=0; factorForceIter<compareFactorForceVec.size(); ++factorForceIter ){
                cout << "start for FactorForce = " << compareFactorForceVec[factorForceIter] << endl;
                DeformationOptimizationShellEnergyInfo<RealType> energyInfo;
                energyInfo.template loadFromFile<ParameterParserType>( "energyInfo", saveDirectoryForceIterVec[factorForceIter] + "/FinalResult" );
                TikzEpotVec[factorForceIter] = energyInfo._potentialEnergy;
                TikzEfreeVec[factorForceIter] = energyInfo._dissipationEnergy;
                TikzEstoredVec[factorForceIter] = energyInfo._storedElasticEnergy;
                TikzEMemVec[factorForceIter] = energyInfo._membraneEnergy;
                TikzEBenVec[factorForceIter] = energyInfo._bendingEnergy;
                TikzL2Vec[factorForceIter] = energyInfo._L2Norm;
                TikzLInfNodesVec[factorForceIter] = energyInfo._LInfNormAtNodes;
                TikzLInfQuadVec[factorForceIter] = energyInfo._LInfNormAtQuadPoints;
                cout << "finished for FactorForce = " << compareFactorForceVec[factorForceIter] << endl;
        }
        cout << endl << "plot curves" << endl;
        tikzPlotterScaledForce.plotCurve ( compareFactorForceVec, TikzEpotVec, "Epot", "potential energy" );
        tikzPlotterScaledForce.plotCurve ( compareFactorForceVec, TikzEfreeVec, "Efree", "free energy" );
        tikzPlotterScaledForce.plotCurve ( compareFactorForceVec, TikzEstoredVec, "EstoredElastic", "stored elastic energy" );
        tikzPlotterScaledForce.plotCurve ( compareFactorForceVec, TikzEMemVec, "Emem", "membrane energy" );
        tikzPlotterScaledForce.plotCurve ( compareFactorForceVec, TikzEBenVec, "Eben", "bending energy" );
        tikzPlotterScaledForce.plotCurve ( compareFactorForceVec, TikzL2Vec, "L2", "L2 norm of disp" );
        tikzPlotterScaledForce.plotCurve ( compareFactorForceVec, TikzLInfNodesVec, "LInfNodes", "Linf norm of disp (at nodes)" );
        tikzPlotterScaledForce.plotCurve ( compareFactorForceVec, TikzLInfQuadVec, "LinfQuads", "Linf norm of disp (at quadpoints)" );

        const string mark = "x";
        tikzPlotterScaledForce.plotLogLogCurve ( compareFactorForceVec, TikzEpotVec, "EpotLogLog", "potential energy LogLog", mark );
        tikzPlotterScaledForce.plotLogLogCurve ( compareFactorForceVec, TikzEfreeVec, "EfreeLogLog", "free energy LogLog", mark );
        tikzPlotterScaledForce.plotLogLogCurve ( compareFactorForceVec, TikzEstoredVec, "EstoredElasticLogLog", "stored elastic energy LogLog", mark );
        tikzPlotterScaledForce.plotLogLogCurve ( compareFactorForceVec, TikzEMemVec, "EmemLogLog", "membrane energy LogLog", mark );
        tikzPlotterScaledForce.plotLogLogCurve ( compareFactorForceVec, TikzEBenVec, "EbenLogLog", "bending energy LogLog", mark );
        tikzPlotterScaledForce.plotLogLogCurve ( compareFactorForceVec, TikzL2Vec, "L2LogLog", "L2 norm of disp LogLog" );
        tikzPlotterScaledForce.plotLogLogCurve ( compareFactorForceVec, TikzLInfNodesVec, "LInfNodesLogLog", "Linf norm of disp (at nodes) LogLog", mark );
        tikzPlotterScaledForce.plotLogLogCurve ( compareFactorForceVec, TikzLInfQuadVec, "LinfQuadsLogLog", "Linf norm of disp (at quadpoints) LogLog", mark );
            
        //document with all informations for fixed thickness
        cout << endl << "plot summary" << endl;
        TikzPlotterShellBucklingComparison<DataTypeContainer> tikzPlotterCompareScaledForces ( saveDirectory );
        std::vector<string> energyPlotFilesVec = {"CompareScaledForce/LInfNodesLogLog.pdf", "CompareScaledForce/EmemLogLog.pdf", "CompareScaledForce/EbenLogLog.pdf" };
        tikzPlotterCompareScaledForces.plotSummaryCompareForces( "summary", energyPlotFilesVec, subDirectoryForceIterVec, compareFactorForceVec  );
        if( plotResults ) runPdflatexWithBashFile( "summary.tex", saveDirectory );
}



template<const bool IsometryConstraint>
void summarizeResultsToCompareThickness( const string saveDirectory ) {   
    
    ParameterParserType parser ( saveDirectory + "/ParameterParser.ini" );
    parser.set ( "saving.saveDirectory", saveDirectory );
    std::vector<RealType> compareThicknessVec; parser.template getFreeSizeVector<RealType,std::vector<RealType>> ( "Compare.thicknessVec", compareThicknessVec );
    
    const string valueToCompare = parser.template get<string> ("Compare.valueToCompare" );
    RealType boundToCompare;
    if( valueToCompare == "LInf" )    boundToCompare = parser.template get<RealType> ("Compare.LInfBound");
    if( valueToCompare == "Estored" ) boundToCompare = parser.template get<RealType> ("Compare.EstoredBound");
    
    std::vector<RealType> LargestValueOfForceWithBoundVec ( compareThicknessVec.size() ), SmallestValueOfForceAboveBoundVec ( compareThicknessVec.size() );
    std::vector<RealType> EstoredVec ( compareThicknessVec.size() );
    for( int thicknessIter=0; thicknessIter<compareThicknessVec.size(); ++thicknessIter ){
        ParameterParserType parserIterTmp ( parser );
        setSaveDirectoryThickness( parserIterTmp, compareThicknessVec[thicknessIter] );
        const string saveDirectoryIter = parserIterTmp.template get<string> ( "saving.saveDirectory" );
            
        ParameterParserType parserIter ( saveDirectoryIter + "/ParameterParser.ini" );
        parserIter.set ( "saving.saveDirectory", saveDirectoryIter );
        std::vector<RealType> compareFactorForceVec; parserIter.template getFreeSizeVector<RealType,std::vector<RealType>> ( "Compare.factorForceVec", compareFactorForceVec );
        
        std::vector<RealType>TikzLInfNodesVec( compareFactorForceVec.size() ), TikzEstoredVec( compareFactorForceVec );
        for( int factorForceIter=0; factorForceIter<compareFactorForceVec.size(); ++factorForceIter ){
            ParameterParserType parserForceIter ( parserIter );
            string factorForceString = to_string( compareFactorForceVec[factorForceIter] );
            factorForceString.erase ( factorForceString.find_last_not_of('0') + 1, std::string::npos );
            const string subDirectoryForceIter = parserForceIter.createSubDirectory( "FactorForce" + factorForceString );
            DeformationOptimizationShellEnergyInfo<RealType> energyInfo;
            energyInfo.template loadFromFile<ParameterParserType>( "energyInfo", subDirectoryForceIter + "/FinalResult" );
            TikzLInfNodesVec[factorForceIter] = energyInfo._LInfNormAtNodes;
            TikzEstoredVec[factorForceIter] = energyInfo._storedElasticEnergy;
        }
            
        //get largest value of force s.t. Linf-norm of disp is <= boundToCompare
        int IndexLargestValueOfForceWithBound;
        for( int factorForceIter=0; factorForceIter<compareFactorForceVec.size(); ++factorForceIter ){
            if( TikzLInfNodesVec[factorForceIter] < boundToCompare ) IndexLargestValueOfForceWithBound = factorForceIter;
        }
            
        LargestValueOfForceWithBoundVec[thicknessIter] = compareFactorForceVec[IndexLargestValueOfForceWithBound];
        SmallestValueOfForceAboveBoundVec[thicknessIter] = compareFactorForceVec[IndexLargestValueOfForceWithBound+1];
        EstoredVec[thicknessIter] = TikzEstoredVec[IndexLargestValueOfForceWithBound];
            
    }// end for thicknessIter
    const string subDirectoryCompareThickness = parser.createSubDirectory( "CompareThickness" );
    TikzPlotterCurve<RealType,std::vector<RealType>> tikzPlotterCompareThickness( subDirectoryCompareThickness );
    const string mark = "x";
    string boundString = to_string( boundToCompare );
    boundString.erase ( boundString.find_last_not_of('0') + 1, std::string::npos );
    tikzPlotterCompareThickness.plotLogLogCurve ( compareThicknessVec, LargestValueOfForceWithBoundVec, "CompareThicknessLeq" + boundString + "Below", "compare thickness", mark );
    tikzPlotterCompareThickness.plotLogLogCurve ( compareThicknessVec, SmallestValueOfForceAboveBoundVec, "CompareThicknessLeq" + boundString + "Above", "compare thickness", mark );
    tikzPlotterCompareThickness.plotLogLogCurve ( compareThicknessVec, EstoredVec, "CompareThicknessLeq" + boundString + "Estored", "compare thickness Estored", mark );
}



template<bool linearEnergyType,  bool IsometryConstraint>
void compareThicknessAndForce( const int ElastEnergy, const int numAdaptiveRefinementSteps, const bool useNice ) {    
    
        ParameterParserType parser;
        
        string subDirectory;
        if( linearEnergyType ){
            if( IsometryConstraint ){
                parser.merge( "../../../../parameters/dkt/OptDeformIsometryLin.ini" );
                subDirectory = parser.createSubDirectory( "IsometryLinElast" );
            }else{
                //TODO
                //parser.merge( "../../../../parameters/dkt/OptDeformLin.ini" );
                //subDirectory = parser.createSubDirectory( "LinElast" );
            }
        }else{
            if( IsometryConstraint ){
                parser.merge( "../../../../parameters/dkt/OptDeformIsometryNonLin.ini" );
                subDirectory = parser.createSubDirectory( "IsometryNonLinElast" );
            }else{
                parser.merge( "../../../../parameters/dkt/OptDeformNonLin.ini" );
                subDirectory = parser.createSubDirectory( "NonLinElast" );
            }
        }
        
        parser.set ( "saving.saveDirectory", subDirectory );
        
        parser.set( "ConstraintProblem.ElastEnergy", ElastEnergy );
        const string ElastEnergyString = getElastEnergyString ( ElastEnergy );

    string numAdaptiveRefinementStepsString = to_string( numAdaptiveRefinementSteps );
    
    const int forceType = parser.template get<int> ( pesopt::strprintf( "Force.Type%d", 1 ).c_str() );
    RealType normForce = 1.;
    if( forceType == 11 ){
        normForce = std::abs( parser.template get<RealType> ( pesopt::strprintf( "Force.Factor%d", 1 ).c_str() ) );
    }else{
        TangentVecType forceVec; parser.template getFixSizeVector<RealType,TangentVecType> ( pesopt::strprintf( "Force.Load%d", 1 ).c_str(), forceVec );
        normForce = forceVec.norm();

    }
    string normForceString = to_string( normForce );
    normForceString.erase ( normForceString.find_last_not_of('0') + 1, std::string::npos );
    
    //
    const RealType expForce = parser.template get<RealType> ( "Force.scaleForceWithThicknessExponent" );
    string scaleForceInThicknessString = to_string( expForce );
    scaleForceInThicknessString.erase ( scaleForceInThicknessString.find_last_not_of('0') + 1, std::string::npos );
    //
    parser.addCounterToSaveDirectory( "../../../../parameters/counter.txt",
                                          "/CompareThickness_" 
                                          + ElastEnergyString 
                                          + "_NormForce" + normForceString
                                          + "_AdRef" + numAdaptiveRefinementStepsString );
    const string saveDirectory = parser.template get<string> ("saving.saveDirectory");
    const bool plotResults = parser.template get<bool> ("saving.plotResults" );
    parser.saveToFile( "ParameterParser.ini" );
    pesopt::AdditionalOutputToFile pAddOut ( pesopt::strprintf ( "%s/log.txt", parser.template get<string> ( "saving.saveDirectory" ).c_str () ) );
    
    std::vector<RealType> compareThicknessVec; parser.template getFreeSizeVector<RealType,std::vector<RealType>> ( "Compare.thicknessVec", compareThicknessVec );
    std::vector<RealType> compareFactorForceVec; parser.template getFreeSizeVector<RealType,std::vector<RealType>> ( "Compare.factorForceVec", compareFactorForceVec );
    
    for( int thicknessIter=0; thicknessIter<compareThicknessVec.size(); ++thicknessIter ){
            ParameterParserType parserIter ( parser );
            string thicknessString = to_string( compareThicknessVec[thicknessIter] );
            thicknessString.erase ( thicknessString.find_last_not_of('0') + 1, std::string::npos );
            const string saveDirectoryIter = parserIter.createSubDirectory( "Thickness" + thicknessString );
            parserIter.set ( "saving.saveDirectory", saveDirectoryIter );
            
            // set thickness in bending factor
            parserIter.set( "Material.factor_bendingEnergy", compareThicknessVec[thicknessIter] * compareThicknessVec[thicknessIter] );  
            parserIter.saveToFile( "ParameterParser.ini" );
            
            std::vector<string> saveDirectoryForceIterVec ( compareFactorForceVec.size() ), subDirectoryForceIterVec ( compareFactorForceVec.size() );
            getSaveDirectoryForceIterVec( parserIter, compareFactorForceVec, saveDirectoryForceIterVec, subDirectoryForceIterVec );
            
    #ifdef PESOPT_WITH_OPENMP   
    #pragma omp parallel for
    #endif
            for( int factorForceIter=0; factorForceIter<compareFactorForceVec.size(); ++factorForceIter ){
                cout << endl << "start for forceIter = " << factorForceIter << endl;
                ParameterParserType parserForceIter ( parserIter );
                string factorForceString = to_string( compareFactorForceVec[factorForceIter] );
                factorForceString.erase ( factorForceString.find_last_not_of('0') + 1, std::string::npos );
                const string saveDirectoryForceIter = parserForceIter.createSubDirectory( "FactorForce" + factorForceString );
                parserForceIter.set ( "saving.saveDirectory", saveDirectoryForceIter );
                
                // set factor_force in parser
                parserForceIter.set( "Force.factor_force", compareFactorForceVec[factorForceIter] );  
                parserForceIter.saveToFile( "ParameterParser.ini" );
                 
                //! SYSTEM COMMAND
                executeProgrammOptDeformForSingleForce<linearEnergyType, IsometryConstraint>( parserForceIter, ElastEnergy, numAdaptiveRefinementSteps );
                
                cout << endl << "finished for forceIter = " << factorForceIter << endl;
            }// end for factorForceIter
            
            
            //summarize for fixed thickness
            summarizeResultsToCompareForces( saveDirectoryIter, saveDirectoryForceIterVec, subDirectoryForceIterVec, compareFactorForceVec, plotResults );
            
            
    }// end for thicknessIter
    
    //plot images afterwards
    if( plotResults == false ) plotResultsToCompareThicknessAndForce<IsometryConstraint>( saveDirectory );
}


template<bool linearEnergyType, bool IsometryConstraint>
void compareLInfBound( const int ElastEnergy, const int numAdaptiveRefinementStepsSingleOpt, const bool useNice = false ) {    
    
        ParameterParserType parser;
        
        string subDirectory;
        if( IsometryConstraint ){
            parser.merge( "../../../../parameters/dkt/OptDeformIsometryLin.ini" );
            subDirectory = parser.createSubDirectory( "Isometry" );
        }
        else{
            if( linearEnergyType ){ 
            }else{ 
                parser.merge( "../../../../parameters/dkt/OptDeformNonLin.ini" );
                subDirectory = parser.createSubDirectory( "NonLinElast" );
            }
        }
        parser.set ( "saving.saveDirectory", subDirectory );
        
        parser.set( "ConstraintProblem.ElastEnergy", ElastEnergy );
        const string ElastEnergyString = getElastEnergyString( ElastEnergy );
    
    //
    parser.addCounterToSaveDirectory( "../../../../parameters/counter.txt",  "/CompareLInfBound_" + ElastEnergyString );
    const string saveDirectory = parser.template get<string> ("saving.saveDirectory");
    const int numBisections = parser.template get<int> ("Compare.numBisections");
    const string valueToCompare = parser.template get<string> ("Compare.valueToCompare" );
    const bool plotResults = parser.template get<bool> ("saving.plotResults" );
    RealType boundToCompare;
    if( valueToCompare == "LInf" )    boundToCompare = parser.template get<RealType> ("Compare.LInfBound");
    if( valueToCompare == "Estored" ) boundToCompare = parser.template get<RealType> ("Compare.EstoredBound");
    parser.saveToFile( "ParameterParser.ini" );
    pesopt::AdditionalOutputToFile pAddOut ( pesopt::strprintf ( "%s/log.txt", parser.template get<string> ( "saving.saveDirectory" ).c_str () ) );
    
    std::vector<RealType> compareThicknessVec; parser.template getFreeSizeVector<RealType,std::vector<RealType>> ( "Compare.thicknessVec", compareThicknessVec );
    std::vector<RealType> BoundVec ( compareThicknessVec.size() );
    
    #ifdef PESOPT_WITH_OPENMP   
    #pragma omp parallel for
    #endif
    for( int thicknessIter=0; thicknessIter<compareThicknessVec.size(); ++thicknessIter ){
            ParameterParserType parserIter ( parser );
            setSaveDirectoryThickness( parserIter, compareThicknessVec[thicknessIter] );
            const string saveDirectoryIter = parserIter.template get<string> ( "saving.saveDirectory" );
        
            // set thickness in bending factor
            parserIter.set( "Material.factor_bendingEnergy", compareThicknessVec[thicknessIter] * compareThicknessVec[thicknessIter] );  

            
            std::vector<RealType> compareFactorForceVec;
            compareFactorForceVec.push_back( 1. ); //TODO possibly different initial values for f_0
            
            bool existLastIterSmallerBound = false, existLastIterLargerBound = false;
            int lastIterSmallerBound, lastIterLargerBound;
            for( int factorForceIter=0; factorForceIter<=numBisections; ++factorForceIter ){
                cout << endl << "start for forceIter = " << factorForceIter << endl;
                ParameterParserType parserForceIter ( parserIter );
                setSaveDirectoryForce( parserForceIter, compareFactorForceVec[factorForceIter] );
                const string saveDirectoryForceIter = parserForceIter.template get<string> ("saving.saveDirectory" );
                
                // set factor_force in parser
                parserForceIter.set( "Force.factor_force", compareFactorForceVec[factorForceIter] );  
                parserForceIter.saveToFile( "ParameterParser.ini" );
                 
                //! SYSTEM COMMAND to execute programm
                executeProgrammOptDeformForSingleForce<linearEnergyType,IsometryConstraint>( parserForceIter, ElastEnergy, numAdaptiveRefinementStepsSingleOpt );
                
                //! get info of displacement
                DeformationOptimizationShellEnergyInfo<RealType> energyInfo;
                energyInfo.template loadFromFile<ParameterParserType>( "energyInfo", saveDirectoryForceIter + "/FinalResult" );
                RealType valueIter;
                if( valueToCompare == "LInf" )    valueIter = energyInfo._LInfNormAtNodes; 
                if( valueToCompare == "Estored" ) valueIter = energyInfo._membraneEnergy + energyInfo._bendingEnergy;
                
                if( factorForceIter < numBisections ){
                    RealType newFactorForce;
                    if(  valueIter < boundToCompare ){
                        lastIterSmallerBound = factorForceIter;
                        existLastIterSmallerBound = true;
                        if( existLastIterLargerBound ){
                            newFactorForce = 0.5 * (compareFactorForceVec[lastIterSmallerBound] + compareFactorForceVec[lastIterLargerBound]);
                        }else newFactorForce = compareFactorForceVec[lastIterSmallerBound] * 2.0;
                    }else{
                        lastIterLargerBound = factorForceIter;
                        existLastIterLargerBound = true;
                        if( existLastIterSmallerBound ){
                            newFactorForce = 0.5 * (compareFactorForceVec[lastIterSmallerBound] + compareFactorForceVec[lastIterLargerBound]);
                        }else newFactorForce = compareFactorForceVec[lastIterLargerBound] * 0.5;
                    }
                    compareFactorForceVec.push_back( newFactorForce );
                }else BoundVec[thicknessIter] = valueIter; 
                cout << endl << "finished for forceIter = " << factorForceIter << endl;
            }// end for factorForceIter
            
            
            
            //summarize forces for fixed thickness
            std::vector<RealType> compareFactorForceVecSorted ( compareFactorForceVec.size() );
            for( int forceIter=0; forceIter<compareFactorForceVec.size(); ++forceIter ) compareFactorForceVecSorted[forceIter] = compareFactorForceVec[forceIter];
            std::sort( compareFactorForceVecSorted.begin(), compareFactorForceVecSorted.end() );
            std::vector<string> saveDirectoryForceIterVecSorted, subDirectoryForceIterVecSorted;
            
            parserIter.setFixSizeVector<std::vector<RealType>> ( "Compare.factorForceVec", compareFactorForceVecSorted, false );
            parserIter.saveToFile( "ParameterParser.ini" );
            
            for( int forceIter=0; forceIter<compareFactorForceVec.size(); ++forceIter ){
                ParameterParserType parserForceIter ( parserIter );
                string factorForceString = to_string( compareFactorForceVecSorted[forceIter] );
                factorForceString.erase ( factorForceString.find_last_not_of('0') + 1, std::string::npos );
                const string saveDirectoryForceIter = parserForceIter.createSubDirectory( "FactorForce" + factorForceString );
                subDirectoryForceIterVecSorted.push_back( "FactorForce" + factorForceString );
                saveDirectoryForceIterVecSorted.push_back( saveDirectoryForceIter );
            }
            summarizeResultsToCompareForces( saveDirectoryIter, saveDirectoryForceIterVecSorted, subDirectoryForceIterVecSorted, compareFactorForceVecSorted, plotResults );
            
    }// end for thicknessIter
    
    // save summary in thickness
    TikzPlotterCurve<RealType,std::vector<RealType>> tikzPlotterCompareThickness( saveDirectory );
    const string mark = "x";
    string boundString = to_string( boundToCompare );
    boundString.erase ( boundString.find_last_not_of('0') + 1, std::string::npos );
    tikzPlotterCompareThickness.plotLogLogCurve ( compareThicknessVec, BoundVec, "CompareThickness_" + valueToCompare + "_bound" + boundString, "compare thickness for " + valueToCompare + " leq " + boundString, mark );
    
    //plot images afterwards
    plotResultsToCompareThicknessAndForce<IsometryConstraint>( saveDirectory );
    
}



//--------------------------------------------------------
int main(int argc, char ** argv) {
    
    string optimizeOrSummarizeResults = "";
    if( argc > 1 ) optimizeOrSummarizeResults = argv[1];
    else{
        cout << endl << endl
            << " Usage of programm : " << endl
            << "  1 -  compareDelAndF" << endl
            << "       compareLInfBound " << endl
            << "       plotDelAndF" << endl
            << "       sumForce" << endl
            << "       sumThickness" << endl 
            << endl;
    
        cout << "in case of compare:" << endl
            << "  2 - LinearEnergy  (options:  linelast, nonlinelast)" << endl
            << "  3 - ElastEnergy   (options:  linear:     1 - Laplace, 3 - KL" << endl
            << "                               nonlinear:  11 - DDGNonlinMembraneAndHinge, 21 - NonlinMembraneLinPlateLaplaceBend, 23 - NonlinMembraneLinKLBend 31 - FullNonlinMemBend" << endl
            << "  4 - isometry constraint (options iso, noiso)" << endl
            << "  5 - numAdaptiveRefinementSteps" << endl
            << "  6 - useNice (optional, default false)" << endl
            << endl << endl;
            
        cout << "in case of plot:" << endl
            << "  2 - EnergyType  (options: linelast, nonlinelast)" << endl
            << "  3 - isometry constraint (options: iso, noiso) )" << endl
            << "  4 - saveDirectory" << endl
            << endl << endl;
            
        cout << "in case of sum:" << endl
            << "  2 - EnergyType  (options: linelast, nonlinelast)" << endl
            << "  3 - isometry constraint (options: iso, noiso) )" << endl
            << "  4 - saveDirectory" << endl
            << "  5 - numAdaptiveRefinementSteps" << endl
            << endl << endl;
      return 23;
    }
    
    
    //=============================================
    if( optimizeOrSummarizeResults == "compareDelAndF" ){
        const string EnergyTypeString = argv[2];
        const int ElastEnergy = stoi(argv[3]);
        const string IsometryConstraintString = argv[4];
        const int numAdaptiveRefinementSteps = stoi(argv[5]);
        bool useNice = false; if( argc > 6 ){  if( stoi(argv[6]) != 0 ) useNice = true;}
        
        if( EnergyTypeString == "linelast" ){
            const bool linearElastEnergy = true;
            if( IsometryConstraintString == "iso" ){
                const bool IsometryConstraint = true;
                compareThicknessAndForce<linearElastEnergy,IsometryConstraint>( ElastEnergy, numAdaptiveRefinementSteps, useNice );
            }else{
                const bool IsometryConstraint = false;
                compareThicknessAndForce<linearElastEnergy,IsometryConstraint>( ElastEnergy, numAdaptiveRefinementSteps, useNice );
            }
        }
        if( EnergyTypeString == "nonlinelast" ){
            const bool linearElastEnergy = false;
            if( IsometryConstraintString == "iso" ){
                const bool IsometryConstraint = true;
                compareThicknessAndForce<linearElastEnergy,IsometryConstraint>( ElastEnergy, numAdaptiveRefinementSteps, useNice );
            }else{
                const bool IsometryConstraint = false;
                compareThicknessAndForce<linearElastEnergy,IsometryConstraint>( ElastEnergy, numAdaptiveRefinementSteps, useNice );
            }
        }
    }
    
    
    //=============================================
    if( optimizeOrSummarizeResults == "compareLInfBound" ){
        const string EnergyTypeString = argv[2];
        const int ElastEnergy = stoi(argv[3]);
        const string IsometryConstraintString = argv[4];
        const int numAdaptiveRefinementSteps = stoi(argv[5]);
        bool useNice = false; if( argc > 6 ){  if( stoi(argv[6]) != 0 ) useNice = true;}
        
        if( EnergyTypeString == "linelast" ){
            const bool linearElastEnergy = true;
            if( IsometryConstraintString == "iso" ){
                const bool IsometryConstraint = true;
                compareLInfBound<linearElastEnergy,IsometryConstraint>(  ElastEnergy, numAdaptiveRefinementSteps, useNice );
            }else{
                const bool IsometryConstraint = false;
                compareLInfBound<linearElastEnergy,IsometryConstraint>(  ElastEnergy, numAdaptiveRefinementSteps, useNice );
            }
        }
        if( EnergyTypeString == "nonlinelast" ){
            const bool linearElastEnergy = false;
            if( IsometryConstraintString == "iso" ){
                const bool IsometryConstraint = true;
                compareLInfBound<linearElastEnergy,IsometryConstraint>(  ElastEnergy, numAdaptiveRefinementSteps, useNice );
            }else{
                const bool IsometryConstraint = false;
                compareLInfBound<linearElastEnergy,IsometryConstraint>(  ElastEnergy, numAdaptiveRefinementSteps, useNice );
            }
        }
    }
    
    
    //=============================================
    if( optimizeOrSummarizeResults == "plotDelAndF" ){
        const string EnergyTypeString = argv[2];
        const string IsometryConstraintString = argv[3];
        const string saveDirectory = argv[4];
    
        if( EnergyTypeString == "linelast" ){
            const bool linearElastEnergy = true;
            if( IsometryConstraintString == "iso" ){
                const bool IsometryConstraint = true;
                plotResultsToCompareThicknessAndForce<IsometryConstraint>( saveDirectory );
            }else{
                const bool IsometryConstraint = false;
                plotResultsToCompareThicknessAndForce<IsometryConstraint>( saveDirectory );
            }
        }
        if( EnergyTypeString == "nonlinelast" ){
            const bool linearElastEnergy = false;
            if( IsometryConstraintString == "iso" ){
                const bool IsometryConstraint = true;
                plotResultsToCompareThicknessAndForce<IsometryConstraint>( saveDirectory );
            }else{
                const bool IsometryConstraint = false;
                plotResultsToCompareThicknessAndForce<IsometryConstraint>( saveDirectory );
            }
        }
    }
    
    if( optimizeOrSummarizeResults == "sumThickness" ){
        const string EnergyTypeString = argv[2];
        const string IsometryConstraintString = argv[3];
        const string saveDirectory = argv[4];
    
        if( EnergyTypeString == "linelast" ){
            const bool linearElastEnergy = true;
            if( IsometryConstraintString == "iso" ){
                const bool IsometryConstraint = true;
                summarizeResultsToCompareThickness<IsometryConstraint>( saveDirectory );
            }else{
                const bool IsometryConstraint = false;
                summarizeResultsToCompareThickness<IsometryConstraint>( saveDirectory );
            }
        }
        if( EnergyTypeString == "nonlinelast" ){
            const bool linearElastEnergy = false;
            if( IsometryConstraintString == "iso" ){
                const bool IsometryConstraint = true;
                summarizeResultsToCompareThickness<IsometryConstraint>( saveDirectory );
            }else{
                const bool IsometryConstraint = false;
                summarizeResultsToCompareThickness<IsometryConstraint>( saveDirectory );
            }
        }
    }
    
    
    if( optimizeOrSummarizeResults == "sumForce" ){
        const string EnergyTypeString = argv[2];
        const string IsometryConstraintString = argv[3];
        const string saveDirectory = argv[4];
        const int numAdaptiveRefinementSteps = stoi( argv[5] );
    
        ParameterParserType parser( saveDirectory + "/ParameterParser.ini" );
        std::vector<RealType> compareFactorForceVec; parser.template getFreeSizeVector<RealType,std::vector<RealType>> ( "Compare.factorForceVec", compareFactorForceVec );
        std::vector<string> saveDirectoryForceIterVec ( compareFactorForceVec.size() ), subDirectoryForceIterVec ( compareFactorForceVec.size() );
        getSaveDirectoryForceIterVec( parser, compareFactorForceVec, saveDirectoryForceIterVec, subDirectoryForceIterVec );
        
        if( EnergyTypeString == "linelast" ){
            const bool linearElastEnergy = true;
            if( IsometryConstraintString == "iso" ){
                const bool IsometryConstraint = true;
                summarizeResultsToCompareForces( saveDirectory, saveDirectoryForceIterVec, subDirectoryForceIterVec, compareFactorForceVec, true );
            }else{
                const bool IsometryConstraint = false;
                summarizeResultsToCompareForces( saveDirectory, saveDirectoryForceIterVec, subDirectoryForceIterVec, compareFactorForceVec, true );
            }
        }
        if( EnergyTypeString == "nonlinelast" ){
            const bool linearElastEnergy = false;
            if( IsometryConstraintString == "iso" ){
                const bool IsometryConstraint = true;
                summarizeResultsToCompareForces( saveDirectory, saveDirectoryForceIterVec, subDirectoryForceIterVec, compareFactorForceVec, true );
            }else{
                const bool IsometryConstraint = false;
                summarizeResultsToCompareForces( saveDirectory, saveDirectoryForceIterVec, subDirectoryForceIterVec, compareFactorForceVec, true );
            }
        }
    }
    

  return 0;
}
