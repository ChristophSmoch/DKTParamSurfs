// #ifndef __PHTEXPLOTTER_H
// #define __PHTEXPLOTTER_H
// 
// #include <pesopt_IO.h>
// #include <pesoptLatexIO.h>
// #include <solverInfo.h>
// 
// 
//     
// template< typename DataTypeContainer >
// class TexPlotterMaterialOptimizationPHMultipleLoad{
// private:
//   
//   typedef typename DataTypeContainer::RealType              RealType;
//   typedef typename DataTypeContainer::PointType             PointType;
//   typedef typename DataTypeContainer::VectorType            VectorType;
//   typedef pesopt::BoostParser ParameterParserType;
//   
//   TikzPlotterHelperClass<RealType> _tikzHelper;
//   const ParameterParserType &_parser;
//   const int _precision;
//   
// public :
//     
//     TexPlotterMaterialOptimizationPHMultipleLoad ( const ParameterParserType & parser, const int precision = 5 ) :
//     _parser ( parser ), _precision ( precision ) {}
// 
//     template<typename MatOptConf>
//     void plotParameters( std::ofstream &out, const MatOptConf &matOpConf ) const{
//        typedef typename MatOptConf::ConfiguratorType ConfiguratorType;
//        
//        out << endl << endl
//            << "\\begin{tcolorbox}[width=\\textwidth,colback={white},title={\\centering \\large Parameters},colbacktitle=yellow,coltitle=black]" << endl
//            << "\\begin{minipage}{0.9\\textwidth}" << endl
//            << "  \\begin{minipage}{0.35\\textwidth}" << endl
//            << "    \\begin{tabular}{ | c | c | }" << endl
//            << "       \\hline" << endl
//            << "       \\multicolumn{2}{|c|}{Mesh} \\\\ \\hline" << endl
//            << "       num vertices & $ " << matOpConf._conf.getMesh().getNumVertices() << " $\\\\" << endl
//            << "                    & $ ( = ";
//            if( ConfiguratorType::dimChartDomain == 2 ) 
//                out <<  matOpConf._conf.getMesh().getNumDofs( 0 ) << " \\times " << matOpConf._conf.getMesh().getNumDofs( 1 ) << " ) $ \\\\ \\hline" << endl;      
//            if( ConfiguratorType::dimChartDomain == 3 ) 
//                out <<  matOpConf._conf.getMesh().getNumDofs( 0 ) << " \\times " << matOpConf._conf.getMesh().getNumDofs( 1 ) << " \\times " << matOpConf._conf.getMesh().getNumDofs( 2 ) << " ) $ \\\\ \\hline" << endl;   
//        out << "       num elements & $ " << matOpConf._conf.getMesh().getNumElements() << " $\\\\ \\hline" << endl
//            << "    \\end{tabular}" << endl
//            << "  \\end{minipage}" << endl
//            << "  %" << endl
//            << "  \\begin{minipage}{0.35\\textwidth}" << endl
//            << "    \\begin{tabular}{ | c | c | }" << endl
//            << "      \\hline" << endl
//            << "      \\multicolumn{2}{|c|}{Parameters} \\\\ \\hline" << endl
//            << "      $\\varepsilon$ (interface) & $ " << matOpConf._epsInterfaceLength << "$ \\\\" << endl
// //            << "                                 & ( $= " << matOpConf._epsFactor << "\\cdot h$ ) \\\\ \\hline" << endl
//            << "                                 & ( $= " << "TODO" << "\\cdot h$ ) \\\\ \\hline" << endl
//            << "      $c_{\\textrm{compl}}$ & $ " << matOpConf._factorComplianceCost  << " $ \\\\ \\hline" << endl
//            << "      $c_{\\textrm{interface}}$ & $ " << matOpConf._factorInterfaceCost  << " $ \\\\ \\hline" << endl
//            << "      $c_{\\textrm{doubleWell}}$ & $ " << matOpConf._factorDoubleWell  << " $ \\\\ \\hline" << endl
//            << "    \\end{tabular}" << endl
//            << "  \\end{minipage}" << endl
//            << "  %" << endl
//            << "  \\begin{minipage}{0.325\\textwidth}" << endl
//            << "    \\begin{tabular}{ | c | c | }" << endl
//            << "      \\hline" << endl
//            << "      \\multicolumn{2}{|c|}{Material} \\\\ \\hline" << endl
//            << "                  & Hard \\\\ \\hline" << endl
//            << "       $E$        & " << matOpConf._HardMaterial.getElastModulus() << "\\\\ \\hline" << endl
//            << "       $\\nu$     & " << matOpConf._HardMaterial.getPoissonRatio() << "\\\\ \\hline" << endl
//            << "       $\\mu$     & " << matOpConf._HardMaterial.getMu() << " \\\\ \\hline" << endl
//            << "       $\\lambda$ & " << matOpConf._HardMaterial.getLambda() << " \\\\ \\hline" << endl
//            << "     \\end{tabular}" << endl
//            << "   \\end{minipage}" << endl
//            << "\\end{minipage}" << endl;
//            
//            
//       //AFFINE DISPLACEMENTS
// 
//         //read loads
//            const int numAffineSymGradDofs = matOpConf._conf.numAffineSymGradDofs;
//            const int numLoads = _parser.template get<int> ( "AffineDisp.numLoads" );
//            std::vector<VectorType> affineDisp;
//            for( int i=1; i <= numLoads; ++i ){
//                VectorType affineDisptmp ( numAffineSymGradDofs );
//                _parser.template getFixSizeVector<RealType,VectorType> ( pesopt::strprintf("AffineDisp.Load%d", i ).c_str(), affineDisptmp );
//                affineDisp.push_back ( affineDisptmp ); 
//            }
//            
//            
//         //! \todo here only for 3D case
//         out << endl << "\\bigskip" << endl << endl
//             << "\\begin{minipage}{\\textwidth}" << endl            
//             << " \\begin{tabular}{ | c | c | c | }" << endl
//             << "  \\hline" << endl
//             << "  \\multicolumn{3}{|c|}{Makroskopic affine displacements $u_\\text{aff}(x) = Ax$," << endl
//             << "                        $\\qquad$" << endl 
//             << "                        $\\xi = \\begin{pmatrix} \\xi_1, & \\ldots, & \\xi_6 \\\\ \\end{pmatrix}" << endl
//             << "                        \\Rightarrow A = \\begin{pmatrix} \\xi_1 & \\xi_4 & \\xi_5 \\\\ \\xi_4 & \\xi_2 & \\xi_6 \\\\ \\xi_5 & \\xi_6 & \\xi_3 \\\\ \\end{pmatrix} $ } \\\\ \\hline" << endl;
//         for( int loadIdx = 0; loadIdx < numLoads; ++loadIdx ){ 
//           out << "    Load" << loadIdx << endl
//               << "    &" << endl
//               << "    $\\xi^{B} = \\begin{pmatrix} ";
//           for( int i=0; i<numAffineSymGradDofs; ++i ) out << (affineDisp[loadIdx])[i] << " & ";
//           out << "\\\\ \\end{pmatrix} $" << endl
//               << "    &" << endl
//               << "    \\\\ \\hline" << endl;
//         }
//         out << " \\end{tabular}" << endl
//             << "\\end{minipage}" << endl;
//         
//         
//         out << endl << "\\bigskip" << endl << endl
//             << "\\begin{minipage}{\\textwidth}" << endl            
//             << " \\begin{tabular}{ | c | }" << endl
//             << "  \\hline" << endl
//             << "$ g(x) = " << matOpConf.getWeightFunctionLoad().description().c_str() << " $ \\\\ \\hline" << endl
//             << " \\end{tabular}" << endl
//             << "\\end{minipage}" << endl;
//             
//             
//             
//             
//         out << "\\end{tcolorbox}" << endl;
// 
//     }
//     
//     
//     void plotConvergence( std::ofstream &out, const SolverInfo<DataTypeContainer> &solverInfo ) const{
//             
//         out << "   \\begin{minipage}{0.35\\textwidth}" << endl
//             << "     \\begin{tabular}{ | c | c | }" << endl
//             << "       \\hline" << endl   
//             << "       \\multicolumn{2}{|c|}{Convergence} \\\\ \\hline" << endl
//             << "       solver status & " << solverInfo.getSolverStatus().c_str() << " \\\\ \\hline" << endl
//             << "       error & $ " << solverInfo.getError() << " $ \\\\ \\hline" << endl
//             << "       num iterations & $ " << solverInfo.getNumIterations() << " $ \\\\ \\hline" << endl
//             << "     \\end{tabular}" << endl
//             << "   \\end{minipage}" << endl;
//     }
//     
//     
//     template<typename MatOptConf>
//     void plotResult( std::ofstream &out, const MatOptConf &matOpConf,
//                      const ParameterParserType & energyInfo,
//                      const string fileName, const string infostring ) const{
// 
//         out << setprecision(_precision);
//                         
//         out << "\\bigskip" << endl << endl 
//             << "\\begin{tcolorbox}[width=\\textwidth,colback={white},title={\\centering \\large "
//             << infostring.c_str() << "},colbacktitle=yellow,coltitle=black]" << endl
//             << "  \\begin{minipage}{\\textwidth}" << endl;
//         out << "     % image " << endl
//             << "     \\begin{minipage}{0.33\\textwidth}" << endl
// //TODO
// //             << "        \\includegraphics[width=\\textwidth]{" << fileName.c_str() << "}" << endl
//             << "     \\end{minipage}" << endl;
//         out << "     % total energy " << endl
//             << "     \\begin{minipage}{0.33\\textwidth}" << endl
//             << "        \\begin{tabular}{ | c | c | }" << endl 
//             << "          \\hline" << endl
//             << "          \\multicolumn{2}{|c|}{Cost Functional} \\\\ \\hline" << endl
//             << "          Compliance & $ " << energyInfo.template get<RealType> ( "Energy.compliance" ) << " $ \\\\ \\hline" << endl
//             << "          Interface & $ " << energyInfo.template get<RealType> ( "Energy.interface" ) << " $ \\\\ \\hline" << endl
//             << "          Total & $ " << matOpConf._factorComplianceCost * energyInfo.template get<RealType> ( "Energy.compliance" ) + matOpConf._factorInterfaceCost *          energyInfo.template get<RealType> ( "Energy.interface" ) << " $ \\\\ \\hline" << endl
//             << "          Volume Hard & $ " << energyInfo.template get<RealType> ("Volume.hard") << " $ \\\\ \\hline" << endl
//             << "          Barycenter & $ ( ";
//     
//       
//       PointType barycenterHard;
//       energyInfo.template getFixSizeVector<RealType,PointType> ( "Barycenter.hard", barycenterHard );
//       for( int i = 0; i < barycenterHard.size(); ++i ){
//         out << std::setprecision(2) << barycenterHard[i];
//         if( i < barycenterHard.size() - 1 ) out << ", ";
//       }
//         out << setprecision(_precision);
//         out << " ) $ \\\\ \\hline" << endl
//             << "        \\end{tabular}" << endl  
//             << "      \\end{minipage}" << endl; 
//         out << "     % compliance in detail" << endl
//             << "     \\begin{minipage}{0.3\\textwidth}" << endl
//             << "       \\begin{tabular}{ | c | c | }" << endl
//             << "         \\hline" << endl
//             << "         & $C_{*} \\xi : \\xi$  \\\\ \\hline" << endl;
//       const int numLoads = energyInfo.template get<int> ( "Energy.numLoads");
//       std::vector<RealType> complianceVec(numLoads);  
//       energyInfo.template getFixSizeVector<RealType,std::vector<RealType>>( "Energy.complianceVec", complianceVec ); 
//       for( int i=0; i<numLoads; ++i ){
//         out << "         Load $" << i << " $ & $ " << complianceVec[i]  << " $ " << " \\\\ \\hline" << endl;
//       }
//         out << "        $ g^{\\rho}(x) = $" << " &  $ " << energyInfo.template get<RealType> ( "Energy.weightCompliance" ) << "$" << " \\\\ \\hline" << endl;
//         out << "       \\end{tabular}" << endl
//             << "     \\end{minipage}" << endl;
//         out << "  \\end{minipage}" << endl;
//         out << "\\end{tcolorbox}" << endl;
// 
//    }
//     
//     
// 
//    template<typename MatOptConf>
//    void plotAll( const MatOptConf &matOpConf,
//                  const ParameterParserType & energyInfoInit, const ParameterParserType & energyInfoSolution,
//                  const SolverInfo<DataTypeContainer> &solverInfo ) const{
//         
//         std::ofstream out ( pesopt::strprintf ( "%s/Results.tex", _parser.template get<string>("saving.saveDirectory").c_str()  ).c_str ()  );
//         
//         _tikzHelper.generateIncludes( out, 21, 29.7, 0.5, 0.5, 0.5, 0.5 );
//         _tikzHelper.generateBeginDocument( out );
// 
//         _tikzHelper.generateComment( out, "Parameters" );
//         this->plotParameters<MatOptConf>( out, matOpConf );
//         
//         _tikzHelper.generateComment( out, "Initial Material" );
//         string fileNameInit = pesopt::strprintf( "InitMaterial_Undeformed.png" ).c_str();
//         this->plotResult<MatOptConf>( out, matOpConf, energyInfoInit, fileNameInit, "initial material" );
//         
//         _tikzHelper.generateComment( out, "Solution Material" );
//         const string fileNameSolution = pesopt::strprintf( "SolMaterial_Undeformed.png" ).c_str();
//         this->plotResult<MatOptConf>( out, matOpConf, energyInfoSolution, fileNameSolution, "solution material" );
//         
//         this->plotConvergence( out, solverInfo );
//         
//          _tikzHelper.generateEndDocument( out );
//     }
//     
//     
//     
//     
//     
//     template<typename MatOptConf>
//     void plotResultSingle( std::ofstream &out, const MatOptConf &matOpConf, const ParameterParserType & energyInfo,
//                     const string fileName, const string infostring, const bool plotWeightFunction = true ) const{
//         
//         out << " \\begin{minipage}{0.75\\textwidth}" << endl
// //TODO
// //             << "        \\includegraphics[width=\\textwidth]{" << fileName.c_str() << "}" << endl
//             << " \\end{minipage}" << endl;
//           
//         out << endl << "  \\bigskip" << endl << endl;
//             
//         out << " \\begin{minipage}{0.75\\textwidth}" << endl
//         << "        \\begin{tabular}{ | c | c | }" << endl 
//         << "          \\hline" << endl
//         << "          \\multicolumn{2}{|c|}{Cost Functional} \\\\ \\hline" << endl
//         << "          Compliance & $ " << energyInfo.template get<RealType> ( "Energy.compliance" ) << " $ \\\\ \\hline" << endl
//         << "          Interface & $ " << energyInfo.template get<RealType> ( "Energy.interface" ) << " $ \\\\ \\hline" << endl
//         << "          Total & $ " << matOpConf._factorComplianceCost * energyInfo.template get<RealType> ( "Energy.compliance" ) + matOpConf._factorInterfaceCost * energyInfo.template get<RealType> ( "Energy.interface" ) << " $ \\\\ \\hline" << endl
//         << "        \\end{tabular}" << endl
//         << "      \\end{minipage}" << endl; 
//                      
//         out << endl << "  \\bigskip" << endl << endl;
//         
//         //Compliance
//        const int numLoads = energyInfo.template get<int> ( "Energy.numLoads");
//        std::vector<RealType> complianceVec(numLoads);  
//        energyInfo.template getFixSizeVector<RealType,std::vector<RealType>>( "Energy.complianceVec", complianceVec ); 
//         out << "     \\begin{minipage}{0.9\\textwidth}" << endl
//             << "       \\begin{tabular}{ | c | c |  }" << endl
//             << "         \\hline" << endl
//             << "         & $C_{*} \\xi : \\xi$  \\\\ \\hline" << endl;
//         for( int i=0; i<numLoads; ++i ){
//         out << "         Load $" << i << " $ & $ " << complianceVec[i]  << " $ "
//                                       << " \\\\ \\hline" << endl;
//         }
//         if( plotWeightFunction ){
//             out << "        $ g^{\\rho}(x) = " << " $ & $ " << energyInfo.template get<RealType> ( "Energy.weightCompliance" ) << "$ \\\\ \\hline" << endl;
//         }
//         out << "       \\end{tabular}" << endl
//             << "     \\end{minipage}" << endl
//             << "     %" << endl;
//    }
//     
//     
//    template<typename MatOptConf>
//    void plotInit( const MatOptConf &matOpConf,
//                   const ParameterParserType & energyInfoInit ) const{
//         
//         std::ofstream out ( pesopt::strprintf ( "%s/ResultsInit.tex", _parser.template get<string>("saving.saveDirectory").c_str()  ).c_str ()  );
//         
//         const int numLoads = _parser.template get<int> ( "AffineDisp.numLoads" );
//         
//         _tikzHelper.generateIncludes( out, 6, 8 + static_cast<RealType> ( numLoads ) / 2., 0.0, 0.0, 0.0, 0.0 );
//         _tikzHelper.generateBeginDocument( out );
//         
//         _tikzHelper.generateComment( out, "Initial Material" );
//         string fileNameInit = pesopt::strprintf( "InitMaterial_Undeformed.png" ).c_str();
//         this->plotResultSingle<MatOptConf>( out, matOpConf, energyInfoInit, fileNameInit, "initial material", false );
//         
//          _tikzHelper.generateEndDocument( out );
//     }
//     
//     
//     
// //     void plotResults_CompareDesigns ( const string fileName, const int numDesigns ) const {
// //         
// //         std::ofstream out ( pesopt::strprintf ( "%s/%s.tex", _parser.template get<string>("saving.saveDirectory").c_str (), fileName.c_str ()  ) );
// //         
// //         const int numLoads = _parser.template get<int> ( "AffineDisp.numLoads" );
// //         
// //         _tikzHelper.generateIncludes( out, 10 + 5 * numDesigns, 16 + numLoads, 0.5, 0.5, 0.5, 0.5 );
// //         out << "\\usepackage{graphicx}" << endl
// //             << "\\usepackage{array}" << endl
// //             << "\\usepackage{colortbl}" << endl;
// //             
// //         _tikzHelper.generateBeginDocument( out );
// // 
// //         
// //         //plot Parameters (use from Design0)
// //         _tikzHelper.generateComment( out, "Parameters" );
// //         out << "\\begin{minipage}{0.9\\textwidth}" << endl
// //             << " \\include{Design-0/Parameter_0}" << endl
// //             << " \\end{minipage}" << endl << endl << endl;
// //         
// //         _tikzHelper.generateComment( out, "Initial Material" );
// //         out << "\\begin{minipage}{\\textwidth}" << endl
// //             << "\\begin{tabular}{ "; for(int i=0; i<numDesigns; ++i ) out << " | c ";
// //         out << "| }" << endl
// //             << "  \\hline" << endl;
// //         for( int i=0; i<numDesigns; ++i ){
// //              if( i > 0 ) out << "  &" << endl;
// //              string fileNamePDF = "Design-" + pesopt::strprintf( "%d", i ) + "/ResultsInit_0.pdf"; 
// //              out << "  \\begin{minipage}{" << 1. / static_cast<RealType> ( numDesigns + 2 ) << "\\textwidth}" << endl
// //                  << "      \\includegraphics[width=\\linewidth, height=10cm]{" << fileNamePDF.c_str() << "}" << endl
// //                  << "  \\end{minipage}" << endl;
// //         }
// //         out << "  \\\\ \\hline" << endl;
// //         
// //         out << "\\end{tabular}" << endl 
// //             << "\\end{minipage}" << endl;
// //             
// //      _tikzHelper.generateEndDocument( out );
// //    }
//    
// //for hom tensor
// //    template<typename MatOptConf>
// //    void generateTexCodeForPaper ( const string fileName, const ParameterParserType & energyInfo ) const {
// //        
// //        std::ofstream out ( pesopt::strprintf ( "%s/%s.tex", _parser.template get<string>("saving.saveDirectory").c_str (), fileName.c_str ()  ) );
// // 
// //        out << setprecision(_precision);
// //        
// //        out << "\\begin{figure}[!htbp]" << endl
// //            << "  \\begin{tabular}{ c  c  c }" << endl;
// //            
// //        const int numLoads = energyInfo.template get<int> ( "Energy.numLoads");
// //        std::vector<RealType> complianceVec(numLoads);  
// //        energyInfo.template getFixSizeVector<RealType,std::vector<RealType>>( "Energy.complianceVec", complianceVec ); 
// //        out << "     \\begin{minipage}{0.25\\textwidth}" << endl
// //            << "        \\begin{tabular}{ | c | c | }" << endl
// //            << "           \\hline" << endl
// //            << "              & B \\\\ \\hline" << endl;
// //        for( int i=0; i<numLoads; ++i ){
// //            out << "         Load $" << i + 1 
// //            << " $ & $ " << complianceVec[i] << " $ "
// //            << " \\\\ \\hline" << endl;
// //        }
// //        out << "         vol & " << energyInfo.template get<RealType>( "Volume.hard" ) << " \\\\ \\hline" << endl;
// //        out << "        \\end{tabular} " << endl
// //            << "      \\end{minipage}" << endl;
// //            
// //        out << "  & " << endl
// //            << "  \\begin{minipage}{0.25\\textwidth} " << endl
// //            << "    {\\includegraphics[width=0.9\\textwidth]{images/3DCompression_Cell.png}}" << endl
// //            << "  \\end{minipage}" << endl;
// //            
// //        out << "  & " << endl
// //            << "  \\begin{minipage}{0.25\\textwidth} " << endl
// //            << "    {\\includegraphics[width=0.9\\textwidth]{images/3DCompression_Block.png}}" << endl
// //            << "  \\end{minipage}" << endl;
// // 
// //        out << " \\\\ " << endl
// //            << "\\end{tabular}" << endl;
// //            
// //        out.close();
// //    }
// 
// };
// 
// 
// #endif //__PHMULTIPLELOADPLOTTER_H
