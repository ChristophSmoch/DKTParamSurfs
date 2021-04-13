#ifndef __SHELLBUCKLINGTIKZPLOTTER_H
#define __SHELLBUCKLINGTIKZPLOTTER_H

#include <pesopt_IO.h>
#include <pesopt_DKT.h>



class TikzPlotterShellBucklingFileInfo {
public :
    string _fileChart, _fileUndeformed, _fileDeformed; 
    string _fileMemStress, _fileBenStress, _fileTotStress;
    string _fileEnergyVsArea;
    
    //for isometries
    string _fileStressWithoutMaterialFactor;
//     string _fileAffinePart;
    
    TikzPlotterShellBucklingFileInfo ( ) {}
    
    void setDeformationFiles( const string chart, const string undeformed, const string deformed ) { _fileChart = chart; _fileUndeformed = undeformed; _fileDeformed = deformed; };
    void setStressFiles( const string mem, const string ben, const string tot ) { _fileMemStress = mem; _fileBenStress = ben; _fileTotStress = tot; };
    void setEnergyVsAreaFile ( const string file ) { _fileEnergyVsArea = file; };
    void setStressForIsometryFiles( const string fileStressWithoutMat ){ _fileStressWithoutMaterialFactor = fileStressWithoutMat; };
//     void setStressForIsometryFiles( const string fileStressWithoutMat, const string fileAffinePart ){ _fileStressWithoutMaterialFactor = fileStressWithoutMat; _fileAffinePart = fileAffinePart; };
};




template< typename DataTypeContainer >
class TikzPlotterShellBuckling{
private:
  
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef pesopt::BoostParser ParameterParserType;
  
  TikzPlotterHelperClass<RealType> _tikzHelper;
  const string _saveDirectory;
  
public :
    
    TikzPlotterShellBuckling ( const ParameterParserType &parser ) : _saveDirectory ( parser.template get<string> ("saving.saveDirectory") ) {}
    
    TikzPlotterShellBuckling ( const string saveDirectory ) : _saveDirectory ( saveDirectory ) {}
    
    void plotDeformation( std::ofstream &out, const string fileNameChart, const string fileNameUndeformed, const string fileNameDeformed, const string infostring ) const{
        out << "\\bigskip" << endl << endl 
            << "\\begin{tcolorbox}[width=\\textwidth,colback={white},title={\\centering \\large " << infostring.c_str() << "},colbacktitle=yellow,coltitle=black]" << endl
            << "  \\resizebox{\\textwidth}{!}{" << endl
            << "  \\begin{minipage}{\\textwidth}" << endl;
        out << "    \\begin{tabular}{ c c c  }" << endl;
                
        out << "      \\begin{minipage}{0.3\\textwidth}" << endl
            << "          \\includegraphics[width=0.9\\textwidth]{" << fileNameChart.c_str() << "}" << endl
            << "      \\end{minipage}" << endl;
        out << "      &" << endl;
        out << "      \\begin{minipage}{0.3\\textwidth}" << endl
            << "          \\includegraphics[width=0.9\\textwidth]{" << fileNameUndeformed.c_str() << "}" << endl
            << "      \\end{minipage}" << endl;
        out << "      &" << endl;
        out << "      \\begin{minipage}{0.3\\textwidth}" << endl
            << "          \\includegraphics[width=0.9\\textwidth]{" << fileNameDeformed.c_str() << "}" << endl
            << "      \\end{minipage}" << endl;
        out << "      \\\\" << endl;
        out << "      chart & undeformed & deformed" << endl;
                    
        out << "  \\end{tabular}" << endl;            
        out << "  \\end{minipage}" << endl
            << " }" << endl
            << "\\end{tcolorbox}" << endl;
   }
     
    void plotStress( std::ofstream &out, const string fileNameMem, const string fileNameBen, const string fileNameTot, const string infostring ) const{
        out << "\\bigskip" << endl << endl 
            << "\\begin{tcolorbox}[width=\\textwidth,colback={white},title={\\centering \\large " << infostring.c_str() << "},colbacktitle=yellow,coltitle=black]" << endl
            << "  \\resizebox{\\textwidth}{!}{" << endl
            << "  \\begin{minipage}{\\textwidth}" << endl;
        out << "    \\begin{tabular}{ c c c  }" << endl;
                
        out << "      \\begin{minipage}{0.3\\textwidth}" << endl
            << "          \\includegraphics[width=0.9\\textwidth]{" << fileNameMem.c_str() << "}" << endl
            << "      \\end{minipage}" << endl;
        out << "      &" << endl;
        out << "      \\begin{minipage}{0.3\\textwidth}" << endl
            << "          \\includegraphics[width=0.9\\textwidth]{" << fileNameBen.c_str() << "}" << endl
            << "      \\end{minipage}" << endl;
        out << "      &" << endl;
        out << "      \\begin{minipage}{0.3\\textwidth}" << endl
            << "          \\includegraphics[width=0.9\\textwidth]{" << fileNameTot.c_str() << "}" << endl
            << "      \\end{minipage}" << endl;
        out << "      \\\\" << endl;
        out << "      membrane & bending & total" << endl;
                    
        out << "  \\end{tabular}" << endl;            
        out << "  \\end{minipage}" << endl
            << " }" << endl
            << "\\end{tcolorbox}" << endl;
   }
   
   void plotStressForIsometry( std::ofstream &out, const string fileNameStress, const string fileNameStressWithoutMat, const string infostring ) const{
        out << "\\bigskip" << endl << endl 
            << "\\begin{tcolorbox}[width=\\textwidth,colback={white},title={\\centering \\large " << infostring.c_str() << "},colbacktitle=yellow,coltitle=black]" << endl
            << "  \\resizebox{\\textwidth}{!}{" << endl
            << "  \\begin{minipage}{\\textwidth}" << endl;
        out << "    \\begin{tabular}{ c c c  }" << endl;
                
        out << "      \\begin{minipage}{0.3\\textwidth}" << endl
            << "          \\includegraphics[width=0.9\\textwidth]{" << fileNameStress.c_str() << "}" << endl
            << "      \\end{minipage}" << endl;
        out << "      &" << endl;
//         out << "      \\begin{minipage}{0.3\\textwidth}" << endl
//             << "          \\includegraphics[width=0.9\\textwidth]{" << fileNameStressWithoutMat.c_str() << "}" << endl
//             << "      \\end{minipage}" << endl;
        out << "      &" << endl;
//         out << "      \\begin{minipage}{0.3\\textwidth}" << endl
//             << "          \\includegraphics[width=0.9\\textwidth]{" << fileNameAffinePart.c_str() << "}" << endl
//             << "      \\end{minipage}" << endl;
        out << "      \\\\" << endl;
//         out << "      stress $\\chi(v) |D^2 u|^2$ & stress $|D^2 u|^2$ & affine part" << endl;
        out << "      stress $\\chi(v) | \\nabla \\nabla_h u|^2$ & & " << endl;
                    
        out << "  \\end{tabular}" << endl;            
        out << "  \\end{minipage}" << endl
            << " }" << endl
            << "\\end{tcolorbox}" << endl;
   }
   
   
    template<bool IsometryConstraint = false>
    void plotEnergyInfo ( std::ofstream &out, 
                          const DeformationOptimizationShellEnergyInfo<RealType> &energyInfo,
                          const string fileNameEnergyVsArea,
                          const string infostring
                        ) const{

        out << "\\bigskip" << endl << endl 
            << "\\begin{tcolorbox}[width=\\textwidth,colback={white},title={\\centering \\large " << infostring.c_str() << "},colbacktitle=yellow,coltitle=black]" << endl
            << "  \\resizebox{\\textwidth}{!}{" << endl
            << "  \\begin{minipage}{\\textwidth}" << endl;
        out << "  \\begin{tabular}{ c c }" << endl;
                
        out << endl
            << "   \\begin{minipage}{0.45\\textwidth}" << endl
            << "      \\begin{tabular}{ | c | c | }" << endl 
            << "          \\hline" << endl
            << "          \\multicolumn{2}{|c|}{Energies} \\\\ \\hline" << endl
            << "          $E_{\\text{pot}}$ & $ " << energyInfo._potentialEnergy << " $ \\\\ \\hline" << endl
            << "          $E_{\\text{stored}}$ & $ " << energyInfo._storedElasticEnergy << " $ \\\\ \\hline" << endl
            << "          $E_{\\text{free}}$ & $ " << energyInfo._dissipationEnergy << " $ \\\\ \\hline" << endl
            << "          $E_{\\text{mem}}$ & $ " << energyInfo._membraneEnergy << " $ \\\\ \\hline" << endl
            << "          $E_{\\text{ben}}$ & $ " << energyInfo._bendingEnergy << " $ \\\\ \\hline" << endl
            << "          residual & $ " << energyInfo._residualConstraint << " $ \\\\ \\hline" << endl
            << "          disp in $L^2$ & $ " << energyInfo._L2Norm << " $ \\\\ \\hline" << endl
            << "          disp in $L^\\infty$ (nodes) & $ " << energyInfo._LInfNormAtNodes << " $ \\\\ \\hline" << endl
            << "          disp in $L^\\infty$ (quad) & $ " << energyInfo._LInfNormAtQuadPoints << " $ \\\\ \\hline" << endl
            << "      \\end{tabular}" << endl
            << "   \\end{minipage}" << endl; 
                
        out << "   &" << endl;
                
        out << endl 
            << "  \\begin{minipage}{0.45\\textwidth}" << endl
            << "    \\includegraphics[width=0.9\\textwidth]{" << fileNameEnergyVsArea.c_str() << "}" << endl
            << "  \\end{minipage}" << endl;
                
                
        out << "\\end{tabular}" << endl;             
        out << "  \\end{minipage}" << endl;
        out << " } " << endl;
        out << "\\end{tcolorbox}" << endl;

   }
    
    
    
   template<bool IsometryConstraint = false>
   void plotAll( const string outFileName,
                 const DeformationOptimizationShellEnergyInfo<RealType> & energyInfo,
                 const TikzPlotterShellBucklingFileInfo &fileInfo ) const{
        
        std::ofstream out ( pesopt::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outFileName.c_str()  ).c_str ()  );
        
        _tikzHelper.generateIncludes( out, 21, 29.7, 0.5, 0.5, 0.5, 0.5 );
        _tikzHelper.generateBeginDocument( out );
        
        _tikzHelper.generateComment( out, "Deformation" );
        this->plotDeformation( out, fileInfo._fileChart, fileInfo._fileUndeformed, fileInfo._fileDeformed, "deformation" );

        _tikzHelper.generateComment( out, "Stress" );
        if( IsometryConstraint )  this->plotStressForIsometry( out, fileInfo._fileBenStress, fileInfo._fileStressWithoutMaterialFactor, "stress" );
        else                      this->plotStress( out, fileInfo._fileMemStress, fileInfo._fileBenStress, fileInfo._fileTotStress, "stress" );
        
        _tikzHelper.generateComment( out, "energies" );
        this->plotEnergyInfo<IsometryConstraint>( out, energyInfo, fileInfo._fileEnergyVsArea, "energies" );
        
         _tikzHelper.generateEndDocument( out );
    }
    
    
    
   void plotConvegenceOfAdaptiveRefinement( const string outFileName, const std::vector<string> &convergenceFiles) const{
        
        std::ofstream out ( pesopt::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outFileName.c_str()  ).c_str ()  );
        
        _tikzHelper.generateIncludes( out, 21, 8, 0.5, 0.5, 0.5, 0.5 );
        _tikzHelper.generateBeginDocument( out );
        
        out << "\\bigskip" << endl << endl 
            << "\\begin{tcolorbox}[width=\\textwidth,colback={white},title={\\centering \\large " << "convergence for adaptive refinement" << "},colbacktitle=yellow,coltitle=black]" << endl
            << "  \\resizebox{\\textwidth}{!}{" << endl;
            
           
        out << "    \\begin{tabular}{ "; 
        for( int i=0; i<convergenceFiles.size(); ++i ) out << "c ";             
        out << " }" << endl;
        
        for(int i=0; i<convergenceFiles.size(); ++i ){
            out << "      \\begin{minipage}{0.3\\textwidth}" << endl
                << "          \\includegraphics[width=\\textwidth]{" << convergenceFiles[i].c_str() << "}" << endl
                << "      \\end{minipage}" << endl;
            if( i < convergenceFiles.size() - 1 ) out << "      &" << endl;
            else                                  out << "      \\\\" << endl;
        }
        
        out << "  \\end{tabular}" << endl;            
        out << " }" << endl
            << "\\end{tcolorbox}" << endl;
        
         _tikzHelper.generateEndDocument( out );
    }

};




template< typename DataTypeContainer >
class TikzPlotterShellBucklingComparison{
private:
  
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef pesopt::BoostParser ParameterParserType;
  
  TikzPlotterHelperClass<RealType> _tikzHelper;
//   const ParameterParserType &_parser;
  const string _saveDirectory;
  
public :
    
  TikzPlotterShellBucklingComparison ( const ParameterParserType &parser ) : _saveDirectory ( parser.template get<string> ("saving.saveDirectory") ) {}
    
  TikzPlotterShellBucklingComparison ( const string saveDirectory ) : _saveDirectory ( saveDirectory ) {}
    
  void plotSummaryEnergies( std::ofstream &out, const std::vector<string> &energyPlotFilesVec  ) const{
       out << "\\textbf{\\huge Compare scaled forces}" << endl;
       for( int i=0; i<energyPlotFilesVec.size(); ++i ){
       out  << "\\begin{figure}[!htbp]" << endl
            << "\\resizebox{0.7\\paperwidth}{!}{" << endl
            << "\\begin{minipage}{\\textwidth} \\includegraphics{" << energyPlotFilesVec[i].c_str() << "} \\end{minipage}" << endl
            << "}" << endl
            << "\\end{figure}" << endl;
       }
  }
    
  void plotSingleForce( std::ofstream &out, const string saveDirectoryForce, const RealType scaleForce ) const{
       
       out << "\\textbf{\\huge $F=" << scaleForce << "$}" << endl
           << "\\begin{figure}[!htbp]" << endl
           << "\\resizebox{0.7\\paperwidth}{!}{" << endl
           << "  \\includegraphics{" << saveDirectoryForce.c_str() << "/FinalResult/AllResults.pdf" << "}" << endl
           << "}" << endl
           << "\\end{figure}" << endl
           << "\\begin{figure}[!htbp]" << endl
           << "\\resizebox{0.7\\paperwidth}{!}{" << endl
           << " \\begin{minipage}{\\textwidth} \\includegraphics{" << saveDirectoryForce.c_str() << "/ConvergenceAdaptiveRef/AllResultsConvAdRef.pdf" << "} \\end{minipage}" << endl
           << "}" << endl
           << "\\end{figure}" << endl;
  }

  void plotSummaryCompareForces( const string outFileName,
                                 const std::vector<string> &energyPlotFilesVec,
                                 const std::vector<string> &saveDirectoryForceIterVec,
                                 const std::vector<RealType> &compareFactorForceVec
                               ) const{
        
        std::ofstream out ( pesopt::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outFileName.c_str()  ).c_str ()  );
        
        _tikzHelper.generateIncludes( out, 21, 29.7, 0.5, 0.5, 0.5, 0.5 );
        _tikzHelper.generateBeginDocument( out );
        
        _tikzHelper.generateComment( out, "Summary" );
        this->plotSummaryEnergies( out, energyPlotFilesVec );

        for( int scaleForceIter=0; scaleForceIter<compareFactorForceVec.size(); ++scaleForceIter ){
            out << endl << "\\newpage" << endl;
            this->plotSingleForce( out, saveDirectoryForceIterVec[scaleForceIter], compareFactorForceVec[scaleForceIter] );
        }
        
         _tikzHelper.generateEndDocument( out );
    }
    
};




#endif //__SHELLBUCKLINGTIKZPLOTTER_H
