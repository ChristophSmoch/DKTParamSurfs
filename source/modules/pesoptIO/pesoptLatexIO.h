#ifndef __PESOPTLATEXIO_H
#define __PESOPTLATEXIO_H


//! run pdflatex for given tex-file
bool runPdflatex ( const std::string& CommandFileName, const std::string& SaveDirectory = "", const bool AddDotTex = false )  {
  std::string systemCommand = "pdflatex ";
  if ( SaveDirectory != "" ) systemCommand += " -output-directory " + SaveDirectory + " ";
  systemCommand += CommandFileName;
  if( AddDotTex ) systemCommand += ".tex";
  cout << systemCommand << endl;
  bool failed = true;
  for ( unsigned short int i = 0; i < 2; ++i ) {
    failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
    if ( failed ) cerr << "runPdflatex: Calling pdflatex returned an error.\n";
  }
  return !failed;
}


bool runPdflatexWithBashFile ( const std::string& CommandFileName, const std::string& SaveDirectory = "", const bool AddDotTex = false )  {
  std::ofstream BashFilePDF ( pesopt::strprintf ( "%s/%s_%s.sh", SaveDirectory.c_str(), "generatePDF", CommandFileName.c_str()  ) );
  BashFilePDF << "cd " << SaveDirectory.c_str() << endl;
  BashFilePDF << "pdflatex " << CommandFileName.c_str();
  if( AddDotTex ) BashFilePDF << ".tex"; 
  BashFilePDF << endl;
  std::string systemCommandPDF;
  systemCommandPDF += "nohup ";
  systemCommandPDF += "bash " + SaveDirectory + "/generatePDF_" + CommandFileName + ".sh";
  systemCommandPDF += " >" + SaveDirectory + "/ConsoleOutputPrintPDF_" + CommandFileName + ".out 2>&1 ";
  //   systemCommandPDF += " &>/dev/null&";
  bool failedPDF; failedPDF = ( system ( systemCommandPDF.c_str() ) != EXIT_SUCCESS );
  if ( failedPDF ){
      cerr << "pdflatex returned an error." << endl;
      cerr << "systemCommandPDF = " << systemCommandPDF << endl;
  }
  return !failedPDF;
}


//! run lualatex for given tex-file
bool runLualatex ( const std::string& CommandFileName, const std::string& SaveDirectory = "", const bool AddDotTex = false )  {
  std::string systemCommand = "lualatex ";
  if ( SaveDirectory != "" ) systemCommand += " -output-directory " + SaveDirectory + " ";
  systemCommand += CommandFileName;
  if( AddDotTex ) systemCommand += ".tex";
  cout << systemCommand << endl;
  bool failed = true;
  for ( unsigned short int i = 0; i < 2; ++i ) {
    failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
    if ( failed ) cerr << "runLualatex: Calling lualatex returned an error.\n";
  }
  return !failed;
}


bool runLualatexWithBashFile ( const std::string& CommandFileName, const std::string& SaveDirectory = "", const bool AddDotTex = false )  {
  std::ofstream BashFilePDF ( pesopt::strprintf ( "%s/%s_%s.sh", SaveDirectory.c_str(), "generatePDF", CommandFileName.c_str()  ) );
  BashFilePDF << "cd " << SaveDirectory.c_str() << endl;
  BashFilePDF << "lualatex " << CommandFileName.c_str();
  if( AddDotTex ) BashFilePDF << ".tex"; 
  BashFilePDF << endl;
  std::string systemCommandPDF;
  systemCommandPDF += "nohup ";
  systemCommandPDF += "bash " + SaveDirectory + "/generatePDF_" + CommandFileName + ".sh";
  systemCommandPDF += " >" + SaveDirectory + "/ConsoleOutputPrintPDF_" + CommandFileName + ".out 2>&1 ";
//   systemCommandPDF += " &>/dev/null&";
  cout << "systemCommandPDF = " << systemCommandPDF << endl;
  bool failedPDF; failedPDF = ( system ( systemCommandPDF.c_str() ) != EXIT_SUCCESS );
  if ( failedPDF ){
      cerr << "lualatex returned an error." << endl;
  }
  return !failedPDF;
}




template< typename RealType >
class TikzPlotterHelperClass{
public :
    
    const std::vector<std::string> colorDictionary 
    = { "Apricot", "Cyan", "Mahogany", "Aquamarine", "DarkOrchid", "Gray", "Red", "Periwinkle", "Melon", "Brown",
        "Goldenrod", "Green", "OrangeRed", "RoyalPurple", "Turquoise", "Tan", "Thistle", "YellowOrange", "SeaGreen", "Sepia"  };
        
        
    void generateArticle ( std::ofstream &out,
         const RealType paperwidth, const RealType paperheight,
         const RealType leftBorder = 0.0, const RealType rightBorder = 0.0, const RealType topBorder = 0.0, const RealType bottomBorder = 0.0) const {

        out << "\\documentclass{article}" << endl;
        out << "\\usepackage[paperwidth=" << paperwidth << "cm,paperheight=" << paperheight << "cm]{geometry}" << endl;
        out << "\\geometry{"
            << "left=" << leftBorder << "cm," 
            << "right=" << rightBorder << "cm,"
            << "top=" << topBorder << "cm,"
            << "bottom=" << bottomBorder << "cm}" << endl;
    }
    
    void generateStandalone ( std::ofstream &out ) const {
        out << "\\documentclass{standalone}" << endl;
    }
        
    void generateIncludes ( std::ofstream &out ) const {
            
         out << "\\usepackage{amsmath}" << endl
            << "\\usepackage[dvipsnames]{xcolor}" << endl
            << "\\usepackage{pdfpages}" << endl
            << "\\usepackage{pgfplots}" << endl
            << "\\usepackage{tcolorbox}" << endl
            << "\\pgfplotsset{compat=newest}" << endl
            << "\\pagestyle{empty}" << endl << endl
            << "\\usetikzlibrary{arrows}" << endl 
            << "\\usetikzlibrary{arrows.meta}" << endl << endl;
        // only for 3d plot
        out << "\\usepackage{tikz}" << endl
            << "\\usepackage{tikz-3dplot}" << endl
            << "\\usetikzlibrary{patterns}" << endl << endl << endl;
    }
    
    void generateBeginDocument ( std::ofstream &out ) const { out << "\\begin{document}" << endl << endl;}
    void generateEndDocument ( std::ofstream &out ) const { out << endl << endl << "\\end{document}" << endl;}
    void generateNewPage ( std::ofstream &out ) const { out << "\\newpage \\clearpage" << endl << endl; };
    
    
   void generateComment ( std::ofstream &out, const std::string &comment ) const {
        out << endl << endl
            << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
            << "%%%%" << comment << endl
            << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
            << endl;
    }
    
    
    void generateBeginTikzPicture ( std::ofstream &out, const RealType resizebox, const int dimension = 2, const std::string &comment = "" ) const {
        if ( comment != "" ) out << "% " << comment << endl;
        out << "%fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff" << endl
            << "\\begin{figure}[!htbp]" << endl
            << "\\resizebox{" << resizebox << "\\paperwidth}{!}{" << endl;
        switch( dimension ){
                case 1 :{
                    out << "\\begin{tikzpicture}[x=\\paperwidth,y=\\paperwidth]" << endl;
                }break;
            
                case 2 :{
                    out << "\\begin{tikzpicture}[x=\\paperwidth,y=\\paperwidth]" << endl;
                }break;
                
                case 3 :{
                    out << "\\tdplotsetmaincoords{70}{120} %rotation in x axis and z axis" << endl
                        << "\\begin{tikzpicture}[scale=10, tdplot_main_coords]" << endl;
                }break;
                
                default :
                        throw std::invalid_argument ( pesopt::strprintf ( "wrong mode in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
         }
    }
        
    void generateEndTikzPicture ( std::ofstream &out ) const {
        out << "\\end{tikzpicture}" << endl
            << "}" << endl
            << "\\end{figure}" << endl
            << "%fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff" << endl;
    }

    void generateBarChartAxis ( std::ofstream &out ) const {
        out << "\\draw (0,0) -- (1,0);" << endl
            << "\\draw (0,0) -- (0,-0.1);" << endl
            << "\\draw (1,0) -- (1,-0.1);" << endl
            << "\\draw (-0,0) -- (-0,1);" << endl << endl;   
    }
    
};








template< typename RealType, typename VectorType >
class TikzPlotterCurve{
private:
  
  TikzPlotterHelperClass<RealType> _tikzHelper;
  const std::string _saveDirectory;
  const std::string _texType;
  
public :
    
    TikzPlotterCurve ( const std::string saveDirectory, const std::string texType = "pdflatex" ) : 
    _saveDirectory ( saveDirectory ), _texType ( texType ) { }
    
    void getMinMaxValues( const VectorType &vec, RealType &minVal, RealType &maxVal ) const {
        minVal = 2.e+19; maxVal = -2.e+19;
        for( int i=0; i<vec.size(); ++i ){
          if( vec[i] > maxVal ) maxVal = vec[i];
          if( vec[i] < minVal ) minVal = vec[i];
        }
    }
    
    void getMinMaxValues( const std::vector<VectorType> &vecs, RealType &minVal, RealType &maxVal ) const {
        minVal = 2.e+19; maxVal = -2.e+19;
        for( int i=0; i<vecs.size(); ++i ){
            for( int j=0; j<vecs[i].size(); ++j ){
                if( vecs[i][j] > maxVal ) maxVal = vecs[i][j];
                if( vecs[i][j] < minVal ) minVal = vecs[i][j];
            }
        }
    }
    
    void generateDataFile ( const VectorType &t, const VectorType &vecs, const std::string &outputFileName,
        const int MaxEntriesDataFile
    ) const {
        
        if( MaxEntriesDataFile < t.size() ){
            //save data file with all entries
            std::ofstream outDataFile_full ( pesopt::strprintf ( "%s/%s_full.dat", _saveDirectory.c_str(), outputFileName.c_str() ).c_str()  );
            outDataFile_full << "x y" << endl;
            for( int j=0; j<t.size(); ++j ) outDataFile_full << std::fixed << std::setprecision( 16 ) << t[j] << " " << vecs[j] << endl;
            outDataFile_full.close();
                    
            //generate data file with MaxEntriesDataFile+1 entries
            std::ofstream outDataFile ( pesopt::strprintf ( "%s/%s.dat", _saveDirectory.c_str(), outputFileName.c_str() ).c_str()  );
            outDataFile << "x y" << endl;
            const int size_full = t.size();
            for( int j=0; j<MaxEntriesDataFile; ++j ){
                int entry = j * size_full / MaxEntriesDataFile; 
                outDataFile << std::fixed << std::setprecision( 16 ) << t[entry] << " " << vecs[entry] << endl;
            }
            //insert last entry
                outDataFile << std::fixed << std::setprecision( 16 ) << t[size_full-1] << " " << vecs[size_full-1] << endl;
            outDataFile.close();
        }else{
            std::ofstream outDataFile ( pesopt::strprintf ( "%s/%s.dat", _saveDirectory.c_str(), outputFileName.c_str() ).c_str()  );
            outDataFile << "x y" << endl;
            for( int j=0; j<t.size(); ++j )
                outDataFile << std::fixed << std::setprecision( 16 ) << t[j] << " " << vecs[j] << endl;
            outDataFile.close();
        }
    }
    
    void generateDataFiles ( const VectorType &t, const std::vector<VectorType> &vecs, const std::string &outputFileName ) const {
        for( int i=0; i<vecs.size(); ++i ){
            std::ofstream outDataFile ( pesopt::strprintf ( "%s/%s_%d.dat", _saveDirectory.c_str(), outputFileName.c_str(), i ).c_str()  );
            outDataFile << "x y" << endl;
            for( int j=0; j<t.size(); ++j ) outDataFile << std::fixed << std::setprecision( 16 ) << t[j] << " " << vecs[i][j] << endl;
            outDataFile.close();
        }
    }
    

    void plotCurveFromDatFile ( const std::string &datFileName, const std::string &pointMeta, const std::string &outputFileName, const std::string &title ) const {
        
//         RealType min_x, max_x; this->getMinMaxValues( t, min_x, max_x );
//         RealType minVal, maxVal; this->getMinMaxValues( vec, minVal, maxVal );
        RealType min_x = -0.5, max_x = 1.0, minVal = -1.0, maxVal = 0.5;
        
        std::string texFileName = pesopt::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outputFileName.c_str()  ).c_str ();
        std::ofstream out ( texFileName.c_str()  );
        
        _tikzHelper.generateArticle( out,  20, 20 );
        _tikzHelper.generateIncludes( out );
        _tikzHelper.generateBeginDocument( out );
        _tikzHelper.generateBeginTikzPicture( out, 0.9 );
        
        //plot
        out << "\\begin{axis}" << endl
            << " [ title=" << title.c_str() << "," << endl
            << "   width=0.8\\paperwidth," << endl
            << "   at={(0.0,0)}," << endl
            << "   xmin=" << min_x  << ", xmax=" << max_x << "," << endl
            << "   ymin=" << minVal << ", ymax=" << maxVal << "," << endl
            << "   legend style={at={(0.1,0.8)},anchor=south west}," << endl
            << "    point meta = \\thisrow{" << pointMeta.c_str() << "}," << endl
            << " ]" << endl;
        out << " \\addplot[mesh,tick] table {" << pesopt::strprintf ( "%s/%s", _saveDirectory.c_str(), datFileName.c_str() ).c_str() << "};" << endl;
        out << "\\end{axis}" << endl;
        // end plot
        
        _tikzHelper.generateEndTikzPicture( out );
        _tikzHelper.generateEndDocument( out );
        
        out.close();
        
        if( _texType == "pdflatex" ) runPdflatexWithBashFile( outputFileName, _saveDirectory, true );
        if( _texType == "lualatex" ) runLualatexWithBashFile( outputFileName, _saveDirectory, true );
    }
    
    
   
   void plotCurveOfFct ( const VectorType &t, const VectorType &vec, 
                         const std::string &outputFileName, const std::string &title,
                         const std::string mark, const bool loglog,
                         const int MaxEntriesDataFile = 100
                       ) const {
        
        this->generateDataFile( t, vec, outputFileName, MaxEntriesDataFile );
        
        RealType min_x, max_x; this->getMinMaxValues( t, min_x, max_x );
        RealType minVal, maxVal; this->getMinMaxValues( vec, minVal, maxVal );
        
        std::string texFileName = pesopt::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outputFileName.c_str()  ).c_str ();
        std::ofstream out ( texFileName.c_str()  );
        
        _tikzHelper.generateArticle( out,  20, 20 );
        _tikzHelper.generateIncludes( out );
        _tikzHelper.generateBeginDocument( out );
        _tikzHelper.generateBeginTikzPicture( out, 0.9 );
        
        //plot
        if( loglog ) out << "\\begin{loglogaxis}" << endl;
        else         out << "\\begin{axis}" << endl;
        out << " [ title=" << title.c_str() << "," << endl
            << "   width=0.8\\paperwidth," << endl
            << "   at={(0.0,0)}," << endl;
        if( !loglog ){
        out << "   xmin=" << min_x  << ", xmax=" << max_x << "," << endl
            << "   ymin=" << minVal << ", ymax=" << maxVal << "," << endl;
        }
        out << "   legend style={at={(0.1,0.8)},anchor=south west}," << endl
            << " ]" << endl;
        out << " \\addplot";
        if( mark == "" ) out << "[color = black]"; 
        else             out << "[color = black, mark = " << mark << "]";
        out << " table {" 
            << pesopt::strprintf ( "%s.dat", outputFileName.c_str() ).c_str() 
            << "};" << endl;
        if( loglog ) out << "\\end{loglogaxis}" << endl;
        else         out << "\\end{axis}" << endl;
        // end plot
        
        _tikzHelper.generateEndTikzPicture( out );
        _tikzHelper.generateEndDocument( out );
        
        out.close();
        
        if( _texType == "pdflatex" ) runPdflatexWithBashFile( outputFileName, _saveDirectory, true );
        if( _texType == "lualatex" ) runLualatexWithBashFile( outputFileName, _saveDirectory, true );
    }
   
   void plotCurve ( const VectorType &t, const VectorType &vec, const std::string &outputFileName, const std::string &title, const std::string mark = "", const int MaxEntriesDataFile = 100 ) const {
       this->plotCurveOfFct( t, vec, outputFileName, title, mark, false, MaxEntriesDataFile );
    }
    
    void plotCurve ( const VectorType &vec, const std::string &outputFileName, const std::string &title, const std::string mark = "", const int MaxEntriesDataFile = 100 ) const {
        VectorType t (vec.size()); for(int i=0;i<t.size();++i) t[i] = i;
        this->plotCurveOfFct( t, vec, outputFileName, title, mark, false, MaxEntriesDataFile );
    }
    
    
    void plotLogLogCurve ( const VectorType &t, const VectorType &vec, const std::string &outputFileName, const std::string &title, const std::string mark = "", const int MaxEntriesDataFile = 100 ) const {
        this->plotCurveOfFct( t, vec, outputFileName, title, mark, true, MaxEntriesDataFile );
    }
    
    void plotLogLogCurve ( const VectorType &vec, const std::string &outputFileName, const std::string &title, const std::string mark = "", const int MaxEntriesDataFile = 100 ) const {
        VectorType t (vec.size()); for(int i=0;i<t.size();++i) t[i] = i;
        this->plotCurveOfFct( t, vec, outputFileName, title, mark, true, MaxEntriesDataFile );
    }
    
    
    void plotMutlipleCurves ( const VectorType &t, const std::vector<VectorType> &vecs, 
                              const std::string &outputFileName, const std::string &title, std::vector<std::string> &describtions) const {
        
        this->generateDataFiles( t, vecs, outputFileName );
        
        RealType min_x, max_x; this->getMinMaxValues( t, min_x, max_x );
        RealType minVal, maxVal; this->getMinMaxValues( vecs, minVal, maxVal );
        
        std::string texFileName = pesopt::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outputFileName.c_str()  ).c_str ();
        std::ofstream out ( texFileName.c_str()  );
        
        _tikzHelper.generateArticle( out,  20, 20 );
        _tikzHelper.generateIncludes( out );
        _tikzHelper.generateBeginDocument( out );
        _tikzHelper.generateBeginTikzPicture( out, 0.9 );
        
        //plot
        out << "\\begin{axis}" << endl
            << " [ title=" << title.c_str() << "," << endl
            << "   width=0.8\\paperwidth," << endl
            << "   at={(0.0,0)}," << endl
            << "   xmin=" << min_x  << ", xmax=" << max_x << "," << endl
            << "   ymin=" << minVal << ", ymax=" << maxVal << "," << endl
            << "   legend style={at={(1.1,0.0)},anchor=south west}," << endl
            << " ]" << endl;
        for( int i=0; i<vecs.size(); ++i ){
            out << " \\addplot[color = " << _tikzHelper.colorDictionary[i].c_str() << "] table {" << pesopt::strprintf ( "%s/%s_%d.dat", _saveDirectory.c_str(), outputFileName.c_str(), i ).c_str() << "};" << endl;
        }
        out << "  \\legend{";
        for( int i=0; i<vecs.size(); ++i ){
             out << describtions[i].c_str();
             if( i < vecs.size() - 1 ) out << ", ";
             else out << "}" << endl;
        }
        out << "\\end{axis}" << endl;
        // end plot
        
        _tikzHelper.generateEndTikzPicture( out );
        _tikzHelper.generateEndDocument( out );
        
        out.close();
        
        if( _texType == "pdflatex" ) runPdflatexWithBashFile( outputFileName, _saveDirectory, true );
        if( _texType == "lualatex" ) runLualatexWithBashFile( outputFileName, _saveDirectory, true );
    }
    
    
    void plotMutlipleCurves ( const std::vector<VectorType> &t, const std::vector<VectorType> &vecs, 
                              const std::string &outputFileName, const std::string &title, std::vector<std::string> &describtions, 
                              bool plotLegend = true, const int MaxEntriesDataFile = 100) const {
        
        
        for( int i=0; i<vecs.size(); ++i ){
            if( t[i].size() > MaxEntriesDataFile ){
                
                std::ofstream outDataFile_full ( pesopt::strprintf ( "%s/%s_full_%d.dat", _saveDirectory.c_str(), outputFileName.c_str(), i ).c_str()  );
                outDataFile_full << "x y" << endl;
                for( int j=0; j<t[i].size(); ++j ) outDataFile_full << std::fixed << std::setprecision( 12 ) << t[i][j] << " " << vecs[i][j] << endl;
                outDataFile_full.close();
                
                std::ofstream outDataFile ( pesopt::strprintf ( "%s/%s_%d.dat", _saveDirectory.c_str(), outputFileName.c_str(), i ).c_str()  );
                outDataFile << "x y" << endl;
                const int size_full = t[i].size();
                for( int j=0; j<MaxEntriesDataFile; ++j ){
                    int entry = j * size_full / MaxEntriesDataFile; 
                    outDataFile << std::fixed << std::setprecision( 12 ) << t[i][entry] << " " << vecs[i][entry] << endl;
                }
                outDataFile.close();
                
            } else {
                std::ofstream outDataFile ( pesopt::strprintf ( "%s/%s_%d.dat", _saveDirectory.c_str(), outputFileName.c_str(), i ).c_str()  );
                outDataFile << "x y" << endl;
                for( int j=0; j<t[i].size(); ++j ) outDataFile << std::fixed << std::setprecision( 12 ) << t[i][j] << " " << vecs[i][j] << endl;
                outDataFile.close();
            }
        }
        RealType min_x, max_x; this->getMinMaxValues( t, min_x, max_x );
        RealType minVal, maxVal; this->getMinMaxValues( vecs, minVal, maxVal );
        
        std::string texFileName = pesopt::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outputFileName.c_str()  ).c_str ();
        std::ofstream out ( texFileName.c_str()  );
        
        _tikzHelper.generateArticle( out,  20, 20 );
        _tikzHelper.generateIncludes( out );
        _tikzHelper.generateBeginDocument( out );
        _tikzHelper.generateBeginTikzPicture( out, 0.9 );
        
        //plot
        out << "\\begin{axis}" << endl
            << " [ title=" << title.c_str() << "," << endl
            << "   width=0.8\\paperwidth," << endl
            << "   at={(0.0,0)}," << endl
            << "   xmin=" << min_x  << ", xmax=" << max_x << "," << endl
            << "   ymin=" << minVal << ", ymax=" << maxVal << "," << endl
            << "   legend style={at={(1.1,0.0)},anchor=south west}," << endl
            << " ]" << endl;
        for( int i=0; i<vecs.size(); ++i ){
            out << " \\addplot[color = " << _tikzHelper.colorDictionary[i].c_str() << "] table {" << pesopt::strprintf ( "%s/%s_%d.dat", _saveDirectory.c_str(), outputFileName.c_str(), i ).c_str() << "};" << endl;
        }
        if( !plotLegend ) out << "%"; 
        out << "  \\legend{";
        for( int i=0; i<vecs.size(); ++i ){
             out << describtions[i].c_str();
             if( i < vecs.size() - 1 ) out << ", ";
             else out << "}" << endl;
        }
        out << "\\end{axis}" << endl;
        // end plot
        
        _tikzHelper.generateEndTikzPicture( out );
        _tikzHelper.generateEndDocument( out );
        
        out.close();
        
        if( _texType == "pdflatex" ) runPdflatexWithBashFile( outputFileName, _saveDirectory, true );
        if( _texType == "lualatex" ) runLualatexWithBashFile( outputFileName, _saveDirectory, true );
    }
    

};





template< typename RealType, typename VectorType >
class TikzPlotterDiagramm{
private:
  
  TikzPlotterHelperClass<RealType> _tikzHelper;
  const std::string _saveDirectory;
  
public :
    
    TikzPlotterDiagramm ( const std::string saveDirectory ) : _saveDirectory ( saveDirectory ) { }
    
    void getMinMaxValues( const VectorType &vec, RealType &minVal, RealType &maxVal ) const {
        minVal = 2.e+19; maxVal = -2.e+19;
        for( int i=0; i<vec.size(); ++i ){
          if( vec[i] > maxVal ) maxVal = vec[i];
          if( vec[i] < minVal ) minVal = vec[i];
        }
    }
    
    void getMinMaxValues( const std::vector<VectorType> &vecs, RealType &minVal, RealType &maxVal ) const {
        minVal = 2.e+19; maxVal = -2.e+19;
        for( int i=0; i<vecs.size(); ++i ){
            for( int j=0; j<vecs[i].size(); ++j ){
                if( vecs[i][j] > maxVal ) maxVal = vecs[i][j];
                if( vecs[i][j] < minVal ) minVal = vecs[i][j];
            }
        }
    }
    
    void generateDataFile ( const VectorType &vec, const std::string &outputFileName ) const {
        std::ofstream outDataFile ( pesopt::strprintf ( "%s/%s.dat", _saveDirectory.c_str(), outputFileName.c_str() ).c_str()  );
        outDataFile << "x y" << endl;
        for( int i=0; i<vec.size(); ++i ) outDataFile << std::fixed << std::setprecision( 12 ) << i << " " << vec[i] << endl;
        outDataFile.close();

    }
    
    void generateDataFile ( const VectorType &t, const VectorType &vecs, const std::string &outputFileName ) const {
        std::ofstream outDataFile ( pesopt::strprintf ( "%s/%s.dat", _saveDirectory.c_str(), outputFileName.c_str() ).c_str()  );
        outDataFile << "x y" << endl;
        for( int j=0; j<t.size(); ++j ) outDataFile << std::fixed << std::setprecision( 12 ) << t[j] << " " << vecs[j] << endl;
        outDataFile.close();
    }
    
    void generateDataFiles ( const VectorType &t, const std::vector<VectorType> &vecs, const std::string &outputFileName ) const {
        for( int i=0; i<vecs.size(); ++i ){
            std::ofstream outDataFile ( pesopt::strprintf ( "%s/%s_%d.dat", _saveDirectory.c_str(), outputFileName.c_str(), i ).c_str()  );
            outDataFile << "x y" << endl;
            for( int j=0; j<t.size(); ++j ) outDataFile << std::fixed << std::setprecision( 12 ) << t[j] << " " << vecs[i][j] << endl;
            outDataFile.close();
        }
    }
    
    
    void plotDiagramm ( const VectorType &vec, const std::string &outputFileName ) const {
        
        const RealType withBar = 0.5 / static_cast<RealType> ( vec.size() );
        
        this->generateDataFile( vec, outputFileName );
        RealType minVal, maxVal; this->getMinMaxValues( vec, minVal, maxVal );
        
        std::string texFileName = pesopt::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outputFileName.c_str()  ).c_str ();
        std::ofstream out ( texFileName.c_str()  );
        
        _tikzHelper.generateArticle( out,  20, 20 );
        _tikzHelper.generateIncludes( out );
        _tikzHelper.generateBeginDocument( out );
        _tikzHelper.generateBeginTikzPicture( out, 0.9 );
        _tikzHelper.generateBarChartAxis( out );
        
        //plot 
        for( int nodeIdx = 0; nodeIdx < vec.size(); ++nodeIdx ){
            RealType weight = vec[nodeIdx];
            RealType step = static_cast<RealType> ( nodeIdx ) / static_cast<RealType> ( vec.size() );
            out << "\\draw[fill=blue] ( " << step << ",0) rectangle ( " << step + withBar << ", " << weight << ");" << endl;
        }
        // end plot
        
        _tikzHelper.generateEndTikzPicture( out );
        _tikzHelper.generateEndDocument( out );
        
        out.close();
        
        runPdflatexWithBashFile( outputFileName, _saveDirectory, true );
    }
    
};










template< typename RealType, typename VectorType >
class TexGeneratorTable{
private:
  
  TikzPlotterHelperClass<RealType> _tikzHelper;
  const std::string _saveDirectory;
  const std::string _texType;
  
public :
    
    TexGeneratorTable ( const std::string saveDirectory, const std::string texType = "pdflatex" ) : 
    _saveDirectory ( saveDirectory ), _texType ( texType ) { }
    

    void plotTable ( const std::vector<VectorType> &vecs, 
                     const std::string &outputFileName, const std::string &title, 
                     std::vector<std::string> &describtions) const {
        

        const int numRows = describtions.size();
        const int numCols = vecs[0].size();
                         
        std::string texFileName = pesopt::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outputFileName.c_str()  ).c_str ();
        std::ofstream out ( texFileName.c_str()  );
        
        _tikzHelper.generateArticle( out,  20, 20 );
        _tikzHelper.generateIncludes( out );
        _tikzHelper.generateBeginDocument( out );
        
        out << "\\begin{tcolorbox}[width=\\textwidth,colback={white},title={\\centering \\large " << title.c_str() << "},colbacktitle=yellow,coltitle=black]" << endl
            << "\\resizebox{\\textwidth}{!}{" << endl;
        
        out << "\\begin{tabular}{ |";
        for( int i=0; i<numRows; ++i ) out << " c | ";
        out << "  }" << endl;
        
        out << "\\hline" << endl;
          
        for( int i=0; i<numRows; ++i ){
            out << describtions[i].c_str();
            if( i <numRows - 1 ) out << " & ";
            else                 out << " \\\\ \\hline" << endl;
        }
        
        
           
        for( int j=0; j<numCols; ++j ){
            for( int i=0; i<numRows; ++i ){
                out << (vecs[i])[j];
                if( i <numRows - 1 ) out << " & ";
                else                 out << " \\\\ \\hline" << endl;
            }
            
        }
     
        out << "\\end{tabular}" << endl;
        
        out << "}" << endl
            << "\\end{tcolorbox}" << endl;
        
        _tikzHelper.generateEndDocument( out );
        
        out.close();
        
        if( _texType == "pdflatex" ) runPdflatexWithBashFile( outputFileName, _saveDirectory, true );
        if( _texType == "lualatex" ) runLualatexWithBashFile( outputFileName, _saveDirectory, true );
    }

};

#endif //__PESOPTLATEXIO_H
