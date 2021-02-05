#ifndef __PESOPTBOOSTIO_H
#define __PESOPTBOOSTIO_H

#ifdef __GNUC__
#pragma GCC system_header
#endif


#include <pesoptCoreIO.h>

#ifdef PESOPT_WITH_BOOST

// #include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/optional/optional.hpp>
#include <boost/program_options.hpp>


namespace pesopt {
    
    
std::string createSubDirectory ( const std::string OldDirectory, const std::string SubDirectory ) {
   std::string NewDirectory = OldDirectory + "/" + SubDirectory;
   boost::filesystem::create_directory ( NewDirectory );
   return NewDirectory;
}
    
void moveFile( const std::string FileName, const std::string OldDirectory, const std::string NewDirectory ){
    //TODO use boost::filesystem::rename ( OldDirectory + "/" + FileName, NewDirectory + "/" + FileName );
   std::string systemCommand;
   systemCommand += "mv ";
   systemCommand += OldDirectory + "/" + FileName + " ";
   systemCommand += NewDirectory + "/" + FileName;
   bool failed;
   failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
   if ( failed ){
       std::cout << "systemCommand move file = " << systemCommand << std::endl;
       std::cerr << "moveFile returned an error." << std::endl;
   }
}

void moveDirectory( const std::string OldDirectory, const std::string NewDirectory ){
    boost::filesystem::rename ( OldDirectory, NewDirectory );
//    std::string systemCommand;
//    systemCommand += "mv ";
//    systemCommand += OldDirectory + " " + NewDirectory;
//    bool failed;
//    failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
//    if ( failed ){
//        std::cout << "systemCommand move directory = " << systemCommand << std::endl;
//        std::cerr << "moveDirectory returned an error." << std::endl;
//    }
}


void copyFile( const std::string FileName, const std::string OldDirectory, const std::string NewDirectory ){
   //TODO use boost::filesystem::copy ( OldDirectory + "/" + FileName, NewDirectory + "/" + FileName );
   string systemCommand;
   systemCommand += "cp ";
   systemCommand += OldDirectory + "/" + FileName + " ";
   systemCommand += NewDirectory + "/" + FileName;
   bool failed;
   failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
   if ( failed ){
       std::cout << "systemCommand copy file = " << systemCommand << std::endl;
       std::cerr << "copyFile returned an error." << std::endl;
   }
}


void copyDirectory( const std::string OldDirectory, const std::string NewDirectory ){
   //TODO use 
//     boost::filesystem::copy ( OldDirectory, NewDirectory );
   std::string systemCommand;
   systemCommand += "cp -r ";
   systemCommand += OldDirectory + " " + NewDirectory;
   bool failed;
   failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
   if ( failed ){
       std::cout << "systemCommand copy directory = " << systemCommand << std::endl;
       std::cerr << "copyDirectory returned an error." << std::endl;
   }
}
    
void copyDirectoryContent( const std::string OldDirectory, const std::string NewDirectory ){
   //TODO use 
//     boost::filesystem::copy ( OldDirectory, NewDirectory );
   std::string systemCommand;
   systemCommand += "cp -r ";
   systemCommand += OldDirectory + "/* " + NewDirectory;
   bool failed;
   failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
   if ( failed ){
       std::cout << "systemCommand copy content of directory = " << systemCommand << std::endl;
       std::cerr << "copyDirectoryContent returned an error." << std::endl;
   }
}
    
void deleteFile( const std::string FileName, const std::string Directory ){
   std::string systemCommand;
   systemCommand += "rm ";
   systemCommand += Directory + "/" + FileName;
   bool failed;
   failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
   if ( failed ){
       std::cout << "systemCommand delete file = " << systemCommand << std::endl;
       std::cerr << "deleteFile returned an error." << std::endl;
   }
}
   
void deleteDirectory( const std::string Directory ){
   std::string systemCommand;
   systemCommand += "rm -r ";
   systemCommand += Directory;
   bool failed;
   failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
   if ( failed ){
       std::cout << "systemCommand delete directory = " << systemCommand << std::endl;
       std::cerr << "deleteDirectory returned an error." << std::endl;
   }
}
    
    
    
    
class BoostParser {

    boost::property_tree::ptree _pt;
    
public:
    
    BoostParser ( ) {}
    
    BoostParser ( const std::string & ParFilename ){ boost::property_tree::ini_parser::read_ini ( ParFilename, _pt );}
    
    BoostParser ( const std::string & ParFilename, const int Counter, const std::string additionalName = "", bool dump = true ) : BoostParser( ParFilename ) {
        this->addCounterToSaveDirectory( Counter, additionalName );
        if( dump ) this->saveToFile( "ParameterParser.ini" );
    }
    
    BoostParser ( const std::string & ParFilename, const std::string CounterFileName, const std::string additionalName = "", bool dump = true ) : BoostParser( ParFilename ) {
        this->addCounterToSaveDirectory( CounterFileName, additionalName );
        if( dump ) this->saveToFile( "ParameterParser.ini" );
    }
    
    BoostParser& operator = ( const BoostParser& parser ) {
        // Beware of self-assignment
        if ( this == &parser ) return *this;
        _pt = parser._pt;
        return *this;
    }
    
    //TODO: only for one level: (A.B.C not allowed),
    //insert leaves (see https://stackoverflow.com/questions/8154107/how-do-i-merge-update-a-boostproperty-treeptree)
    void merge( const std::string & ParFilename ){
        boost::property_tree::ptree ptUpdates;
        boost::property_tree::ini_parser::read_ini ( ParFilename, ptUpdates );
        BOOST_FOREACH( auto& update, ptUpdates ){
            _pt.put_child( update.first, update.second );
        }
    }    
    
    template<typename VariableType>
    VariableType get ( const std::string variableName ) const { 
        VariableType val;
        try {
            val = _pt.get<VariableType> ( variableName );
        } catch (...){
            throw std::invalid_argument ( pesopt::strprintf ( "error: cannot find %s in File %s at line %d.", variableName.c_str(),  __FILE__, __LINE__ ).c_str() );
        }
        return val;
    }
  
  template<typename VariableType, typename VecType>
  void getFreeSizeVector(const std::string& variableName, VecType &vec) const{
      try {
        std::string s = _pt.template get<std::string>(variableName);   
        std::stringstream ss(s);
        std::string item;
        while(std::getline(ss, item, ',')) {
            vec.push_back(boost::lexical_cast<VariableType>(item));
        }
      } catch (...){
            throw std::invalid_argument ( pesopt::strprintf ( "error: cannot find %s in File %s at line %d.", variableName.c_str(),  __FILE__, __LINE__ ).c_str() );
      }
  }
  
  template<typename VariableType, typename VecType>
  void getFixSizeVector(const std::string& variableName, VecType &vec) const{
      try {
        std::string s = _pt.template get<std::string>(variableName);   
        std::stringstream ss(s);
        std::string item;
        int counter = 0;
        while(std::getline(ss, item, ',')) {
            vec[counter] = boost::lexical_cast<VariableType>(item);
            counter++;
        }
      } catch (...){
            throw std::invalid_argument ( pesopt::strprintf ( "error: cannot find %s in File %s at line %d.", variableName.c_str(),  __FILE__, __LINE__ ).c_str() );
      }
  }
  
  template<typename VariableType, typename MatType>
  void getFixSizeMatrix(const std::string& variableName, MatType &mat) const{
      try {
        const int numRows = mat.rows();
        const int numCols = mat.cols();
        std::string s = _pt.template get<std::string>(variableName);   
        std::stringstream ss(s);
        std::string item;
        int counter = 0;
        int row=-1, col;
        while(std::getline(ss, item, ',')) {
            if( (counter % numCols) == 0 ){
              row++;
              col = 0;   
            }else{
              col++;   
            }
            mat(row,col) = boost::lexical_cast<VariableType>(item);
            counter++;
        }
      } catch (...){
            throw std::invalid_argument ( pesopt::strprintf ( "error: cannot find %s in File %s at line %d.", variableName.c_str(),  __FILE__, __LINE__ ).c_str() );
      }
  }
  
  //checks if variable does exist
  bool hasVariable ( const std::string variableName ) const {
    boost::optional< const boost::property_tree::ptree& > child = _pt.get_child_optional( variableName );
    if( child ) return true;
    else return false;
   }
    
    
    
  template<typename VariableType>
  void set ( const std::string variableName, const VariableType & variable ) { _pt.put( variableName, variable );}
    
  template<typename VecType>
  void setFixSizeVector(const std::string& variableName, const VecType &vec, const bool useInteger ) {
      std::string s = "";
      if( useInteger ){
        for( int i=0; i<vec.size(); ++i ){
            if( i < vec.size() - 1 ) s += pesopt::strprintf ( "%d,", vec[i] ).c_str();
            else s += pesopt::strprintf ( "%d", vec[i] ).c_str();
        }
      }else{
        for( int i=0; i<vec.size(); ++i ){
            if( i < vec.size() - 1 ) s += pesopt::strprintf ( "%f,", vec[i] ).c_str();
            else s += pesopt::strprintf ( "%f", vec[i] ).c_str();
        } 
      }
      _pt.put( variableName, s );
  }
    
  void setDirectoryName( const std::string additionalName ){
        std::string newFileName = _pt.get < std::string > ( "saving.saveDirectory" );
        newFileName += additionalName;
        boost::filesystem::create_directory ( newFileName );
        _pt.put ( "saving.saveDirectory", newFileName );
  }
  
  std::string createSubDirectory ( const std::string additionalName ) const{
        std::string newFileName = _pt.get < std::string > ( "saving.saveDirectory" );
        newFileName += "/" + additionalName;
        boost::filesystem::create_directory ( newFileName );
        return newFileName;
  }
    
  void addCounterToSaveDirectory ( const int Counter, const std::string additionalName = "" ) {
        std::string newFileName = _pt.get < std::string > ( "saving.saveDirectory" );
        newFileName += additionalName;
        newFileName += "_Counter";
        newFileName += std::to_string(Counter);
        boost::filesystem::create_directory ( newFileName );
        _pt.put ( "saving.saveDirectory", newFileName );
  }
    
  int addCounterToSaveDirectory ( const std::string CounterFileName, const std::string additionalName = "" ) {
        if ( !boost::filesystem::exists( CounterFileName ) ) {
            std::ofstream out ( CounterFileName.c_str() );
            out << 0 << std::endl;
            out.close ( );
        }
        std::fstream counterFile;
        counterFile.open ( CounterFileName.c_str () );
        int counter = 0;
        if ( counterFile.is_open () ) {
            std::string temp;
            std::getline ( counterFile, temp );
            counter = atoi ( temp.c_str () );
            counterFile.seekg ( std::ios::beg );
            counterFile << ++counter;
        }
        else throw std::invalid_argument ( "cannot open counter file for writing in File" +  std::string(__FILE__) + "at line " + std::to_string( __LINE__ ) );
        counterFile.close ();
        
        this->addCounterToSaveDirectory( counter, additionalName );

        return counter;
  }
  
  void saveToFile ( const std::string& FileName ) const {
      std::string saveName = _pt.get < std::string > ( "saving.saveDirectory" );
      saveName += "/" + FileName;
      boost::property_tree::write_ini ( saveName, _pt );
  }
    
};


} // namespace pesopt

#endif  // PESOPT_WITH_BOOST

    
#endif // __PESOPTBOOSTIO_H
