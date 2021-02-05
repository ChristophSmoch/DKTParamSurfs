#ifndef __PESOPTCOREIO_H
#define __PESOPTCOREIO_H

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC system_header
#endif


// C standard libraries
#include <cmath>
#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <initializer_list>

// STL
#include <algorithm>
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <typeinfo>
#include <unistd.h> //for sleep function
#include <unordered_set>

using namespace std;

#include <stdexcept>
#include <stdint.h>

#include <functional>

// #include <any.h>


#ifdef PESOPT_WITH_BOOST
#include <boost/filesystem.hpp>
#endif

namespace pesopt {
    
    
//! table of ansi color codes
namespace color {
const string reset       = "\033[0;0m";
const string invert      = "\033[0;7m";
const string black       = "\033[0;30m";
const string red         = "\033[0;31m";
const string green       = "\033[0;32m";
const string brown       = "\033[0;33m";
const string blue        = "\033[0;34m";
const string purple      = "\033[0;35m";
const string cyan        = "\033[0;36m";
const string light_grey  = "\033[0;37m";
const string dark_grey   = "\033[1;30m";
const string light_red   = "\033[1;31m";
const string light_green = "\033[1;32m";
const string yellow      = "\033[1;33m";
const string light_blue  = "\033[1;34m";
const string pink        = "\033[1;35m";
const string light_cyan  = "\033[1;36m";
const string white       = "\033[1;37m";
const string beep        = "\007";
const string error       = beep + red;
const string ok          = green;
}
    
    
    
//! Give back formatted string, analogously to sprintf, but save the long way 'round with char arrays.
std::string strprintf(const char * format, ...) {
  // declare variable argument list
  va_list az;
  // copy from my input into this list(second argument is nothing really used, but the last named argument of myself)
  va_start(az, format);
  // give this argument list variable to vscprintf instead of my own arguments (that is in what functions like vsprintf differ from functions like sprintf)
  const int sizeNeeded = vsnprintf ( NULL, 0, format, az ) + 1;
  // restore stack into clean state:
  va_end(az);

  char *buffer = new char[sizeNeeded];

  va_start(az, format);
  vsprintf (buffer, format, az);
  va_end(az);

  // automatic return type conversion into string:
  std::string ret = buffer;
  delete[] buffer;
  return ret;
}



std::string getFileNameExtension ( const std::string & fileName ) {
    #ifdef PESOPT_WITH_BOOST
    return boost::filesystem::extension(fileName);
    #else
    return fileName.substr(fileName.find_last_of(".") + 1);
    #endif
}
    
 //Output
template<typename VectorType>
void printVector( const VectorType &vec, const int breakNumber = 5, const std::string name = "", const unsigned precision = 5 ){

  cout << "------------------------------------------------------------------------------------------" << endl;
  cout << "-----------" << name.c_str() << "-----------------------" << endl;

  int iter = 0;
  cout.precision(precision);
  for ( int i = 0; i < vec.size(); ++i ){
      if( std::abs(vec[i]) < 1.e-16 ) cout << "0" << " ; \t "; 
      else cout << vec[i] << " ; \t ";
      if( iter == breakNumber - 1 ){
        cout << endl;
        iter = 0;
      }else iter++;
  }
  cout << "------------------------------------------------------------------------------------------" << endl;
  cout.precision(5);
}


template<typename VectorType>
void printErrorVector( const VectorType &vec, const int breakNumber = 5, const std::string name = "", const unsigned precision = 5, 
                       const double tolerance = 1.e-14 ){

  cout << "------------------------------------------------------------------------------------------" << endl;
  cout << "-----------" << name.c_str() << "-----------------------" << endl;

  int iter = 0;
  cout.precision(precision);
  for ( int i = 0; i < vec.size(); ++i ){
      if( std::abs(vec[i]) < tolerance ) cout << pesopt::color::green <<  vec[i] <<  pesopt::color::reset << " ; \t "; 
      else                               cout << pesopt::color::red << vec[i] << pesopt::color::reset << " ; \t ";
      if( iter == breakNumber - 1 ){
        cout << endl;
        iter = 0;
      }else iter++;
  }
  cout << "------------------------------------------------------------------------------------------" << endl;
  cout.precision(5);
}

template<typename VectorType>
void printVectorToFile( const VectorType &vec, const string fileName, const unsigned precision = 5 ){
  std::ofstream out ( fileName.c_str() );
  out.precision(precision);
  for ( int i = 0; i < vec.size(); ++i ){
      if( std::abs(vec[i]) < 1.e-16 ) out << "0" << endl; 
      else out << vec[i] << endl;
  }
  out.close();
}

template<typename MatrixType>
void printDenseMatrixToFile( const MatrixType &mat, const string fileName, const unsigned precision = 5 ){
  std::ofstream out ( fileName.c_str() );
  out.precision(precision);
  for ( int row = 0; row < mat.rows(); ++row ){
   for ( int col = 0; col < mat.cols(); ++col ){
      if( std::abs( mat(row,col) ) < 1.e-16 ) out << "0"; 
      else out << mat(row,col);
      out << " \t ";
   }
   out << endl;
  }
  out.close();
}

// TODO check if size of vec = size of file
template<typename VectorType>
void loadVectorFromFile(  VectorType &vec, const string fileName ){
  std::ifstream File; 
  File.open ( fileName.c_str() );
  if ( File.fail() ) 
      throw std::invalid_argument( pesopt::strprintf ( "Cannot open file with name %s. In File %s at line %d.", fileName.c_str(), __FILE__, __LINE__ ).c_str() );
  int iter = 0;
  while ( !File.eof() && iter < vec.size() ) {
      char line[256];
      File.getline ( line, 256 );  
      char * substring = strtok ( line, " " );
      vec[iter] = atof( substring );
      iter++;
  }
  File.close();
}

// TODO this requires 2x to copy vector (since result is first stored in vec)
// TODO replace double by general type
template<typename VectorType>
VectorType loadVectorFromFile( const string fileName ){
  
  std::ifstream is(fileName.c_str());
  std::istream_iterator<double> startIter(is), endIter;
  std::vector<double> vec(startIter, endIter);
  
  VectorType out( vec.size() );
  for ( int i = 0; i < vec.size(); ++i ) {
    out[i] = vec[i];   
  }
  
  return out;
  
}



template<typename SparseMatrixType>
void printSparseMatrixToFile( const SparseMatrixType &mat, const string fileName, const unsigned precision = 5 ){
  std::ofstream out ( fileName.c_str() );
  out.precision(precision);
  for (int k=0; k<mat.outerSize(); ++k)
     for (typename SparseMatrixType::InnerIterator it(mat,k); it; ++it){
         out << it.row() << " " << it.col() << " " << it.value() << endl;
  }
  out.close();
}


template<typename TripletType, typename SparseMatrixType>
void loadSparseMatrixFromFile(  SparseMatrixType &mat, const string fileName ){
  std::ifstream File; 
  File.open ( fileName.c_str() );
  if ( File.fail() ) 
      throw std::invalid_argument( pesopt::strprintf ( "Cannot open file with name %s. In File %s at line %d.", fileName.c_str(), __FILE__, __LINE__ ).c_str() );
  int iter = 0;
  std::vector<TripletType> tripletList;
  while ( !File.eof() ) {
      int row; File >> row;
      int col; File >> col;
      double val; File >> val;
      tripletList.push_back( TripletType(row, col, val ) );
  }
  File.close();
  mat.setFromTriplets( tripletList.begin(), tripletList.end() );
}







void consoleOutputStartProgramm( const std::string message) {
    cout << endl << endl << endl  
    << "================================================================================" << endl 
    << " Start Programm: " << message << endl 
    << "================================================================================" << endl << endl << endl;
}

void consoleOutput( const std::string message) {
    cout << endl << endl << endl  
    << "================================================================================" << endl 
    << message << endl 
    << "================================================================================" << endl << endl << endl;
}   
    
void consoleOutputItem( const std::string message) {
    cout << "* " << message << endl;
} 
    
    
void consoleOutputFinishProgramm( const std::string message) {
    cout << endl << endl << endl  
    << "================================================================================" << endl 
    << " Finished Programm " << message << endl 
    << "================================================================================" << endl << endl << endl;
}
    
    



//! \brief Output stream buffer that sends its content to multiple other stream buffers.
class MultiStreambuf : public streambuf {
public:
  MultiStreambuf() { setp(0, 0); }

  size_t addStreambuf(streambuf * buffer) {
    _streambufPtrList.push_back(buffer);
    return _streambufPtrList.size();
  }
  size_t removeStreambuf(streambuf * buffer) {
    _streambufPtrList.remove(buffer);
    return _streambufPtrList.size();
  }

protected:
  typedef std::char_traits<char> traits_type;
  typedef traits_type::int_type  int_type;

  typedef list<streambuf*> StreambufList;
  StreambufList _streambufPtrList;

  int_type overflow( int_type c) {
    if (!traits_type::eq_int_type(c, traits_type::eof())) {
        StreambufList::iterator iter = _streambufPtrList.begin();
        for (; iter != _streambufPtrList.end(); ++iter)
        (*iter)->sputc(c);
    }
    return traits_type::not_eof(c);
  }
  //! synchronizes all buffers. Returns for failure if at least one participation buffer has failed and returned -1. 
  int_type sync() {
    int ret = 0;
    StreambufList::iterator iter = _streambufPtrList.begin();
    for (; iter != _streambufPtrList.end(); ++iter)
        // if ret has already set to value "-1" (failed), do not overwrite this, but keep.
        if ( ret == -1 )
            (*iter)->pubsync();
        // otherwise, give *iter a chance to indicate failure.
        else
            ret = (*iter)->pubsync();
        return ret;
  }
};

//! \brief Print cout and cerr to a file.
class AdditionalOutputToFile {
public:
  AdditionalOutputToFile ( string Filename ) : 
  _filename ( Filename ), _filestream ( Filename.c_str() ), _previousCoutStreambuf ( cout.rdbuf() ), _previousCerrStreambuf ( cerr.rdbuf() ), _previousClogStreambuf ( clog.rdbuf() ) {
    _multiCoutStreambuf.addStreambuf ( _previousCoutStreambuf ); _multiCerrStreambuf.addStreambuf ( _previousCerrStreambuf ); _multiClogStreambuf.addStreambuf ( _previousClogStreambuf );
    _multiCoutStreambuf.addStreambuf ( _filestream.rdbuf() ); _multiCerrStreambuf.addStreambuf ( _filestream.rdbuf() ); _multiClogStreambuf.addStreambuf ( _filestream.rdbuf() );
    cout.rdbuf ( &_multiCoutStreambuf ); cerr.rdbuf ( &_multiCerrStreambuf ); clog.rdbuf ( &_multiClogStreambuf );
  }
  
  ~AdditionalOutputToFile () {
    cout.flush(); cerr.flush(); clog.flush();
    cout.rdbuf ( _previousCoutStreambuf ); cerr.rdbuf ( _previousCerrStreambuf ); clog.rdbuf ( _previousCerrStreambuf );
    clog << "All output has been written to file " << _filename << endl;
  }

protected:
  string      _filename;
  ofstream    _filestream;
  streambuf * _previousCoutStreambuf; streambuf * _previousCerrStreambuf; streambuf * _previousClogStreambuf;
  MultiStreambuf _multiCoutStreambuf; MultiStreambuf _multiCerrStreambuf; MultiStreambuf _multiClogStreambuf;
};



// to read csv files
class CSVRow{
    public:
        // TODO for C++17 use string_view
        std::string operator[](std::size_t index) const {
            return std::string(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
        }
        std::size_t size() const{
            return m_data.size() - 1;
        }
        void readNextRow(std::istream& str){
            std::getline(str, m_line);

            m_data.clear();
            m_data.emplace_back(-1);
            std::string::size_type pos = 0;
            while((pos = m_line.find(',', pos)) != std::string::npos)
            {
                m_data.emplace_back(pos);
                ++pos;
            }
            // This checks for a trailing comma with no data after it.
            pos   = m_line.size();
            m_data.emplace_back(pos);
        }
    private:
        std::string         m_line;
        std::vector<int>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data){
    data.readNextRow(str);
    return str;
}   

} // namespace pesopt

    
#endif // __PESOPTCOREIO_H
