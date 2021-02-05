#ifndef __SOLVERINFO_H
#define __SOLVERINFO_H

#include <pesopt_IO.h>


template<typename DataTypeContainer>
class SolverInfo{
    
protected:
    typedef typename DataTypeContainer::RealType  RealType;
    
    RealType _error;
    int _numIterations;
    string _solverStatus;
    RealType _runningTime;

public:
  SolverInfo( ) {}
     
  void setError( const RealType error ) { _error = error; }
  RealType getError ( ) const { return _error; }
  void setNumIterations( const int num ) { _numIterations = num; }
  int getNumIterations ( ) const { return _numIterations; }
  void setSolverStatus( const string status ) { _solverStatus = status; }
  string getSolverStatus ( ) const { return _solverStatus; }
  void setRunningTime ( const RealType time ) { _runningTime = time; }
  RealType getRunningTime ( ) const { return _runningTime; }
  
  
  void save( const string saveDirectory,  const string fileName ) const {
      
      pesopt::BoostParser infoParser;
      infoParser.set( "saving.saveDirectory", saveDirectory );
      
      infoParser.set( "SolverInfo.error", _error );
      infoParser.set( "SolverInfo.numIterations", _numIterations );
      infoParser.set( "SolverInfo.solverStatus", _solverStatus );
      infoParser.set( "SolverInfo.runningTime", _runningTime );
      
      infoParser.saveToFile( fileName );
  }
  
};  



#endif
