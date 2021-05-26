#ifndef __DERIVATIVETESTER_H
#define __DERIVATIVETESTER_H

#include <pesopt_IO.h>
 
#include <energyDefines.h>

#include <pesoptLatexIO.h>

// Compare DE(testPoint)(testDirection) with diffQuotient 1/h ( E(testPoint + h testDirection ) - E(testPoint) )
template<typename DataTypeContainer, typename EnergyOp = pesopt::NonlinearEnergyOp<DataTypeContainer> >
class DerivativeTester {

protected:

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType MatrixType;
    
  const EnergyOp &_energyOp;
  const RealType _stepSize;
  const VectorType _testPoint;

  const int _outputLevel;
  const RealType _tolerance;
  
public:
  DerivativeTester ( const EnergyOp &E, const RealType stepSize, const VectorType & testPoint, const int outputLevel = 1, const RealType tolerance = 1.e-5  ) 
  : _energyOp ( E ), _stepSize ( stepSize ), _testPoint ( testPoint ), _outputLevel ( outputLevel ), _tolerance ( tolerance ) {}        
   
  RealType testFirstDerivative_SingleDirection ( const VectorType & testDirection ) const {
    
    // evalualte J(m)
    RealType energy;
    _energyOp.evaluateEnergy( _testPoint, energy );
    
    //evaluate J(m + h testDirection )
    RealType energyShifted;
    VectorType shiftedPoint = _testPoint + _stepSize * testDirection;    
    _energyOp.evaluateEnergy( shiftedPoint, energyShifted );
    RealType diffQuotient = ( energyShifted - energy ) / (_stepSize ) ;
    if(_outputLevel > 1 ) {
     cout << "E(p+ht) = " << energyShifted << endl;
     cout << "E(p)    = " << energy << endl;
     cout << "approxDerivative = " << diffQuotient << endl;
    }
    return diffQuotient;
  }
  
  void testFirstDerivative_AllDirections( ) const {
      
      pesopt::consoleOutput( "TEST First Derivative" );
      
      unsigned numDirections = _testPoint.size();
      
      //evaluate DJ(m)
      VectorType derivative ( numDirections );
      _energyOp.evaluateJacobian( _testPoint, derivative );
      
      // diffQuotient
      VectorType approxDerivative ( numDirections );
      
      VectorType testDirection( numDirections );
      cout << "num directions to test = " << numDirections << endl;
      pesopt::printVector<VectorType> ( _testPoint, 10, "test point" );
      for( unsigned i=0; i<numDirections; ++i ){
        if( _outputLevel > 2 )  cout << "direction = " << i << endl;
        testDirection.setZero(); testDirection[i] = 1.;
        approxDerivative[i] = testFirstDerivative_SingleDirection ( testDirection );
      }
      
      VectorType diff = derivative - approxDerivative;
      if( _outputLevel > 0 ){
          pesopt::printVector<VectorType> ( derivative, 10, "gradient" );
          pesopt::printVector<VectorType> ( approxDerivative, 10, "diffQuotient" );
          pesopt::printVector<VectorType> ( diff, 10, "grad - diffQuotient" );
      }
      int numErrors = 0;
      for( int i=0; i<numDirections; ++i ){
       if( std::abs( diff[i] ) > _tolerance ){
           numErrors++;
           cout << std::fixed << std::setprecision( 12 ) << "error in direction " << i << ": " << std::abs( diff[i] ) << endl 
                                                         << " deriv  = " << derivative[i] << endl 
                                                         << " approx = " << approxDerivative[i] << endl;
           if( std::abs( derivative[i] ) > 1.e-15 ){
              cout << std::fixed << std::setprecision( 12 ) << " relativ error in direction " << i << ": " << std::abs( diff[i] / derivative[i] ) << endl;
           }
       }
      }
      cout << endl << "num errors in first derivative = " << numErrors << endl << endl;
  }
  
    
  RealType testSecondDerivative_SingleDirection ( const VectorType & testDirection1, const VectorType & testDirection2 ) const {
      
    // evalualte DJ(m)
    VectorType derivative ( testDirection1.size() );
    _energyOp.evaluateJacobian( _testPoint, derivative );
    
    //evaluate DJ(m + h testDirection2 )
    VectorType derivativeShifted ( testDirection2.size() );
    VectorType shiftedPoint = _testPoint + _stepSize * testDirection2;    
    _energyOp.evaluateJacobian( shiftedPoint, derivativeShifted );
    RealType diffQuotient = ( derivativeShifted.dot(testDirection1) - derivative.dot(testDirection1) ) / (_stepSize ) ;
    
    return diffQuotient;
  }
  
  void testSecondDerivative_AllDirections( ) const {
      
      pesopt::consoleOutput( "TEST Second Derivative" );
      
      unsigned numDirections = _testPoint.size();
      //evaluate D2J(m)
      MatrixType hessian ( numDirections, numDirections ); 
      Eigen::MatrixXd approxHessian (numDirections, numDirections );
      _energyOp.evaluateHessian( _testPoint, hessian );
      VectorType testDirection1( numDirections ), testDirection2( numDirections );
      for( unsigned i=0; i<numDirections; ++i ){
          testDirection1.setZero();
          testDirection1[i] = 1.;
          for( unsigned j=0; j<numDirections; ++j ){
              testDirection2.setZero();
              testDirection2[j] = 1.;
              approxHessian(i,j) = testSecondDerivative_SingleDirection ( testDirection1, testDirection2 );
          }
      }
      
      Eigen::MatrixXd hessianDense =  Eigen::MatrixXd( hessian );
      Eigen::MatrixXd diff = hessianDense - approxHessian;
      if( _outputLevel > 0 ){
          cout << "hessian = " << endl << Eigen::MatrixXd( hessian ) << endl << endl;
          cout << "diffQuotient = " << endl << approxHessian << endl;
          cout << endl << endl << endl << "hessian - diffQuotient = " << endl << diff << endl;
      }
      int numErrors = 0;
      for( unsigned i=0; i<numDirections; ++i ){
          for( unsigned j=0; j<numDirections; ++j ){
            if( std::abs( diff(i,j) ) > 1.e-4 ){ 
                numErrors++;
                cout << pesopt::color::red << "error in direction (" << i << " , " << j << " ) :\t" << "hessian = " << hessianDense(i,j) << "\t approx = " << approxHessian(i,j) << pesopt::color::reset << endl; 
            }else{ 
              if( std::abs( hessianDense(i,j) ) > 1.e-16  ) 
                cout << pesopt::color::green << "correct in direction (" << i << " , " << j << " ) :\t" << "hessian = " << hessianDense(i,j) << "\t approx = " << approxHessian(i,j) << pesopt::color::reset << endl;
            }
          }
      }
      cout << endl << "num errors in first derivative = " << numErrors << endl << endl;
      
      cout << endl << "check if hessian is symmetric: " << endl;
      MatrixType hessianT, diffT; hessianT = MatrixType( hessian.transpose());
      diffT = hessian - hessianT;
      cout << "norm( hessian - hessianT ) = " << diffT.squaredNorm() << endl << endl;
   }
   
   
   void testSecondDerivative_AllDirections( const string &saveDirectory, const string &fileName ) const {
      
      unsigned numDirections = _testPoint.size();
      //evaluate D2J(m)
      MatrixType hessian ( numDirections, numDirections ); 
      Eigen::MatrixXd approxHessian (numDirections, numDirections );
      _energyOp.evaluateHessian( _testPoint, hessian );
      VectorType testDirection1( numDirections ), testDirection2( numDirections );
      for( unsigned i=0; i<numDirections; ++i ){
          testDirection1.setZero();
          testDirection1[i] = 1.;
          for( unsigned j=0; j<numDirections; ++j ){
              testDirection2.setZero();
              testDirection2[j] = 1.;
              approxHessian(i,j) = testSecondDerivative_SingleDirection ( testDirection1, testDirection2 );
          }
      }
      
      Eigen::MatrixXd hessianDense =  Eigen::MatrixXd( hessian );
      Eigen::MatrixXd diff = hessianDense - approxHessian;
      
      //plot with latex
      std::ofstream out ( pesopt::strprintf( "%s/%s.tex", saveDirectory.c_str(), fileName.c_str()  ) );
      std::ofstream outDiff ( pesopt::strprintf( "%s/%s_diff.tex", saveDirectory.c_str(), fileName.c_str()  ) );
      // out.precision(precision);
      
      TikzPlotterHelperClass<RealType> tikzHelper;
        
      tikzHelper.generateStandalone( out );
      tikzHelper.generateIncludes( out );
      tikzHelper.generateBeginDocument( out );
      out << "\\setcounter{MaxMatrixCols}{" << numDirections << "}" << endl;
      out << "\\begin{align}" << endl
          << "\\begin{pmatrix}" << endl;
          
          
      tikzHelper.generateStandalone( out );
      tikzHelper.generateIncludes( outDiff );
      tikzHelper.generateBeginDocument( outDiff );   
      outDiff << "\\setcounter{MaxMatrixCols}{" << numDirections << "}" << endl;
      outDiff << "\\begin{align}" << endl
              << "\\begin{pmatrix}" << endl;
          
     for( unsigned i=0; i<numDirections; ++i ){
          for( unsigned j=0; j<numDirections; ++j ){
            if( std::abs( diff(i,j) ) > 1.e-4 ){
                out << "{\\color{red} " << hessianDense(i,j) << " } ";
                outDiff << "{\\color{red} " << std::abs( diff(i,j) ) << " } ";
            }else{
              if( std::abs( hessianDense(i,j) ) > 1.e-16  ){
                out << "{\\color{green} " << hessianDense(i,j) << " } ";
                outDiff << "{\\color{green} " << std::abs( diff(i,j) ) << " } ";
              }else{
                out << "{\\color{black} " << hessianDense(i,j) << " } ";
                outDiff << "{\\color{black} " << std::abs( diff(i,j) ) << " } ";
              }
            }
            if (j < numDirections - 1){
                out << " & ";
                outDiff << " & ";
            }
          }
          out << " \\\\" << endl;
          outDiff << " \\\\" << endl;
      }
          
      out << "\\end{pmatrix}" << endl;
      out << "\\end{align}" << endl;
      tikzHelper.generateEndDocument( out );
      
      
      outDiff << "\\end{pmatrix}" << endl;
      outDiff << "\\end{align}" << endl;
      tikzHelper.generateEndDocument( outDiff );
      
      std::ofstream BashFilePDF ( pesopt::strprintf ( "%s/%s.sh", saveDirectory.c_str (), "generatePDF"  ) );
      BashFilePDF << "cd " << saveDirectory.c_str() << endl;
      BashFilePDF << "lualatex " << fileName.c_str() << ".tex" << endl;
      BashFilePDF << "lualatex " << fileName.c_str() << "_diff.tex" << endl;
      string systemCommand = "bash " + saveDirectory + "/generatePDF.sh";
      cout << "systemCommand = " << systemCommand << endl;
      bool failed;
      failed= ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
      if ( failed ) cerr << "programm returned an error." << endl;
      
   }
                      
};

#endif //__DERIVATIVETESTER_H
