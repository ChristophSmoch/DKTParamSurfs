#ifndef __PESOPTIMAGEIO_H
#define __PESOPTIMAGEIO_H
 

class ImageHandler {
    
public :
    ImageHandler( ) {}
    
    void trimImage( const string inputFileName, const string outputFileName ) const{
      
      string systemCommand;
      systemCommand += "convert ";
      systemCommand += inputFileName;
      systemCommand += " -trim ";
      systemCommand += outputFileName;
      bool failed;
      failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
      if ( failed ){
          cerr << "systemCommandTrimImage = " << systemCommand << endl;
          cerr << "trimImage returned an error." << endl;
      }
    }
    
    void trimAndResizeImage( const string inputFileName, const string outputFileName, 
                             const int width, const int height,
                             const bool extent = true, //if extent = false, image might be smaller then width x height
                             const string background = "white"
                           ) const{
                               
      string systemCommand;
      systemCommand += "convert ";
      systemCommand += inputFileName;
      systemCommand += " -trim ";
      systemCommand += " -resize " + std::to_string(width) + "x" + std::to_string(height) + " ";
      systemCommand += " -gravity center ";
      systemCommand += " -background " + background + " ";
      if( extent ) systemCommand += " -extent " + std::to_string(width) + "x" + std::to_string(height) + " ";
      systemCommand += outputFileName;
      bool failed; failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
      if ( failed ){
          cout << "systemCommandTrimAndResizeImage = " << systemCommand << endl;
          cerr << "trimAndResizeImage returned an error." << endl;
      }
    }

    
    void composeImages( const string inputFileName1, const string inputFileName2, const string outputFileName ) const{
        
      string systemCommand;
      systemCommand += "convert ";
      systemCommand += inputFileName1;
      systemCommand += " ";
      systemCommand += inputFileName2;
      systemCommand += " -background black -compose difference -maximum ";
      systemCommand += outputFileName;
      bool failed;
      failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
      if ( failed ){
          cerr << "systemCommandComposeImages = " << systemCommand << endl;
          cerr << "composeImages returned an error." << endl;
      }
    }
    
    
    void convertImage( const string inputFileName, const string outputFileName ) const{
      
      string systemCommand;
      systemCommand += "convert ";
      systemCommand += inputFileName;
      systemCommand += " ";
      systemCommand += outputFileName;
      bool failed;
      failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
      if ( failed ){
          cerr << "systemCommandConvertImage = " << systemCommand << endl;
          cerr << "convertImage returned an error." << endl;
      }
    }
    
    
    template <typename VectorType>
    void save2DArrayAsPGM( const int sizeX, const int sizeY,
                           const VectorType &val, 
                           const string outputFileName,
                           const int maxIntValue = 255
                         ) const{
        std::ofstream out ( outputFileName.c_str() );
        out << "P2" << endl; //Gray value type
        out << sizeX << " " << sizeY << endl;
        out << maxIntValue << endl; //0=black, maxIntValue=white
        
        auto max = val.maxCoeff();
        
        for( int i=0; i<val.size(); ++i )
            out << static_cast<int> ( val[i] * maxIntValue / max ) << " ";
        out.close();
    }

};







#endif  //  __PESOPTIMAGEIO_H
