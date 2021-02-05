#ifndef __QUOCFEDEFINES_H
#define __QUOCFEDEFINES_H

#include <feDefines.h>



struct Quoc1DDataTypeContainer 
: public FEDataTypeContainerBase {
public :
  typedef Eigen::Matrix<int, 1, 1>          IntVecChart;
  typedef Eigen::Matrix<RealType, 1, 1 >    RealVecChart;
  typedef Eigen::Matrix<RealType, 1, 1 >    PointType;
  typedef Eigen::Matrix<RealType, 1, 1 >    DerivativeVectorValuedType;
};

struct Quoc2DDataTypeContainer 
: public FEDataTypeContainerBase {
public :
  typedef Eigen::Vector2i                   IntVecChart;
  typedef Eigen::Matrix<RealType, 1, 1 >    RealVecChart1D;
  typedef Eigen::Vector2d                   RealVecChart;
  typedef RealType                          RealVecChartBoundary;
  typedef Eigen::Vector2d                   PointType;
  typedef Matrix22                          DerivativeVectorValuedType;
};


struct Quoc3DDataTypeContainer 
: public FEDataTypeContainerBase {
public :
  typedef Eigen::Vector3i                   IntVecChart;
  typedef Eigen::Matrix<RealType, 1, 1 >    RealVecChart1D;
  typedef Eigen::Vector3d                   RealVecChart;
  typedef Eigen::Vector2d                   RealVecChartBoundary;
  typedef Eigen::Vector3d                   PointType;
  typedef Matrix33                          DerivativeVectorValuedType;
  typedef Eigen::Matrix<RealType, 6, 6 >    VoigtMatrixType;
  typedef Eigen::Matrix<RealType,6,1>       VoigtVecType; 
};




//TODO different for 2d, 3d?
enum QuocBoundaryType {
      NOBOUNDARY = 0,
      ALL = 100,
      //2D and 3D
      LEFT = 1,
      RIGHT = 2,
      TOP = 3,
      BOTTOM = 4,
      //only 3D
      FRONT = 5,
      BACK = 6
};



#endif //__QUOCMESHDEFINES_H
