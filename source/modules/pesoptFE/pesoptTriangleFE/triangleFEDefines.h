#ifndef __TRIANGLEFEDEFINES_H
#define __TRIANGLEFEDEFINES_H

#include <pesopt_IO.h>
#include <feDefines.h>

    
struct DataTypeContainerTriangleFE 
: public FEDataTypeContainerBase {
    
public :
  typedef Eigen::Vector2d                                   RealVecChart;
  typedef Eigen::Vector2d                                   TangentVecType;
  typedef Eigen::Vector2d                                   PointType;
  typedef Eigen::Vector2i                                   IntVecChart;
  typedef Matrix22                                          DerivativeVectorValuedType;
};



#endif //__TRIANGLEFEDEFINES_H
