#ifndef __DKTFEDEFINES_H
#define __DKTFEDEFINES_H

#include <pesopt_IO.h>
#include <pesopt_Eigen.h>
    
    
struct DataTypeContainerDKTFE {
public :
  typedef double          RealType;
  typedef Eigen::Vector2d DomVecType;
  typedef Eigen::Vector3d TangentVecType;
  typedef Eigen::Vector3d Point3DType;
  typedef Eigen::Vector3d PointType;
  typedef Eigen::Vector3i Indices3DType;
  typedef Eigen::Vector2i Vec2i; //for example used for edges
  
  typedef Eigen::Matrix<RealType, 2, 2> Matrix22;
  typedef Eigen::Matrix<RealType, 3, 2> Matrix32;
  typedef Eigen::Matrix<RealType, 3, 3> Matrix33;
  
  typedef pesopt::Tensor222<RealType, Matrix22> Tensor222Type;
  typedef pesopt::Tensor322<RealType, Matrix22, TangentVecType> Tensor322Type;
  typedef pesopt::Tensor332<RealType, Matrix32> Tensor332Type;
  typedef pesopt::Tensor333<RealType, Matrix33> Tensor333Type;
  
  typedef Eigen::VectorXd  VectorType;
  
  typedef std::vector<bool>  MaskType;
  
  typedef Eigen::MatrixXd FullMatrixType;
  typedef Eigen::SparseMatrix<RealType, 0, int> SparseMatrixType;
 
  typedef Eigen::Triplet<RealType> TripletType;

};


#endif //__DKTFEDEFINES_H
