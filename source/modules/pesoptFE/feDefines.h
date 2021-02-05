#ifndef __FEDEFINES_H
#define __FEDEFINES_H

#include <pesopt_IO.h>
 
//Include Eigen library
#ifdef PESOPT_WITH_EIGEN
#include <pesopt_Eigen.h>
#endif   // PESOPT_WITH_EIGEN


struct FEDataTypeContainerBase {
public :
  typedef double                                        RealType;
  
  typedef Eigen::Matrix<RealType, 2, 2>                 Matrix22;
  typedef Eigen::Matrix<RealType, 3, 2>                 Matrix32;
  typedef Eigen::Matrix<RealType, 3, 3>                 Matrix33;
 
  typedef Eigen::VectorXd                               VectorType;
  typedef std::vector<bool>                             MaskType;
  typedef Eigen::MatrixXd                               FullMatrixType;
  
  //! NOTE default is ColMajor, best performance for bicgstab is RowMajor
  typedef Eigen::SparseMatrix<RealType, 0, int>         SparseMatrixType;
  //typedef Eigen::SparseMatrix<RealType, Eigen::RowMajor, int> SparseMatrixType;
  typedef Eigen::Triplet<RealType>                      TripletType;
  //typedef pesopt::BoostParser                           ParameterParserType;

};


enum PESOPT_FEMeshType { QUOCMESH, TRIANGLE };


#endif
