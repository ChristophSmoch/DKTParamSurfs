#ifndef __INCLUDESEIGEN_H
#define __INCLUDESEIGEN_H

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
// #include <Eigen/SparseQR>
#include <Eigen/LU>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>


#ifdef PESOPT_WITH_SUITESPARSE
#include <Eigen/UmfPackSupport>
#include <Eigen/CholmodSupport>
#endif // PESOPT_WITH_SUITESPARSE

// #include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/IterativeSolvers>


//  FullMatrixType matDense = FullMatrixType(matSparse);
//  Eigen::FullPivLU<FullMatrixType> lu(matDense);
//  FullMatrixType A_null_space = lu.kernel();
//  cout << "Null space:" << endl << A_null_space.transpose() << endl;

#endif
