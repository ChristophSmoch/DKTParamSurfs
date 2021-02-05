 // This file is part of Eigen, a lightweight C++ template library
 // for linear algebra.
 //
 // Copyright (C) 2011-2014 Gael Guennebaud <gael.guennebaud@inria.fr>
 // Copyright (C) 2012 Désiré Nuentsa-Wakam <desire.nuentsa_wakam@inria.fr>
 //
 // This Source Code Form is subject to the terms of the Mozilla
 // Public License v. 2.0. If a copy of the MPL was not distributed
 // with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
 
 
 #include <fstream>
 
 #ifndef EIGEN_BICGSTAB_NEW_H
 #define EIGEN_BICGSTAB_NEW_H
 
 namespace Eigen { 
 
 namespace internal {
 
 template<typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
 bool bicgstab_new(const MatrixType& mat, const Rhs& rhs, Dest& x,
                   const Preconditioner& precond, Index& iters,
                   typename Dest::RealScalar& tol_error,
                   const bool output, std::ofstream & outstream )
 {
   typedef typename Dest::RealScalar RealScalar;
   typedef typename Dest::Scalar Scalar;
   typedef Matrix<Scalar,Dynamic,1> VectorType;
   RealScalar tol = tol_error;
   Index maxIters = iters;
 
   Index n = mat.cols();
   VectorType r  = rhs - mat * x;
   VectorType r0 = r;
   
   RealScalar r0_sqnorm = r0.squaredNorm();
   RealScalar rhs_sqnorm = rhs.squaredNorm();

   if( output ){
    outstream << std::endl << std::endl << "in bicgstab_new" << std::endl;
    outstream << "init.norm = " << x.norm() << std::endl;
    outstream << "rhs_sqnorm = " << rhs_sqnorm << std::endl;
    outstream << "tol = " << tol << std::endl;
    outstream << "tol * tol = " << tol * tol << std::endl;
   }
   
   //this is the modification of Eigen bicgstab
   if(rhs_sqnorm < tol * tol){
     x.setZero();
     return true;
   }
   
   Scalar rho    = 1;
   Scalar alpha  = 1;
   Scalar w      = 1;
   
   VectorType v = VectorType::Zero(n), p = VectorType::Zero(n);
   VectorType y(n),  z(n);
   VectorType kt(n), ks(n);
 
   VectorType s(n), t(n);
 
   RealScalar tol2 = tol*tol;
   if( rhs_sqnorm > 1. ) tol2 *= rhs_sqnorm;
   RealScalar eps2 = NumTraits<Scalar>::epsilon()*NumTraits<Scalar>::epsilon();
   Index i = 0;
   Index restarts = 0;

   Scalar minResidual = r.norm();
   Index minResidualIndex = 0;
   VectorType minSolution ( n );
   
   if( output ){
      outstream << "r.squaredNorm = " << r.squaredNorm() << std::endl 
            << "tol2 = " << tol2 << std::endl << std::endl;
   }

   while ( r.squaredNorm() > tol2 && i<maxIters ){
       
    if( output ) outstream << "iter = " << i << std::endl;
       
     Scalar rho_old = rho;
 
     rho = r0.dot(r);
     if( output ) outstream << "<r0,r> = " << std::abs(rho) << std::endl;    
     if (std::abs(rho) < eps2*r0_sqnorm){
       if( output ){
        outstream << " new residual vector to orhtogonal: <r0,r> = " << std::abs(rho) << std::endl;    
       }
       // The new residual vector became too orthogonal to the arbitrarily chosen direction r0
       // Let's restart with a new r0:
       r  = rhs - mat * x;
       r0 = r;
       rho = r0_sqnorm = r.squaredNorm();
       if(restarts++ == 0) i = 0;
     }
     Scalar beta = (rho/rho_old) * (alpha / w);
     p = r + beta * (p - w * v);
     
     y = precond.solve(p);
     
     v.noalias() = mat * y;
 
     alpha = rho / r0.dot(v);
     s = r - alpha * v;
 
     z = precond.solve(s);
     t.noalias() = mat * z;
 
     RealScalar tmp = t.squaredNorm();
     if(tmp>RealScalar(0)) w = t.dot(s) / tmp;
     else w = Scalar(0);
     x += alpha * y + w * z;
     r = s - w * t;
     
     Scalar currentResidual = r.norm();
     if( currentResidual < minResidual ){
       minResidual = currentResidual;
       minResidualIndex = i;
       minSolution = x;
     }
     
     if( output ){
        outstream << "tmp = " << tmp << std::endl;
        outstream << "residual = " << r.norm() << std::endl;
        outstream << "residual sqr = " << r.squaredNorm() << std::endl;
     }

     ++i;
   }   
   
   if( output ){
    outstream << std::endl << std::endl
         << "minResidual = " << minResidual << std::endl
         << "minResidualIndex = " << minResidualIndex << std::endl;
   }

   x = minSolution;
   
   if( rhs_sqnorm > 1. ) tol_error = std::sqrt(minResidual * minResidual/rhs_sqnorm);
   else tol_error = minResidual;
   iters = minResidualIndex;
   
//    if( rhs_sqnorm > 1. ) tol_error = std::sqrt(r.squaredNorm()/rhs_sqnorm);
//    else tol_error = std::sqrt(r.squaredNorm());
//    iters = i;
   
   return true; 
 }
 
 }
 
 template< typename _MatrixType,
           typename _Preconditioner = DiagonalPreconditioner<typename _MatrixType::Scalar> >
 class BiCGSTAB_NEW;
 
 namespace internal {
 
 template< typename _MatrixType, typename _Preconditioner>
 struct traits<BiCGSTAB_NEW<_MatrixType,_Preconditioner> >
 {
   typedef _MatrixType MatrixType;
   typedef _Preconditioner Preconditioner;
 };
 
 }
 
 template< typename _MatrixType, typename _Preconditioner>
 class BiCGSTAB_NEW : public IterativeSolverBase<BiCGSTAB_NEW<_MatrixType,_Preconditioner> >
 {
   typedef IterativeSolverBase<BiCGSTAB_NEW> Base;
   using Base::matrix;
   using Base::m_error;
   using Base::m_iterations;
   using Base::m_info;
   using Base::m_isInitialized;
   
//    int _outputLevel;
   bool _output;
   std::ofstream *_outstream;
   
 public:
   typedef _MatrixType MatrixType;
   typedef typename MatrixType::Scalar Scalar;
   typedef typename MatrixType::RealScalar RealScalar;
   typedef _Preconditioner Preconditioner;
 
 public:
 
   BiCGSTAB_NEW() : Base(), _output(false), _outstream(NULL) {}
 
   template<typename MatrixDerived>
   explicit BiCGSTAB_NEW(const EigenBase<MatrixDerived>& A) : Base(A.derived()), _output(false), _outstream(NULL) {}
 
   ~BiCGSTAB_NEW() {}
 
//    int getOutputLevel( ) const { return _outputLevel;}
//    void setOutputLevel( int level ) {_outputLevel = level;}
 
   void writeOutput( const bool output ) { _output = output; }
   void setOutputStream( std::ofstream &out, const int precision = 32 ) { _outstream = &out; _output = true; _outstream->precision( precision );  }
 
   template<typename Rhs,typename Dest>
   void _solve_with_guess_impl(const Rhs& b, Dest& x) const
   {    
     bool failed = false;
     for(Index j=0; j<b.cols(); ++j)
     {
       m_iterations = Base::maxIterations();
       m_error = Base::m_tolerance;
       
       typename Dest::ColXpr xj(x,j);
       if(!internal::bicgstab_new(matrix(), b.col(j), xj, Base::m_preconditioner, m_iterations, m_error, _output, *_outstream)) failed = true;
     }
     //use acceptTol = 100 * m_tolerance
     Scalar acceptTol = Base::m_tolerance * 100.;
     m_info = failed ? NumericalIssue
            : (m_error <= acceptTol ) ? Success
            : NoConvergence;
     m_isInitialized = true;
   }
 
   using Base::_solve_impl;
   template<typename Rhs,typename Dest>
   void _solve_impl(const MatrixBase<Rhs>& b, Dest& x) const
   {
     x.resize(this->rows(),b.cols());
     x.setZero();
     _solve_with_guess_impl(b,x);
   }
 
 protected:
 
 };
 
 } // end namespace Eigen
 
 #endif // EIGEN_BICGSTAB_NEW_H
