// /*
//  * ahmedPGMRes.h
//  *  Contains a linear solver using the preconditioned GMRES method of the AHMED library.
//  */
// 
// #ifndef AHMEDPGMRES_H_
// #define AHMEDPGMRES_H_
// 
// #include <iostream>
// #include <sstream>
// #include <string>
// #include <vector>
// 
// #ifdef USE_EXTERNAL_AHMED
// #include <ahmedIncludes.h>
// 
// 
// //! \note AHMED CRS matrices are quadratic compressed row format matrices.
// template<typename ConfiguratorType, typename PrecondType>
// class AHMEDPGMRes :
// public AHMED::CRSMatrix<typename ConfiguratorType::RealType> {
// 
//   typedef typename ConfiguratorType::RealType RealType;
//   typedef typename ConfiguratorType::VectorType VectorType;
// 
//   void createLUDecomposition ( unsigned int rankMax, RealType eps, AHMED::CRSMatrix<RealType> crsMat, AHMED::BlockClusterTree_sptr &blTree, bool preCoarse ) {
//       
// //    AHMED::displayInfo ( *blTree );
//     const AHMED::SvdParam param ( eps, rankMax, AHMED::relative, AHMED::Euclidean );
// 
//     AHMED::CRS_perm ( crsMat, _perm );
//     // AHMED::CRS_perm ( *this, _perm );
// 
//     const AHMED::SvdParam paramCRS ( param.eps * 1E-4, param.maxRank, param.typeOfAcc, param.typeOfNrm );
// 
//     AHMED::GeHmat<PrecondType> workingCopy ( blTree );
// 
//     workingCopy.add ( crsMat, paramCRS );
//     // workingCopy.add ( *this, paramCRS );
// 
//     _L = new AHMED::LtHmat<PrecondType> ( blTree );
//     _U = new AHMED::UtHmat<PrecondType> ( blTree );
// 
// #define AHMED_PRESERVATION
// #ifdef AHMED_PRESERVATION
//     // create haar basis
//     VectorType presVec ( 2 * crsMat.numRows () );
// 
//     // generate basis
//     VectorType basis ( crsMat.numRows () );
//     for( int i=0; i<basis.size(); ++i ) basis[i] = 1.;
// 
//     // constants to be set
//     // para1: dimension
//     // para2: higher precision for not admissible blocks (approximation accuracy)
//     // para3: threshold for vectors to be preserved
//     // para4: number of unkowns
//     AHMED::contConst constants ( static_cast<unsigned int>( ConfiguratorType::DimOfWorld ), 1.0, 1e-6 * eps, crsMat.numRows () );
// 
//     // para1: temp vector
//     // para2: leading dimension of temp vector
//     // para3: basis vector to be preserved
//     // para4: constants set contConst container
//     // para5: blockcluster tree
//     // para6: cluster Tree
//     // para7: function to initialize vectors to be preserved
//     AHMED::contBasis<PrecondType> haar ( presVec.data(), crsMat.numRows (), basis.data(), constants, blTree->root (), blTree->root ()->rowCluster (), AHMED::initVec );
// 
//     if ( !workingCopy.genLUDecomp_noCopy ( *_L, *_U, param, preCoarse, &haar ) ) {
// #else
//     if ( !workingCopy.genLUDecomp_noCopy ( *_L, *_U, param, preCoarse ) ) {
// #endif
// //       throw Exception ( "LU decomposition not successful", __FILE__, __LINE__, __FUNCTION__ );
//         //TODO 
//         cout << "exception" << endl;
//     }
//   }
// 
// public:
//   //! \brief Constructor creating the block cluster tree and the LU decomposition.
//  // Eigen::SparseMatrix: 
//     //     m.valuePtr();      // Pointer to the values
//     //     m.innerIndextr();  // Pointer to the indices.
//     //     m.outerIndexPtr(); // Pointer to the beginning of each inner vector
//   AHMEDPGMRes ( const SparseMatrixType &mat, unsigned int bmin, RealType eta, unsigned int rankMax, RealType eps,
//                 unsigned int numSolverSteps, unsigned int solverRestart, RealType solverTol, bool useAlgorithmicClustering, bool preCoarse = false )
//   : AHMED::CRSMatrix<RealType> ( mat.rows(), mat.innerIndextr(), mat.outerIndexPtr(),
//       mat.valuePtr() ),
//       _perm ( mat.getNumRows () ), _bmin ( bmin ), _eta ( eta ),
//       _solverTol ( solverTol ), _numSolverSteps ( numSolverSteps ), _solverRestart ( solverRestart ),
//       _log ( NULL ), _logDeleteFlag ( false ) {
//     AHMED::CRSMatrix<RealType> crsMat ( mat.rows(), 
//                                         mat.innerIndextr(), 
//                                         mat.outerIndexPtr(),
//                                         mat.valuePtr() );
// 
//     AHMED::BlockClusterTree_sptr blTree;
// //     if ( useAlgorithmicClustering ) {
// // #ifdef USE_EXTERNAL_METIS
//       blTree = AHMED::genBlockClusterTree_alg ( crsMat.noRows, crsMat.iA, crsMat.jA, _bmin, _eta, _perm );
// // #else
// //      //throw Exception ( "External Metis is required for algorithmic clustering", __FILE__, __LINE__, __FUNCTION__ );
// //       cout << "External Metis is required for algorithmic clustering" << endl;
// // #endif
// //     }
// //     else {
// //       AHMED_SPTR<QuocMeshClustering<ConfiguratorType> > clustering ( new QuocMeshClustering<ConfiguratorType> ( _mesh ) );
// //       auto clTree = AHMED::genClusterTree2d_bbx ( clustering, mat.getNumRows(), _bmin, _perm );
// //       blTree = AHMED::genBlockClusterTree ( clTree, clTree, _eta );
// //     }
// 
//     createLUDecomposition ( rankMax, eps, crsMat, blTree, preCoarse );
//   }
// 
//   void setLogging ( std::ofstream &log, bool deleteFlag, std::string logPrefix = "" ) {
//     if ( _logDeleteFlag ) delete _log;
//     _log = &log;
//     _logDeleteFlag = deleteFlag;
//     _logPrefix = logPrefix;
//   }
// 
//   void solve ( const VectorType &arg, VectorType &dest ) const {
//     auto info = GMRes ( *this, _solverTol, _numSolverSteps, _solverRestart, arg.data(), dest.data() );
// 
//     std::stringstream logMessage;
//     logMessage << "steps " << info.second << " res " << info.first;
//     log ( logMessage.str ().c_str () );
// 
//     if ( info.first > _solverTol ) std::cerr << "Solver iteration did not converge." << std::endl;
//   }
// 
//   void precond_apply ( RealType *x ) const override {
//     std::vector<RealType> tmp ( x, x + this->numRows () );
// 
//     AHMED::permuteVec ( _perm, x, tmp.data() );
//     pApply ( tmp.data(), std::is_same<RealType, PrecondType> () );
//     depermuteVec ( _perm, tmp.data(), x );
//   }
// 
// protected:
//   AHMED::Perm _perm;
//   unsigned int _bmin;
//   RealType _eta;
//   RealType _solverTol;
//   unsigned int _numSolverSteps;
//   unsigned int _solverRestart;
// 
//   AHMED::LtHmat<PrecondType> *_L; //!< lower triangular H-matrix of LU decomposition
//   AHMED::UtHmat<PrecondType> *_U; //!< upper triangular H-matrix of LU decomposition
// 
//   std::ofstream *_log;                         //!< File for logging of solver data.
//   bool _logDeleteFlag;                         //!< Delete flag for _log.
//   std::string _logPrefix;                      //!< Prefix for each log line.
// 
//   /// wrapper function for precond_apply with value_type == precond_type
//   void pApply ( RealType *x, const std::true_type ) const {  AHMED::HLU_solve ( *_L, *_U, x );}
// 
//   /// wrapper function for precond_apply with value_type != precond_type
//   void pApply ( RealType *x, const std::false_type ) const {
//     std::vector<PrecondType> tmp ( this->numRows () );
//     tmp.assign ( x, x + this->numRows () );
//     HLU_solve ( *_L, *_U, tmp.data () );
//     std::copy_n ( tmp.data (), this->numRows (), x );
//   }
// 
//   inline void log ( const char* s ) const {
//     if ( _log != NULL ) {
//       *( _log ) << _logPrefix << ' ' << s << std::endl;
//     }
//   }
// 
// private:
//   AHMEDPGMRes ()
// : _perm ( 0 ), _bmin ( 0 ), _eta ( 0.0 ), _solverTol ( 0.0 ), _numSolverSteps ( 0 ), _solverRestart ( 0 ),
//   _log ( NULL ), _logDeleteFlag ( false ) {}
// };
// 
// 
// 
// //! \brief Linear solver using the preconditioned GMRES method of the AHMED library
// //! \note The MultiVector version uses the vector version of this class.
// //! \note AHMED CRS matrices are quadratic compressed row format matrices.
// // template<typename ConfiguratorType, typename PrecondType>
// // class AHMEDPGMRes_MultiVector<ConfiguratorType, PrecondType> {
// //   typedef typename ConfiguratorType::RealType RealType;
// // 
// // protected:
// //   AHMEDPGMRes<ConfiguratorType, PrecondType> _vecOp;
// // public:
// //   AHMEDPGMRes_MultiVector ( const RPOpType &rp,
// //       const CSRMatrix<RealType, unsigned int> &mat, unsigned int bmin, RealType eta, unsigned int rankMax, RealType eps,
// //       unsigned int numSolverSteps, unsigned int solverRestart, RealType solverTol, bool useAlgorithmicClustering, bool preCoarse = false )
// // : _vecOp ( rp, mat, bmin, eta, rankMax, eps, numSolverSteps, solverRestart, solverTol, useAlgorithmicClustering, preCoarse ) {}
// // 
// //   void setLogging ( std::ofstream &log, bool deleteFlag, std::string logPrefix = "" ) {
// //     _vecOp.setLogging ( log, deleteFlag, logPrefix );
// //   }
// // 
// //   void solve ( const std::vector<VectorType> &arg, std::vector<VectorType> &dest ) const {
// //       
// //     VectorType rhs ( arg.getTotalSize () );
// //     VectorType sol ( dest.getTotalSize () );
// // 
// //     rhs.copyUnblockedFrom ( arg );
// //     _vecOp.apply ( rhs, sol );
// //     dest.copySplitFrom ( sol );
// //   }
// // };
// 
// #endif
// 
// 
// // #else
// // template<typename ConfiguratorType>
// // class QuocMeshClustering {
// // public:
// //   template<typename... Types>
// //   QuocMeshClustering ( Types... ) {
// //     throw Exception ( "External AHMED required", __FILE__, __LINE__, __FUNCTION__ );
// //   }
// // };
// // 
// // template<typename ConfiguratorType, typename VectorType, typename PrecondType>
// // class AHMEDPGMRes : public InverseOp<VectorType> {
// //   typedef typename ConfiguratorType::RealType RealType;
// // 
// // public:
// //   template<typename RPOpType>
// //   AHMEDPGMRes ( const RPOpType &,
// //       const CSRMatrix<RealType, unsigned int> &, unsigned int , RealType , unsigned int , RealType ,
// //       unsigned int , unsigned int , RealType , bool , bool preCoarse = false ) {
// //     throw Exception ( "External AHMED required", __FILE__, __LINE__, __FUNCTION__ );
// //   }
// // 
// //   virtual void applyAdd ( const VectorType &, VectorType & ) const {
// //     throw Exception ( "External AHMED required", __FILE__, __LINE__, __FUNCTION__ );
// //   }
// // 
// //   void setLogging ( std::ofstream &, bool , std::string logPrefix = "" ) {
// //     throw Exception ( "External AHMED required", __FILE__, __LINE__, __FUNCTION__ );
// //   }
// // };
// // 
// // }
// 
// 
// 
// 
// //      CSRMatrix<RealType, unsigned int> csrMat;
// //      d2lcs.apply ( arg, csrMat );
// //
// //      const unsigned int bmin = 25;
// //      const RealType eta = 0.8;
// //      const unsigned int rankMax = 1000;
// //      const RealType eps = 1e-2; const RealType eps = 1e-4;
// //      AHMEDPGMRes<RealType, Vector<RealType>, RealType> ahmedpgmres ( csrMat, bmin, eta, rankMax, eps, 1e5, 50, 1e-16 );
// //      AHMEDPGMRes<ConfiguratorType, MultiVector<RealType>, float> ( _rp, *( this->_pMatDF ), bmin, eta, rankMax, eps, 1e3, 50, 1e-16, _useAlgorithmicClustering );
// //      Vector<RealType> rhs1 ( csrMat.getNumCols () );
// //      rhs1.setAll ( 1.0 );
// //      Vector<RealType> sol1 ( rhs1, STRUCT_COPY );
// //      ahmedpgmres.apply ( rhs1, sol1 );
// //
// //
// //      Vector<RealType> check ( sol1, STRUCT_COPY );
// //      csrMat.apply ( sol1, check );
// //      for ( int i = 0; i < check.size (); ++i ) {
// //        if ( check[i] - rhs1[i] > 1e-6 )
// //          cout << "Difference at " << i << ": " << check [i] << " != " << rhs1[i] << endl;
// //      }
// //
// //      exit ( 0 );
// 
// 
//     
// 
//     
// 
// #endif /* AHMEDPGMRES_H_ */
