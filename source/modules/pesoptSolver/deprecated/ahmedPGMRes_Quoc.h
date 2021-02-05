// /*
//  * ahmedPGMRes.h
//  *
//  *  Created on: May 14, 2014
//  *      Author: Toelkes
//  *
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
// #include <adaptiveElement.h>
// #include <adaptiveQuocRPOps.h>
// #include <multiVector.h>
// #include <pointerClasses.h>
// #include <qmException.h>
// #include <solver.h>
// #include <sparseMatrices.h>
// #include <vec.h>
// 
// #ifdef USE_EXTERNAL_AHMED
// #include <ahmedIncludes.h>
// 
// namespace pesopt {
// 
// //! \brief Helper class for geometric clustering
// //! \author Toelkes
// template<typename ConfiguratorType>
// class AdaptiveQuocClustering : public AHMED::GeoInfoForClustering {
//   typedef typename ConfiguratorType::RealType RealType;
//   typedef typename ConfiguratorType::InitType GridType;
//   typedef typename ConfiguratorType::VecType  VecType;
//   typedef qc::AdaptiveQuocRPOp<ConfiguratorType> RPOpType;
// 
// protected:
//   mutable std::unordered_map<unsigned int, qc::AdaptiveHashKey<>> _reverseMap;
// 
// public:
//   AdaptiveQuocClustering ( const RPOpType &rp ) {
//     for ( typename RPOpType::DofIterator it = rp.begin (); it != rp.end (); ++it ) {
//       _reverseMap[rp.getRestrictedGlobalNodeIndex ( it->_full )] = *it;
//     }
//   }
// 
//   double readCoordinate ( const unsigned i, const unsigned j ) const {
//     VecType x;
//     qc::AdaptiveHashKeyConversion<ConfiguratorType::DimOfWorld>::toCartesianCoords ( _reverseMap[i], x );
//     return x[j];
//   }
// 
//   double squaredRadius ( const unsigned i ) const {
//     return pesopt::Sqr ( 1 << _reverseMap[i]._block[0] );
//   }
// };
// 
// template<typename ConfiguratorType, typename VectorType, typename PrecondType>
// class AHMEDPGMRes {};
// 
// //! \brief Linear solver using the preconditioned GMRES method of the AHMED library
// //! \author Toelkes
// //! \note AHMED CRS matrices are quadratic compressed row format matrices.
// template<typename ConfiguratorType, typename PrecondType>
// class AHMEDPGMRes<ConfiguratorType, pesopt::Vector<typename ConfiguratorType::RealType>, PrecondType > : public pesopt::InverseOp<pesopt::Vector<typename ConfiguratorType::RealType> >,
// public AHMED::CRSMatrix<typename ConfiguratorType::RealType> {
//   typedef pesopt::DeleteFlagPointer<AHMED::LtHmat<PrecondType> > LTHMATPOINTER;
//   typedef pesopt::DeleteFlagPointer<AHMED::UtHmat<PrecondType> > UTHMATPOINTER;
// 
//   typedef typename ConfiguratorType::RealType RealType;
//   typedef typename ConfiguratorType::InitType GridType;
//   typedef qc::AdaptiveQuocRPOp<ConfiguratorType> RPOpType;
// 
//   //! \brief Creates the LU decomposition of a matrix using a block cluster tree.
//   void createLUDecomposition ( unsigned int rankMax, RealType eps, AHMED::CRSMatrix<RealType> crsMat,
//       AHMED::BlockClusterTree_sptr &blTree, bool preCoarse ) {
// //    AHMED::displayInfo ( *blTree );
// 
//     const AHMED::SvdParam param ( eps, rankMax, AHMED::relative, AHMED::Euclidean );
// 
//     AHMED::CRS_perm ( crsMat, _perm );
//     // AHMED::CRS_perm ( *this, _perm );
// 
//     const AHMED::SvdParam paramCRS ( param.eps * 1E-4, param.maxRank,
//         param.typeOfAcc, param.typeOfNrm );
// 
//     AHMED::GeHmat<PrecondType> workingCopy ( blTree );
// 
//     workingCopy.add ( crsMat, paramCRS );
//     // workingCopy.add ( *this, paramCRS );
// 
//     _L.reset ( new AHMED::LtHmat<PrecondType> ( blTree ), true );
//     _U.reset ( new AHMED::UtHmat<PrecondType> ( blTree ), true );
// 
// #define AHMED_PRESERVATION
// #ifdef AHMED_PRESERVATION
//     // create haar basis
//     pesopt::Vector<PrecondType> presVec ( 2 * crsMat.numRows () );
// 
//     // generate basis
//     pesopt::Vector<PrecondType> basis ( crsMat.numRows () );
//     basis.setAll ( 1.0 );
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
//     AHMED::contBasis<PrecondType> haar ( presVec.getData (), crsMat.numRows (), basis.getData (), constants,
//         blTree->root (), blTree->root ()->rowCluster (), AHMED::initVec );
// 
//     if ( !workingCopy.genLUDecomp_noCopy ( *_L, *_U, param, preCoarse, &haar ) ) {
// #else
//     if ( !workingCopy.genLUDecomp_noCopy ( *_L, *_U, param, preCoarse ) ) {
// #endif
//       throw pesopt::Exception ( "LU decomposition not successful", __FILE__, __LINE__, __FUNCTION__ );
//     }
//   }
// 
// public:
//   //! \brief Constructor creating the block cluster tree and the LU decomposition.
//   AHMEDPGMRes ( const RPOpType &rp,
//       const pesopt::CSRMatrix<RealType, unsigned int> &mat, unsigned int bmin, RealType eta, unsigned int rankMax, RealType eps,
//       unsigned int numSolverSteps, unsigned int solverRestart, RealType solverTol, bool useAlgorithmicClustering, bool preCoarse = false )
//   : AHMED::CRSMatrix<RealType> ( mat.getNumRows (), mat.getRowPointerReference ().getData (), mat.getColumnIndexReference ().getData (),
//       mat.getValueReference ().getData () ),
//       _perm ( mat.getNumRows () ), _bmin ( bmin ), _eta ( eta ),
//       _solverTol ( solverTol ), _numSolverSteps ( numSolverSteps ), _solverRestart ( solverRestart ),
//       _log ( NULL ), _logDeleteFlag ( false ) {
//     AHMED::CRSMatrix<RealType> crsMat ( mat.getNumRows (), mat.getRowPointerReference ().getData (), mat.getColumnIndexReference ().getData (),
//         mat.getValueReference ().getData () );
// 
//     AHMED::BlockClusterTree_sptr blTree;
//     if ( useAlgorithmicClustering ) {
// #ifdef USE_EXTERNAL_METIS
//       blTree = AHMED::genBlockClusterTree_alg ( crsMat.noRows, crsMat.iA, crsMat.jA, _bmin, _eta, _perm );
// #else
//       throw pesopt::Exception ( "External Metis is required for algorithmic clustering", __FILE__, __LINE__, __FUNCTION__ );
// #endif
//     }
//     else {
//       AHMED_SPTR<AdaptiveQuocClustering<ConfiguratorType> > clustering ( new AdaptiveQuocClustering<ConfiguratorType> ( rp ) );
//       auto clTree = AHMED::genClusterTree2d_bbx ( clustering, mat.getNumRows(), _bmin, _perm );
//       blTree = AHMED::genBlockClusterTree ( clTree, clTree, _eta );
//     }
// 
//     createLUDecomposition ( rankMax, eps, crsMat, blTree, preCoarse );
//   }
// 
//   void setLogging ( std::ofstream &log, bool deleteFlag, std::string logPrefix = "" ) {
//     if ( _logDeleteFlag )
//       delete _log;
// 
//     _log = &log;
//     _logDeleteFlag = deleteFlag;
//     _logPrefix = logPrefix;
//   }
// 
//   void applyAdd ( const pesopt::Vector<RealType> &, pesopt::Vector<RealType> & ) const {
//     throw pesopt::UnimplementedCodeException ( "Only apply is implemented", __FILE__, __LINE__ );
//   }
// 
//   void apply ( const pesopt::Vector<RealType> &arg, pesopt::Vector<RealType> &dest ) const {
//     auto info = GMRes ( *this, _solverTol, _numSolverSteps, _solverRestart, arg.getData (), dest.getData () );
// //    auto info = CG ( *this, _solverTol, _numSolverSteps, arg.getData (), dest.getData () );
// 
//     std::stringstream logMessage;
//     logMessage << "steps " << info.second << " res " << info.first;
//     log ( logMessage.str ().c_str () );
// 
//     if ( info.first > _solverTol )
//       std::cerr << "Solver iteration did not converge." << std::endl;
//   }
// 
//   void precond_apply ( RealType *x ) const override {
//     std::vector<RealType> tmp ( x, x + this->numRows () );
// 
//     AHMED::permuteVec ( _perm, x, tmp.data () );
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
//   LTHMATPOINTER _L;   //!< lower triangular H-matrix of LU decomposition
//   UTHMATPOINTER _U;   //!< upper triangular H-matrix of LU decomposition
// 
//   std::ofstream *_log;                         //!< File for logging of solver data.
//   bool _logDeleteFlag;                         //!< Delete flag for _log.
//   std::string _logPrefix;                      //!< Prefix for each log line.
// 
//   /// wrapper function for precond_apply with value_type == precond_type
//   void pApply ( RealType *x, const std::true_type ) const {
//     AHMED::HLU_solve ( *_L, *_U, x );
//   }
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
// //! \brief Linear solver using the preconditioned GMRES method of the AHMED library
// //! \author Toelkes
// //! \note The MultiVector version uses the vector version of this class.
// //! \note AHMED CRS matrices are quadratic compressed row format matrices.
// template<typename ConfiguratorType, typename PrecondType>
// class AHMEDPGMRes<ConfiguratorType, pesopt::MultiVector<typename ConfiguratorType::RealType>, PrecondType> : public pesopt::InverseOp<pesopt::MultiVector<typename ConfiguratorType::RealType> > {
//   typedef typename ConfiguratorType::RealType RealType;
//   typedef typename ConfiguratorType::InitType GridType;
//   typedef qc::AdaptiveQuocRPOp<ConfiguratorType> RPOpType;
// 
// protected:
//   pesopt::AHMEDPGMRes<ConfiguratorType, pesopt::Vector<RealType>, PrecondType> _vecOp;
// public:
//   AHMEDPGMRes ( const RPOpType &rp,
//       const pesopt::CSRMatrix<RealType, unsigned int> &mat, unsigned int bmin, RealType eta, unsigned int rankMax, RealType eps,
//       unsigned int numSolverSteps, unsigned int solverRestart, RealType solverTol, bool useAlgorithmicClustering, bool preCoarse = false )
// : _vecOp ( rp, mat, bmin, eta, rankMax, eps, numSolverSteps, solverRestart, solverTol, useAlgorithmicClustering, preCoarse ) {}
// 
//   void setLogging ( std::ofstream &log, bool deleteFlag, std::string logPrefix = "" ) {
//     _vecOp.setLogging ( log, deleteFlag, logPrefix );
//   }
// 
//   void solve ( const pesopt::MultiVector<RealType> &arg, pesopt::MultiVector<RealType> &dest ) const {
//     pesopt::Vector<RealType> rhs ( arg.getTotalSize () );
//     pesopt::Vector<RealType> sol ( dest.getTotalSize () );
// 
//     rhs.copyUnblockedFrom ( arg );
//     _vecOp.apply ( rhs, sol );
//     dest.copySplitFrom ( sol );
//   }
// };
// 
// }
// 
// #else
// namespace pesopt {
// template<typename ConfiguratorType>
// class AdaptiveQuocClustering {
// public:
//   template<typename... Types>
//   AdaptiveQuocClustering ( Types... ) {
//     throw pesopt::Exception ( "External AHMED required", __FILE__, __LINE__, __FUNCTION__ );
//   }
// };
// 
// template<typename ConfiguratorType, typename VectorType, typename PrecondType>
// class AHMEDPGMRes : public pesopt::InverseOp<VectorType> {
//   typedef typename ConfiguratorType::RealType RealType;
// 
// public:
//   template<typename RPOpType>
//   AHMEDPGMRes ( const RPOpType &,
//       const pesopt::CSRMatrix<RealType, unsigned int> &, unsigned int , RealType , unsigned int , RealType ,
//       unsigned int , unsigned int , RealType , bool , bool preCoarse = false ) {
//     throw pesopt::Exception ( "External AHMED required", __FILE__, __LINE__, __FUNCTION__ );
//   }
// 
//   virtual void applyAdd ( const VectorType &, VectorType & ) const {
//     throw pesopt::Exception ( "External AHMED required", __FILE__, __LINE__, __FUNCTION__ );
//   }
// 
//   void setLogging ( std::ofstream &, bool , std::string logPrefix = "" ) {
//     throw pesopt::Exception ( "External AHMED required", __FILE__, __LINE__, __FUNCTION__ );
//   }
// };
// 
// }
// #endif
// 
// 
// 
// #endif /* AHMEDPGMRES_H_ */
