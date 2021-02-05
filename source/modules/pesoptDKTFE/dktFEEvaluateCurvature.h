#ifndef __DKTFEEVALUATECURVATURE_H
#define __DKTFEEVALUATECURVATURE_H

#include <pesopt_IO.h>
#include <dktFEHandler.h>

template< typename ConfiguratorType >
class CurvatureEvaluator {
    
protected:
    typedef typename ConfiguratorType::DTContainer DataTypeContainer;
    typedef typename ConfiguratorType::InitType MeshType;
    typedef typename DataTypeContainer::RealType RealType;
    typedef pesopt::BoostParser ParameterParserType;
    typedef typename DataTypeContainer::PointType PointType;
    typedef typename ConfiguratorType::ElementType ElementType;
    typedef typename DataTypeContainer::VectorType VectorType;
    typedef typename DataTypeContainer::Matrix22 Matrix22;
    typedef typename DataTypeContainer::Matrix33 Matrix33;
    typedef typename DataTypeContainer::DomVecType DomVecType;
    typedef typename DataTypeContainer::TangentVecType TangentVecType;
    
    const ShellHandler<ConfiguratorType,FirstAndSecondOrder> & _shellHandler;
    const string _saveDirectory;
    const ConfiguratorType &_conf;
    
    const DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> _xAFullCached;
    
public:
    CurvatureEvaluator ( const ShellHandler<ConfiguratorType,FirstAndSecondOrder> & shellHandler, const string saveDirectory = "" ) :
    _shellHandler ( shellHandler ), _saveDirectory ( saveDirectory ), _conf( shellHandler._conf ),
    _xAFullCached ( _conf, shellHandler.getChartToUndeformedShell(), 3 ) {}
    
    
    void plotMeshUndeformed( ) const {
        std::vector<PointType> pointVec; pointVec.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        for ( int elementIdx = 0; elementIdx < _conf.getInitializer().getNumTriangs(); ++elementIdx){
            const ElementType& El ( _conf.getInitializer().getTriang( elementIdx ) );
            for ( int localQuadPointIndex = 0; localQuadPointIndex < _conf.maxNumQuadPoints (); ++localQuadPointIndex){
                pointVec.push_back( _xAFullCached.getCoords(elementIdx,localQuadPointIndex) );
            }
        }
    }
    
    
    void plotCurvatureUndeformed( ) const {
        std::vector<PointType> pointVec; pointVec.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        std::vector<RealType> curvatureVec; curvatureVec.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        std::vector<RealType> WeingartenVec; WeingartenVec.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        for ( int elementIdx = 0; elementIdx < _conf.getInitializer().getNumTriangs(); ++elementIdx){
            const ElementType& El ( _conf.getInitializer().getTriang( elementIdx ) );
            for ( int localQuadPointIndex = 0; localQuadPointIndex < _conf.maxNumQuadPoints (); ++localQuadPointIndex){
                pointVec.push_back( _xAFullCached.getCoords(elementIdx,localQuadPointIndex) );
                curvatureVec.push_back ( (_xAFullCached.getSecondFF(elementIdx,localQuadPointIndex)).norm() );
                WeingartenVec.push_back ( (_xAFullCached.getFirstFFInv(elementIdx,localQuadPointIndex) *  _xAFullCached.getSecondFF(elementIdx,localQuadPointIndex)).norm() );
            }
        }
    }
    
    void plotCurvatureDeformed( const VectorType &disp ) const {
        
        VectorType deform ( disp + _shellHandler.getChartToUndeformedShell() );
        DiscreteVectorFunctionStorage<ConfiguratorType,FirstAndSecondOrder> xBFullCached ( _conf, deform, 3 );
        
        std::vector<PointType> pointVec; pointVec.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        std::vector<RealType> curvatureVec; curvatureVec.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        std::vector<RealType> WeingartenVec; WeingartenVec.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        std::vector<RealType> WeingartenVecWrtxA; WeingartenVecWrtxA.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        std::vector<TangentVecType> EigenVecsWeingarten1, EigenVecsWeingarten2, EigenVecsWeingarten3; 
        EigenVecsWeingarten1.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        EigenVecsWeingarten2.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        EigenVecsWeingarten3.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        std::vector<TangentVecType> EigenVecsHB1, EigenVecsHB2; 
        EigenVecsHB1.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        EigenVecsHB2.reserve( _conf.getInitializer().getNumTriangs() * _conf.maxNumQuadPoints () );
        for ( int elementIdx = 0; elementIdx < _conf.getInitializer().getNumTriangs(); ++elementIdx){
            const ElementType& El ( _conf.getInitializer().getTriang( elementIdx ) );
            for ( int localQuadPointIndex = 0; localQuadPointIndex < _conf.maxNumQuadPoints (); ++localQuadPointIndex){
                pointVec.push_back( xBFullCached.getCoords(elementIdx,localQuadPointIndex) );
                curvatureVec.push_back ( (xBFullCached .getSecondFF(elementIdx,localQuadPointIndex)).norm() );
                Matrix22 Weingarten = xBFullCached .getFirstFFInv(elementIdx,localQuadPointIndex) *  xBFullCached.getSecondFF(elementIdx,localQuadPointIndex);
                WeingartenVec.push_back ( (Weingarten).norm() );
                Matrix22 WeingartenWrtxA = _xAFullCached.getFirstFFInv(elementIdx,localQuadPointIndex) *  xBFullCached .getSecondFF(elementIdx,localQuadPointIndex);
                WeingartenVecWrtxA.push_back ( WeingartenWrtxA.norm() );

                Eigen::EigenSolver<Matrix22> es(xBFullCached .getSecondFF(elementIdx,localQuadPointIndex));
                cout << "eigenvalues of h_B = " << endl << es.eigenvalues().real() << endl;
                DomVecType eigenVecHB1 = es.eigenvectors().real().col(0), eigenVecHB2 = es.eigenvectors().real().col(1);
                EigenVecsHB1.push_back( xBFullCached.getGradient(elementIdx,localQuadPointIndex) * eigenVecHB1 );
                EigenVecsHB2.push_back( xBFullCached.getGradient(elementIdx,localQuadPointIndex) * eigenVecHB2 );
                
                Matrix33 Weingarten33 = xBFullCached.getGradient(elementIdx,localQuadPointIndex) * Weingarten * xBFullCached.getFirstFFInv(elementIdx,localQuadPointIndex) * xBFullCached.getGradient(elementIdx,localQuadPointIndex).transpose();
                Eigen::EigenSolver<Matrix33> es33(Weingarten33);
                //cout << "eigenvalues of Weingartenmap33 = " << endl << es33.eigenvalues().real() << endl;
                RealType eigenVal1 = es33.eigenvalues().real()(0), eigenVal2 = es33.eigenvalues().real()(1), eigenVal3 = es33.eigenvalues().real()(2);
                int eigenValIndex1, eigenValIndex2, eigenValIndex3;
                if( (eigenVal1 >= eigenVal2) && (eigenVal1 >= eigenVal3) ){
                    eigenValIndex1 = 0;
                    if( eigenVal2 >= eigenVal3 ){
                      eigenValIndex2 = 1; eigenValIndex3 = 2;   
                    }else{
                      eigenValIndex2 = 2; eigenValIndex3 = 1;  
                    }
                }else{
                    if( eigenVal2 >= eigenVal3 ){
                      eigenValIndex2 = 0;
                      if( eigenVal3 >= eigenVal1 ){
                          eigenValIndex3 = 1; eigenValIndex1 = 2; 
                      }else{
                          eigenValIndex3 = 2; eigenValIndex1 = 1; 
                      }  
                    }else{
                      eigenValIndex3 = 0;  
                       if( eigenVal2 >= eigenVal1 ){
                          eigenValIndex2 = 1; eigenValIndex1 = 2; 
                      }else{
                          eigenValIndex2 = 2; eigenValIndex1 = 1; 
                      } 
                    }
                }
                //TangentVecType eigenVec1 = es33.eigenvectors().real().col(0), eigenVec2 = es33.eigenvectors().real().col(1), eigenVec3 = es33.eigenvectors().real().col(2);
                TangentVecType eigenVec1 = es33.eigenvectors().real().col(eigenValIndex1), eigenVec2 = es33.eigenvectors().real().col(eigenValIndex2), eigenVec3 = es33.eigenvectors().real().col(eigenValIndex3);
                EigenVecsWeingarten1.push_back( eigenVec1 );
                EigenVecsWeingarten2.push_back( eigenVec2 );
                EigenVecsWeingarten3.push_back( eigenVec3 );
                
            }
        }
    
    }

};



#endif //__EVALUATECURVATURE_H
