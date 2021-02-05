#ifndef __QUOCFEQUADRATURE_H
#define __QUOCFEQUADRATURE_H


// template <typename RealType, int _NumQuadPoints>
// class FEQuadrature {
// public:
//   static const int numQuadPoints = _NumQuadPoints;
// 
// protected:
//   RealVecChart _points [ numQuadPoints ];
//   RealType   _weights[ numQuadPoints ];
//   
// public:
//   FEQuadrature() {}
// 
//   inline const RealVecChart &getRefCoord ( int QuadPoint ) const { return _points[ QuadPoint ];}
//   inline RealType getWeight ( int QuadPoint ) const { return _weights[ QuadPoint ];}
//    
// };
// 
// 
// 
// template <typename RealType>
// class SimpsonQuadrature1D : public FEQuadrature<RealType,3> {
// public:
//   SimpsonQuadrature1D( ) {
//     this->_points[0][0] = 0.0;
//     this->_points[1][0] = 0.5;
//     this->_points[2][0] = 1.0;
//     this->_weights[0] = 1./6.; this->_weights[1] = 2./3.; this->_weights[2] = 1./6.
//   }
// };
// 
// 
// 
// 
// template <typename RealType, typename Quadrature1DType>
// class TensorProductQuadrature2D 
//   : public FEQuadrature<RealType,Quadrature1DType::numQuadPoints * Quadrature1DType::numQuadPoints> {
//       
// public:
//   typedef RealType RealType;
// 
//   static const int numQuadPoints = Quadrature1DType::numQuadPoints * Quadrature1DType::numQuadPoints;
// 
//   TensorProductQuadrature2D( ) {
// 
//     Quadrature1DType quad1D;
// 
//     for ( int k = 0, i = 0; i < Quadrature1DType::numQuadPoints; i++ ) {
//       for ( int j = 0; j < Quadrature1DType::numQuadPoints; j++, k++ ) {
//         this->_points[k][0] = quad1D.getRefCoord ( i )[0];
//         this->_points[k][1] = quad1D.getRefCoord ( j )[0];
//         this->_weights[k] = quad1D.getWeight ( i ) * quad1D.getWeight ( j );
//       }
//     }
//   }
// };
// 
// template <typename RealType, typename Quadrature1DType>
// class TensorProductQuadrature3D 
//   : public FEQuadrature<RealType,Quadrature1DType::numQuadPoints * Quadrature1DType::numQuadPoints * Quadrature1DType::numQuadPoints> {
//       
// public:
//   typedef RealType RealType;
// 
//   static const int numQuadPoints = Quadrature1DType::numQuadPoints * Quadrature1DType::numQuadPoints * Quadrature1DType::numQuadPoints;
// 
//   TensorProductQuadrature3D( ) {
// 
//     Quadrature1DType quad1D;
// 
//     for ( int k = 0, i = 0; i < Quadrature1DType::numQuadPoints; i++ ) {
//       for ( int j = 0; j < Quadrature1DType::numQuadPoints; j++ ) {
//         for ( int l = 0; l < Quadrature1DType::numQuadPoints; l++, k++ ) {
//           this->_points[k][0] = quad1D.getRefCoord ( i )[0];
//           this->_points[k][1] = quad1D.getRefCoord ( j )[0];
//           this->_points[k][2] = quad1D.getRefCoord ( l )[0];
//           this->_weights[k] = quad1D.getWeight ( i ) * quad1D.getWeight ( j ) * quad1D.getWeight ( l );
//         }
//       }
//     }
//   }
// };



//! ==================================================================================================================
//! ==================================================================================================================
//! ==================================================================================================================
//                           Simpson quadrature rules
//! ==================================================================================================================
//! ==================================================================================================================
//!  =================================================================================================================

template <typename RealType, typename RealVecChart>
class SimpsonQuadrature1D {
public:
  static const int numQuadPoints = 3;
  
protected:
  RealVecChart _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
  
public:
  SimpsonQuadrature1D(){
        _points[0] = RealVecChart(0.0);
        _points[1] = RealVecChart(0.5);
        _points[2] = RealVecChart(1.0);
        _weights[0] = 1./6.; _weights[1] = 2./3.; _weights[2] = 1./6.;
  }
  
  inline const RealVecChart &getRefCoord ( int QuadPoint ) const { return _points[ QuadPoint ];}
  inline RealType getWeight ( int QuadPoint ) const { return _weights[ QuadPoint ];}
};


template <typename RealType, typename RealVecChart, typename RealVecChart1D>
class SimpsonQuadrature2D {
public:
  typedef SimpsonQuadrature1D<RealType, RealVecChart1D> QuadRuleType1D;
  static const int numQuadPoints1D = QuadRuleType1D::numQuadPoints;
  static const int numQuadPoints =  numQuadPoints1D * numQuadPoints1D;
    
protected:
  RealVecChart _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
    
public:
    
  SimpsonQuadrature2D( ) {

    QuadRuleType1D quad1D;

    for ( int k = 0, i = 0; i < numQuadPoints1D; i++ )
      for ( int j = 0; j < numQuadPoints1D; j++, k++ ) {
        _points[k][0] = quad1D.getRefCoord( i )[0];
        _points[k][1] = quad1D.getRefCoord ( j )[0];
        _weights[k] = quad1D.getWeight ( i ) * quad1D.getWeight ( j );
      }
  }

  inline const RealVecChart &getRefCoord ( int QuadPoint ) const { return _points[ QuadPoint ];}
  inline RealType getWeight ( int QuadPoint ) const { return _weights[ QuadPoint ];}
};





template <typename RealType, typename RealVecChart, typename RealVecChart1D>
class SimpsonQuadrature3D {
public:
  typedef SimpsonQuadrature1D<RealType, RealVecChart1D> QuadRuleType1D;  
  static const int numQuadPoints1D = QuadRuleType1D::numQuadPoints;
  static const int numQuadPoints =  numQuadPoints1D * numQuadPoints1D * numQuadPoints1D;

protected:
  RealVecChart _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
  
public:
  
  SimpsonQuadrature3D( ) {

    QuadRuleType1D quad1D;

    for ( int k = 0, i = 0; i < numQuadPoints1D; i++ ) {
      for ( int j = 0; j < numQuadPoints1D; j++ ) {
        for ( int l = 0; l < numQuadPoints1D; l++, k++ ) {
          _points[k][0] = quad1D.getRefCoord ( i )[0];
          _points[k][1] = quad1D.getRefCoord ( j )[0];
          _points[k][2] = quad1D.getRefCoord ( l )[0];
          _weights[k] = quad1D.getWeight ( i ) * quad1D.getWeight ( j ) * quad1D.getWeight ( l );
        }
      }
    }
  }

  inline const RealVecChart &getRefCoord ( int QuadPoint ) const { return _points[ QuadPoint ];}
  inline RealType getWeight ( int QuadPoint ) const { return _weights[ QuadPoint ];}
};












//! ==================================================================================================================
//! ==================================================================================================================
//! ==================================================================================================================
//                           Gauss quadrature rules
//! ==================================================================================================================
//! ==================================================================================================================
//!  =================================================================================================================


template <typename RealType, typename RealVecChart, int _Order>
class GaussQuadrature1D {};

template <typename RealType,  typename RealVecChart>
class GaussQuadrature1D<RealType, RealVecChart, 1> {
public:
  static const int numQuadPoints = 1;

protected:
  RealVecChart _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
  
public:
  
  GaussQuadrature1D( ) {
      _points[0] = RealVecChart(0.5);
      _weights[0] = 1.;
  }

  inline const RealVecChart& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _weights[QuadPoint];
  }
};


template <typename RealType,  typename RealVecChart>
class GaussQuadrature1D<RealType, RealVecChart, 3> {
public:
  static const int numQuadPoints = 2;

protected:
  RealVecChart _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
  
public:
  
  GaussQuadrature1D( ) {
    _points[0] = RealVecChart(0.2113248654051871);
    _points[1] = RealVecChart(0.7886751345948129);
    _weights[0] = 0.5; _weights[1] = 0.5;
  }

  inline const RealVecChart& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _weights[QuadPoint];
  }
};


template <typename RealType,  typename RealVecChart>
class GaussQuadrature1D<RealType, RealVecChart, 5> {
public:
  static const int numQuadPoints = 3;

protected:
  RealVecChart _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
  
public:
  
  GaussQuadrature1D( ) {
    _points[0] = RealVecChart(0.1127016653792583);
    _points[1] = RealVecChart(0.5);
    _points[2] = RealVecChart(0.8872983346207417);
    _weights[0] = 0.2777777777777778; _weights[1] = 0.4444444444444444; _weights[2] = 0.2777777777777778;
  }

  inline const RealVecChart& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _weights[QuadPoint];
  }
};


template <typename RealType,  typename RealVecChart>
class GaussQuadrature1D<RealType, RealVecChart, 7> {
public:
  static const int numQuadPoints = 4;

protected:
  RealVecChart _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
  
public:
  
  GaussQuadrature1D( ) {
    _points[0] = RealVecChart(0.6943184420297371e-1);
    _points[1] = RealVecChart(0.3300094782075719);
    _points[2] = RealVecChart(0.6699905217924281);
    _points[3] = RealVecChart(0.9305681557970263);
    
    _weights[0] = 0.1739274225687269; 
    _weights[1] = 0.3260725774312731;
    _weights[2] = 0.3260725774312731;
    _weights[3] = 0.1739274225687269;
  }

  inline const RealVecChart& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _weights[QuadPoint];
  }
};


template <typename RealType,  typename RealVecChart>
class GaussQuadrature1D<RealType, RealVecChart, 9> {
public:
  static const int numQuadPoints = 5;

protected:
  RealVecChart _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
  
public:
  
  GaussQuadrature1D( ) {
    
    _points[0] = RealVecChart(0.4691007703066800e-1);
    _points[1] = RealVecChart(0.2307653449471585);
    _points[2] = RealVecChart(0.5);
    _points[3] = RealVecChart(0.7692346550528415);
    _points[4] = RealVecChart(0.9530899229693320);
    
    _weights[0] = 0.1184634425280945;
    _weights[1] = 0.2393143352496832;
    _weights[2] = 0.2844444444444444;
    _weights[3] = 0.2393143352496832;
    _weights[4] = 0.1184634425280945;
  }

  inline const RealVecChart& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _weights[QuadPoint];
  }
};




template <typename RealType,  typename RealVecChart>
class GaussQuadrature1D<RealType, RealVecChart, 11> {
public:
  static const int numQuadPoints = 6;

protected:
  RealVecChart _points [ numQuadPoints ];
//   RealType   _weights[ numQuadPoints ];
  
public:
  
  GaussQuadrature1D( ) {
    _points[0] = RealVecChart(0.3376524289842399e-1);
    _points[1] = RealVecChart(0.1693953067668677);
    _points[2] = RealVecChart(0.3806904069584015);
    _points[3] = RealVecChart(0.6193095930415985);
    _points[4] = RealVecChart(0.8306046932331323);
    _points[5] = RealVecChart(0.9662347571015760);
  }

  inline const RealVecChart& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    static const RealType _weights [numQuadPoints] = {
                                                       .8566224618958517e-1, .1803807865240693, .2339569672863455,
                                                       .2339569672863455, .1803807865240693, .8566224618958517e-1
                                                     };
    return _weights [QuadPoint];
  }
};


template <typename RealType,  typename RealVecChart>
class GaussQuadrature1D<RealType, RealVecChart, 13> {
public:
  static const int numQuadPoints = 7;

protected:
  RealVecChart _points [ numQuadPoints ];
//   RealType   _weights[ numQuadPoints ];
  
public:
  
  GaussQuadrature1D( ) {
    _points[0] = RealVecChart(0.2544604382862074e-1);
    _points[1] = RealVecChart(0.1292344072003028);
    _points[2] = RealVecChart(0.2970774243113014);
    _points[3] = RealVecChart(0.5);
    _points[4] = RealVecChart(0.7029225756886986);
    _points[5] = RealVecChart(0.8707655927996972);
    _points[6] = RealVecChart(0.9745539561713793);
  }

  inline const RealVecChart& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    static const RealType _weights [numQuadPoints] = {
                                                       .6474248308443485e-1, .1398526957446383, .1909150252525595, .2089795918367347,
                                                       .1909150252525595, .1398526957446383, .6474248308443485e-1
                                                     };
    return _weights [QuadPoint];
  }
};



template <typename RealType,  typename RealVecChart>
class GaussQuadrature1D<RealType, RealVecChart, 15> {
public:
  static const int numQuadPoints = 8;

protected:
  RealVecChart _points [ numQuadPoints ];
//   RealType   _weights[ numQuadPoints ];
  
public:
  
  GaussQuadrature1D( ) {
    _points[0] = RealVecChart(0.1985507175123188e-1);
    _points[1] = RealVecChart(0.1016667612931866);
    _points[2] = RealVecChart(0.2372337950418355);
    _points[3] = RealVecChart(0.4082826787521751);
    _points[4] = RealVecChart(0.5917173212478249);
    _points[5] = RealVecChart(0.7627662049581645);
    _points[6] = RealVecChart(0.8983332387068134);
    _points[7] = RealVecChart(0.9801449282487681);
  }

  inline const RealVecChart& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    static const RealType _weights [numQuadPoints] = {
                                                       .5061426814518813e-1, .1111905172266872, .1568533229389436, .1813418916891810,
                                                       .1813418916891810, .1568533229389436, .1111905172266872, .5061426814518813e-1
                                                     };
    return _weights [QuadPoint];
  }
};


template <typename RealType,  typename RealVecChart>
class GaussQuadrature1D<RealType, RealVecChart, 17> {
public:
  static const int numQuadPoints = 9;

protected:
  RealVecChart _points [ numQuadPoints ];
//   RealType   _weights[ numQuadPoints ];
  
public:
  
  GaussQuadrature1D( ) {
    _points[0] = RealVecChart(0.1591988024618696e-1);
    _points[1] = RealVecChart(0.8198444633668210e-1);
    _points[2] = RealVecChart(0.1933142836497048);
    _points[3] = RealVecChart(0.3378732882980955);
    _points[4] = RealVecChart(0.5);
    _points[5] = RealVecChart(0.6621267117019045);
    _points[6] = RealVecChart(0.8066857163502952);
    _points[7] = RealVecChart(0.9180155536633179);
    _points[8] = RealVecChart(0.9840801197538130);
  }

  inline const RealVecChart& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    static const RealType _weights [numQuadPoints] = {
                                                       .4063719418078721e-1, .9032408034742870e-1, .1303053482014677,
                                                       .1561735385200014, .1651196775006299, .1561735385200014,
                                                       .1303053482014677, .9032408034742870e-1, .4063719418078721e-1
                                                     };
    return _weights [QuadPoint];
  }
};



template <typename RealType,  typename RealVecChart>
class GaussQuadrature1D<RealType, RealVecChart, 19> {
public:
  static const int numQuadPoints = 10;

protected:
  RealVecChart _points [ numQuadPoints ];
//   RealType   _weights[ numQuadPoints ];
  
public:
  
  GaussQuadrature1D( ) {
    _points[0] = RealVecChart(0.1304673574141414e-1);
    _points[1] = RealVecChart(0.6746831665550774e-1);
    _points[2] = RealVecChart(0.1602952158504878);
    _points[3] = RealVecChart(0.2833023029353764);
    _points[4] = RealVecChart(0.4255628305091844);
    _points[5] = RealVecChart(0.5744371694908156);
    _points[6] = RealVecChart(0.7166976970646236);
    _points[7] = RealVecChart(0.8397047841495122);
    _points[8] = RealVecChart(0.9325316833444923);
    _points[9] = RealVecChart(0.9869532642585859);
  }

  inline const RealVecChart& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    static const RealType _weights [numQuadPoints] = {
                                                       .3333567215434407e-1, .7472567457529030e-1, .1095431812579910, .1346333596549982,
                                                       .1477621123573764, .1477621123573764, .1346333596549982,
                                                       .1095431812579910, .7472567457529030e-1, .3333567215434407e-1
                                                     };
    return _weights [QuadPoint];
  }
};



// template <typename RealType, int _Order>
// class GaussQuadrature {
// public:
//   typedef RealType RealType;
//   static const qc::Dimension Dim = qc::QC_1D;
//   static const int Order = _Order;
// 
//   enum { numQuadPoints =
//            GaussQuadrature<RealType, qc::QC_1D, _Order>::numQuadPoints *
//          GaussQuadrature<RealType, qc::QC_1D, _Order>::numQuadPoints };
// 
//   GaussQuadrature( ) {
// 
//     GaussQuadrature<RealType, qc::QC_1D, _Order> quad1D;
// 
//     for ( int k = 0, i = 0; i < GaussQuadrature<RealType, qc::QC_1D, _Order>::numQuadPoints; i++ ) {
//       for ( int j = 0; j < GaussQuadrature<RealType, qc::QC_1D, _Order>::numQuadPoints; j++, k++ ) {
//         _points[k][0] = quad1D.getRefCoord ( i )[0];
//         _weights[k] = quad1D.getWeight ( i ) * quad1D.getWeight ( j );
//       }
//     }
//   }
// 
//   inline const Vec<1, RealType> &getRefCoord ( int QuadPoint ) const {
//     return _points[ QuadPoint ];
//   }
// 
//   inline RealType getWeight ( int QuadPoint ) const {
//     return _weights[ QuadPoint ];
//   }
// 
// protected:
//   Vec<1, RealType> _points [ numQuadPoints ];
//   RealType         _weights[ numQuadPoints ];
// };


template <typename RealType, typename RealVecChart, typename RealVecChart1D, int _Order>
class GaussQuadrature2D {
public:
  typedef GaussQuadrature1D<RealType, RealVecChart1D, _Order> QuadRuleType1D;  
  static const int numQuadPoints1D = QuadRuleType1D::numQuadPoints;
  static const int numQuadPoints =  numQuadPoints1D * numQuadPoints1D;
    
protected:
  RealVecChart _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
    
public:
    
  GaussQuadrature2D( ) {

    QuadRuleType1D quad1D;

    for ( int k = 0, i = 0; i < numQuadPoints1D; i++ )
      for ( int j = 0; j < numQuadPoints1D; j++, k++ ) {
        _points[k][0] = quad1D.getRefCoord( i )[0];
        _points[k][1] = quad1D.getRefCoord ( j )[0];
        _weights[k] = quad1D.getWeight ( i ) * quad1D.getWeight ( j );
      }
  }

  inline const RealVecChart &getRefCoord ( int QuadPoint ) const { return _points[ QuadPoint ];}
  inline RealType getWeight ( int QuadPoint ) const { return _weights[ QuadPoint ];}
};





template <typename RealType, typename RealVecChart, typename RealVecChart1D, int _Order>
class GaussQuadrature3D {
public:
  typedef GaussQuadrature1D<RealType, RealVecChart1D, _Order> QuadRuleType1D;  
  static const int numQuadPoints1D = QuadRuleType1D::numQuadPoints;
  static const int numQuadPoints =  numQuadPoints1D * numQuadPoints1D * numQuadPoints1D;

protected:
  RealVecChart _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
  
public:
  
  GaussQuadrature3D( ) {

    QuadRuleType1D quad1D;

    for ( int k = 0, i = 0; i < numQuadPoints1D; i++ ) {
      for ( int j = 0; j < numQuadPoints1D; j++ ) {
        for ( int l = 0; l < numQuadPoints1D; l++, k++ ) {
          _points[k][0] = quad1D.getRefCoord ( i )[0];
          _points[k][1] = quad1D.getRefCoord ( j )[0];
          _points[k][2] = quad1D.getRefCoord ( l )[0];
          _weights[k] = quad1D.getWeight ( i ) * quad1D.getWeight ( j ) * quad1D.getWeight ( l );
        }
      }
    }
  }

  inline const RealVecChart &getRefCoord ( int QuadPoint ) const { return _points[ QuadPoint ];}
  inline RealType getWeight ( int QuadPoint ) const { return _weights[ QuadPoint ];}
};


#endif
