#ifndef __TRIANGLEQUADRATURE_H
#define __TRIANGLEQUADRATURE_H
    
    
template <typename RealType, typename RealVecChart>
class CenterQuadrature {
public:
  static const int numQuadPoints = 1;

  inline const RealVecChart& getRefCoord ( int /*quadpoint*/ ) const  {
    static const RealVecChart c ( 1. / 3., 1. / 3. );
    return c;
  }

  inline RealType getWeight ( int /*quadpoint*/ ) const { return 1.;}
};

template <typename RealType, typename RealVecChart>
class TriQuadrature {
public:
   static const int numQuadPoints = 3;

  inline const RealVecChart& getRefCoord ( int quadpoint ) const  {
    static const RealVecChart c[3] = { RealVecChart ( 1. / 6., 1. / 6. ),
                                     RealVecChart ( 1. / 6., 2. / 3. ),
                                     RealVecChart ( 2. / 3., 1. / 6. )
                                    };
    return c[quadpoint];
  }

  inline RealType getWeight ( int /*quadpoint*/ ) const { return 1. / 3.;}
};


template <typename RealType, typename RealVecChart>
class EdgeQuadrature {
public:
   static const int numQuadPoints = 3;

  inline const RealVecChart& getRefCoord ( int quadpoint ) const  {
    static const RealVecChart c[3] = {    RealVecChart ( 1./2., 1./2. ),
                                        RealVecChart ( 0., 1./2. ),
                                        RealVecChart ( 1./2., 0. )
                                   };
    return c[quadpoint];
  }

  inline RealType getWeight ( int /*quadpoint*/ ) const { return 1./3.;}
};


// template <typename RealType>
// class CornerQuadrature {
// public:
//    static const int numQuadPoints = 3;
// 
//   inline const RealVecChart& getRefCoord ( int quadpoint ) const  {
//     static const RealVecChart c[3] = {
//                                               RealVecChart ( 0., 0. ),
//                                               RealVecChart ( 0., 1. ),
//                                               RealVecChart ( 1., 0. )
//                                             };
//     return c[quadpoint];
//   }
// 
//   inline RealType getWeight ( int ) const {
//     return 1./3.;
//   }
// };


//exact for polynomials of degree 3
template <typename RealType, typename RealVecChart>
class QuadriQuadrature {
public:
   static const int numQuadPoints = 4;

  inline const RealVecChart& getRefCoord ( int quadpoint ) const  {
    static const RealVecChart c[4] = {
                                              RealVecChart ( 1. / 3., 1. / 3. ),
                                              RealVecChart ( 3. / 5., 1. / 5. ),
                                              RealVecChart ( 1. / 5., 1. / 5. ),
                                              RealVecChart ( 1. / 5., 3. / 5. )
                                            };
    return c[quadpoint];
  }

  inline RealType getWeight ( int quadpoint) const {
    static const RealType w[4] = { 
                    RealType (-27./48. ),
                    RealType ( 25./48. ),
                    RealType ( 25./48. ),
                    RealType ( 25./48. )
                  };
    return w[quadpoint];
  }
};


//exact for polynomials of degree 6?
template <typename RealType, typename RealVecChart>
class DuodecQuadrature {
                                   
public:
   static const int numQuadPoints = 12;

  inline const RealVecChart& getRefCoord ( int quadpoint ) const  {
       static const RealVecChart c[12] = { RealVecChart ( 0.24928674517091 , 0.24928674517091 ),
                                         RealVecChart ( 0.24928674517091 , 0.50142650965818 ),
                                         RealVecChart ( 0.50142650965818 , 0.24928674517091 ),
                                         RealVecChart ( 0.06308901449150 , 0.06308901449150 ),
                                         RealVecChart ( 0.06308901449150 , 0.87382197101700 ),
                                         RealVecChart ( 0.87382197101700 , 0.06308901449150 ),
                                         RealVecChart ( 0.31035245103378 , 0.63650249912140 ),
                                         RealVecChart ( 0.63650249912140 , 0.05314504984482 ),
                                         RealVecChart ( 0.05314504984482 , 0.31035245103378 ),
                                         RealVecChart ( 0.63650249912140 , 0.31035245103378 ),
                                         RealVecChart ( 0.31035245103378 , 0.05314504984482 ),
                                         RealVecChart ( 0.05314504984482 , 0.63650249912140 )
                                       };
   return c[quadpoint];
  }
  inline RealType getWeight ( int quadpoint) const {
      static const RealType w[12] = { 
                     0.11678627572638 ,
                     0.11678627572638 ,
                     0.11678627572638 ,
                     0.05084490637021 ,
                     0.05084490637021 ,
                     0.05084490637021 ,
                     0.08285107561837 ,
                     0.08285107561837 ,
                     0.08285107561837 ,
                     0.08285107561837 ,
                     0.08285107561837 ,
                     0.08285107561837 
                  };    
    return w[quadpoint];
  }
};



// // TredecQuadrature, exact for polynomials of degree 6
// template <typename RealType>
// class TredecQuadrature {
// public:
//    static const int numQuadPoints = 13;
// 
//   inline const RealVecChart& getRefCoord ( int quadpoint ) const  {
//     static const RealVecChart c[13] = {
//                                               RealVecChart ( 0.33333333333333 , 0.33333333333333 ),
//                                               RealVecChart ( 0.26034596607904 , 0.26034596607904 ),
//                                               RealVecChart ( 0.26034596607904 , 0.47930806784192 ),
//                                               RealVecChart ( 0.47930806784192 , 0.26034596607904 ),
//                                               RealVecChart ( 0.06513010290222 , 0.06513010290222 ),
//                                               RealVecChart ( 0.06513010290222 , 0.86973979419557 ),
//                                               RealVecChart ( 0.86973979419557 , 0.06513010290222 ),
//                                               RealVecChart ( 0.31286549600487 , 0.63844418856981 ),
//                                               RealVecChart ( 0.63844418856981 , 0.04869031542532 ),
//                                               RealVecChart ( 0.04869031542532 , 0.31286549600487 ),
//                                               RealVecChart ( 0.63844418856981 , 0.31286549600487 ),
//                                               RealVecChart ( 0.31286549600487 , 0.04869031542532 ),
//                                               RealVecChart ( 0.04869031542532 , 0.63844418856981 )
//                                             };
//     return c[quadpoint];
//   }
// 
//   inline RealType getWeight ( int quadpoint) const {
//     static const RealType w[13] = { 
//                     RealType ( -0.14957004446768 ),
//                     RealType ( 0.17561525743321 ),
//                     RealType ( 0.17561525743321 ),
//                     RealType ( 0.17561525743321 ),
//                     RealType ( 0.05334723560884 ),
//                     RealType ( 0.05334723560884 ),
//                     RealType ( 0.05334723560884 ),
//                     RealType ( 0.07711376089026 ),
//                     RealType ( 0.07711376089026 ),
//                     RealType ( 0.07711376089026 ),
//                     RealType ( 0.07711376089026 ),
//                     RealType ( 0.07711376089026 ),
//                     RealType ( 0.07711376089026 )
//                   };
//     return w[quadpoint];
//   }
// };


//  exact for polynomials of degree 8
// template <typename RealType>
// class SexdecQuadrature {
// public:
//    static const int numQuadPoints = 16;
// 
//   inline const RealVecChart& getRefCoord ( int quadpoint ) const  {
//     static const RealVecChart c[16] = {
//     RealVecChart ( 0.333333333333333333333333333333333 , 0.333333333333333333333333333333333 ),
//     RealVecChart ( 0.459292588292723156028815514494169 , 0.459292588292723156028815514494169 ),
//     RealVecChart ( 0.459292588292723156028815514494169 , 0.081414823414553687942368971011661 ),
//     RealVecChart ( 0.081414823414553687942368971011661 , 0.459292588292723156028815514494169 ),
//     RealVecChart ( 0.170569307751760206622293501491464 , 0.170569307751760206622293501491464 ),
//     RealVecChart ( 0.170569307751760206622293501491464 , 0.658861384496479586755412997017071 ),
//     RealVecChart ( 0.658861384496479586755412997017071 , 0.170569307751760206622293501491464 ),
//     RealVecChart ( 0.050547228317030975458423550596598 , 0.050547228317030975458423550596598 ),
//     RealVecChart ( 0.050547228317030975458423550596598 , 0.898905543365938049083152898806802 ),
//     RealVecChart ( 0.898905543365938049083152898806802 , 0.050547228317030975458423550596598 ),
//     RealVecChart ( 0.263112829634638113421785786284643 , 0.728492392955404281241000379176062 ),
//     RealVecChart ( 0.728492392955404281241000379176062 , 0.008394777409957605337213834539294 ),
//     RealVecChart ( 0.008394777409957605337213834539294 , 0.263112829634638113421785786284643 ),
//     RealVecChart ( 0.728492392955404281241000379176062 , 0.263112829634638113421785786284643 ),
//     RealVecChart ( 0.263112829634638113421785786284643 , 0.008394777409957605337213834539294 ),
//     RealVecChart ( 0.008394777409957605337213834539294 , 0.728492392955404281241000379176062 )
//     };
//     return c[quadpoint];
//   }
// 
//   inline RealType getWeight ( int quadpoint) const {
//     static const RealType w[16] = { 
//                     RealType ( 0.144315607677787168251091110489064 ),
//                     RealType ( 0.095091634267284624793896104388584 ),
//                     RealType ( 0.095091634267284624793896104388584 ),
//                     RealType ( 0.095091634267284624793896104388584 ),
//                     RealType ( 0.103217370534718250281791550292129 ),
//                     RealType ( 0.103217370534718250281791550292129 ),
//                     RealType ( 0.103217370534718250281791550292129 ),
//                     RealType ( 0.032458497623198080310925928341780 ),
//                     RealType ( 0.032458497623198080310925928341780 ),
//                     RealType ( 0.032458497623198080310925928341780 ),
//                     RealType ( 0.027230314174434994264844690073908 ),
//                     RealType ( 0.027230314174434994264844690073908 ),
//                     RealType ( 0.027230314174434994264844690073908 ),
//                     RealType ( 0.027230314174434994264844690073908 ),
//                     RealType ( 0.027230314174434994264844690073908 ),
//                     RealType ( 0.027230314174434994264844690073908 )
//                   };
//     return w[quadpoint];
//   }
// };


//  exact for polynomials of degree 10
template <typename RealType,typename RealVecChart>
class GaussDegree10Quadrature {
public:
   static const int numQuadPoints = 25;

  inline const RealVecChart& getRefCoord ( int quadpoint ) const  {
    static const RealVecChart c[25] = {
        RealVecChart ( 0.333333333333333 , 0.333333333333333 ),
        RealVecChart ( 0.028844733232685 , 0.485577633383657 ),
        RealVecChart ( 0.485577633383657 , 0.028844733232685 ),
        RealVecChart ( 0.485577633383657 , 0.485577633383657 ),
        RealVecChart ( 0.781036849029926 , 0.109481575485037 ),
        RealVecChart ( 0.109481575485037 , 0.781036849029926 ),
        RealVecChart ( 0.109481575485037 , 0.109481575485037 ),
        RealVecChart ( 0.141707219414880 , 0.307939838764121 ),
        RealVecChart ( 0.141707219414880 , 0.550352941820999 ),
        RealVecChart ( 0.307939838764121 , 0.141707219414880 ),
        RealVecChart ( 0.307939838764121 , 0.550352941820999 ),
        RealVecChart ( 0.550352941820999 , 0.141707219414880 ),
        RealVecChart ( 0.550352941820999 , 0.307939838764121 ),
        RealVecChart ( 0.025003534762686 , 0.246672560639903 ),
        RealVecChart ( 0.025003534762686 , 0.728323904597411 ),
        RealVecChart ( 0.246672560639903 , 0.025003534762686 ),
        RealVecChart ( 0.246672560639903 , 0.728323904597411 ),
        RealVecChart ( 0.728323904597411 , 0.025003534762686 ),
        RealVecChart ( 0.728323904597411 , 0.246672560639903 ),
        RealVecChart ( 0.009540815400299 , 0.066803251012200 ),
        RealVecChart ( 0.009540815400299 , 0.923655933587500 ),
        RealVecChart ( 0.066803251012200 , 0.009540815400299 ),
        RealVecChart ( 0.066803251012200 , 0.923655933587500 ),
        RealVecChart ( 0.923655933587500 , 0.009540815400299 ),
        RealVecChart ( 0.923655933587500 , 0.066803251012200 )
    };
    return c[quadpoint];
  }

  inline RealType getWeight ( int quadpoint) const {
    static const RealType w[25] = { 
                    RealType ( 0.090817990382754 ),
                    RealType ( 0.036725957756467 ),
                    RealType ( 0.036725957756467 ),
                    RealType ( 0.036725957756467 ),
                    RealType ( 0.045321059435528 ),
                    RealType ( 0.045321059435528 ),
                    RealType ( 0.045321059435528 ),
                    RealType ( 0.072757916845420 ),
                    RealType ( 0.072757916845420 ),
                    RealType ( 0.072757916845420 ),
                    RealType ( 0.072757916845420 ),
                    RealType ( 0.072757916845420 ),
                    RealType ( 0.072757916845420 ),
                    RealType ( 0.028327242531057 ),
                    RealType ( 0.028327242531057 ),
                    RealType ( 0.028327242531057 ),
                    RealType ( 0.028327242531057 ),
                    RealType ( 0.028327242531057 ),
                    RealType ( 0.028327242531057 ),
                    RealType ( 0.009421666963733 ),
                    RealType ( 0.009421666963733 ),
                    RealType ( 0.009421666963733 ),
                    RealType ( 0.009421666963733 ),
                    RealType ( 0.009421666963733 ),
                    RealType ( 0.009421666963733 )
                  };
    return w[quadpoint];
  }
};



// template <typename RealType, typename RealVecChart, int N>
// class TestQuadrature{
//     
// protected:
//     std::vector<RealVecChart> _c;
//     
// public:
//   enum { numQuadPoints = (N+1) * (N+2) / 2 };
// 
//   
//   TestQuadrature() {
//         for( int x=0; x<= N; ++x )
//             for( int y=0; y<=N; ++y ){
//                 RealType xCoord = static_cast<RealType> ( x ) / static_cast<RealType> ( N );
//                 RealType yCoord = static_cast<RealType> ( y ) / static_cast<RealType> ( N );
//                 if( (xCoord + yCoord) <= 1. ){
//                     _c.push_back( RealVecChart ( xCoord, yCoord ) );
//                 }
//             }
//   }
//   
//   inline const RealVecChart& getRefCoord ( int quadpoint ) const  {
//     return _c[quadpoint];
//   }
// 
//   inline RealType getWeight ( int /*quadpoint*/ ) const { return 1./13.;}
// };






// template <typename RealType, typename RealVecChart>
// class TestQuadratureNodesEdgesBary {
// public:
//   enum { numQuadPoints = 7 };
// 
//   inline const RealVecChart& getRefCoord ( int quadpoint ) const  {
// 
//     static const RealVecChart c[7] = {    RealVecChart ( 0., 0. ),
//                                         RealVecChart ( 1., 0. ),
//                                         RealVecChart ( 0., 1. ),
//                                         RealVecChart ( 1./2., 1./2. ),
//                                         RealVecChart ( 0., 1./2. ),
//                                         RealVecChart ( 1./2., 0. ),
//                                         RealVecChart ( 1./3., 1./3. )
//                                    };
//     return c[quadpoint];
//   }
// 
//   inline RealType getWeight ( int /*quadpoint*/ ) const { return 1./7.;}
// };



#endif
