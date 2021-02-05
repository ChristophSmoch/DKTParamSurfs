#ifndef __PESOPTTENSOR_H
#define __PESOPTTENSOR_H

#include <iostream>
#include  <ostream>


namespace pesopt {

//! Tensor of order n+1 if VecType is of order n
template <class VecType, int dim>
class Tensor {
public:
  Tensor () {}

  Tensor ( const Tensor& tensor ) {
    for ( int i = 0; i < dim; ++i )  this->_vec [i] = tensor._vec [i];
  }

  Tensor& operator = ( const Tensor& tensor ) {
    // Beware of self-assignment
    if ( this == &tensor ) return *this;
    for ( int i = 0; i < dim; ++i ) this->_vec[i] = tensor._vec[i];
    return *this;
  }

 ~Tensor () {}

  const VecType& operator [] ( const int i ) const {return this->_vec [i];}
  VecType& operator [] ( const int i ) { return this->_vec [i];}

  std::ostream& print ( std::ostream& os ) const {
    for ( int i = 0; i < dim; ++i ) os << this->_vec [i] << std::endl;
    return os;
  }

  void setZero () { for ( int i = 0; i < dim; ++i )  this->_vec[i].setZero ();}
  
  template < typename RealType >
  Tensor &operator*= ( const RealType Alpha ) {
    for ( int i = 0; i < dim; ++i )  this->_vec[i] *= Alpha;
    return *this; 
  }
  
protected:
  VecType _vec[dim];
};


template <class DataType, class Matrix22>
class Tensor222 : public Tensor<Matrix22, 2 > {
public:
  Tensor222 () : Tensor<Matrix22, 2 > () {}

  Tensor222 ( const Tensor222<DataType, Matrix22>& tensor ) : Tensor<Matrix22, 2 > ( tensor ) {}

  Tensor222<DataType,Matrix22>& operator = ( const Tensor222<DataType,Matrix22>& tensor ) {
    // Beware of self-assignment
    if ( this == &tensor ) return *this;
    Tensor<Matrix22, 2 >::operator = ( tensor );
    return *this;
  }

  const DataType& operator () ( const int i, const int j, const int k ) const {  return this->_vec [i](j,k);}
  DataType& operator () ( const int i, const int j, const int k ) { return this->_vec[i](j,k);}
  DataType get ( const int i, const int j, const int k ) const { return this->_vec[i](j,k);}
  void set ( const int i, const int j, const int k, const DataType value ) {this->_vec[i](j,k) = value;}
  void add ( const int i, const int j, const int k, const DataType value ) {this->_vec[i](j,k) += value;}
};


template <class DataType, class Matrix22, class Vec3Type>
class Tensor322 : public Tensor<Matrix22, 3 > {
public:
  Tensor322 () : Tensor<Matrix22, 3 > () {}

  Tensor322 ( const Tensor322<DataType,Matrix22,Vec3Type>& tensor ) : Tensor<Matrix22, 3 > ( tensor ) {}

  Tensor322<DataType,Matrix22,Vec3Type>& operator = ( const Tensor322<DataType,Matrix22,Vec3Type>& tensor ) {
    // Beware of self-assignment
    if ( this == &tensor ) return *this;
    Tensor<Matrix22, 3 >::operator = ( tensor );
    return *this;
  }

  const DataType& operator () ( const int i, const int j, const int k ) const {return this->_vec [i](j,k);}
  DataType& operator () ( const int i, const int j, const int k ) {return this->_vec [i](j,k);}
  DataType get ( const int i, const int j, const int k ) const {return this->_vec [i](j,k);}
  void set ( const int i, const int j, const int k, const DataType value ) {this->_vec [i](j,k) = value;}

  void add ( const int i, const int j, const int k, const DataType value ) {this->_vec [i](j,k) += value;}
  
  Vec3Type getVector(const int j, const int k) const{
    Vec3Type vector(this->_vec[0](j,k), this->_vec[1](j,k), this->_vec[2](j,k));
    return vector;
  }
  
  void getVector(Vec3Type& vector, const int j, const int k) const{ vector = Vec3Type(this->_vec[0](j,k), this->_vec[1](j,k), this->_vec[2](j,k));}
  void setVector(const Vec3Type& vector , const int j, const int k){
    for (int i =0;i<3;i++) this->_vec[i](j,k) = vector[i];
  }
  
  DataType normSqr ( ) const{
    DataType result = 0.0;
    for( int i=0; i<3; ++i )
      for( int j=0; j<2; ++j )
        for( int k = 0; k <2; ++k )
          result += pesopt::Sqr( this->_vec[i](j,k) );
    return result;
  }
  
  DataType ddprod ( const pesopt::Tensor322<DataType,Matrix22,Vec3Type> &T ) const {
    DataType aux = 0.;
    for( int i=0; i<3; ++i )
      for( int j=0; j<2; ++j )
        for( int k = 0; k <2; ++k )
          aux += this->get(i,j,k) * T.get(i,j,k);
    return aux;
  }
};




template <class DataType, class Matrix32>
class Tensor332 : public Tensor<Matrix32, 3 > {
public:
  Tensor332 () : Tensor<Matrix32, 3 > () {}

  Tensor332 ( const Tensor332<DataType,Matrix32>& tensor ) : Tensor<Matrix32, 3 > ( tensor ) {}

  Tensor332<DataType,Matrix32>& operator = ( const Tensor332<DataType,Matrix32>& tensor ) {
    // Beware of self-assignment
    if ( this == &tensor ) return *this;
    Tensor<Matrix32, 3 >::operator = ( tensor );
    return *this;
  }

  const DataType& operator () ( const int i, const int j, const int k ) const {return this->_vec[i](j,k);}
  DataType& operator () ( const int i, const int j, const int k ) {return this->_vec[i](j,k);}
  DataType get ( const int i, const int j, const int k ) const {return this->_vec[i](j,k);}
  void set ( const int i, const int j, const int k, const DataType value ) {this->_vec[i](j,k) = value;}

  void add ( const int i, const int j, const int k, const DataType value ) {this->_vec[i](j,k) += value;}
  
  
  DataType normSqr ( ) const{
    DataType result = 0.0;
    for( int i=0; i<3; ++i )
      for( int j=0; j<3; ++j )
        for( int k = 0; k <2; ++k )
          result += pesopt::Sqr( this->_vec[i](j,k) );
    return result;
  }
  
  DataType ddprod ( const pesopt::Tensor332<DataType,Matrix32> &T ) const {
    DataType aux = 0.;
    for( int i=0; i<3; ++i )
      for( int j=0; j<3; ++j )
        for( int k = 0; k <2; ++k )
          aux += this->get(i,j,k) * T.get(i,j,k);
    return aux;
  }
  
};

template <class DataType, class Matrix33>
class Tensor333 : public Tensor<Matrix33, 3 > {
public:
  Tensor333 () : Tensor<Matrix33, 3 > () {}

  Tensor333 ( const Tensor333<DataType,Matrix33>& tensor ) : Tensor<Matrix33, 3 > ( tensor ) {}

  Tensor333<DataType,Matrix33>& operator = ( const Tensor333<DataType,Matrix33>& tensor ) {
    // Beware of self-assignment
    if ( this == &tensor ) return *this;
    Tensor<Matrix33, 3 >::operator = ( tensor );
    return *this;
  }

  const DataType& operator () ( const int i, const int j, const int k ) const {return this->_vec[i](j,k);}
  DataType& operator () ( const int i, const int j, const int k ) {return this->_vec[i](j,k);}
  DataType get ( const int i, const int j, const int k ) const {return this->_vec[i](j,k);}
  void set ( const int i, const int j, const int k, const DataType value ) {this->_vec[i](j,k) = value;}

  void add ( const int i, const int j, const int k, const DataType value ) {this->_vec[i](j,k) += value;}
  
  
  DataType normSqr ( ) const{
    DataType result = 0.0;
    for( int i=0; i<3; ++i )
      for( int j=0; j<3; ++j )
        for( int k = 0; k <3; ++k )
          result += pesopt::Sqr( this->_vec[i](j,k) );
    return result;
  }
  
  DataType ddprod ( const pesopt::Tensor333<DataType,Matrix33> &T ) const {
    DataType aux = 0.;
    for( int i=0; i<3; ++i )
      for( int j=0; j<3; ++j )
        for( int k = 0; k <3; ++k )
          aux += this->get(i,j,k) * T.get(i,j,k);
    return aux;
  }
  
};

} // namespace pesopt

#endif // __TENSOR_H
