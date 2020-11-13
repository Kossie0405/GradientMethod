#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <valarray>
//#include <iomanip>

using namespace std;

template <class Tp>
class Vector {
  // Associations
  // Attributes
  private:
    valarray<Tp> v;
    size_t n;
  // Operations
  public:
    Vector<Tp> (size_t m=1);
    // Copy constructor
    Vector<Tp> ( const Vector<Tp>&  );

    // Assignment operators
    Tp& operator() ( size_t i );
    Tp  operator[] ( size_t i )const;

    size_t dim(){return n;};
    void resize(size_t m){n=m;v.resize(n);};
    void null(void);

    //unary operators
    Vector<Tp> operator+ (  )const{ return *this;};
    Vector<Tp> operator- (  )const;


    Vector<Tp> operator+ (const Vector<Tp>& w );
    Vector<Tp> operator* (const double& c );

    Vector<Tp> operator- (const Vector<Tp>& u );

    void operator= (const Vector<Tp>& w );

    bool operator== (const Vector<Tp>& w );
    bool operator!= (const Vector<Tp>& w );

    friend double Norm( const Vector<Tp>& v){
      return sqrt(I(v,v));
    };

    friend double I( const Vector<Tp>& v1,const Vector<Tp>& v2 ){
      double i_p(0);
      for(size_t i=0; i< v1.n ;i++){
	i_p=i_p+(v1.v)[i]*(v2.v)[i] ;
      }
      return i_p;
    };

    friend        istream&  operator>>( istream& str, Vector<Tp>& w){
                  for(size_t i=0; i< w.n ;i++){
                    cin  >> (w.v)[i] ;
                  }
                  return str;
    };
   friend       ostream&  operator<<( ostream& str, Vector<Tp>& w){
                  cout << "col[ " ;
                  for(size_t i=0;i< w.n ;i++){
                    cout << (w.v)[i] << "  "   ;
                  }
                  cout << " ]" << endl  ;
                  return str;
   };
   
};



#include "vector.cpp"
#endif
