/*

v=[v(0),...,v(n-1)] \in R^{n}
*/

template <class Tp> 
Vector<Tp>::Vector (size_t m  ){  //ベクトル宣言。引数は要素数。
  n=m;
  Tp zero;
  v.resize(n);
  for(size_t i=0;i<n;i++){
    v[i]=zero;
  }
}
template <class Tp>
Vector<Tp>::Vector (const Vector<Tp>& w  ){  //ベクトルを複製
  n=w.n;
  v.resize(n);
  for(size_t i=0;i<n;i++){
    v[i]=w.v[i];
  }
}

//指定した要素を返す
template <class Tp>
Tp& Vector<Tp>::operator() ( size_t i ){
  return v[i];
}
template <class Tp>
Tp Vector<Tp>::operator[] ( size_t i )const{
  return v[i];
}

//一行演算子 -
template <class Tp>
Vector<Tp> Vector<Tp>::operator-()const{
  Vector<Tp> ans(n);
  for(size_t i=0;i<n;i++){
    ans.v[i]=-v[i];
  }
  return ans;
}

//減算
template <class Tp>
Vector<Tp> Vector<Tp>::operator- (const Vector<Tp>& w ){
  Vector<Tp> ans(n);
  for(size_t i=0;i<n;i++){
    ans.v[i]=v[i]-w.v[i];
  }
  return ans;
}

//加算
template <class Tp>
Vector<Tp> Vector<Tp>::operator+ (const Vector<Tp>& w ){
  Vector<Tp> ans(n);
  for(size_t i=0;i<n;i++){
    ans.v[i]=v[i]+w.v[i];
  }
  return ans;
}

//スカラー倍
template <class Tp>
Vector<Tp> Vector<Tp>::operator* (const double& c ){
  Vector<Tp> ans(n);
  for(size_t i=0;i<n;i++){
    ans.v[i]=v[i]*c;
  }
  return ans;
}

//代入
template <class Tp>
void Vector<Tp>::operator= (const Vector<Tp>& w ){
  if( n != w.n ){
    n=w.n;
    v.resize(n);
  }
  for(size_t i=0;i<n;i++){
    v[i]=w.v[i];
  }
}

//同値判定
template <class Tp>
bool Vector<Tp>::operator== (const Vector<Tp>& w ){
  for(size_t i=0;i<n;i++){
    if(v[i]!=w.v[i]){
      return false;
    }
  }
  return true;
}

//不同値判定
template <class Tp>
bool Vector<Tp>::operator!= (const Vector<Tp>& w ){
  return !((*this)==w);
}

//0クリア？
template <class Tp>
void Vector<Tp>::null(void){
  Tp zero;
  for(size_t i=0;i<n;i++){
    v[i]=zero;
  }
}
/*
template <class Tp>
double Norm( const Vector<Tp>& v){
      return sqrt(I(v,v));
}
template <class Tp>
double I( const Vector<Tp>& v1,const Vector<Tp>& v2 ){
      T i_p(0);
      for(size_t i=0; i< n ;i++){
	i_p=i_p+(v1.v)[i]*(v2.v)[i] ;
      }
      return i_p;
}
*/
