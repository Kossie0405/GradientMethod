#include <iostream>
#include <cmath>
#include <fstream>
#include <octave/Matrix.h>
#include "vector.h"
using namespace std;
typedef ColumnVector X;
typedef ColumnVector U;
typedef Vector<X> X_T;
typedef Vector<U> U_T;
typedef Matrix Mat;
typedef double R;
typedef int Z;
typedef size_t N;
const N n(4);
const N m(1);
const R eAlpha(1.0e-9);
const R eU(1.0e-3);
const R eX(1.0e-3);
const R eJ(1.0);

//物理パラメータ
const R g(9.8);
const R l(0.4);

/*
U_T sumVectorSet(U_T u1, U_T u2, N k){
    U_T ans(k+1);
    N i;
    for(i = 0; i <= k; i++){
        ans(i) = u1(i)+u2(i);
    }
    return ans;
}
*/

X calDphiDX(X x, X x1){
    X dx(n);
    dx = x - x1;
    Mat S(n, n, 0.0);
    N i;
    for(i = 0; i < n; i++){
        S(i, i) = 1.0;
    }
    S(1, 1) = 10.0;

    S *= 1.0e+4;

    X ans(n);
    ans = S.transpose()*dx;

    return ans;
}

R calPhi(X x, X x1){
    X dx(n);
    dx = x - x1;
    Mat S(n, n, 0.0);
    N i;
    for(i = 0; i < n; i++){
        S(i, i) = 1.0;
    }
    S(1, 1) = 10.0;

    S *= 1.0e+4;

    R ans;

    ans = (dx.transpose() * S * dx)/2.0;

    return ans;
}

R calL(X x, X x1, U u){
    X dx(n);
    dx = x - x1;
    Mat Q(n, n, 0.0);
    N i;
    for(i = 0; i < n; i++){
        Q(i, i) = 1.0;
    }
    R r(1.0), ans;

    ans = (dx.transpose()*Q*dx)/2.0 + (u.transpose()*r*u)/2.0;

    return ans;
}

X calDLDX(X x, X x1){
    X dx(n);
    dx = x - x1;
    Mat Q(n, n, 0.0);
    N i;
    for(i = 0; i < n; i++){
        Q(i, i) = 1.0;
    }
    X ans(n);

    ans = Q.transpose()*dx;

    return ans;
}

R calDLDU(U u){
    R r(1.0);
    return r*u(0);
}

R integralL(X_T x_t, U_T u_t, X x1, R dt, N k){
    R ans(0.0);
    N i;
    ans += calL(x_t(0), x1, u_t(0)) + calL(x_t(k), x1, u_t(k));
    for(i = 1; i < k; i++){
        if(i%2){
            ans += 4.0*calL(x_t(i), x1, u_t(i));
        }else{
            ans += 2.0*calL(x_t(i), x1, u_t(i));
        }
    }
    ans *= (dt/3.0);
    return ans;
}

X calX(X x, U u){
    X ans(n);
    ans(0) = x(2);
    ans(1) = x(3);
    ans(2) = u(0);
    ans(3) = -g/l*sin(x(1)) - u(0)/l*cos(x(1));

    return ans;
}

X calXD(X x, U u, R dt){
    X ans(n);
    X dx1(n), dx2(n), dx3(n), dx4(n);
    const R w1(1.0 / 6.0), w2(1.0 / 3.0), w3(w2), w4(w1);

    dx1 = calX(x, u) * dt;
    dx2 = calX(x + dx1 * 0.5, u) * dt;
    dx3 = calX(x + dx2 * 0.5, u) * dt;
    dx4 = calX(x + dx3, u) * dt;

    ans = x + (w1 * dx1 + w2 * dx2 + w3 * dx3 + w4 * dx4);
    
    return ans;
}

X_T calIVP_X(X_T x_t, U_T u_t, X x0, R dt, N k){
    N i;
    x_t(0) = x0;
    for(i = 0; i < k; i++){
        x_t(i+1) = calXD(x_t(i), u_t(i), dt);
    }
    return x_t;
}

X calRamda(X ramda, X x, X x1, U u){
    X ans(n);
    ans = calDLDX(x, x1);

    ans(1) += ramda(3)*(-g/l*cos(x(1))+u(0)/l*sin(x(1)));
    ans(2) += ramda(0);
    ans(3) += ramda(1);

    ans *= -1.0;

    return ans;
}

X calRamdaD(X ramda, X x, U u, X x1, R dt){
    X ans(n);
    X dx1(n), dx2(n), dx3(n), dx4(n);
    const R w1(1.0 / 6.0), w2(1.0 / 3.0), w3(w2), w4(w1);

    dx1 = calRamda(ramda, x, x1, u) * dt;
    dx2 = calRamda(ramda + dx1 * 0.5, x, x1, u) * dt;
    dx3 = calRamda(ramda + dx2 * 0.5, x, x1, u) * dt;
    dx4 = calRamda(ramda + dx3, x, x1, u) * dt;

    ans = ramda + (w1 * dx1 + w2 * dx2 + w3 * dx3 + w4 * dx4);

    return ans;
}

X_T calFVP_R(X_T ramda_t, X_T x_t, U_T u_t, X x1, R dt, N k){
    N i;
    for(i = k; i > 0; i--){
        ramda_t(i-1) = calRamdaD(ramda_t(i), x_t(i), u_t(i), x1, -dt);
    }
    return ramda_t;
}

R calJ(X_T x_t, U_T u_t, X x1, X x0, R dt, N k){
    R ans(0.0);
    N i;
    x_t = calIVP_X(x_t, u_t, x0, dt, k);

    ans += calPhi(x_t(k), x1);
    ans += integralL(x_t, u_t, x1, dt, k);

    return ans;
}

R calGradJ(X x, X ramda, U u){
    return (calDLDU(u)+ramda(2)-ramda(3)/l*cos(x(1)));
}

U_T calDK(X_T x_t, X_T ramda_t, U_T u_t, N k){
    U_T dk(k+1);
    U dk_temp(m);
    N i;
    for(i = 0; i <= k; i++){
        dk_temp(0) = -calGradJ(x_t(i), ramda_t(i), u_t(i));
        dk(i) = dk_temp;
    }
    return dk;
}


R NormJ(X_T x_t, X_T ramda_t, U_T u_t, N k){
    N i;
    R ans(0.0);
    for(i = 0; i < k; i++){
        ans += pow(calGradJ(x_t(i), ramda_t(i), u_t(i)), 2);
    }
    ans = sqrt(ans);
    cout << "Nolm of gradJ: " << ans << endl;
    return ans;
}

R NormU(U_T u1, U_T u2, N k){
    U u1_k(m), u2_k(m);
    N i;
    R ans(0.0);
    for(i = 0; i <= k; i++){
        u1_k = u1(i);
        u2_k = u2(i);
        ans += pow(u1_k(0)-u2_k(0), 2);
    }
    ans = sqrt(ans);
    cout << "Norm of du :" << ans << endl;
    return ans;
}

R error(X x, X x1){
    X dx(n);
    dx = x - x1;
    N i;
    R ans(0.0);
    for(i = 0; i < n;i++){
        ans += pow(dx(i), 2);
    }
    ans = sqrt(ans);
    return ans;
}

R searchAlpha(X_T x_t, X_T ramda_t, U_T u_t, U_T dk, X x1, X x0, R dt, N k){
    R ansAlpha(0.0), alphaMax, alphaMin(0.0), a1, a2, h;
    const R r(0.618);
    N count(0);

    if(!(error(x_t(k), x1) < eX)){
        h = 1.0e-6;
    }else{
        h = 1.0e-9;
    }

    //囲い込み
    while(!(calJ(x_t, (u_t+dk*(ansAlpha+h)), x1, x0, dt, k) >= calJ(x_t, (u_t+(dk*ansAlpha)), x1, x0, dt, k))){
        ansAlpha+=h;
        count++;
    }

    //alphaの存在区間の定義
    if(count){
        alphaMin = ansAlpha - h;
    }
    alphaMax = ansAlpha + h;

    //a1, a2の初期設定
    a1 = alphaMin + (1.0 - r)*(alphaMax - alphaMin);
    a2 = alphaMin + r * (alphaMax - alphaMin);

    //黄金分割法
    while(!(abs(alphaMax-alphaMin) < eAlpha)){
        if(calJ(x_t, (u_t+(dk*a1)), x1, x0, dt, k) < calJ(x_t, (u_t+(dk*a2)), x1, x0, dt, k)){
            alphaMax = a2;
            a2 = a1;
            a1 = alphaMin + (1.0 - r)*(alphaMax - alphaMin);
        }else{
            alphaMin = a1;
            a1 = a2;
            a2 = alphaMin + r * (alphaMax - alphaMin);
        }
    }

    ansAlpha = (alphaMax + alphaMin)/2.0;

    return ansAlpha;
}

/*
R searchAlpha(X_T x_t, X_T ramda_t, U_T u_t, U_T dk, X x1, X x0, R dt, N k){
    N i;
    R ansAlpha(1.0);
    R A(1.0e-4), B(0.5), C(0.0);
    X dk_temp(n);

    for(i = 0; i <= k; i++){
        dk_temp = dk(i);
        C += A * calGradJ(ramda_t(i), x_t(i), u_t(i)) * dk_temp(0);
    }

    while(!(((calJ(x_t, (u_t+dk*(ansAlpha)), x1, x0, dt, k) - calJ(x_t, u_t, x1, x0, dt, k)) <= (ansAlpha * C)))){
        ansAlpha *= B;
    }
    return ansAlpha;
}
*/

void dispResult(X_T x_t, U_T u_t, X x0, R dt, N k){
    N i;
    x_t = calIVP_X(x_t, u_t, x0, dt, k);
    for(i = 0; i < k+1; i++){
        cout << dt*i << " " << x_t(i).transpose() << " " << u_t(i).transpose() << endl; 
    }
}

void writeU(U_T u_t, N k){
    ofstream out_u("u.dat");
    N i;
    for(i = 0; i <= k; i++){
        out_u << u_t(i).transpose() << endl;
    }
    out_u.close();
}

Z main(void){
    X x0(n), x1(n);
    R dt(0.001), t0(0.0), t1(1.587), alpha(0.0);
    N k((N)((t1-t0)/dt)), k0(0), k1(k0 + k);      //kは時間ステップ
    N count(0), i;
    X_T x_t(k + 1), ramda_t(k + 1);
    U_T u_t(k + 1), dk(k + 1), u_Next(k + 1);
    U udk_initial(m);
    udk_initial(0) = 0.0;

    for(i = 0; i <= k; i++){
        u_t(i) = dk(i) = u_Next(i) = udk_initial;
    }

    //U_T unext_t(k + 1);
    R nolm;

    x0(0) = 0.0;
    x0(1) = 0.0;
    x0(2) = 0.0;
    x0(3) = 0.0;

    x1(0) = 0.0;
    x1(1) = M_PI;
    x1(2) = 0.0;
    x1(3) = 0.0;

    do{
        count++;
        x_t(0) = x0;
        u_t = u_Next;


        cout << count << " start" << endl;

        //step 1
        cout << "step 1" << endl;

        x_t = calIVP_X(x_t, u_t, x0, dt, k);
        
        ramda_t(k) = calDphiDX(x_t(k), x1);
        ramda_t = calFVP_R(ramda_t, x_t, u_t, x1, dt, k);
        
        dk = calDK(x_t, ramda_t, u_t, k);

        //step 2
        cout << "step 2" << endl;
        cout << "Final value of x: " << x_t(k).transpose() << endl;

        alpha = searchAlpha(x_t, ramda_t, u_t, dk, x1, x0, dt, k);

        cout << "step 3" << endl;
        cout << "alpha = " << alpha << endl;
        //step 3
        u_Next = u_t+dk*alpha;
        
        if(!(count%10)){
            writeU(u_t, k);
        }
    }while(!(NormU(u_t, u_Next, k) < eU));

    dispResult(x_t, u_t, x0, dt, k);

    return 0;
}