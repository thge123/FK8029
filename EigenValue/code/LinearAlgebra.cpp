#include "HEADER.h"

Matrix<double,Dynamic,1> zeros(int N){

    /* Returns N-dimensional zero vector. */
    
    Matrix<double,Dynamic,1> V;
    V.resize(N,1);
    for (int i=0; i<N; i++) V(i) = (0,0);
    return V;
}

Matrix<double,Dynamic,Dynamic> zeros(int M, int N){

    /* Returns MxN sized zero matrix. */
    
    Matrix<double,Dynamic,Dynamic> V;
    V.resize(M,N);
    for (int i=0; i<M; i++){
        for (int j=0; j<N; j++){
             V(i,j) = (0,0);
        }
    }
    return V;
}

void PowerIter(MatrixXd &A, VectorXd &x0){

    // Power iteration with matrix A using initial guess x0

    int N = x0.rows();
    double norm;
    for (int i=0; i<1000; i++){
        
        x0 = A*x0;
        x0.normalize();
    }
}
        
        

Matrix<complex<double>,Dynamic,1> solve(MatrixXcd &A, VectorXcd &b){

    // Solve using partial LU
    Matrix<complex<double>,Dynamic,1> x;
    x.resize(b.rows(),1);
    x = A.partialPivLu().solve(b);
    return x;
}
