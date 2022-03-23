#include "HEADER.h"

Matrix<complex<double>,Dynamic,1> zeros(int N){

    /* Returns N-dimensional zero vector. */
    
    Matrix<complex<double>,Dynamic,1> V;
    V.resize(N,1);
    for (int i=0; i<N; i++) V(i) = (0,0);
    return V;
}

Matrix<complex<double>,Dynamic,Dynamic> zeros(int M, int N){

    /* Returns MxN sized zero matrix. */
    
    Matrix<complex<double>,Dynamic,Dynamic> V;
    V.resize(M,N);
    for (int i=0; i<M; i++){
        for (int j=0; j<N; j++){
             V(i,j) = (0,0);
        }
    }
    return V;
}

Matrix<complex<double>,Dynamic,1> solve(MatrixXcd &A, VectorXcd &b){

    Matrix<complex<double>,Dynamic,1> x;
    x.resize(b.rows(),1);
    x = A.colPivHouseholderQr().solve(b);
    return x;
}
