#include "HEADER.h"

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

Matrix<double,Dynamic,1> zeros(int N){

    /* Returns MxN sized zero matrix. */
    
    Matrix<double,Dynamic,1> V;
    V.resize(N);
    for (int j=0; j<N; j++){
         V(j) = 0;
    }
    return V;
}
