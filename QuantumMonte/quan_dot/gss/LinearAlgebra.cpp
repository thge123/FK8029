#include "HEADER.h"

Matrix<double,Dynamic,1> zeros(int N){

    /* Returns N-dimensional zero vector. */
    
    Matrix<double,Dynamic,1> V;
    V.resize(N,1);
    for (int i=0; i<N; i++) V(i) = 0;
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

VectorXd invPowerIter(SparseMatrix<double> &B, 
                      double shift,
                      double *lam,
                      double *err){

    // Power iteration with matrix A using initial guess x0

    int N = B.rows();
    VectorXd x0;
    x0 = VectorXd::Random(N);

    MatrixXd I(N,N);
    I = MatrixXd::Identity(N,N);
    SparseMatrix<double> A = B-shift*I;

    SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    VectorXd y(N),Delta(N);
    double MAX;
    int k = 1;
    while (1) {

        y = x0.normalized();
        x0 = solver.solve(y);
        if (k%10 == 0){
            *lam = x0.dot(B*x0)/x0.dot(x0);
            Delta = (B*x0-(*lam)*x0);
            MAX = 0;
            for (int i=0; i<N; i++){
                if (abs(Delta(i)) > MAX)
                    MAX = abs(Delta(i));
            }
        if (MAX < 1e-6){
            *err = MAX;
            break;
        }
        if (k>1000){
            cout << "Eigenvector may not have been found for shift = " << shift << endl;
            *err = MAX;
            break;
            }
        }
        k++;
    }
    return x0;
}
        
Vector4d copy(Vector4d *X){
    
    int N = (*X).rows();
    Vector4d Y(N);
    for (int i=0; i<N; i++) Y(i) = (*X)(i);
    return Y;
}        

//Matrix<complex<double>,Dynamic,1> solve(MatrixXcd &A, VectorXcd &b){
//
//    // Solve using partial LU
//    Matrix<complex<double>,Dynamic,1> x;
//    x.resize(b.rows(),1);
//    x = A.partialPivLu().solve(b);
//    return x;
//}
