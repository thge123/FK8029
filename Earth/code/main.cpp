#include "HEADER.h"

VectorXd get_eps(int N){

    VectorXd eps(N+1);
    for (int i=0; i<N+1; i++)
        eps(i) = 0.78;
    eps[0] = 1;
    return eps;
}

MatrixXd get_A(int N){
    
    MatrixXd A(N+1,N+1);
    VectorXd eps(N+1);
    A   = zeros(N+1,N+1);
    eps = get_eps(N+1);

    for (int i=1; i<N+1; i++){
        
        // Lower layers
        for (int j=0; j<i; j++){
            A(i,j) = -eps(j)*eps(i);
            for (int k=j+1; k<i; k++)
                A(i,j) *= (1-eps(k));
        }

        A(i,i) = 2*eps(i);

        // Upper layers
        for (int j=i+1; j<N+1; j++){
            A(i,j) = -eps(j)*eps(i);
            for (int k=i+1; k<j; k++)
                A(i,j) *= (1-eps(k));
        }
    }

    // Last eq. Total out = Total in
    for (int j=0; j<N+1; j++){
        A(0,j) = eps(j);
        for (int k=j+1; k<N+1; k++)
            A(0,j) *= (1-eps(k));
    }

    return A;
}

VectorXd get_b(int N){

    VectorXd b(N+1);
    b = zeros(N+1);

    // (S0/4)*(1-alpha) in W/m^2
    b(0) = 239.05; 

    return b;

}

int main(){

    // Number of layers
    int N = 2;

    MatrixXd A;
    VectorXd b;
    VectorXd x;

    A = get_A(N);
    b = get_b(N);

    // Solution sigma*T^4
    x = A.colPivHouseholderQr().solve(b);

    double T;
    double sigma = 5.670374e-8;
    for (int i=0; i<x.rows(); i++){
        T = x[i]/sigma;
        T = pow(T,0.25);
        cout << T << endl;
    }
    cout << A << endl;
}


