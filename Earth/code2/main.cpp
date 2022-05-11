#include "HEADER.h"

VectorXd get_eps(int N){

    VectorXd eps(N+1);
    double sigma = 0.5;
    double h     = (double) 10/N;
    double H     = 10.4;

    eps(0) = 1;
    for (int i=1; i<N+1; i++){
        eps(i) = sigma;
        eps(i)*= 0.5*exp(-i*h/H);
    }

    return eps;
}

MatrixXd get_A(int N, VectorXd *eps){
    
    MatrixXd A(N+1,N+1);
    A   = zeros(N+1,N+1);

    for (int i=1; i<N+1; i++){
        
        // Lower layers
        for (int j=0; j<i; j++){
            A(i,j) = -(*eps)(j)*(*eps)(i);
            for (int k=j+1; k<i; k++)
                A(i,j) *= exp(-(*eps)(k));
        }

        A(i,i) = 2*(*eps)(i);

        // Upper layers
        for (int j=i+1; j<N+1; j++){
            A(i,j) = -(*eps)(j)*(*eps)(i);
            for (int k=i+1; k<j; k++)
                A(i,j) *= exp(-(*eps)(k));
        }
    }

    // Last eq. Total out = Total in
    for (int j=0; j<N+1; j++){
        A(0,j) = (*eps)(j);
        for (int k=j+1; k<N+1; k++)
            A(0,j) *= exp(-(*eps)(k));
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

void OneLayer(){

    ofstream OutStream;
    OutStream.open("OneLayer.dat");
    
    // Number of layers
    int N = 1;

    MatrixXd A;
    VectorXd b;
    VectorXd x;
    VectorXd eps(2);
    b = get_b(N);

    double e;
    double T;
    double sigma = 5.670374e-8;
    eps(0) = 1;

    for (int i=0; i<1001; i++){
        e = (double) i/1000;
        eps(1) = e;
        A = get_A(N,&eps);

        // Solution sigma*T^4
        x = A.colPivHouseholderQr().solve(b);

        T = x[0]/sigma;
        T = pow(T,0.25);
        OutStream << e << ";" << T << ";";
        T = x[1]/sigma;
        T = pow(T,0.25);
        OutStream << T << "\n";
    }
    OutStream.close();
}
    
void TwoLayer(){

    ofstream OutStream;
    OutStream.open("TwoLayer.dat");
    
    // Number of layers
    int N = 2;

    MatrixXd A;
    VectorXd b;
    VectorXd x;
    VectorXd eps(3);
    b = get_b(N);
    eps(0) = 1;

    double e1,e2;
    double T;
    double sigma = 5.670374e-8;

    double n = 50;
    for (int i=0; i<n+1; i++){
        e1 = (double) i/n;
        eps(1) = e1;
        for (int i=0; i<n+1; i++){
            e2 = (double) i/n;
            eps(2) = e2;
            A = get_A(N,&eps);

            // Solution sigma*T^4
            x = A.colPivHouseholderQr().solve(b);

            T = x(0)/sigma;
            T = pow(T,0.25);
            OutStream << e1 << ";" << e2 << ";" << T << ";";
            T = x(1)/sigma;
            T = pow(T,0.25);
            OutStream << T << ";";
            T = x(2)/sigma;
            T = pow(T,0.25);
            OutStream << T << "\n";
        }
    }

    OutStream.close();
}

void US_STD_ATM(){

    ofstream OutStream;
    OutStream.open("US_STD_ATM.dat");

    double rho[20];
    double z[21];
    rho[0]  = 1.225;      z[0]  = 0;
    rho[1]  = 1.112;      z[1]  = 1000;
    rho[2]  = 1.007;      z[2]  = 2000;
    rho[3]  = 0.9093;     z[3]  = 3000;       
    rho[4]  = 0.8194;     z[4]  = 4000;      
    rho[5]  = 0.7364;     z[5]  = 5000;      
    rho[6]  = 0.6601;     z[6]  = 6000;     
    rho[7]  = 0.5900;     z[7]  = 7000;     
    rho[8]  = 0.5258;     z[8]  = 8000;    
    rho[9]  = 0.4671;     z[9]  = 9000;    
    rho[10] = 0.4135;     z[10] = 10000;   
    rho[11] = 0.1948;     z[11] = 15000;          
    rho[12] = 0.08891;    z[12] = 20000;
    rho[13] = 0.04008;    z[13] = 25000;
    rho[14] = 0.01841;    z[14] = 30000;
    rho[15] = 0.003996;   z[15] = 40000;
    rho[16] = 0.001027;   z[16] = 50000;
    rho[17] = 0.0003097;  z[17] = 60000;
    rho[18] = 0.00008283; z[18] = 70000;  
    rho[19] = 0.00001846; z[19] = 80000; 
                          z[20] = 90000;
    
    // Number of layers
    int N = 19;

    MatrixXd A;
    VectorXd b;
    VectorXd x;
    VectorXd eps(N+1);
    b = get_b(N);

    for (int i=0; i<N+1; i++){
        eps(i) = rho[i]/rho[0];
        eps(i)*= (z[i+1]-z[i])/(z[1]-z[0]);
    }

    double T;
    double sigma = 5.670374e-8;

    A = get_A(N,&eps);
    x = A.colPivHouseholderQr().solve(b);

    for (int i=0; i<x.rows(); i++){
        T = x(i)/sigma;
        T = pow(T,0.25);
        cout << T << endl;
    }
    cout << eps << endl;

    OutStream.close();
}
    

int main(){

    US_STD_ATM();

}

