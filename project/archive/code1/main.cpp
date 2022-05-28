#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double EOS(double p, double K, double n){

    return 1e-3;
    double E;
    E = K*pow(p,n/(n+1));
    return E;
}

int main(){

    ofstream pmStream;
    pmStream.open("pm.dat");

    int    N = 10000;
    cout << "Write N: ";
    cin  >> N;
    double X = 100;
    double h = X/N;

    double x[N];
    double m[N];
    double m_Newt[N];
    double p[N];
    double p_Newt[N];
    double E[N];
    double E_Newt[N];

    for (int i=0; i<N; i++){
        x[i] = h+i*h;
        m[i] = 0;
        m_Newt[i] = 0;
        E[i] = 0;
        E_Newt[i] = 0;
        p[i] = 0;
        p_Newt[i] = 0;
    }
    
    cout << "p0: ";
    cin  >> p[0];
    p_Newt[0] = p[0];

    double K;
    double n = 0.5;
    cout << "K: ";
    cin  >> K;
    
    E[0] = EOS(p[0],K,n);
    E_Newt[0] = EOS(p_Newt[0],K,n);
    pmStream << x[0] << ";"
             << p[0] << ";"
             << p_Newt[0] << ";"
             << m[0] << ";"
             << m_Newt[0] << "\n";

    double f1,f2;
    double f3,f4;
    for (int j=0; j<N-1; j++){

        // Heun's method for GR
        f1 = -(m[j]+x[j]*x[j]*x[j]*p[j]);
        f1*= (E[j]+p[j])/(x[j]*(x[j]-2*m[j]));
        p[j+1] = p[j] + h*f1;
        E[j+1] = EOS(p[j+1],K,n);

        f2 = x[j]*x[j]*E[j];
        m[j+1] = m[j] + h*f2;

        f3 = -(m[j+1]+x[j+1]*x[j+1]*x[j+1]*p[j+1]);
        f3*= (E[j+1]+p[j+1])/(x[j+1]*(x[j+1]-2*m[j+1]));
        p[j+1] = p[j] + h*(f1+f3)/2;
        E[j+1] = EOS(p[j+1],K,n);
        
        f4 = x[j+1]*x[j+1]*E[j+1];
        m[j+1] = m[j] + h*(f2+f4)/2;

        // Heun's method for Newt
        f1 = -m[j];
        f1*= E[j]/(x[j]*x[j]);
        p_Newt[j+1] = p_Newt[j] + h*f1;
        E_Newt[j+1] = EOS(p_Newt[j+1],K,n);

        f2 = x[j]*x[j]*E_Newt[j];
        m_Newt[j+1] = m_Newt[j] + h*f2;

        f3 = -m_Newt[j+1];
        f3*= E_Newt[j+1]/(x[j+1]*x[j+1]);
        p_Newt[j+1] = p_Newt[j] + h*(f1+f3)/2;
        E_Newt[j+1] = EOS(p_Newt[j+1],K,n);
        
        f4 = x[j+1]*x[j+1]*E_Newt[j+1];
        m_Newt[j+1] = m_Newt[j] + h*(f2+f4)/2;

        // Write to file
        pmStream << x[j+1] << ";"
                 << p[j+1] << ";"
                 << p_Newt[j+1] << ";"
                 << m[j+1] << ";"
                 << m_Newt[j+1] << "\n";
        
        // EOS not defined for P<0? Either way end of star
        if (p[j+1] < 0)
            break;
    }

    pmStream.close();
    
}



    

    
