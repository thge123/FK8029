#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

// Star struct
struct Star{

    double h;   // Step size
    double pc;  // Central pressure
    double E;   // EOS constant

};

double P(double n){

    double a = 2.10019e-1;
    double xF = a*pow(3*M_PI*M_PI*n,0.3333);
    double yF = sqrt(1+xF*xF);
    double p = 2*xF*xF*xF*yF/3;
    p-= xF*yF;
    p+= log(xF+yF);
    return p;
}

double E(double n){

    double a = 2.10019e-1;
    double xF = a*pow(3*M_PI*M_PI*n,0.3333);
    double yF = sqrt(1+xF*xF);
    double e  = xF*yF;
    e*= (xF*xF+yF*yF);
    e-= log(xF+yF);
    return e;
}


double bisection(double (*f)(double x),
                  double a, double b,
                  int *FLAG, double shift=0){

    int max_iters = 10000;
    double TOL = 1e-8;
    int iters = 0;
    double fa,fc,c;
    while (iters<max_iters){
        c = (a+b)/2;
        fc = f(c)-shift;
        fa = f(a)-shift;
        if (fc*fc < TOL*TOL){
            return c;
        }
        else{
            if (fc*fa < 0)
                b = c;
            else
                a = c;
        }
        iters ++;
    }
    *FLAG = 1;
    return c;
}

double EOS(double p){

    return 1e-3;
}


void TOV_Heun(Star *star){

    ofstream GRpm;
    GRpm.open("GRpm.dat");
    GRpm.precision(16);

    double max_iters = 1000000;
    double iters = 0;
    double x1,x2;
    double p1,p2;
    double m1,m2;
    double E1,E2;
    double h = star -> h;
    double E = star -> E;

    x1 = 1e-8;
    x2 = x1+h;
    p1 = star -> pc;
    E1 = E;
    m1 = E1*x1*x1*x1/3;

    GRpm << x1 << ";"
         << p1 << ";"
         << m1 << ";"
         << E1 << "\n";

    double f1,f2;
    double f3,f4;
    while (p1 > 0 && iters<max_iters){

        // Heun's method
        f1 = -(m1+x1*x1*x1*p1);
        f1*= (E1+p1)/(x1*(x1-2*m1));
        p2 = p1 + h*f1;
        E2 = E;

        f2 = x1*x1*E1;
        m2 = m1 + h*f2;

        f3 = -(m2+x2*x2*x2*p2);
        f3*= (E2+p2)/(x2*(x2-2*m2));
        p2 = p1 + h*(f1+f3)/2;
        E2 = E;
        
        f4 = x2*x2*E2;
        m2 = m1 + h*(f2+f4)/2;

        // Write to file
        GRpm << x2 << ";"
             << p2 << ";"
             << m2 << ";"
             << E2 << "\n";

        // Rename variables
        p1 = p2;
        m1 = m2;
        E1 = E2;
        x1 = x2;
        x2 = x2+h;
        iters++; 
    }

    GRpm.close();
    
    // Append last point (pc,x,m) to xm file
    fstream GRxm;
    GRxm.precision(16);
    GRxm.open("GRxm.dat",ios::app);
    GRxm << star -> pc << ";" 
         << x1 << ";" 
         << m1 << endl;
    
    GRxm.close();
    
}

void NEWT_Heun(Star *star){

    ofstream NEWTpm;
    NEWTpm.open("NEWTpm.dat");
    NEWTpm.precision(16);

    double max_iters = 1000000;
    double iters = 0;
    double x1,x2;
    double p1,p2;
    double m1,m2;
    double E1,E2;
    double h = star -> h;
    double E = star -> E;

    x1 = 1e-8;
    x2 = x1+h;
    p1 = star -> pc;
    E1 = E;
    m1 = E1*x1*x1*x1/3;

    NEWTpm << x1 << ";"
           << p1 << ";"
           << m1 << ";"
           << E1 << "\n";

    double f1,f2;
    double f3,f4;
    while (p1 > 0 && iters<max_iters){

        // Heun's method
        f1 = -m1;
        f1*= E1/(x1*x1);
        p2 = p1 + h*f1;
        E2 = E;

        f2 = x1*x1*E1;
        m2 = m1 + h*f2;

        f3 = -m2;
        f3*= E2/(x2*x2);
        p2 = p1 + h*(f1+f3)/2;
        E2 = E;
        
        f4 = x2*x2*E2;
        m2 = m1 + h*(f2+f4)/2;

        // Write to file
        NEWTpm << x2 << ";"
               << p2 << ";"
               << m2 << ";"
               << E2 << "\n";

        // Rename variables
        p1 = p2;
        m1 = m2;
        E1 = E2;
        x1 = x2;
        x2 = x2+h;
        iters++; 
    }

    NEWTpm.close();
    
    // Append last point (pc,x,m) to xm file
    fstream NEWTxm;
    NEWTxm.precision(16);
    NEWTxm.open("NEWTxm.dat",ios::app);
    NEWTxm << star -> pc << ";" 
           << x1 << ";" 
           << m1 << endl;
    
    NEWTxm.close();
    
}

int main(){

    Star star;

    cout << "Write step interval: ";
    cin  >> star.h;
    cout << "Write EOS constant: ";
    cin  >> star.E;
    
    double N;
    double pc1;
    double pc2;

    cout << "Write N: ";
    cin  >> N;
    cout << "Write starting pressure: ";
    cin  >> pc1;
    pc2 = pc1;
    
    if (N>1){
        cout << "Write end pressure: ";
        cin  >> pc2;
        for (int i=0; i<N; i++){
            star.pc = pc1 + i*(pc2-pc1)/(N-1);
            cout << star.pc << endl;
            TOV_Heun(&star);
            NEWT_Heun(&star);
        }
    } else{
        star.pc = pc1;
        cout << star.pc << endl;
        TOV_Heun(&star);
        NEWT_Heun(&star);
    }
}



