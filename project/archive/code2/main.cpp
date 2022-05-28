#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

// Star struct
struct Star{

    double h;   // Step size
    double b;   // End of bisection interval
    double pc;  // Central pressure

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

double EOS(double p, double b, int *FLAG){

    /* a = start of interval 
     * b = end of interval */

    double a=0;
    double n = bisection(P,a,b,FLAG,p);
    return E(n);
}


void TOV_Heun(Star *star){

    ofstream pmStream;
    pmStream.open("pm.dat");
    pmStream.precision(16);

    double max_iters = 1000000;
    double iters = 0;
    double x1,x2;
    double p1,p2;
    double m1,m2;
    double E1,E2;
    double h = star -> h;
    double b = star -> b;
    int FLAG=0;

    x1 = 1e-8;
    x2 = x1+h;
    p1 = star -> pc;
    E1 = EOS(p1,b,&FLAG);
    m1 = E1*x1*x1*x1/3;

    pmStream << x1 << ";"
             << p1 << ";"
             << m1 << ";"
             << E1 << "\n";
    //cout     << x1 << ";"
    //         << p1 << ";"
    //         << m1 << ";"
    //         << E1 << "\n";

    double f1,f2;
    double f3,f4;
    while (p1 > 0 && iters<max_iters){

        // Heun's method
        f1 = -(m1+x1*x1*x1*p1);
        f1*= (E1+p1)/(x1*(x1-2*m1));
        p2 = p1 + h*f1;
        E2 = EOS(p2,b,&FLAG);

        f2 = x1*x1*E1;
        m2 = m1 + h*f2;

        f3 = -(m2+x2*x2*x2*p2);
        f3*= (E2+p2)/(x2*(x2-2*m2));
        p2 = p1 + h*(f1+f3)/2;
        E2 = EOS(p2,b,&FLAG);
        
        f4 = x2*x2*E2;
        m2 = m1 + h*(f2+f4)/2;

        // Write to file
        pmStream << x2 << ";"
                 << p2 << ";"
                 << m2 << ";"
                 << E2 << "\n";

        if (FLAG == 1){
            break;
        }

        // Rename variables
        p1 = p2;
        m1 = m2;
        E1 = E2;
        x1 = x2;
        x2 = x2+h;
        iters++; 
    }

    pmStream.close();
    
    // Append last point (pc,x,m) to xm file
    fstream xmStream;
    xmStream.precision(16);
    xmStream.open("xm.dat",ios::app);
    xmStream << star -> pc << ";" 
             << x1 << ";" 
             << m1 << endl;
    
    xmStream.close();
    
}

int main(){

    Star star;
    cout << "Write step interval: ";
    cin  >> star.h;
    cout << "Write end for bisection interval: ";
    cin  >> star.b;
    
    double N;
    double pc_1;
    double pc_2;
    cout << "Write N: ";
    cin  >> N;
    cout << "Write pc1: ";
    cin  >> pc_1;
    cout << "Write pc2: ";
    cin  >> pc_2;
    for (int i=0; i<N-1; i++){
        star.pc = pc_1 + i*(pc_2-pc_1)/(N-1);
        cout << star.pc << endl;
        TOV_Heun(&star);
    }

}



















    
