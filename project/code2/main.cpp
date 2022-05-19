#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double EOS(double p){

    return 1e-3;
}

double bissection(double (*f)(double x),double a, double b){

    int max_iters = 10000;
    double TOL = 1e-6;
    int iters = 0;
    double c;
    while (iters<max_iters){
        c = (a+b)/2;
        if (f(c)*f(c) < TOL*TOL) 
            return c;
        else{
            if (f(c)*f(a) < 0)
                b = c;
            else
                a = c;
        }
    }
    cout << "Bissection error" << sqrt(f(c)*f(c)) << endl;
    return c;
}

void TOV_Heun(){

    ofstream pmStream;
    pmStream.open("pm.dat");

    double max_iters = 10000;
    double iters = 0;
    double h;
    double x1,x2;
    double p1,p2;
    double m1,m2;
    double E1,E2;

    cout << "p0: ";
    cin  >> p1;
    cout << "h: ";
    cin  >> h;
    x1 = h;
    x2 = x1+h;
    E1 = EOS(p1);

    pmStream << x1 << ";"
             << p1 << ";"
             << m1 << "\n";


    double f1,f2;
    double f3,f4;
    while (p1 > 0 && iters<max_iters){

        // Heun's method
        f1 = -(m1+x1*x1*x1*p1);
        f1*= (E1+p1)/(x1*(x1-2*m1));
        p2 = p1 + h*f1;
        E2 = EOS(p2);

        f2 = x1*x1*E1;
        m2 = m1 + h*f2;

        f3 = -(m2+x2*x2*x2*p2);
        f3*= (E2+p2)/(x2*(x2-2*m2));
        p2 = p1 + h*(f1+f3)/2;
        E2 = EOS(p2);
        
        f4 = x2*x2*E2;
        m2 = m1 + h*(f2+f4)/2;

        // Write to file
        pmStream << x2 << ";"
                 << p2 << ";"
                 << m2 << "\n";

        // Rename variables
        p1 = p2;
        m1 = m2;
        E1 = E2;
        x1 = x2;
        x2 = x2+h;
        iters++; 
    }

    pmStream.close();
}
        

int main(){

    TOV_Heun();


}

    

    


    

    
