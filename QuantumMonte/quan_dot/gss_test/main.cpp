#include <iostream>
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;

struct Test{

    Vector4d E;
    Vector4f alphas;
};

double E(double x){
    
    return -2/x + 2/(x*x);
}

struct Test get_test(){

    struct Test test;
    (test.alphas)(0) = 0;
    (test.alphas)(3) = 3;
    return test;
}

void DummyMonteCarlo(struct Test *test,int j){
    
    double x = (test -> alphas)(j);
    (test -> E)(j) = E(x);
}
    

int main(){

    double phi = 1.618033988749;

    struct Test test = get_test();
    double a = (test.alphas)(0);
    double b = (test.alphas)(3);
    double c = b-(b-a)/phi;
    double d = a+(b-a)/phi;
    (test.alphas)(1) = c;
    (test.alphas)(2) = d;

    DummyMonteCarlo(&test,1);
    DummyMonteCarlo(&test,2);
    
    
    double Ec = (test.E)(1);
    double Ed = (test.E)(2);
    
    int k = 0;
    while (fabs(a-b)>1e-6 && k<1000){

        if (Ec<Ed){
            b = d;
            (test.alphas)(3) = d;
        } else {
            a = c;
            (test.alphas)(0) = c;
        }

        c = b-(b-a)/phi;
        d = a+(b-a)/phi;
        (test.alphas)(1) = c;
        (test.alphas)(2) = d;

        DummyMonteCarlo(&test, 1);
        DummyMonteCarlo(&test, 2);
        
        Ec = (test.E)(1);
        Ed = (test.E)(2);
        k++;
    }
    cout << (a+b)/2 << endl;
}


