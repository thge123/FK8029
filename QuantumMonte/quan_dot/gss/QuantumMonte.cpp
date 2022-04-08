#include "HEADER.h"

struct Params get_parameters(){
        
    struct Params params;
    
    cout << "Write number of particles: ";
    cin  >> params.N;
    
    cout << "Write Delta: ";
    cin  >> params.Delta;

    cout << "Write iters: ";
    cin  >> params.iters;

    params.alphas     = zeros(4);
    params.AvgElocals = zeros(4);
    params.VarElocals = zeros(4);


    double a;
    cout << "Write first alpha: ";
    cin  >> a;
    (params.alphas)(0) = a;

    double b;
    cout << "Write second alpha: ";
    cin  >> b;
    (params.alphas)(3) = b;

    return params;
}

double r1(Vector4d *X){
    
    double out = pow((*X)(0),2.0) + pow((*X)(1),2.0);
    return sqrt(out);
}

double r2(Vector4d *X){
    
    double out = pow((*X)(2),2.0) + pow((*X)(3),2.0);
    return sqrt(out);
}

double r12(Vector4d *X){
    
    // X = (x1,y1,x2,y2)

    double out  = pow((*X)(0)-(*X)(2),2.0);  // x1-x2
    out        += pow((*X)(1)-(*X)(3),2.0);  // y1-y2
    return sqrt(out);
}


double psi(Vector4d *X, struct Params *params, int j){
    
    double prod = 1;
    double alpha  = (params -> alphas)(j);
    double lambda = (params -> lambda);

    prod *= exp(-(r1(X)*r1(X)+r2(X)*r2(X))/2);
    prod *= exp(lambda*r12(X)/(1+alpha*r12(X)));

    return prod;
}

double Hpsi(Vector4d *X, struct Params *params, int j){

    double alpha  = (params -> alphas)(j);
    double lambda = (params -> lambda);
    double beta = lambda/(1+alpha*r12(X));
    double out=0;
    double Psi = psi(X,params,j);
    double R = r12(X);
    double x12 = (*X)(0)-(*X)(2);
    double y12 = (*X)(1)-(*X)(3);

    // **** Second derivative **** 
    // Square part
    out += pow(beta*(x12/R) -(*X)(0),2.0);
    out += pow(beta*(y12/R) -(*X)(1),2.0);
    out += pow(beta*(-x12/R)-(*X)(2),2.0);
    out += pow(beta*(-y12/R)-(*X)(3),2.0);

    out += -2*(1+2*lambda*alpha/pow(1+alpha*R,3.0));
    out += -(1/r1(X))*(r1(X)-lambda/pow(1+alpha*R,2.0));
    out += -(1/r2(X))*(r1(X)-lambda/pow(1+alpha*R,2.0));

    out *= -Psi/2;
    // **** End of second derivative ****

    // **** Potential part ****
    // Harmonic oscillator
    out += ((*X)(0)/2+(*X)(1)/2+(*X)(2)/2+(*X)(3)/2)*Psi;
    out += lambda/r12(X);

    return out;
}

double ELocal(Vector4d *X, struct Params *params, int j){

    double out = Hpsi(X,params,j);
    return out/psi(X,params,j);
}

void perturb(Vector4d *X, struct Params *params){

    int      j = (rand() % (params -> N));
    double dx1 = (params -> Delta);
    double dx2 = (params -> Delta);
    dx1 *= 2*(rand() % 1000)/1000.0 - 1;
    dx2 *= 2*(rand() % 1000)/1000.0 - 1;
    (*X)(2*j)   += dx1;
    (*X)(2*j+1) += dx2;
}

Vector4d init_config(struct Params *params){

    Vector4d X;
    X(0) = -1;
    X(1) =  0;
    X(2) =  1;
    X(3) =  0;
    return X;
}


void sample(Vector4d *X, struct Params *params, int j, int k){

    /* j is number of iteration. 
     * k is which instance of a single )(k) = (j-1)*(params -> Avparameter. */

    double Elocal = ELocal(X,params,k);

    (params -> VarElocals)(k) += pow((params -> AvgElocals)(k),2.0);

    (params -> AvgElocals)(k) = (j-1)*(params -> AvgElocals)(k);
    (params -> AvgElocals)(k) += Elocal;
    (params -> AvgElocals)(k) /= j;

    (params -> VarElocals)(k) = (j-1)*(params -> VarElocals)(k);
    (params -> VarElocals)(k) += Elocal*Elocal;
    (params -> VarElocals)(k) /= j;
    (params -> VarElocals)(k) -= pow((params -> AvgElocals)(k),2.0);

    cout << "Avg. E = " << (params -> AvgElocals)(k) << endl;
    cout << "Var. E = " << (params -> VarElocals)(k) << endl;
    cout << endl;
}

void gss_MonteCarlo(Vector4d *X, struct Params *params){

    /* a -> 0
     * b -> 3
     * c -> 1
     * d -> 2 */
    
    int filenumber;
    cout << "Write filenumber: ";
    cin  >> filenumber;
    struct Files files = get_files(filenumber);
    ofstream OutStreamA;
    OutStreamA.open(files.filenameA);
    OutStreamA.precision(16);

    double phi = 1.618033988749;
    double a = (params -> alphas)(0);
    double b = (params -> alphas)(3);
    double c = b-(b-a)/phi;
    double d = a+(b-a)/phi;
    (params -> alphas)(1) = c;
    (params -> alphas)(2) = d;

    MonteCarlo(X, params, 1);
    MonteCarlo(X, params, 2);

    double Ec = (params -> AvgElocals)(1);
    double Ed = (params -> AvgElocals)(2);
    
    int k = 0;
    while (fabs(a-b)>1e-6 && k<10000){

        if (Ec<Ed){
            b = d;
            (params -> alphas)(3) = d;
        } else {
            a = c;
            (params -> alphas)(0) = c;
        }

        c = b-(b-a)/phi;
        d = a+(b-a)/phi;
        (params -> alphas)(1) = c;
        (params -> alphas)(2) = d;

        MonteCarlo(X, params, 1);
        MonteCarlo(X, params, 2);
        
        Ec = (params -> AvgElocals)(1);
        Ed = (params -> AvgElocals)(2);
        k++;
        OutStreamA << c << ";" << Ec << ";" << (params -> VarElocals)(1) << endl;
        OutStreamA << d << ";" << Ed << ";" << (params -> VarElocals)(2) << endl;
    }

    OutStreamA.close();
}

void MonteCarlo(Vector4d *X, struct Params *params,int k){
    
    int iters = (params -> iters);
    int N = (params -> N);
    int totAcc=0;
    double p,random;

    (params -> AvgElocals)(k) = 0;
    (params -> VarElocals)(k) = 0;

    
    /* Copy X */
    Vector4d oldX;
    oldX = copy(X);
    cout << "Made it" << endl;

    int equilib = 5000;
    for (int i=1; i<iters; i++){
        perturb(X,params);
        p = psi(X,params,k)/psi(&oldX,params,k);
        p = p*p;
        if (p>=1){
            oldX = copy(X);
            totAcc++;
        } else {
            random = (rand() % 1000)/1000.0;
            if (random<p){
                oldX = copy(X); 
                totAcc++;
            } else
                *X = copy(&oldX);
        }
        if (i>equilib)
            sample(X,params,i-equilib,k);
    }
    cout << "Acceptence ratio: " << (double) totAcc/iters << endl;
}

