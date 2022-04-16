#include "HEADER.h"

struct Params get_parameters(){
        
    struct Params params;
    
    cout << "Write Delta: ";
    cin  >> params.Delta;

    cout << "Write iters: ";
    cin  >> params.iters;

    cout << "Write lambda: ";
    cin  >> params.lambda;
    

    for (int i=0; i<4; i++){
        params.alphas(i) = 0.0;
        params.AvgElocals(i) = 0.0;
        params.VarElocals(i) = 0.0;
    }

    cout << "Write first alpha: ";
    cin  >> (params.alphas)(0);

    cout << "Write second alpha: ";
    cin  >> (params.alphas)(3);

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
    
    double alpha  = (params -> alphas)(j);
    double lambda = (params -> lambda);
    double R1  = r1(X);
    double R2  = r2(X);
    double R12 = r12(X);

    double prod = 1.0;
    prod *= exp(-(R1*R1+R2*R2)/2.0);
    prod *= exp(lambda*R12/(1.0+alpha*R12));

    return prod;
}

double Hpsi(Vector4d *X, struct Params *params, int j){

    double alpha  = (params -> alphas)(j);
    double lambda = (params -> lambda);
    double R = r12(X);                
    double gamma = 1+alpha*R;
    double beta = lambda/gamma;
    double Psi = psi(X,params,j);
    double x1  = (*X)(0);
    double y1  = (*X)(1);
    double x2  = (*X)(2);
    double y2  = (*X)(3);
    double x12 = x1-x2;
    double y12 = y1-y2;

    double out=0;
    // **** Second derivative **** 
    // Square part
    out += pow(beta*(x12/R) - x1,2.0);
    out += pow(beta*(y12/R) - y1,2.0);
    out += pow(beta*(-x12/R)- x2,2.0);
    out += pow(beta*(-y12/R)- y2,2.0);

    out +=  2*lambda/(gamma*gamma*R); 
    out += -4*lambda*alpha/(gamma*gamma*gamma);
    out += -2;
    out += -2;

    out *= -Psi/2;
    // **** End of second derivative ****

    // **** Potential part ****
    // Harmonic oscillator
    out += (x1*x1/2+y1*y1/2+x2*x2/2+y2*y2/2)*Psi;
    out += lambda*Psi/R;

    return out;
}

double ELocal(Vector4d *X, struct Params *params, int j){

    double out = Hpsi(X,params,j);
    return out/psi(X,params,j);
}

Vector4d dX(double h, int j){
    
    Vector4d X = zeros(4);
    X(j) = h;
    return X;
}

double ELocal2(Vector4d *X, struct Params *params, int j){
    double alpha  = (params -> alphas)(j);
    double lambda = (params -> lambda);
    double R = r12(X);
    double Psi = psi(X,params,j);
    double x1  = (*X)(0);
    double y1  = (*X)(1);
    double x2  = (*X)(2);
    double y2  = (*X)(3);

    double h=0.00001;

    Vector4d Dummy = *X + dX(-h,0);
    double f2x1 = psi(&Dummy,params,j);
    f2x1 += -2*Psi;
    Dummy = *X + dX(h,0);
    f2x1 += psi(&Dummy,params,j);
    f2x1 /= h*h;

    Dummy = *X+dX(-h,1);
    double f2y1 = psi(&Dummy,params,j);
    f2y1 += -2*Psi;
    Dummy = *X+dX(h,1);
    f2y1 += psi(&Dummy,params,j);
    f2y1 /= h*h;
    
    Dummy = *X+dX(-h,2);
    double f2x2 = psi(&Dummy,params,j);
    f2x2 += -2*Psi;
    Dummy = *X+dX(h,2);
    f2x2 += psi(&Dummy,params,j);
    f2x2 /= h*h;

    Dummy = *X+dX(-h,3);
    double f2y2 = psi(&Dummy,params,j);
    f2y2 += -2*Psi;
    Dummy = *X+dX(h,3);
    f2y2 += psi(&Dummy,params,j);
    f2y2 /= h*h;

    double out = -0.5*(f2x1+f2y1+f2x2+f2y2)/Psi;
    out += (x1*x1/2+y1*y1/2+x2*x2/2+y2*y2/2);
    out += lambda/R;

    return out;
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

    double Elocal = ELocal2(X,params,k);

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
    ofstream OutStreamB;
    OutStreamB.open(files.filenameB);

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
            OutStreamB << d << ";" << Ed << ";" << (params -> VarElocals)(2) << endl;
        } else {
            a = c;
            (params -> alphas)(0) = c;
            OutStreamB << c << ";" << Ec << ";" << (params -> VarElocals)(1) << endl;
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
    }

    OutStreamB.close();
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
        if (i>equilib){
            sample(X,params,i-equilib,k);
        }
    }
    cout << "Acceptence ratio: " << (double) totAcc/iters << endl;
}

