#include "HEADER.h"

struct Params get_parameters(){
        
    struct Params params;
    
    cout << "Write number of particles: ";
    cin  >> params.N;
    
    cout << "Write Delta: ";
    cin  >> params.Delta;

    cout << "Write iters: ";
    cin  >> params.iters;



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

double trial(VectorXd *X){
    
    double prod = 1;
    for (int j=0; j<(*X).rows(); j++)
        prod *= exp(-(*X)(j)*(*X)(j));
    return prod;
}

VectorXd init_function(struct Params *params){

    int    N = (params -> N);
    VectorXd R(N);
    
    for (int i=0; i<N; i++)
        R(i) = 0;

    return R;
}

void perturb(VectorXd *X, struct Params *params){

    int j = (rand() % (params -> N));
    double dx = (params -> Delta);
    dx *= 2*(rand() % 1000)/1000.0 - 1;
    (*X)(j) += dx;
}

void sample(struct Params *params, VectorXd *X, int j, int k){ 

    /* j is number of iteration. 
     * k is which instance of a single )(k) = (j-1)*(params -> Avparameter. */

    double Elocal = (*X)(0)*(*X)(0);
    Elocal *= (0.5-2*pow((params -> alphas)(k),2.0));
    Elocal += (params -> alphas)(k);

    (params -> VarElocals)(k) += pow((params -> AvgElocals)(k),2.0);

    (params -> AvgElocals)(k) = (j-1)*(params -> AvgElocals)(k);
    (params -> AvgElocals)(k) += Elocal;
    (params -> AvgElocals)(k) /= j;

    (params -> VarElocals)(k) = (j-1)*(params -> VarElocals)(k);
    (params -> VarElocals)(k) += Elocal*Elocal;
    (params -> VarElocals)(k) /= j;
    (params -> VarElocals)(k) -= pow((params -> AvgElocals)(k),2.0);
}

void gss_MonteCarlo(VectorXd *X, struct Params *params){

    /* a -> 0
     * b -> 3
     * c -> 1
     * d -> 2 */

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
    while (fabs(a-b)>1e-3 && k<100){

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
    }
}

void MonteCarlo(VectorXd *X, struct Params *params,int k){
    
    int iters = (params -> iters);
    int N = (params -> N);
    int totAcc=0;
    double p,random;
    
    /* Copy X */
    VectorXd oldX(N);
    oldX = copy(X);
    
    for (int i=1; i<iters; i++){
        perturb(X,params);
        p = trial(X)/trial(&oldX);
        p = p*p;
        if (p>=1){
            oldX = copy(X);
            totAcc++;
        } else {
            random = (rand() % 1000)/1000.0;
            if (p < random){
                oldX = copy(X); 
                totAcc++;
            } else
                *X = copy(&oldX);
        }
        sample(params,X,i,k);
    }
    cout << "Acceptence ratio: " << (double) totAcc/iters << endl;
}

