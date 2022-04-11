#include "HEADER.h"

struct Params get_parameters(){
        
    struct Params params;
    
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

double trial(VectorXd *X, struct Params *params, int j){
    
    double prod = 1;
    double a = (params -> alphas)(j);
    for (int j=0; j<(*X).rows(); j++)
        prod *= exp(-a*(*X)(j)*(*X)(j));
    return prod;
}

void perturb(VectorXd *X, struct Params *params){

    int j = (rand() % (params -> N));
    double dx = (params -> Delta);
    dx *= 2*(rand() % 1000)/1000.0 - 1;
    (*X)(j) += dx;
}

VectorXd init_function(struct Params *params){

    VectorXd X;
    X = zeros(params -> N);
    return X;
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

    cout << "Avg. E = " << (params -> AvgElocals)(k) << endl;
    cout << "Var. E = " << (params -> VarElocals)(k) << endl;
    cout << endl;
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
    }

}

void MonteCarlo(VectorXd *X, struct Params *params,int k){
    
    int filenumber;
    cout << "Write filenumber: ";
    cin  >> filenumber;
    struct Files files = get_files(filenumber);
    ofstream OutStreamA;
    OutStreamA.open(files.filenameA);
    OutStreamA.precision(16);

    int iters = (params -> iters);
    int N = (params -> N);
    int totAcc=0;
    double p,random;

    (params -> AvgElocals)(k) = 0;
    (params -> VarElocals)(k) = 0;
    
    /* Copy X */
    VectorXd oldX(N);
    oldX = copy(X);

    int equilib = 5000;
    for (int i=1; i<iters; i++){
        perturb(X,params);
        p = trial(X,params,k)/trial(&oldX,params,k);
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
            sample(params,X,i-equilib,k);
            write_vector(*X,OutStreamA);
        }
    }
//    cout << "Acceptence ratio: " << (double) totAcc/iters << endl;
}

