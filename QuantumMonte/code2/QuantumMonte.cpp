#include "HEADER.h"

struct Params get_parameters(){
        
    struct Params params;
    
    cout << "Write number of particles: ";
    cin  >> params.N;
    
    cout << "Write Delta: ";
    cin  >> params.Delta;

    cout << "Write iters: ";
    cin  >> params.iters;

    cout << "Write first alpha: ";
    cin  >> params.a;

    cout << "Write second alpha: ";
    cin  >> params.b;

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

    cout << "Avg. local energy: " << (params -> AvgElocals)(0) << endl;
    cout << "Var. local energy: " << (params -> VarElocals)(0) << endl;
    cout << endl;
    
}

void gcs(struct Params *params, int k){

    double phi = 1.618033988749;
}

    
void MonteCarlo(VectorXd *X, struct Params *params){
    
    int iters = (params -> iters);
    int N = (params -> N);
    double p,random;
    
    int totAcc=0;
    /* Copy X */
    VectorXd oldX(N);
    oldX = copy(X);
    
    for (int k=1; k<(params -> alphas).rows(); k++){

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
    gcs(params);
    }
}

