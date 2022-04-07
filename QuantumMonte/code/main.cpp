#include "HEADER.h"

int main(){
    
    srand(time(0));
    struct Params params = get_parameters();
    VectorXd X;
    X = init_function(&params);
    gss_MonteCarlo(&X,&params);
    cout << "a = " << (params.alphas)(0) << endl;
    cout << "Ea = " << (params.AvgElocals)(0) << endl;
    cout << "b = " << (params.alphas)(3) << endl;
    cout << "Eb = " << (params.AvgElocals)(3) << endl;
    
}
