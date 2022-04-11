#include "HEADER.h"

int main(){
    
    srand(time(0));
    struct Params params = get_parameters();
    VectorXd X;
    X = init_function(&params);
    gss_MonteCarlo(&X,&params);
}
