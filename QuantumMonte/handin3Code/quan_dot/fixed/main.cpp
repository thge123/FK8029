#include "HEADER.h"

int main(){
    
    srand(time(0));
    struct Params params = get_parameters();
    Vector4d X;
    X = init_config(&params);
    MonteCarlo(&X,&params,0);
}
