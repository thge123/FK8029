#include "HEADER.h" 

struct Params get_parameters(){

    struct Params params;
    
    int N;  
    cout << "Write N: ";
    cin  >> N;
    params.N = N;
    if (N<10) {
        cout << "N too small. Aborting .. " << endl;
        exit(1);
    }

    cout << "Write h: ";
    cin  >> params.h;

    return params;

}
    

MatrixXd EvenMatrix(struct Params *params){

    double h = params -> h;
    int    N = params -> N;

    MatrixXd A = zeros(N,N);

    return A;
    

}
