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
    if (params.h<=0) {
        cout << "h not valid. Aborting .. " << endl;
        exit(1);
    }
    
    cout << "Write alpha factor: ";
    cin  >> params.alpha;

    return params;

}

SparseMatrix<double> EvenMatrix(struct Params *params){


    double h     = params -> h;
    int    N     = params -> N;
    double alpha = params -> alpha;

    vector<triplet> tripletList;

    tripletList.push_back(triplet(0,0,  -80.0/3 - 30.0));
    tripletList.push_back(triplet(0,1,   72.0));
    tripletList.push_back(triplet(0,2,  -18.0));
    tripletList.push_back(triplet(0,3,   8.0/3));
    
    tripletList.push_back(triplet(1,0,   10.0/3+16.0));
    tripletList.push_back(triplet(1,1,  -36.0));
    tripletList.push_back(triplet(1,2,   18.0));
    tripletList.push_back(triplet(1,3,  -4.0/3));


    for (int i=2; i<N-3; i++){

        tripletList.push_back(triplet(i,i-2,  -1.0));
        tripletList.push_back(triplet(i,i-1,  16.0));
        tripletList.push_back(triplet(i,i,   -30.0));
        tripletList.push_back(triplet(i,i+1,  16.0));
        tripletList.push_back(triplet(i,i+2,  -1.0));

    }

    tripletList.push_back(triplet(N-3,N-5,  -1.0));
    tripletList.push_back(triplet(N-3,N-4,   16.0));
    tripletList.push_back(triplet(N-3,N-3,  -30.0));
    tripletList.push_back(triplet(N-3,N-2,   16.0));
    
    tripletList.push_back(triplet(N-2,N-6,   1.0));
    tripletList.push_back(triplet(N-2,N-5,  -6.0));
    tripletList.push_back(triplet(N-2,N-4,   14.0));
    tripletList.push_back(triplet(N-2,N-3,  -4.0));
    tripletList.push_back(triplet(N-2,N-2,  -15.0));

    // Potential
    double x;
    for (int i=0; i<N-1; i++){

        x = i*h;
        if (alpha>0.5)
            tripletList.push_back(triplet(i,i,  - (x*x+alpha*exp(-x*x) - (1+log(alpha)))*(12*h*h)));    // modif. har. osc.
        else
            //tripletList.push_back(triplet(i,i,  + alpha*(12*h*h)));    // sq. well
            tripletList.push_back(triplet(i,i,  - x*x*(12*h*h)));    // har osc. 
            //tripletList.push_back(triplet(i,i,  - x*(12*h*h)));    // |x| pot
            //tripletList.push_back(triplet(i,i,  - pow(x,8.0)*(12*h*h)));    // x^(2n) potential

    }

    SparseMatrix<double> A(N-1,N-1);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    A = -(1/(12*h*h)) * A;
    

    return A;
    

}

SparseMatrix<double> OddMatrix(struct Params *params){


    double h = params -> h;
    int    N = params -> N;
    double alpha = params -> alpha;

    vector<triplet> tripletList;

    tripletList.push_back(triplet(0,0,  -15.0));
    tripletList.push_back(triplet(0,1,  -4.0));
    tripletList.push_back(triplet(0,2,   14.0));
    tripletList.push_back(triplet(0,3,  -6.0));
    tripletList.push_back(triplet(0,4,   1.0));
    
    tripletList.push_back(triplet(1,0,   16.0));
    tripletList.push_back(triplet(1,1,  -30.0));
    tripletList.push_back(triplet(1,2,   16.0));
    tripletList.push_back(triplet(1,3,  -1.0));

    for (int i=2; i<N-3; i++){

        tripletList.push_back(triplet(i,i-2,  -1.0));
        tripletList.push_back(triplet(i,i-1,  16.0));
        tripletList.push_back(triplet(i,i,   -30.0));
        tripletList.push_back(triplet(i,i+1,  16.0));
        tripletList.push_back(triplet(i,i+2,  -1.0));

    }

    tripletList.push_back(triplet(N-3,N-5,  -1.0));
    tripletList.push_back(triplet(N-3,N-4,   16.0));
    tripletList.push_back(triplet(N-3,N-3,  -30.0));
    tripletList.push_back(triplet(N-3,N-2,   16.0));
    
    tripletList.push_back(triplet(N-2,N-6,   1.0));
    tripletList.push_back(triplet(N-2,N-5,  -6.0));
    tripletList.push_back(triplet(N-2,N-4,   14.0));
    tripletList.push_back(triplet(N-2,N-3,  -4.0));
    tripletList.push_back(triplet(N-2,N-2,  -15.0));

    // Potential
    double x;
    for (int i=0; i<N-1; i++){

        x = (i+1)*h;
        // tripletList.push_back(triplet(i,i,  - (x*x+alpha*(exp(-x*x)-1))*(12*h*h)));    // harmonic_osc
        if (alpha>0.5)
            tripletList.push_back(triplet(i,i,  - (x*x+alpha*exp(-x*x) - (1+log(alpha)))*(12*h*h)));    // modif. har. osc.
        else
            //tripletList.push_back(triplet(i,i,  + alpha*(12*h*h)));    // sq. well
            //tripletList.push_back(triplet(i,i,  - x*(12*h*h)));    // |x| pot
            tripletList.push_back(triplet(i,i,  - x*x*(12*h*h)));    // har osc. 
            //tripletList.push_back(triplet(i,i,  - pow(x,8.0)*(12*h*h)));    // har osc. 

    }

    SparseMatrix<double> A(N-1,N-1);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    A = -(1/(12*h*h)) * A;
    

    return A;
    

}
