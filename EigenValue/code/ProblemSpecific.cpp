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

    return params;

}

SparseMatrix<double> EvenMatrix(struct Params *params){


    double h = params -> h;
    int    N = params -> N;

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
    tripletList.push_back(triplet(N-2,N-3,  -4));
    tripletList.push_back(triplet(N-2,N-2,  -15.0));

    // Potential
    double x;
    for (int i=0; i<N-1; i++){

        x = i*h;
        tripletList.push_back(triplet(i,i,  - x*x*(12*h*h)));

    }

    SparseMatrix<double> A(N-1,N-1);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    A = -(1/(12*h*h)) * A;
    

    return A;
    

}

SparseMatrix<double> OddMatrix(struct Params *params){


    double h = params -> h;
    int    N = params -> N;

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
        tripletList.push_back(triplet(i,i,  - x*x*(12*h*h)));

    }

    SparseMatrix<double> A(N-1,N-1);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    A = -(1/(12*h*h)) * A;
    

    return A;
    

}
