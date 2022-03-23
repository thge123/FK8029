#include "HEADER.h"

void write_vector(VectorXcd &X, ofstream &OutStream){

    int N = X.rows();
    for (int i=0; i<N; i++) OutStream << X(i) << ';';
    OutStream << endl;
}

void write_vector(VectorXd &X, ofstream &OutStream){

    int N = X.rows();
    for (int i=0; i<N; i++) OutStream << X(i) << ';';
    OutStream << endl;
}

void check_fail(ofstream &OutStream){

    if (OutStream.fail()){
        cout << "Failed to open output file" << endl;
        exit(1);
    }
}
