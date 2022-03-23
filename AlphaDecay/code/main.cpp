#include "HEADER.h"

int main(){
    
    // File to write to
    ofstream OUTstreamX,OUTstreamV,OUTstreamR;
    OUTstreamX.open("X.dat");
    OUTstreamV.open("V.dat");
    OUTstreamR.open("r.dat");
    check_fail(OUTstreamX);
    check_fail(OUTstreamV);
    check_fail(OUTstreamR);
    
    // Problem spceific parameters
    Parameters params;
    params = get_parameters();

    // Solution of matrix problem
    Matrix<complex<double>,Dynamic,Dynamic> A = get_A(&params);
    Matrix<complex<double>,Dynamic,1>       b = get_b(&params);
    Matrix<complex<double>,Dynamic,1>       x = solve(A,b);

    // Write vector to file
    write_vector(x,OUTstreamX);
    write_vector(params.V,OUTstreamV);
    write_vector(params.r,OUTstreamR);
    OUTstreamX.close();
    OUTstreamV.close();
    OUTstreamR.close();
    
}

