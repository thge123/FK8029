#include "HEADER.h"

int main(){
    
    // Problem spceific parameters
    Parameters params;
    params = get_parameters();

    // File to write to
    ofstream OUTstreamX,OUTstreamV,OUTstreamR, OUTstreamO;
    struct Files files = get_files(params.N);
    OUTstreamX.open(files.filenameX); check_fail(OUTstreamX);
    OUTstreamV.open(files.filenameV); check_fail(OUTstreamV);
    OUTstreamR.open(files.filenameR); check_fail(OUTstreamR);
    OUTstreamO.open(files.filenameO); check_fail(OUTstreamO);
    
    // Solution of matrix problem
    Matrix<complex<double>,Dynamic,Dynamic> A = get_A(&params);
    Matrix<complex<double>,Dynamic,1>       b = get_b(&params);
    Matrix<complex<double>,Dynamic,1>       x = solve(A,b);

    // Write to files
    write_vector(x,OUTstreamX);
    write_vector(params.V,OUTstreamV);
    write_vector(params.r,OUTstreamR);
    write_vector(params.omega,OUTstreamO);
    OUTstreamX.close();
    OUTstreamV.close();
    OUTstreamR.close();
    OUTstreamO.close();
    
}

