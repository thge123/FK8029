#include "HEADER.h"

using namespace std;
using namespace Eigen;

int main() {

    /****** File handling ******/
    int filenumber;
    cout << "Write filenumber: ";
    cin  >> filenumber;
    struct Files files   = get_files(filenumber);

    ofstream OutStreamA;    // Eigenvects (have to create dir A)
    ofstream OutStreamB;    // Energy and other things (have to create dir B)
    OutStreamA.open(files.filenameA); check_fail(OutStreamA);
    OutStreamB.open(files.filenameB); check_fail(OutStreamB);

    /****** Problem ******/

    struct Params params = get_parameters();
    
    SparseMatrix<double> A  = EvenMatrix(&params);
    SparseMatrix<double> B  = OddMatrix(&params);

    VectorXd x0(params.N-1);    // Eigenvects
    VectorXd E(4);              // Energy and other things
    E(1) = params.N;            // Number of steps
    E(2) = params.h;            // Step size
    
    double shift = 0;
    for (int i=0; i<20; i++){
    
        cout << "Shift: " << shift << endl;
        if (i%2==0){
            x0 = invPowerIter(A,shift,&(params.E),&(params.err));
            cout << "Eigenvalue found: " << params.E 
                 << "  with error: " << params.err << endl;
            //shift = params.E+0.1;
            shift = params.E + rand()%1000/1000.0;
        }
        if (i%2==1){
            x0 = invPowerIter(B,shift,&(params.E),&(params.err));   
            cout << "Eigenvalue found: " << params.E 
                 << "  with error: " << params.err << endl;
            //shift = params.E+0.1;
            shift = params.E + rand()%1000/1000.0;
        }
        E(0) = params.E;    // Energy
        E(3) = params.err;  // Error in eigenvalue eq.
        write_vector(x0, OutStreamA);
        write_vector(E , OutStreamB);
            
    }


}
