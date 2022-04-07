// HEADER.h
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "Eigen/Sparse"
#include <time.h>
using namespace std;
using namespace Eigen;
typedef Triplet<double> triplet;


// *********************** LinearAlgebra.cpp *******************

Matrix<double,Dynamic,Dynamic> zeros(int M, int N);
Matrix<double,Dynamic,1>       zeros(int M);
Matrix<complex<double>,Dynamic,1>       solve(MatrixXcd &A, VectorXcd &b);
VectorXd copy(VectorXd *X);
VectorXd invPowerIter(SparseMatrix<double>&,double shift,double *lam,double *err);



// ********************* IO.cpp *********************

struct Files{

    char *filenameA; 
    char *filenameB; 
    char *filenameC;
    char *filenameD; 
};

void get_filenames(struct Files *);
void write_vector(VectorXcd &X, ofstream &OutStream);
void write_vector(VectorXd  &X, ofstream &OutStream);
void check_fail(ofstream &OutStream);
char* get_filename(char name, int filenumber);
void update_filename(char *filename, int filenumber);
struct Files get_files(int);


// ******************** ProblemSpecific.cpp *******************

struct Params{

    int    N;         // Number of coords.
    double Delta  = 0.1; 
    int    iters  = 100;
    Vector4d AvgElocals;
    Vector4d VarElocals;
    Vector4d alphas;

};

struct Params get_parameters();
VectorXd init_function(struct Params*);
void perturb(VectorXd *, struct Params *);
void MonteCarlo(VectorXd *, struct Params*,int);
void sample (VectorXd*, struct Params*,int);
void gss_MonteCarlo(VectorXd *X, struct Params *params);


// ***************** END OF HEADER ************************
