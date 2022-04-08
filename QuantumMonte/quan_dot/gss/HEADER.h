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
Vector4d copy(Vector4d *X);
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
void write_vector(Vector4d  &X, ofstream &OutStream);
void check_fail(ofstream &OutStream);
char* get_filename(char name, int filenumber);
void update_filename(char *filename, int filenumber);
struct Files get_files(int);


// ******************** ProblemSpecific.cpp *******************

struct Params{

    int    N;         // Number of coords.
    double Delta  = 0.1; 
    int    iters  = 100;
    double lambda = 1;
    Vector4d AvgElocals;
    Vector4d VarElocals;
    Vector4d alphas;

};

struct Params get_parameters();
void perturb(Vector4d*, struct Params*);
Vector4d init_config(struct Params*);
double r12(Vector4d*);
double r1(Vector4d*);
double r2(Vector4d*);
double psi(Vector4d*, struct Params*,int);
double Hpsi(Vector4d*, struct Params*,int);
double ELocal(Vector4d*, struct Params*,int);
void MonteCarlo(Vector4d*, struct Params*,int);
void sample (Vector4d*,struct Params*,int);
void gss_MonteCarlo(Vector4d*,struct Params*);


// ***************** END OF HEADER ************************
