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

    int    N      = 2;     // Number of particles 
    double Delta  = 0.1;   // Maximum displacement
    int    iters  = 100;   // MC iterations
    double lambda = 1;     // Coulomb coupling
    Vector4d AvgElocals;   // Samples of avg. local energies
    Vector4d VarElocals;   // Samples of var. local energies
    Vector4d alphas;       // Trial function params.

};

struct Params get_parameters();
void perturb(Vector4d*, struct Params*);
Vector4d init_config(struct Params*);
double r12(Vector4d*);
double r1(Vector4d*);
double r2(Vector4d*);
double psi(Vector4d*, struct Params*,int);
double Hpsi(Vector4d*, struct Params*,int);
Vector4d dX(double h, int j);
double ELocal(Vector4d*, struct Params*,int);
double ELocal2(Vector4d*, struct Params*,int);
void MonteCarlo(Vector4d*, struct Params*,int);
void sample (Vector4d*,struct Params*,int);
void gss_MonteCarlo(Vector4d*,struct Params*);


// ***************** END OF HEADER ************************
