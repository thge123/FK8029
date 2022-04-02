// HEADER.h
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "Eigen/Sparse"
#include <complex>
#include <vector>
using namespace std;
using namespace Eigen;
typedef Triplet<double> triplet;


// *********************** LinearAlgebra.cpp *******************

Matrix<double,Dynamic,Dynamic> zeros(int M, int N);
Matrix<double,Dynamic,1>       zeros(int M);
Matrix<complex<double>,Dynamic,1>       solve(MatrixXcd &A, VectorXcd &b);

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

    double h;      // step size
    int    N;      // 
    double E;
    double err;
    double alpha=10;

};

struct Params get_parameters();
SparseMatrix<double> EvenMatrix(struct Params *);
SparseMatrix<double> OddMatrix(struct Params *);




// ***************** END OF HEADER ************************
