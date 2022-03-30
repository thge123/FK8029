// HEADER.h
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "Eigen/Dense"
#include <complex>
using namespace std;
using namespace Eigen;


// *********************** LinearAlgebra.cpp *******************

Matrix<complex<double>,Dynamic,Dynamic> zeros(int M, int N);
Matrix<complex<double>,Dynamic,1>       zeros(int M);
Matrix<complex<double>,Dynamic,1>       solve(MatrixXcd &A, VectorXcd &b);


// ********************* IO.cpp *********************

struct Files{

    char *filenameX;                   // filename for solution
    char *filenameV;                   // filename for potential
    char *filenameO;                   // filename for omegas in exp(omega*x)
    char *filenameR;                   // filename for gridpoints
};

void get_filenames(struct Files *);
void write_vector(VectorXcd &X, ofstream &OutStream);
void write_vector(VectorXd  &X, ofstream &OutStream);   // use template?
void check_fail(ofstream &OutStream);
char* get_filename(char name, int filenumber);
void update_filename(char *filename, int filenumber);
struct Files get_files(int);

// ******************** ProblemSpecific.cpp *************************

struct Parameters{
    
    int N;                             // Number of gridpoints
    int number_of_unknowns;
    double E;                          // Total energy
//    double alpha = 4e-4;               
    double V0    = 134;                // Potential in well
    double alpha = 0.1913;             // 2m'*(1fm)^2*1MeV/hbar^2
    int Z;                             // Atomic number
    int A;                             // Mass number
    int M;                             // Number of intervals after barrier
    Matrix<double,Dynamic,1> r;        // Gridpoints
    Matrix<double,Dynamic,1> V;        // Potential
    Matrix<double,Dynamic,1> omega;    // omegas in exp(omega*x)
};

Matrix<complex<double>,Dynamic,Dynamic> get_A(struct Parameters *);
Matrix<complex<double>,Dynamic,1>       get_b(struct Parameters *);
struct Parameters get_parameters();


// ***************** END OF HEADER ************************
