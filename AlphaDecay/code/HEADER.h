// HEADER.h
#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include <complex>
#define BARRIERCONST 0.0179    // e^2/(2pi*eps_0*R0*V0), R0 = 1.2fm and V0 = 134MeV
using namespace std;
using namespace Eigen;


// *********************** LinearAlgebra.cpp *******************

Matrix<complex<double>,Dynamic,Dynamic> zeros(int M, int N);
Matrix<complex<double>,Dynamic,1>       zeros(int M);
Matrix<complex<double>,Dynamic,1>       solve(MatrixXcd &A, VectorXcd &b);


// ********************* IO.cpp *********************

void write_vector(VectorXcd &X, ofstream &OutStream);
void write_vector(VectorXd  &X, ofstream &OutStream);   // use template?
void check_fail(ofstream &OutStream);


// ******************** ProblemSpecific.cpp *************************

struct Parameters{
    
    int N;                        // Number of gridpoints
    double dr;                    // Space interval
    int number_of_unknowns;
    double E;                     // Total energy
    double alpha;                 // Problem spec. constant
    Matrix<double,Dynamic,1> r;   // Gridpoints
    Matrix<double,Dynamic,1> V;   // Potential
};

Matrix<complex<double>,Dynamic,Dynamic> get_A(struct Parameters *);
Matrix<complex<double>,Dynamic,1>       get_b(struct Parameters *);
struct Parameters get_parameters();


// ***************** END OF HEADER ************************
