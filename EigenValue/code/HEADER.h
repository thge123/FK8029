// HEADER.h
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "Eigen/Dense"
#include <complex>
using namespace std;
using namespace Eigen;


// *********************** LinearAlgebra.cpp *******************

Matrix<double,Dynamic,Dynamic> zeros(int M, int N);
Matrix<double,Dynamic,1>       zeros(int M);
Matrix<complex<double>,Dynamic,1>       solve(MatrixXcd &A, VectorXcd &b);

void PowerIter(MatrixXd&,VectorXd&);



// ********************* IO.cpp *********************

struct Files{

    char *filename1; 
    char *filename2; 
    char *filename3;
    char *filename4; 
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
    double L;      // End of interval

};

struct Params get_parameters();
MatrixXd EvenMatrix(struct Params *);






// ***************** END OF HEADER ************************
