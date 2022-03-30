#include <iostream>
#include "HEADER.h"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

int main() {

    struct Params params = get_parameters();
    
    MatrixXd A  = EvenMatrix(&params);

    cout << A << endl;

}
