#include "HEADER.h"

struct Parameters get_parameters(){
    
    struct Parameters parameters;

    parameters.V0    = 134;
    parameters.alpha = 0.1913;

    int N;
    cout << "Write number of barrier segments N: ";
    cin  >> N;
    parameters.N = N;
    parameters.number_of_unknowns = 2*N+3;

    int M;    // Number of intervals after barrier - 1
    cout << "Write number of intervals after barrier: ";
    cin >> M; M-=1;
    parameters.M = M;
    parameters.number_of_unknowns += 2*M;

    double E;
    cout << "Write energy E in MeV: ";
    cin  >> E;
    parameters.E = E;

    int Z;
    cout << "Write atom number Z of parent nucleus: ";
    cin  >> Z;
    parameters.Z = Z;

    int A;
    cout << "Write mass number A of parent nucleus: ";
    cin  >> A;
    parameters.A = A;

    const double K = 2.879;              // K = e^2/(2*pi*eps0*1MeV*1fm)
    const double end_barrier = K*(Z-2)/E;   

    //const double start_barrier = 1.2*(pow(4.0,0.3333) + pow(A-4.0,0.3333));
    const double start_barrier = 1.5*(pow(4,0.3333)+pow(A-4.0,0.3333));

    /* Write grid points in vector r. 
     * The distance is normalized with
     * the length of the lowest potential
     * V0, hence r(1) = 1. 
     * ¦----¦               
     *      ¦----¦          barrier with N=3
     *           ¦----¦     last segnemt is <E
     *                ¦----¦ 
     * ¦----¦----¦----¦----¦-*/
    
    parameters.r.resize(N+2+M);
    parameters.V.resize(N+2+M);
    parameters.omega.resize(N+2+M);
    parameters.r(0) = 0;

    double dr = (end_barrier-start_barrier)/N;
    for (int j=0; j<N+1+M; j++)
        // grid for barrier segments
        parameters.r(j+1) = j*dr + start_barrier;

    /* Write values of the potential in
     * vector V. The values are normalized
     * such that V(0) = -1. */
    parameters.V(0) = -parameters.V0 + (Z-2)*K/parameters.r(1);

    for (int i=0; i<N+1+M; i++) {
        parameters.V(i+1) = (Z-2)*K/parameters.r(i+1);
        parameters.omega(i) = 0;    // initializing
    }

    double d = parameters.V(N+M) - parameters.V(N+M+1);
    for (int i=0; i<N+1+M; i++) 
        /* Shift the potential slightly so that 
         * E not equal to V anywhere. */
        parameters.V(i+1) -= d/2;

    parameters.V(N+M+1) = 0;

    return parameters;
}

MatrixXcd get_A(struct Parameters *paramsPTR){
    
    /* Use signature -++++++------- to see when to
     * use trigonometric functions and when to use 
     * exponential functions, i.e. - for E>V and +
     * for E<V. - gives trigs and + gives exps. The 
     * vector we want to solve is 
     *      ¦   A0   ¦   ---> This one we set = 1
     *      ¦   B0   ¦
     *      ¦   A1   ¦
     *      ¦   B1   ¦ 
     *      ¦   A2   ¦
     *      ¦   B2   ¦ 
     *      ¦   ..   ¦
     *      ¦ A(N-1) ¦ 
     *      ¦ B(N-1) ¦   ---> This one we don't need
     * which are 2N-2 unknowns. */

    int    Q     = paramsPTR -> number_of_unknowns;
    int    N     = paramsPTR -> N;
    int    M     = paramsPTR -> M;
    double alpha = paramsPTR -> alpha;
    double E     = paramsPTR -> E;

    
    /* Construction of the matrix A. It is 
     * (2N-1)x(2N-1) where N is grid points. */
    Matrix<complex<double>,Dynamic,Dynamic> A;
    A = zeros(Q,Q);

    // First equation
    A(0,0) = 1;
    
    double V = paramsPTR -> V(0);    
    if (V>=E){
        cout << "The energy must be higher than V0. Terminating ..." << endl;
        exit(1);
    }

    double r = paramsPTR -> r(1);    // =  1
    double omega = sqrt(alpha*(E-V));
    paramsPTR -> omega(0) = omega;

    // Continuity at first bndry.
    A(1,0) = +exp(+1i*omega*r);
    A(1,1) = +exp(-1i*omega*r);

    // Continuity of deriv. at first bndry.
    A(2,0) = +1i*omega*exp(+1i*omega*r);
    A(2,1) = -1i*omega*exp(-1i*omega*r);

    int P = N+3+M;    // Number of gridpoints
    for (int i=1; i<P-2; i++){
        
        V = paramsPTR -> V(i);
        r = paramsPTR -> r(i);

        if (V>E) {
            // Exponential functions

            omega = sqrt(alpha*(V-E));
            paramsPTR -> omega(i) = omega;

            // Continuity of function to the left
            A(2*i-1,2*i+0) = -exp(+omega*r);
            A(2*i-1,2*i+1) = -exp(-omega*r);

            // Continuity of deriv. to the left
            A(2*i-0,2*i+0) = -omega*exp(+omega*r);
            A(2*i-0,2*i+1) = +omega*exp(-omega*r);

            r = paramsPTR -> r(i+1);

            // Continuity of function to the right
            A(2*i+1,2*i+0) = +exp(+omega*r);
            A(2*i+1,2*i+1) = +exp(-omega*r);

            // Continuity of deriv to the right
            A(2*i+2,2*i+0) = +omega*exp(+omega*r);
            A(2*i+2,2*i+1) = -omega*exp(-omega*r);
            
        } else {
        
        if (V<E) {
            // Trigonometric functions

            omega = sqrt(alpha*(E-V));
            paramsPTR -> omega(i) = omega;

            // Continuity of function to the left
            A(2*i-1,2*i+0) = -exp(+1i*omega*r);
            A(2*i-1,2*i+1) = -exp(-1i*omega*r);

            // Continuity of deriv. to the left
            A(2*i-0,2*i+0) = -1i*omega*exp(+1i*omega*r);
            A(2*i-0,2*i+1) = +1i*omega*exp(-1i*omega*r);

            r = paramsPTR -> r(i+1);

            // Continuity of function to the right
            A(2*i+1,2*i+0) = +exp(+1i*omega*r);
            A(2*i+1,2*i+1) = +exp(-1i*omega*r);

            // Continuity of deriv. to the right
            A(2*i+2,2*i+0) = +1i*omega*exp(+1i*omega*r);
            A(2*i+2,2*i+1) = -1i*omega*exp(-1i*omega*r);

        } else { 
            /* Linear function. Haven't taken this
             * into account. */ 
            cout << "E=V at some point. Terminating..." << endl;    
            exit(1);
            }
        }
    }

    // The last bndry.
    r = paramsPTR -> r(P-2);
    V = paramsPTR -> V(P-2);

    if (V>E) {
        // Exponential function

        omega = sqrt(alpha*(V-E));
        paramsPTR -> omega(P-2) = omega;

        // Continuity of function to the left
        A(Q-2,Q-1) = -exp(-omega*r);

        // Continuity of deriv. to the left 
        A(Q-1,Q-1) = +omega*exp(-omega*r);
    } else {

    if (V<E) {
        // Trigonometric function

        omega = sqrt(alpha*(E-V));
        paramsPTR -> omega(P-2) = omega;

        // Continuity of function to the left
        A(Q-2,Q-1) = -exp(1i*omega*r);

        // Continuity of deriv. to the left 
        A(Q-1,Q-1) = -1i*omega*exp(1i*omega*r);
    
    } else { 
        /* Linear function. Haven't taken this
         * into account. */ 
        cout << "E=V at some point. Terminating..." << endl;    
        exit(1);
        }
    }
    
    return A;
}

VectorXcd get_b(struct Parameters *paramsPTR){
    
    VectorXcd b;
    int N = paramsPTR -> number_of_unknowns;
    b = zeros(N);
    b(0) = 1;    // b = (1,0,0,...)
    return b;
}
