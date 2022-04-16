from PLOTTING import *
from numpy import array,zeros,exp,cos,pi
from numpy.linalg import solve
from scipy.interpolate import BSpline

def sigma(x):
    return 1.0
    alpha = 10
    return 4*alpha**3*exp(-2*alpha*x)/3

def numsol(N):

    K = 4                                    # order of B splines
    X = [j/N for j in range(N+1)]            # grid
    print(len(X))
    t = array(K*[X[0]] + X + K*[X[-1]])      # knots 
    n = len(t)-K-1                           # = number of B splines
    print(n)

    A = zeros((n,n))
    b = zeros(n)

    X = [0.0]
    for k in range(1,n-1):
        phi = (2*k-1)*pi/(2*(n-2))
        X.append(0.5+0.5*cos(phi))
    X.append(1.0)
    
    # Construct matrix eq.
    A[0,0]     = 1
    A[n-1,n-1] = 1
    c = zeros(n)
    for j in range(1,n-1):
        c[j] = 1
        spl  = BSpline(t,c,K,extrapolate=False)
        spl  = spl.derivative(2)
        for k in range(1,n-1):
            A[k,j] = spl(X[k])
        c[j] = 0
        b[j] = -3*sigma(X[j])*X[j]
    c = solve(A,b)
    return BSpline(t,c,K,extrapolate=False)

def main():
    N = 10
    K = 4
    X = [j/N for j in range(N+1)]     # grid
    t = 4*[X[0]] + X + 4*[X[-1]]      # knots 
    n = len(t)-K-1                    # = number of B splines

    spl = numsol(N)
    X = array([j/100 for j in range(101)])
    plt.plot(X,spl(X)+X)
    
    plt.show()

main()
        
    



    
