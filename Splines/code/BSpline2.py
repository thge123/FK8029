from Distributions import *
from PLOTTING import *
from numpy import array,zeros,exp,cos,pi,matmul
from numpy.linalg import solve
from scipy.interpolate import BSpline
from numpy import max as MAX
from numpy import abs as ABS

def numsol(N,i,r):

    K = 4                                    # order of B splines
    X = [0.0]
    for k in range(2,N-1):
        X.append(k/(N-1))
    #for k in range(1,N-1):
    #    phi = (2*k-1)*pi/(2*(N-2))
    #    X.append(0.5+0.5*cos(phi))
    X.append(1.0)
    X.sort()
    t = array(K*[X[0]] + X + K*[X[-1]])      # knots 
    print("Knots: ", t)
    n = len(t)-K-1                           # = number of B splines
    print("Number of splines: ", n)

    X = array([j/1000 for j in range(1001)])
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    c = zeros(n)
    for k in range(n):
        c[k] = 1
        ax.plot(X,BSpline(t,c,K,extrapolate=False)(X))
        c = zeros(n)

    A = zeros((n,n))
    b = zeros(n)

    X = [0.0]
    X.append(0.05)
    for k in range(2,n-1):
        X.append(k/(n-1))
        #phi = (2*k-1)*pi/(2*(n-2))
        #X.append(0.5+0.5*cos(phi))
    X.append(1.0)
    X.sort()
    print('Collocation points: ', X)

    sigma = [sigma1,sigma2,sigma3][i]
    
    # Construct matrix eq.
    A[0,0]     = 1.0
    A[n-1,n-1] = 1.0
    c = zeros(n)
    for j in range(1,n-1):
        c[j] = 1
        spl  = BSpline(t,c,K,extrapolate=False)
        spl  = spl.derivative(2)
        for k in range(1,n-1):
            A[k,j] = spl(X[k])
        c[j] = 0.0
        b[j] = -sigma(X[j],r)*X[j]
    c = solve(A,b)
    print("b: ", b)
    print("Ac: ", matmul(A,c))
    print("c: ", c)

    X = array([j/1000 for j in range(1001)])
    ax.plot(X,BSpline(t,c,K,extrapolate=False)(X))
    return BSpline(t,c,K,extrapolate=False)
