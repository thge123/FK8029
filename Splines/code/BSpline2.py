from Distributions import *
from PLOTTING import *
from numpy import array,zeros,exp,cos,pi,matmul
from numpy.linalg import solve
from numpy.random import rand
from scipy.interpolate import BSpline
from numpy import max as MAX
from numpy import abs as ABS

def numsol(N,i,r):

    K = 3                                    # order of B splines
    X = [0.0]
    for k in range(1,N-1):
        #phi = (2*k-1)*pi/(2*(N-2))
        #X.append(0.5+0.5*cos(phi))
        X.append(k/(N-1))
    X.append(1.0)
    X.sort()
    t = array(K*[X[0]] + X + K*[X[-1]])      # knots 
    n = len(t)-K-1                           # = number of B splines
    print("Knots: ", t)
    print("Number of splines: ", n)

    X = array([j/1000 for j in range(1001)])
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    c = zeros(n)
    for k in range(n):
        c[k] = 1
        ax.plot(X,BSpline(t,c,K,extrapolate=False)(X),lw=0.75)
        c = zeros(n)

    A = zeros((n,n))
    b = zeros(n)

    X = [0.0]
    for k in range(1,n-1):
        X.append(k/(n-1))
    X.append(1.0)
    X.sort()
    X = array(X)
    print('Collocation points: ', X)

    sigma = [sigma1,sigma2,sigma3][i]
    sigma_exact = [sigma1_exact,sigma2_exact,sigma3_exact][i]
    
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
        c = zeros(n)
        b[j] = -sigma(X[j],r)*X[j]
    c = solve(A,b)
    print("Error Ac-b: ", matmul(A,c)-b)
    print("Solution c: ", c)

    spl = BSpline(t,c,K,extrapolate=False)
    err = ABS(spl(X)-sigma_exact(X,r))
    print("Error in collocation points: ",err)
    print("Maximum error: ",MAX(err))
    ax.scatter(X,spl(X),facecolor='none',edgecolor='purple',s=200,label='Collocation points')

    X = array([j/1000 for j in range(1001)])
    ax.plot(X,sigma_exact(X,r),c='k',lw=5,label='Exact solution')
    ax.plot(X,spl(X),ls='--',c='r',label='B-spline approximation')
    ax.legend(loc='upper right',framealpha=0)
    ax.set_xlabel(r'$\xi$')
    ax.set_ylabel(r'$g(\xi)$')
    return spl
