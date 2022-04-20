from Distributions import *
from PLOTTING import *
from numpy import array,zeros,exp,cos,pi,matmul,sin
from numpy.linalg import solve
from numpy.random import rand
from scipy.interpolate import BSpline
from numpy import max as MAX
from numpy import abs as ABS

def F(k,N,a,b):
    return a + (b-a)*k/(N+1)

def numsol(N,i,r):

    K = 6                                    # degree of B splines
    # Inner points
    X = [F(j,N,0,1) for j in range(1,N+1)]
    X.sort()
    #print("Inner points: ", X)

    # Knot points
    t = array((K+1)*[0.0] + X + (K+1)*[1.0])

    # Number of B-splines
    n = len(t)-K-1                           # = number of B splines
    #print("Knots: ", t)
    #print("Number of splines: ", n)

    # Plot the B-splines.
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    X = array([j/1000 for j in range(1001)])
    c = zeros(n)
    for k in range(n):
        c[k] = 1
        ax.plot(X,BSpline(t,c,K,extrapolate=False)(X),lw=0.75)
        c[k] = 0


    # Create collocation points.
    h = 1/(n-6)
    X = [0,h/10,h/5]
    X += [j*h for j in range(n-6)]
    X += [1-h/5,1-h/10,1.0]
    X = array(X)
    #print('Collocation points: ', X)

    # Get the distributions and the correct solution
    # to the distributions.
    sigma = [sigma1,sigma2,sigma3][i]
    sigma_exact = [sigma1_exact,sigma2_exact,sigma3_exact][i]
    
    # Construct matrix eq. and solve
    A = zeros((n,n)); b = zeros(n)
    A[0,0] = 1.0; A[n-1,n-1] = 1.0
    for j in range(n):
        # columns 
        c[j] = 1
        spl  = BSpline(t,c,K,extrapolate=False)
        spl  = spl.derivative(2)
        for k in range(1,n-1):
            # rows
            A[k,j] = spl(X[k])
        c[j] = 0
    for j in range(1,n-1):
        b[j] = -sigma(X[j],r)*X[j]
    c = solve(A,b)
    #print("Error Ac-b: ", matmul(A,c)-b)
    #print("Solution c: ", c)

    # The solution spline
    spl = BSpline(t,c,K,extrapolate=False)  
    ax.scatter(X,spl(X),facecolor='none',edgecolor='purple',s=200,label='Collocation points')
    
    # Evaluate error in collocation points
    err = ABS(spl(X)-sigma_exact(X,r))  
    #print("Error in collocation points: ",err)
    #print("Maximum error: ",MAX(err))

    # Plot the correct solution
    X = array([j/1000 for j in range(1001)])
    ax.plot(X,sigma_exact(X,r),c='k',lw=5,label='Exact solution')

    # Plot the approximated solution
    ax.plot(X,spl(X),ls='--',c='r',label='B-spline approximation')
    ax.legend(loc='upper right',framealpha=0)
    ax.set_xlabel(r'$\xi$')
    ax.set_ylabel(r'$g(\xi)$')

    plt.show()
    return spl


