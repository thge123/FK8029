from scipy.interpolate import BSpline
import matplotlib.pyplot as plt
from numpy import zeros,copy,matmul,eye,dot,array,nan,isnan,exp,sqrt
from numpy.linalg import norm,solve
from numpy.random import rand

def InvPower(A,B,x0,shift=0):
    ''' Solve the generalized eigenvalue
        problem Ax = lam*Bx where lam is
        an eigenvalue.

        Input:
        ----------
        A:     nxn numpy array.
        B:     nxn numpy array.
        x0:    nx1 numpy array. Initial guess.
        shift: Approximation of eigenvalue.

        Output:
        ----------
        x:   Approximated eigenvector.
        lam: Approximated eigenvalue.
    '''

    x = copy(x0)
    for i in range(100):
        x = solve(A-shift*B,matmul(B,x))
        x = x/norm(x)
    invB_Ax = solve(B,matmul(A,x))
    lam     = dot(x,invB_Ax)/dot(x,x)

    err = matmul(A,x)-lam*matmul(B,x)
   
    print("Inverse power method error: ", err.max())
    return lam,x

def Simpson(f,a,b,N=100):
    ''' Numerical integration using
        Simpson's rule.
    '''
   
    h = (b-a)/N
    f0 = 0
    Sum = 0
    Sum += f0+4*f(a+h/2)+f(a+N*h)
    for j in range(1,N):
        Sum += 2*f(a+j*h)
        Sum += 4*f(a+j*h+h/2)
    return Sum*h/6
