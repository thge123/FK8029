from numpy import array,zeros,exp,cos,pi,matmul
from numpy import copy,dot
from numpy.linalg import solve,norm
from numpy.random import rand
from scipy.special import roots_legendre


def InvPower(A,B,x0,shift=0,N=100):
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
    for i in range(N):
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

def Quadrature(f,a,b,N=20,deg=8):

    h = (b-a)/N
    Sum = 0
    roots  = roots_legendre(deg)
    weights = roots[1]
    roots  = roots[0]

    beg = a 
    end = a+h
    for i in range(N):
        inner_sum = 0
        for j in range(len(roots)):
            x = (roots[j]-(beg+end)/(beg-end))*((end-beg)/2)
            inner_sum += weights[j]*f(x)
        inner_sum *= (end-beg)/2
        Sum += inner_sum
        beg = beg + h
        end = beg + h

    return Sum


