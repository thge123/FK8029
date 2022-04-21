from PLOTTING import *
from scipy.interpolate import BSpline
import matplotlib.pyplot as plt
from numpy import zeros,copy,matmul,eye,dot,array,nan,isnan,exp,sqrt
from numpy import max as Max
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
   
    print("Inverse power method error: ", Max(err))
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

def getA(t,K,ell=0,beta=1):
    ''' The "hamiltonian matrix". The matrix elements are
        given by 
            A_ij = -beta^2*( -int(B'_ik*B'_jk) 
                             -l(l+1)*int(B_ik*B_jk/x^2)
                             +(2/beta)*int(B_ik*B_jk/x) )
    '''
   
    n = len(t)-K-1
    A = zeros((n-2,n-2))

    for i in range(1,n-1):
        c_i    = zeros(n)
        c_i[i] = 1
        spl_i  = BSpline(t,c_i,K,extrapolate=False)
        dspl_i = spl_i.derivative()
        if i+K+1>n-1:
            m = n-1
        else:
            m = i+K+1
        for j in range(i,m):
            # B_(>i+K) is non-overlapping with B_i

            c_j    = zeros(n)
            c_j[j] = 1
            spl_j  = BSpline(t,c_j,K,extrapolate=False)
            dspl_j = spl_j.derivative()

            f = lambda x: spl_i(x)*spl_j(x)/x
            A[i-1,j-1]  = 2*Simpson(f,t[i],t[i+K+1])/beta

            f = lambda x: spl_i(x)*spl_j(x)/x**2
            A[i-1,j-1] += -ell*(ell+1)*Simpson(f,t[j],t[i+K+1])


            f = lambda x: dspl_i(x)*dspl_j(x)
            A[i-1,j-1] += -Simpson(f,t[i],t[i+K+1])

            xx = array([j/100 for j in range(1001)])
            #plt.plot(xx,spl_i(xx))
            #plt.plot(xx,spl_j(xx))
            #plt.plot(xx,f(xx))
            #plt.scatter(t[j],0)
            #plt.scatter(t[i+K+1],0)
            #plt.show()
    
            A[i-1,j-1] *= -100*beta**2
            A[j-1,i-1]  = A[i-1,j-1]

    return A


def getB(t,K):
    n = len(t)-K-1
    B = zeros((n-2,n-2))
    for i in range(1,n-1):
        c_i    = zeros(n)
        c_i[i] = 1
        spl_i  = BSpline(t,c_i,K,extrapolate=False)
        if i+K+1>n-1:
            m = n-1
        else:
            m = i+K+1
        for j in range(i,m):
            # B_(>i+K) is non-overlapping with B_i
            c_j    = zeros(n)
            c_j[j] = 1
            spl_j  = BSpline(t,c_j,K,extrapolate=False)
            spl_ij = lambda x: spl_i(x)*spl_j(x)
            B[i-1,j-1] = Simpson(spl_ij,t[i],t[i+K+1])
            B[j-1,i-1] = B[i-1,j-1]
    return B

def numsol(E,Beta,Ell):

    # Number of inner points
    N = 50

    # Construct inner points
    #X = [j/N for j in range(1,N)]
    X = [1-(j/N)**0.5 for j in range(1,N)]
    X.sort()

    # Degree K and knots t.
    K = 4
    t = array((K+1)*[0] + X + (K+1)*[1])

    # Total number of splines
    n = len(t)-K-1
    print("Number of splines: ", n)

    # We wish to solve eigenvalue problem
    #           Ac = lam*Bc.
    # The method uses power iteration.

    A = getA(t,K,ell=Ell,beta=Beta)
    B = getB(t,K)

    x0 = rand(len(A))
    eig = InvPower(A,B,x0,shift=E)
    print("Energy: ", eig[0])
    sol = array([0] + [j for j in eig[1]] + [0])
    u_spl = BSpline(t,sol,K,extrapolate=False)
    u_splSqrd = lambda x: u_spl(x)**2

    xx = array([j/1000 for j in range(1001)])
    plt.plot(xx,u_spl(xx))
    plt.scatter(t,[0 for i in t],c='k')
    plt.show()

    A   = Simpson(u_splSqrd,0,1,N=100)  # normalizing
    u_splSqrd = lambda x: u_spl(x)**2/A

    return (lambda x: u_spl(x)/sqrt(A))


def main():

    mu_me_ratio = 1
    Z = 1
    D = 20

    Beta = 1/(Z*mu_me_ratio*D)
    Ell  = 0
    E    = -100

    u = numsol(E,Beta,Ell)

    xx = array([j/1000 for j in range(1,1001)])
    plt.plot(xx,u(xx))


#    ### Testing
#    
#    # Ground state u
#    u_exc = lambda x: x*exp(-D*x)
#    u_excSqrd = lambda x: u_exc(x)**2
#    B   = Simpson(u_excSqrd,0,1,N=100)  # normalizing
#    u_exc = lambda x: x*exp(-D*x)/sqrt(B)
#    plt.plot(xx,u_exc(xx))
#    plt.show()

main()
