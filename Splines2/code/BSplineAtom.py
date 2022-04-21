from numpy import array,zeros,pi
from Numerics import *
from scipy.interpolate import BSpline

Exact = {'u00': lambda r,D: 2*D**1.5*r*exp(-D*r),
         'u10': lambda r,D: D**1.5*r*(2-r*D)*exp(-D*r/2)/2**1.5,
         'u11': lambda r,D: D**2.5*r**2*exp(-D*r/2)/(2*6**0.5),}

def Coulomb(x):
    return -1/x

def UniformSphere(x,R0):
    if x<R0:
        return -(3-(x/R0)**2)/(2*R0)
    else:
        return -1/x

def getA(t,K,R0,ell,beta):

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

            f = lambda x: spl_i(x)*spl_j(x)*UniformSphere(x,R0)
            A[i-1,j-1]  = -2*Simpson(f,t[i],t[i+K+1])/beta

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

def numsol(dct):

    E    = dct['E']
    Z    = dct['Z']
    ell  = dct['ell']
    R0   = dct['R0']
    N    = dct['N']
    mume = dct['mass_ratio']
    D    = dct['D']
    beta = 1/(Z*mume*D)
    
    # Construct inner points
    #X = [j/N for j in range(1,N)]
    X = [1-(j/N)**0.4 for j in range(1,N)]
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

    A = getA(t,K,R0,ell,beta)
    B = getB(t,K)

    x0 = rand(len(A))
    eig = InvPower(A,B,x0,shift=E)
    print("Energy: ", eig[0])
    sol = array([0] + [j for j in eig[1]] + [0])
    u_spl = BSpline(t,sol,K,extrapolate=False)

    u_splSqrd = lambda x: u_spl(x)**2
    A   = Simpson(u_splSqrd,0,1,N=100)  # normalizing

    xx = array([j/1000 for j in range(1001)])
    plt.plot(xx,u_spl(xx))
    plt.scatter(t,[0 for i in t],c='k')
    plt.show()

    return (lambda x: u_spl(x)/sqrt(A))
