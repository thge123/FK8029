from PLOTTING import *
from Numerics import *
from numpy import array,zeros,exp,cos,pi,matmul,sqrt
from numpy.linalg import solve
from numpy.random import rand
from scipy.interpolate import BSpline

def knots1(X,a):
    return [1-j**a for j in X]

def knots2(X,D):
    return [exp(-D*j) for j in X]

def knots3(X,a):
    return [2**(j**a)-1 for j in X]

def numsol_Poisson(Poisson,PlotSplines=False):

    # "Source" function Poisson's eq.
    sigma = Poisson['sigma']
    N     = Poisson['N']
    
    # Order of B-splines
    K = 3

    # Knot points
    X = [j/(N+1) for j in range(1,N+1)]
    X.sort()
    t = array((K+1)*[0] + X + (K+1)*[1])     

    # Number of splines
    n = len(t)-K-1
    print("Number of splines: ", n)

    # Plot the boundary value solution along with B-splines
    if PlotSplines:
        xx = array([j/1000 for j in range(1001)])
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        for k in range(n):
            c = zeros(n)
            c[k] = 1
            ax.plot(xx,BSpline(t,c,K,extrapolate=False)(xx),lw=0.75)
        ax.scatter(t,[0 for j in t])

    # Create collocation points
    # Use knot points plus K-1 other ones
    X = [0] + X + [1]
    X += [0.1*X[1],1-0.1*X[1]]
    X.sort()
    X = array(X)
    
    # Construct matrix eq. and solve
    A = zeros((n,n))
    b = zeros(n)
    A[0,0]     = 1.0
    A[n-1,n-1] = 1.0
    for j in range(1,n-1):
        c    = zeros(n)
        c[j] = 1
        spl  = BSpline(t,c,K,extrapolate=False)
        spl  = spl.derivative(2)
        for k in range(1,n-1):
            A[k,j] = spl(X[k])
        b[j] = -sigma(X[j])*X[j]
    c = solve(A,b)
    err = matmul(A,c)-b
    print("Error Ac-b: ", max(abs(err.min()),err.max()))
    #print("Solution c: ", c)
    
    # The B-spline solution
    spl = BSpline(t,c,K,extrapolate=False)    # The B-spline solution

    if PlotSplines: 
        ax.scatter(X,spl(X),facecolor='none',
                   edgecolor='purple',s=200,
                   label='Collocation points')
        X = array([j/1000 for j in range(1001)])
        ax.plot(X,spl(X),ls='--',c='r',label='B-spline approximation')
        ax.legend(loc='upper right',framealpha=0)
        ax.set_xlabel(r'$\xi$')
        ax.set_ylabel(r'$g(\xi)$')

        plt.show()

    return spl

def getA(t,K,R0,ell,beta,V,quadN=20,quadDeg=8):

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

            f = lambda x: spl_i(x)*spl_j(x)*V(x)
            A[i-1,j-1]  = -2*Quadrature(f,t[i],t[i+K+1],N=quadN,deg=quadDeg)/beta
            #A[i-1,j-1]  = -2*Simpson(f,t[i],t[i+K+1])/beta

            f = lambda x: spl_i(x)*spl_j(x)/x**2
            A[i-1,j-1] += -ell*(ell+1)*Quadrature(f,t[j],t[i+K+1],N=quadN,deg=quadDeg)
            #A[i-1,j-1] += -ell*(ell+1)*Simpson(f,t[j],t[i+K+1])


            f = lambda x: dspl_i(x)*dspl_j(x)
            A[i-1,j-1] += -Quadrature(f,t[i],t[i+K+1],N=quadN,deg=quadDeg)
            #A[i-1,j-1] += -Simpson(f,t[i],t[i+K+1])

            xx = array([j/100 for j in range(1001)])
    
            A[i-1,j-1] *= -100*beta**2
            A[j-1,i-1]  = A[i-1,j-1]

    return A


def getB(t,K,quadN=20,quadDeg=8):
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
            B[i-1,j-1] = Quadrature(spl_ij,t[i],t[i+K+1],N=quadN,deg=quadDeg)
            #B[i-1,j-1] = Simpson(spl_ij,t[i],t[i+K+1])
            B[j-1,i-1] = B[i-1,j-1]

    return B

def numsol_radial(dct,PlotSplines=False):

    E       = dct['E']
    Z       = dct['Z']
    ell     = dct['ell']
    R0      = dct['R0']
    N       = dct['N']
    mume    = dct['mass_ratio']
    D       = dct['D']
    beta    = 1/(Z*mume*D)
    quadN   = dct['quadN']
    quadDeg = dct['quadDeg']
    V       = dct['pot']
    
    # Construct inner points
    X = [j/N for j in range(1,N)]
    #X = knots1(X,0.7)
    #X = knots2(X,60*D)
    X = knots3(X,2)
    X.sort()

    # Degree K and knots t.
    K = 3
    t = array((K+1)*[0] + X + (K+1)*[1])

    # Total number of splines
    n = len(t)-K-1
    print("Number of splines: ", n)

    # We wish to solve eigenvalue problem
    #           Ac = lam*Bc.
    # The method uses power iteration.

    A = getA(t,K,R0,ell,beta,V,quadN=quadN,quadDeg=quadDeg)
    B = getB(t,K,quadN,quadDeg=quadDeg)

    x0 = rand(len(A))
    eig = InvPower(A,B,x0,shift=E)
    print("Energy: ", eig[0])
    if 'eigs' in dct:
        dct['eigs'].append(eig[0])
    else:
        dct['eigs'] = [eig[0]]
    sol = array([0] + [j for j in eig[1]] + [0])
    u_spl = BSpline(t,sol,K,extrapolate=False)

    # normalizing
    u_splSqrd = lambda x: u_spl(x)**2
    A   = Quadrature(u_splSqrd,0,1,N=quadN,deg=quadDeg)
    #A   = Simpson(u_splSqrd,0,1)  # normalizing
    u   = lambda x: u_spl(x)/sqrt(A)

    if PlotSplines:
        ### plot splines
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        xx = array([j/1000 for j in range(1001)])
        ax.set_ylabel(r'$u(\xi)$')
        ax.set_xlabel(r'$\xi$')

        for i in range(n):
            c_i    = zeros(n)
            c_i[i] = 1
            spl_i  = BSpline(t,c_i,K,extrapolate=False)
            ax.plot(xx,spl_i(xx),lw=1)

        ax.plot(xx,u(xx),lw=4,c='k',label='Numerical')
        ax.scatter(t,[0 for i in t],c='k')
        ax.legend(framealpha=0)
        plt.show()
        ###

    return u

def numsol():

    
    V = lambda x: -1/x
    Atom = {'Z'         : 2,     # Atomic number
            'ell'       : 0,     # Angular momentum
            'E'         : -99,   # Guess at eigenvalue
            'N'         : 20,    # Number of inner points
            'D'         : 12,    # Scaling of interval
            'mass_ratio': 1,     # red. mass/mass of elect.
            'R0'        : 0.0,   # Uniform shell cut-off
            'quadN'     : 20,    # Quadrature nodes 
            'quadDeg'   : 5,     # Quadrature degree
            'pot'       : V}     # Potential 

    Poisson = {'N': 30}

    xx = array([j/500 for j in range(501)])
    rr = Atom['D']*xx

    for i in range(2):

        u = numsol_radial(Atom,PlotSplines=True)
        rho = lambda x: u(x)**2/x**2
        plt.plot(xx,rho(xx)*xx**2)
        plt.show()

        Poisson['sigma'] = lambda x: rho(x)

        g   = numsol_Poisson(Poisson,PlotSplines=True)
        Lam = lambda x: g(x)+2*x/(4*pi)
        
        xx = array([j/100 for j in range(101)])
        plt.plot(xx,Lam(xx))
        plt.show()

        V = lambda x: -1/x+4*pi*Lam(x)/(Atom['Z']*x)
        Atom['pot'] = V
    print(Atom['eigs'])

    return 0

