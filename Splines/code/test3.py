from Distributions import *
from PLOTTING import *
from numpy import array,zeros,exp,cos,sin,pi,nan,sqrt
from numpy.linalg import solve
'''
def sigma1(x,r):
    return 3.0

def sigma1_exact(X,r):
            g''(x) + sigma*x = 0
            g(0) = g(1) = 0
   
    return -X**3/2+0.5*X

def sigma2(x,r):
    return 3*(x>r)/(1-r**3)

def sigma2_exact(X,r):
    
    A1 = ( -1.5*(r**2)/(1-r**3) + 0.5/(1-r) )/(1+r/(1-r))
    A2 = (0.5-A1*r)/(1-r)
    B2 = 1/(2*(1-r**3)) - A2
    Y  = (X<r)*A1*X + (X>=r)*(-X**3/(2*(1-r**3))+A2*X+B2)
    return Y

def sigma3(x,alpha):
    return 4*alpha**3*exp(-2*alpha*x)

def sigma3_exact(X,alpha):
    A = (alpha+1)*exp(-2*alpha)-1
    B = 1
    Y = -(alpha*X+1)*exp(-2*alpha*X)+A*X+B
    return Y

def sigma4(x,alpha):
    return 4*sqrt(alpha*pi)*exp(-alpha*(x-0.5)**2)

def sigma4_exact(x,alpha):
    return nan

'''

def numsol(N,i,r):

    sigma = [sigma1,sigma2,sigma3,sigma4,sigma5][i]

    h = 1/(N-1)
    A = zeros((N-2,N-2))
    b = zeros(N-2)
    
    A[0,0] = -2.0
    A[0,1] = 1.0
    b[0]   = -h**3*sigma(h,r)
    print(h)
    for i in range(1,N-3):
        A[i,i-1] = 1.0
        A[i,i]   = -2.0
        A[i,i+1] = 1.0
        b[i]     = -(i+1)*h**3*sigma((i+1)*h,r)
        print((i+1)*h)
    A[N-3,N-4]   = 1.0
    A[N-3,N-3]   = -2.0
    b[N-3]       = -(N-2)*h**3*sigma((N-1)*h,r)
    print((N-2)*h)

    X = array([j*h for j in range(N)])
    Y = solve(A,b)
    Y = array([0] + [j for j in Y] + [0])
    return X,Y

def main():

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    i = int(input("Write sigma [1,2,3,4,5]: ")) - 1 
    #N = int(input("Write N: "))
    exact_sol = [sigma1_exact,sigma2_exact,sigma3_exact,sigma4_exact,sigma5_exact][i]

    if i == 1:
        r = float(input("Write r: "))
    elif i in [2,3,4]:
        r = float(input("Write alpha: "))
    else:
        r = 1

    for N in [200,250]:
        # Potential times distance
        X,uu  = numsol(N,i,r)
        uu    = uu+X
        u     = exact_sol(X,r) + X 

        # Potentials
        V   = u[1:]/X[1:]
        VV  = uu[1:]/X[1:]

        # Electric fields
        E  = array([0] + [-(V[j+1]-V[j])/(X[1]-X[0]) for j in range(len(V)-1)])
        EE = array([0] + [-(VV[j+1]-VV[j])/(X[1]-X[0]) for j in range(len(V)-1)])

        ax.plot(X,u)
        ax.plot(X[1:],V)
        ax.plot(X[:-1],E)

        ax.plot(X,uu,label='Potential times distance')
        ax.plot(X[1:],VV,label='Potential')
        ax.plot(X[:-1],EE,label='Electric field')
    #ax.legend()
    plt.show()

main()

        

