from PLOTTING import *
from numpy import array,zeros,dot,matmul,copy,eye,exp
from numpy.linalg import solve,norm
from numpy import abs as Abs
from numpy import max as Max
from numpy.random import rand

psi = {'0': lambda t: exp(-t**2/2),
       '1': lambda t: t*exp(-t**2/2),
       '2': lambda t: (2*t**2-1)*exp(-t**2/2)
      }

def RQI(A,x0,iters):

    I = eye(len(A))
    for i in range(iters):
        y  = x0/norm(x0)
        lam = dot(y,matmul(A,y))
        x0 = solve(A-lam*I,y)

    return x0/norm(x0)

#def inverse_power_iteration(A,x0,iters):
#
#    for i in range(iters):
#        y  = x0/norm(x0)
#        x0 = solve(A,y)
#
#    return x0/norm(x0)

def inverse_power_iteration(B,x0,s):

    k = 1
    A = B-s*eye(len(B))
    while True:
        y  = x0/norm(x0)
        x0 = solve(A,y)
        if k%10 == 0:
            lam = dot(x0,matmul(B,x0))/dot(x0,x0)
            print(lam)
            Delta = Abs(matmul(B,x0)-lam*x0)
            if Max(Delta) < 1e-6:
                break
            elif k>1000:
                print('No eigenvector found')
                break
        k+=1

    return x0/norm(x0)

def OddMatrix(params):

    N = params['N']
    h = params['h']

    A = zeros((N-2,N-2))
    
    # Skewed fourth order approx. 
    # of second derivative
    A[0,0] = -15
    A[0,1] = -4
    A[0,2] = 14
    A[0,3] = -6
    A[0,4] = 1

    A[1,0]   = 16
    A[1,1]   = -30
    A[1,2] = 16
    A[1,3] = -1

    # Central fourth order approx.
    # of second derivative.
    for i in range(2,N-4):
        A[i,i-2] = -1
        A[i,i-1] = 16
        A[i,i]   = -30
        A[i,i+1] = 16
        A[i,i+2] = -1

    # Same as above except zeros at 
    # right bndry.
    A[N-4,N-6] = -1
    A[N-4,N-5] = 16
    A[N-4,N-4] = -30
    A[N-4,N-3] = 16

    # Same skewed approx. as above.
    A[N-3,N-7] = 1
    A[N-3,N-6] = -6
    A[N-3,N-5] = 14
    A[N-3,N-4] = -4
    A[N-3,N-3] = -15

    A = -A/(12*h**2)

    # Potential
    for i in range(N-2):
        x = (i+1)*h
        A[i,i] += test_pot(x)
    
    return A

def EvenMatrix(params):
    N = params['N']
    h = params['h']

    A = zeros((N-1,N-1))

    A[0,0] = -80/3 - 30
    A[0,1] = 72
    A[0,2] = -18
    A[0,3] = 8/3
    
    A[1,0] = 10/3 + 16
    A[1,1] = -36
    A[1,2] = 18
    A[1,3] = -4/3

    for i in range(2,N-3):
        A[i,i-2] = -1
        A[i,i-1] = 16
        A[i,i]   = -30
        A[i,i+1] = 16
        A[i,i+2] = -1

    # Same as above except zeros at 
    # right bndry.
    A[N-3,N-5] = -1
    A[N-3,N-4] = 16
    A[N-3,N-3] = -30
    A[N-3,N-2] = 16

    # Same skewed approx. as above.
    A[N-2,N-6] = 1
    A[N-2,N-5] = -6
    A[N-2,N-4] = 14
    A[N-2,N-3] = -4
    A[N-2,N-2] = -15

    A = -A/(12*h**2)

    # Potential
    for i in range(N-1):
        x = ((i)*h)
        A[i,i] += test_pot(x)

    return A

def har_osc(x):
    return x**2

def test_pot(x):
    return x**2

def initial_guess(params):

    N  = params['N']
    if params['parity'] == 'odd':
        x0 = rand(N-2)
    else:
        x0 = rand(N-1)
    x0[0] = 1

    return x0

def normalized_density(X,y):

    I = 0
    y_out = y**2
    for i in range(len(X)-1):
        y_avg = (y_out[i+1]+y_out[i])/2
        I += y_avg*(X[i+1]-X[i])
    return y_out/(2*I)
        
def plot(params,y,density=False):
    
    h = params['h']
    N = params['N']

    x = array([j*h for j in range(N)])
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    if params['parity'] == 'odd':
        y = array([0] + [j for j in y] + [0])
        if density:
            y = normalized_density(x,y)
            ax.plot(x,y,c='k')
            ax.plot(-x,y,c='k')
        else:
            ax.plot(x,y,c='k')
            ax.plot(-x,-y,c='k')
    else:
        y = array([j for j in y] + [0])
        if density:
            y = normalized_density(x,y)
            ax.plot(x,y,c='k')
            ax.plot(-x,y,c='k')
        else:
            ax.plot(x,y,c='k')
            ax.plot(-x,y,c='k')
    pot = [test_pot(j) for j in x]
    pot = [j/max(pot) for j in pot]
    ax.plot(x,pot,ls='--',c='k')
    ax.plot(-x,pot,ls='--',c='k')

    plt.show()

def Energy(A,x0):

    eig = dot(x0,matmul(A,x0))/dot(x0,x0)
    #print(matmul(A,x0)-eig*x0)
    return eig

def shift(A,c):
    B = copy(A)
    for i in range(len(A)):
        B[i,i] -= c
    return B

def main():

    ''' The eigenvalues: 
        
        E =   hf/2  (even)  --> shift = 1
        E =  3hf/2  (odd)   --> shift = 3
        E =  5hf/2  (even)  --> shift = 5
        E =  7hf/2  (odd)   --> shift = 7
        W =  9hf/2  (even)  --> shift = 9
        E = 11hf/2  (odd)   --> shift = 11 etc. 

    '''

    params = {'N': 500, 'h': 0.05, 'parity': 'odd'}

    A = EvenMatrix(params)
    B = OddMatrix(params)

    
    x0 = initial_guess(params)
    x0 = inverse_power_iteration(B,x0,1)
    print('i=',3,'  ', 'E=',Energy(B,x0))

    #if params['parity'] == 'even':
    #    A = EvenMatrix(params)
    #else:
    #    A = OddMatrix(params)

#    for i in range(0,30):
#        if i%2 == 0:
#            params['parity'] = 'even'
#            x0 = initial_guess(params)
#            x0 = inverse_power_iteration(A,x0,i+rand())
#            print('i=',i,'  ', 'E=',Energy(A,x0))
#        else:
#            params['parity'] = 'odd'
#            x0 = initial_guess(params)
#            x0 = inverse_power_iteration(B,x0,i+rand())
#            print('i=',i,'  ', 'E=',Energy(B,x0))

#    x0 = initial_guess(params)
#    shiftA = shift(A,3)
#    x0 = inverse_power_iteration(A,x0,3.2)
#    print(Energy(A,x0))
    plot(params,x0,density=False)

main()
