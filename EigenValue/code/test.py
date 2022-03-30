from PLOTTING import *
from numpy import array,zeros,dot,matmul
from numpy.linalg import solve,norm

def inverse_power_iteration(A,x0,iters):

    for i in range(iters):
        y  = x0/norm(x0)
        x0 = solve(A,y)

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

    # Central fourth order approx.
    # of second derivative.
    for i in range(1,N-4):
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
        A[i,i] += ((i+1)*h)**2
    
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
        A[i,i] += ((i)*h)**2

    return A

def initial_guess(params):

    N  = params['N']
    if params['parity'] == 'odd':
        x0 = zeros(N-2)
    else:
        x0 = zeros(N-1)
    x0[0] = 1

    return x0

def plot(params,y):
    
    h = params['h']
    N = params['N']

    x = array([j*h for j in range(N)])
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    if params['parity'] == 'odd':
        y = array([0] + [j for j in y] + [0])
        ax.plot(x,y)
        ax.plot(-x,-y)
    else:
        y = array([j for j in y] + [0])
        ax.plot(x,y)
        ax.plot(-x,y)

    
    plt.show()

def Energy(A,x0,params):

    eig = dot(x0,matmul(A,x0))/dot(x0,x0)
    h = params['h']
    return eig/2

def main():

    params = {'N': 400, 'h': 0.01, 'parity': 'even'}

    A = EvenMatrix(params)
    B = OddMatrix(params)
    x0 = initial_guess(params)
    x0 = inverse_power_iteration(A,x0,5)

    print(Energy(A,x0,params))
    
    plot(params,x0)

main()
