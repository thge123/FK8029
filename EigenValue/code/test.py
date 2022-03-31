from PLOTTING import *
from numpy import array,zeros,dot,matmul,copy,eye
from numpy.linalg import solve,norm

def RQI(A,x0,iters):

    I = eye(len(A))
    for i in range(iters):
        y  = x0/norm(x0)
        lam = dot(y,matmul(A,y))
        x0 = solve(A-lam*I,y)

    return x0/norm(x0)

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
    y0 = max(y)
    if params['parity'] == 'odd':
        y = array([0] + [j/y0 for j in y] + [0])
        ax.plot(x,y,c='k')
        ax.plot(-x,-y,c='k')
    else:
        y = array([j/y0 for j in y] + [0])
        ax.plot(x,y,c='k')
        ax.plot(-x,y,c='k')

    
    plt.show()

def Energy(A,x0,params):

    eig = dot(x0,matmul(A,x0))/dot(x0,x0)
    return eig/2

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

    params = {'N': 300, 'h': 0.05, 'parity': 'odd'}

    if params['parity'] == 'even':
        A = EvenMatrix(params)
    else:
        A = OddMatrix(params)

    x0 = initial_guess(params)
    shiftA = shift(A,3)
    x0 = inverse_power_iteration(shiftA,x0,20)

    print(Energy(A,x0,params))
    
    plot(params,x0)

main()
