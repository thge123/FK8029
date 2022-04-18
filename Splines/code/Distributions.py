from numpy import array,zeros,exp,cos,pi,nan

def sigma1(x,r):
    return 3.0

def sigma1_exact(X,r):
    ''' Takes in X numpy array and dummy r.
        Gives exact sol. to bndry. problem
            g''(x) + sigma*x = 0
            g(0) = g(1) = 0
    '''
    return -X**3/2+0.5*X

def sigma2(x,r):
    return 3*(x>r)/(1-r**3)

def sigma2_exact(X,r):
    ''' Takes in X numpy array and r=R1/R2.
        Gives exact sol. to bndry. problem
            g''(x) + sigma2*x = 0
            g(0) = g(1) = 0
    '''
    A1 = ( -1.5*(r**2)/(1-r**3) + 0.5/(1-r) )/(1+r/(1-r))
    A2 = (0.5-A1*r)/(1-r)
    B2 = 1/(2*(1-r**3)) - A2
    Y  = (X<r)*A1*X + (X>=r)*(-X**3/(2*(1-r**3))+A2*X+B2)
    return Y

def sigma3(x,alpha):
    ''' Ground state of hydrogen '''
    return 4*alpha**3*exp(-2*alpha*x)

def sigma3_exact(X,alpha):
    A = (alpha+1)*exp(-2*alpha)-1
    B = 1
    Y = -(alpha*X+1)*exp(-2*alpha*X)+A*X+B
    return Y

def sigma4(x,alpha):
    ''' 200 state of hydrogen '''
    return (2-alpha*x)**2*alpha**3*exp(-alpha*x)/8

def sigma4_exact(X,alpha):
    return nan

def sigma5(x,alpha):
    ''' 300 state of hydrogen '''
    return 4*alpha**3*(27-18*alpha*x+2*(alpha*x)**2)**2*exp(-2*alpha*x/3)/(3*81**2)

def sigma5_exact(X,alpha):
    return nan
