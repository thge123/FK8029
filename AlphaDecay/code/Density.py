from PLOTTING import *
from XY import *
from numpy import exp

def main():
    
    N = 21
    E = float(input("Write E: "))
    r = XY()
    r.get_X('r.dat')
    r.get_Y('omega.dat')
    
    Coeffs = XY()
    Coeffs.get_X('X.dat',Complex=True)
    Coeffs.get_Y('V.dat')

    x = []
    y = []

    for i in range(len(r.X)-1):
        xj = [r.X[i] + j*(r.X[i+1]-r.X[i])/(N-1) for j in range(N-1)]
        yj = []
        A = Coeffs.X[2*i][0]   + (Coeffs.X[2*i][1])*1j
        print(A)
        B = Coeffs.X[2*i+1][0] + (Coeffs.X[2*i+1][1])*1j

        V = Coeffs.Y[i]
        w = r.Y[i]
        if V-E > 0:
            f = lambda t: A*exp(w*t)+B*exp(-w*t)
        elif V-E < 0:
            f = lambda t: A*exp(1j*w*t)+B*exp(-1j*w*t)

        for i in xj:
            yj.append(abs(f(i))**2)
        x = x+xj
        y = y+yj

    xj = [r.X[-1] + j*(r.X[-1]-r.X[-2])/(N-1) for j in range(N)]
    yj = []
    F = Coeffs.X[-1][0] + Coeffs.X[-1][1]*1j
    V = Coeffs.Y[-1]
    w = r.Y[-1]
    if V-E > 0:
        f = lambda t: F*exp(-w*t)
    elif V-E < 0:
        f = lambda t: F*exp(1j*w*t)
    for i in xj:
        yj.append(abs(f(i))**2)
    x = x+xj
    y = y+yj
    print(x)
    print(y)

    plt.plot(x,y)
    plt.show()

main()
