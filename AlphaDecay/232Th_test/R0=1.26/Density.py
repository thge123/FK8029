from PLOTTING import *
import FILES
from numpy import exp,cos,sin,log
from os import listdir

def get_data(number):
    numberstr = FILES.number2string(number)

    r      = FILES.get_data('R/R{}.dat'.format(numberstr))
    omega  = FILES.get_data('O/O{}.dat'.format(numberstr))
    Coeffs = FILES.get_data('X/X{}.dat'.format(numberstr),Complex = True)
    V      = FILES.get_data('V/V{}.dat'.format(numberstr))
    
    return r,omega,Coeffs,V

def plot_density(ax,E,filenumber):

    N = 1000
    r,omega,Coeffs,V = get_data(filenumber)

    x = []
    y = []

    for i in range(len(r)-1):
        xj = [r[i] + j*(r[i+1]-r[i])/(N-1) for j in range(N-1)]
        yj = []
        A = Coeffs[2*i][0]   + (Coeffs[2*i][1])*1j
        B = Coeffs[2*i+1][0] + (Coeffs[2*i+1][1])*1j
        w = omega[i]

        if V[i]-E > 0:
            f = lambda t: A*exp(w*t)+B*exp(-w*t)
        elif V[i]-E < 0:
            f = lambda t: A*exp(1j*w*t)+B*exp(-1j*w*t)
        for i in xj:
            yj.append(abs(f(i))**2)
        x = x+xj
        y = y+yj

    xj = [r[-1] + j*(r[-1]-r[-2])/(N-1) for j in range(N)]
    yj = []
    F = Coeffs[-1][0] + (Coeffs[-1][1])*1j
    V = Coeffs[-1]
    w = r[-1]
    if V[-1]-E > 0:
        f = lambda t: F*exp(-w*t)
    elif V[-1]-E < 0:
        f = lambda t: F*exp(1j*w*t)
    for i in xj:
        yj.append(abs(f(i))**2)
    x = x+xj
    y = y+yj

    ax.plot(x,y,lw=2)

    

def main():
    
    E = float(input("Write E: "))
    fig = plt.figure()
    ax = fig.add_subplot(111)

    filenumber = int(input("Write filenumber: "))
    #plot_density(ax,E,filenumber)
    #plot_density(ax,E,filenumber)
    plot_density(ax,E,1)
    plot_density(ax,E,50)
    plt.show()


