from PLOTTING import *
import FILES
from numpy import exp,pi
psi = {1: lambda t: (1/pi)**0.25*exp(-t**2/2),
       2: lambda t: t*exp(-t**2/2),
       3: lambda t: (2*t**2-1)*exp(-t**2/2),
       4: lambda t: (2*t**3-3*t)*exp(-t**2/2)
      }

def main():

    x = []
    for j in range(1,10000):
        X = FILES.get_data('A/A001.dat',j)
        x.append(X[0])
    
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.hist(x,density=True,bins=100)

    plt.show()

    

main()



