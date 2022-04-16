from PLOTTING import *
import FILES
from numpy import exp,pi,array

psi = {1: lambda t: (1/pi)**0.25*exp(-t**2/2),
       2: lambda t: t*exp(-t**2/2),
       3: lambda t: (2*t**2-1)*exp(-t**2/2),
       4: lambda t: (2*t**3-3*t)*exp(-t**2/2)
      }

def main():

    x = []
    for j in range(1,80000):
        X = FILES.get_data('A/A001.dat',j)
        x.append(X[0])
    x = array(x)
    
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.hist(x,density=True,bins=100,histtype='step',lw=1,color='k')
    
    x = array([min(x) + j*(max(x)-min(x))/100 for j in range(101)])
    ax.plot(x,psi[1](x)**2,label=r'$\pi^{-1/2}e^{-\left(x\sqrt{m\omega/\hbar}\right)^2}$',lw=2,ls='dashdot',c='purple')
    ax.legend(framealpha=0,fontsize=32)
    ax.set_xlabel(r'$x\sqrt{m\omega/\hbar}$')
    ax.set_ylabel('Density')
    ax.set_xticks([-2,-1,0,1,2])
    ax.set_yticks([0,0.2,0.4,0.6,0.8,1])
    
    

    plt.show()

    

main()



