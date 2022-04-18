from PLOTTING import *
from BSpline2 import *

def main():

    
    N = int(input("Write number of splines: "))-3
    X   = array([j/1000 for j in range(1001)])

    r = float(input("Write alpha: "))

    fig = plt.figure()
    ax  = fig.add_subplot(211)
    axx = fig.add_subplot(212)
    axs = [ax,axx] 


    legends = {0:r'$\varphi/(q/4\pi\epsilon_0 ({}a_0))$'.format(int(r)),
               1:r'$E/(q/4\pi\epsilon_0 ({}a_0)^2)$'.format(int(r))}

    for i in [3,4]:

        spl = numsol(N,i,r)
        # Potential times distance
        uu  = spl(X) + X

        # Potentials
        VV  = uu[1:]/X[1:]

        # Electric fields
        EE = array([0] + [-(VV[j+1]-VV[j])/(X[1]-X[0]) for j in range(len(VV)-1)])

        #ax.plot(X,u,label='Potential times distance')
        axs[i-3].plot(X[1:],VV,lw=2,c='b',ls='--',label=legends[0])
        axs[i-3].plot(X[:-1],EE,c='k',lw=2,ls='dashdot',label=legends[1])

    for i in axs:
        i.grid(linestyle='--')

    axs[0].legend(framealpha=0,loc='upper right')
    axs[1].set_xlabel(r'$r/({}a_0)$'.format(int(r)))

    axs[0].text(s="2s state",x=0.18,y=20)
    axs[1].text(s="3s state",x=0.18,y=5.5)

    
    plt.show()

main()
