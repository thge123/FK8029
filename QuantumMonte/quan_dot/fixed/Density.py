from PLOTTING import *
import FILES
from numpy import array,exp,log,sqrt
from numpy import abs as ABS
from seaborn import kdeplot
from numpy.random import normal

def data(filenumber):
    data = []
    r1 = []
    r2 = []
    filename   = 'A/A'+FILES.number2string(filenumber)+'.dat'
    data = FILES.get_data(filename,fileline=0)
    for i in data:
        r1.append([i[0],i[1]])
        r2.append([i[2],i[3]])
    r1 = array(r1)
    r2 = array(r2)
    return r1,r2

def main():

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    fig2 = plt.figure()

    ax2 = fig2.add_subplot(221)
    ax3 = fig2.add_subplot(222)
    ax4 = fig2.add_subplot(223)
    ax5 = fig2.add_subplot(224)
    axs = [ax2,ax3,ax4,ax5]
    
    adct = {1:0,8:1,20:2,30:3}

    for filenumber in [1,2,8,20,30]:
#    for filenumber in [999,1]:
        r1,r2 = data(filenumber)
        x1 = array([j[0] for j in r1])
        y1 = array([j[1] for j in r1])
        x2 = array([j[0] for j in r2])
        y2 = array([j[1] for j in r2])
        r1 = x1**2+y1**2
        r2 = x2**2+y2**2
        r12 = (x1-x2)**2+(y1-y2)**2
        ax.hist(r12,bins=100,histtype='step',density=True,lw=2)
        if filenumber in adct:
            axs[adct[filenumber]].hist2d(x1,y1,bins=100,density=True,cmap=plt.cm.Greys)
            axs[adct[filenumber]].set_aspect('equal')
            axs[adct[filenumber]].set_xlim(-3,3)
            axs[adct[filenumber]].set_xticks([-2,0,2])
            axs[adct[filenumber]].set_ylim(-3,3)
            axs[adct[filenumber]].set_yticks([-2,0,2])
            axs[adct[filenumber]].set_title(r'$\lambda = {}$'.format(filenumber),fontsize=22)

            
    axs[2].set_xlabel(r'$x\sqrt{m\omega/\hbar}$')
    axs[2].set_ylabel(r'$y\sqrt{m\omega/\hbar}$')

    ax.set_xlabel(r'$r_{12}^2$')
    ax.set_ylabel(r'Density')

    ax.text(x=3,y=0.2,s=r'$\lambda=1$')
    ax.text(x=5,y=0.15,s=r'$\lambda=2$')
    ax.text(x=7,y=0.12,s=r'$\lambda=8$')
    ax.text(x=13,y=0.09,s=r'$\lambda=20$')
    ax.text(x=19,y=0.07,s=r'$\lambda=30$')
    plt.show()

main()

