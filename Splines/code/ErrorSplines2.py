from PLOTTING import *
from BSpline2 import *

def main():

    
    fig = plt.figure()
    ax1  = fig.add_subplot(211)
    ax2  = fig.add_subplot(212)
    axs  = [ax1,ax2]


    n1   = [24,29,34,39,44,49,54,59,64,69,74,79,84]
    err1 = [0.024,0.021,0.024,0.0189,0.014,0.011,0.009,0.0076,0.0064,0.0055,0.0047,0.0041,0.0036]
    n2   = [24,29,34,39,44,49,54,59,64,69,74,79,84]
    err2 = [0.000189,0.0069,0.000107,0.0054,6.45e-5,0.0042,4.48e-5,0.0036,3.19e-5,0.003,2.44e-5,0.0027,2.07e-5]
    ax2.scatter(n1,err1,edgecolor='k',facecolor='none',s=100,label='Hydrogen charge distribution')
    ax1.scatter(n2,err2,edgecolor='k',facecolor='none',s=100,label='Shell charge distribution')

    for i in axs:
        i.set_ylabel('Error')
        i.grid(ls='dashdot',axis='y')
    axs[0].set_yscale('log')
    axs[1].set_yscale('log')
    axs[1].set_xlabel('Number of splines')
    #axs[0].set_ylim(0.6e-5,1.4e-2)
    axs[1].set_ylim(0.6e-3,1.4e-1)
    axs[0].set_yticks([1e-4,1e-3,1e-2])
    axs[1].set_yticks([1e-3,1e-2,1e-1])
    axs[0].text(s='Shell charge distribution',x=22.7,y=2.2e-5)
    axs[1].text(s='1s state charge distribution',x=22.7,y=0.8e-3)
    
#    axs[1].set_ylim(-0.08,2.60)
#    axs[2].set_ylim(-0.08,2.60)
#    axs[2].set_yticks([0,1,2])

    plt.show()

main()
