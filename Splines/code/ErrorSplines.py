from PLOTTING import *
from BSpline import *

def main():

    
    fig = plt.figure()
    ax1  = fig.add_subplot(311)
    ax2  = fig.add_subplot(312)
    ax3  = fig.add_subplot(313)
    axs  = [ax1,ax2,ax3]
    #i = int(input("Write sigma [1,2,3]: ")) - 1 
    #N = int(input("Write N: "))
    X   = array([j/100 for j in range(101)])
    Ns = [1+j for j in range(50)]
    labels = {0: 'Sphere',1: 'Shell', 2: 'Hydrogen'}
    for i in range(3):
        err = []
        exact_sol = [sigma1_exact,sigma2_exact,sigma3_exact][i]

        if i == 1:
            r = float(input("Write r: "))
        elif i == 2:
            r = float(input("Write alpha: "))
        else:
            r = 1

        Y = exact_sol(X,r)
        for N in Ns:
            print("N: ", N)
            spl = numsol(N,i,r)
            err.append(MAX(ABS(spl(X)-Y)))
            #ax.plot(X,spl(X)+X)

        axs[i].scatter([N+3 for N in Ns],err,s=100,edgecolor='k',facecolor='none')
        axs[i].set_ylabel('Error')
    axs[0].text(s=labels[0],x=15,y=1.3e-9)
    axs[1].text(s=labels[1],x=15,y=2)
    axs[2].text(s=labels[2],x=15,y=2)
    axs[2].set_xlabel('Number of splines')
    axs[1].set_yticks([0,1,2])
    axs[1].set_ylim(-0.08,2.60)
    axs[2].set_ylim(-0.08,2.60)
    axs[2].set_yticks([0,1,2])

    plt.show()

main()
