from PLOTTING import *
from FILES import *
from numpy import pi,sqrt,log

def main():

    fig_p = plt.figure()
    ax_p  = fig_p.add_subplot(111)
    fig_m = plt.figure()
    ax_m  = fig_m.add_subplot(111)
    fig_E = plt.figure()
    ax_E  = fig_E.add_subplot(111)

    for filename in ['GRpm.dat','NEWTpm.dat']:
        X = get_data('GRpm.dat')
        x,p,m,E = [],[],[],[]
        
        for i in X[:-1:100]:
            x.append(i[0])
            p.append(i[1])
            m.append(i[2])
            E.append(i[3])

        ax_p.plot(x,p,c='k')
        ax_p.set_xlabel('$x$')
        ax_p.set_ylabel('$p$')

        ax_m.plot(x,m,c='k')
        ax_m.set_xlabel('$x$')
        ax_m.set_ylabel('$m$')

        ax_E.plot(p,E,c='k')
        ax_E.set_xlabel('$p$')
        ax_E.set_ylabel('$\mathcal{E}$')
        

        print("x = ", x[-1])
        print("M = ", m[-1])

    plt.show()

def MRP2():
    
    fig = plt.figure()
    ax2  = fig.add_subplot(121)
    ax1 = fig.add_subplot(122)

    for filename in ['GRxm.dat','NEWTxm.dat']:
        X = get_data(filename)
        pc,x0,m = [],[],[]
        
        for i in X:
            pc.append(i[0])
            x0.append(i[1])
            m.append(i[2])
            print(i)

        pc = [j*8.874 for j in pc]  # meV/fm^3
        m =  [j  for j in m]   # Solar masses
        x0 = [j*0.003897  for j in x0]  # solar radii

        ax1.scatter(x0,m,s=1)
        ax2.scatter(pc,m,s=1)

    #ax2.legend(framealpha=0,loc='lower center')
    ax1.set_xlabel(r'$R$ [R$_\odot$]')
    ax2.set_ylabel(r'$M$ [M$_\odot$]')
    ax2.set_xlabel(r'$P_c$ [meV/fm$^3$]')
    ax2.set_xscale('log')
    plt.show()

    
main()
MRP2()

