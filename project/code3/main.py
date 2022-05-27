from PLOTTING import *
from FILES import *
from numpy import pi,sqrt,log

def main():

    fig_p = plt.figure()
    ax_p  = fig_p.add_subplot(111)
    fig_m = plt.figure()
    ax_m  = fig_m.add_subplot(111)

    for filename in ['GRpm.dat','NEWTpm.dat']:
        X = get_data(filename)
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

        pc = [j*1.285 for j in pc]  # GeV/fm^3
        m =  [j*4.63  for j in m]   # Solar masses
        x0 = [j*6.84  for j in x0]  # km

        ax1.scatter(x0,m,s=1)
        ax2.scatter(pc,m,s=1)

    #ax2.legend(framealpha=0,loc='lower center')
    ax1.set_xlabel(r'$R$ [km]')
    ax2.set_ylabel(r'$M$ [M$_\odot$]')
    ax2.set_xlabel(r'$P_c$ [GeV/fm$^3$]')
    ax2.set_xscale('log')
    plt.show()

    
MRP2()


