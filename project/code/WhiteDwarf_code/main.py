from PLOTTING import *
from FILES import *
from numpy import pi,sqrt,log

def main():

    fig_pmr = plt.figure()
    ax_pr  = fig_pmr.add_subplot(111)
    ax_mr = ax_pr.twinx()
    fig_E = plt.figure()
    ax_E  = fig_E.add_subplot(111)

    lss = {'GRpm.dat':'-','NEWTpm.dat':'--'}
    for filename in ['GRpm.dat','NEWTpm.dat']:
        X = get_data(filename)
        x,p,m,E = [],[],[],[]
        
        for i in X[:-1]:
            x.append(i[0])
            p.append(i[1])
            m.append(i[2])
            E.append(i[3])

        p = [j*8.874 for j in p]  # meV/fm^3
        m = [j for j in m]   # Solar masses
        x = [j*2711  for j in x]  # km
        E = [j*16.294 for j in E]  # eV/fm^3

        ax_pr.plot(x,p,ls=lss[filename],lw=4)
        ax_mr.plot(x,m,ls=lss[filename],lw=4)
        ax_E.plot(p,E,lw=4,c='k')


    ax_E.set_xlabel('$P$ [mev/fm$^3$]')
    ax_E.set_ylabel('$\mathcal{E}$ [ev/fm$^3$]')
    ax_pr.set_xlabel(r'$r$ [km]')
    ax_pr.set_ylabel('$P$ [mev/fm$^3$]')
    ax_mr.set_xlabel(r'$r$ [km]')
    ax_mr.set_ylabel(r'$M$ [M$_\odot$]')
        

    plt.show()

def MRP2():
    
    fig1 = plt.figure()
    fig2 = plt.figure()
    ax1  = fig1.add_subplot(111)
    ax2  = fig2.add_subplot(111)
    labels = {'GRxm.dat': (500,1.3,5000,1.46,'Relativistic'),
              'NEWTxm.dat': (500,1.46,5000,1.3,'Newtonian')}

    for filename in ['GRxm.dat','NEWTxm.dat']:
        X = get_data(filename)
        pc,x0,m = [],[],[]
        
        for i in X:
            pc.append(i[0])
            x0.append(i[1])
            m.append(i[2])
            print(i)

        pc = [j*8.874 for j in pc]  # GeV/fm^3
        m =  [j for j in m]   # Solar masses
        x0 = [j*2711  for j in x0]  # km

        ax1.scatter(x0,m,s=10,label=labels[filename])
        ax2.scatter(pc,m,s=10,label=labels[filename])

        ax1.text(x=labels[filename][0],
                 y=labels[filename][1],
                 s=labels[filename][4])
        ax2.text(x=labels[filename][2],
                 y=labels[filename][3],
                 s=labels[filename][4])

    ax1.set_xlabel(r'$R$ [km]')
    ax1.set_ylabel(r'$M(R)$ [M$_\odot$]')
    ax2.set_xlabel(r'$P_c$ [meV/fm$^3$]')
    ax2.set_ylabel(r'$M(R)$ [M$_\odot$]')
    ax2.set_xscale('log')
    ax1.set_ylim(0.5,1.55)
    ax2.set_ylim(0.5,1.55)


    plt.show()

    
main()
MRP2()
