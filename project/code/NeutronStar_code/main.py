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

        p = [j*1.285 for j in p]  # GeV/fm^3
        x = [j*6.84  for j in x]  # km
        m = [j*4.63  for j in m]  # km
        E = [j*1.285 for j in E]  # GeV/fm^3

        ax_pr.plot(x,p,ls=lss[filename],lw=4)
        ax_mr.plot(x,m,ls=lss[filename],lw=4)
        ax_E.plot(p,E,lw=4,c='k')


    ax_E.set_xlabel('$P$ [Gev/fm$^3$]')
    ax_E.set_ylabel('$\mathcal{E}$ [Gev/fm$^3$]')
    ax_pr.set_xlabel(r'$r$ [km]')
    ax_pr.set_ylabel('$P$ [Gev/fm$^3$]')
    ax_mr.set_xlabel(r'$r$ [km]')
    ax_mr.set_ylabel(r'$M$ [M$_\odot$]')
        

    plt.show()

def MRP2():
    
    fig1 = plt.figure()
    fig2 = plt.figure()
    ax1  = fig1.add_subplot(111)
    ax2  = fig2.add_subplot(111)
    labels = {'GRxm.dat': (5,0.4,0.4,0.5,'Relativistic'),
              'NEWTxm.dat': (5,1.3,0.4,1.2,'Newtonian')}

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
    ax2.set_xlabel(r'$P_c$ [GeV/fm$^3$]')
    ax2.set_ylabel(r'$M(R)$ [M$_\odot$]')
    ax2.set_xscale('log')


    plt.show()

    
main()
MRP2()

