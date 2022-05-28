from PLOTTING import *
from FILES import *
from scipy.integrate import odeint
from numpy import pi,sqrt,log

def main():

    X = get_data('pm.dat')
    x,p,m,E = [],[],[],[]
    
    for i in X[:-1]:
        x.append(i[0])
        p.append(i[1])
        m.append(i[2])
        E.append(i[3])

    fig_p = plt.figure()
    ax_p  = fig_p.add_subplot(111)
    ax_p.scatter(x,p,s=0.1,c='k')
    ax_p.set_xlabel('$x$')
    ax_p.set_ylabel('$p$')

    fig_m = plt.figure()
    ax_m  = fig_m.add_subplot(111)
    ax_m.scatter(x,m,s=1,c='k')
    ax_m.set_xlabel('$x$')
    ax_m.set_ylabel('$m$')

    fig_E = plt.figure()
    ax_E  = fig_E.add_subplot(111)
    ax_E.scatter(p,E,s=1,c='k')
    ax_E.set_xlabel('$p$')
    ax_E.set_ylabel('$\epsilon$')

    print("x = ", x[-1])
    print("M = ", m[-1])

    plt.show()

def MR():
    
    X = get_data('xm.dat')
    pc,x0,m = [],[],[]
    
    for i in X:
        pc.append(i[0])
        x0.append(i[1])
        m.append(i[2])
        print(i)

    pc = [j*0.11239 for j in pc]  # GeV/fm^3
    m =  [j*4.63  for j in m]   # Solar masses
    x0 = [j*6.84  for j in x0]  # km

   # x0_m = zip(x0,m)
   # x0_m = sorted(x0_m)

   # pc_m = zip(pc,m)
   # pc_m = sorted(pc_m) 
    
    fig = plt.figure()
    ax1  = fig.add_subplot(111)
    ax1.scatter(x0,m,
             color='k',label='$R$',s=1)
    ax1.set_xlabel(r'$R$ [km]')
    ax1.set_ylabel(r'$M$ [M$_\odot$]')
    
    #ax2 = fig.add_subplot(212)
    #ax2.plot([j[0] for j in pc_m],
    #         [j[1] for j in pc_m],color='k')
    ##ax2.scatter([j[0] for j in pc_m],
    ##            [j[1] for j in pc_m],color='k')
    #ax2.set_xlabel(r'$P_c$ [GeV/fm$^3$]')
    #ax2.set_ylabel(r'$M$ [M$_\odot$]')


    plt.show()

def MRP():
    
    X = get_data('xm.dat')
    pc,x0,m = [],[],[]
    
    for i in X:
        pc.append(i[0])
        x0.append(i[1])
        m.append(i[2])
        print(i)

    pc = [j*1.285 for j in pc]  # GeV/fm^3
    m =  [j*4.63  for j in m]   # Solar masses
    x0 = [j*6.84  for j in x0]  # km

    x0_m = zip(x0,m)
    x0_m = sorted(x0_m)

    pc_m = zip(pc,m)
    pc_m = sorted(pc_m) 
    
    fig = plt.figure()
    ax1  = fig.add_subplot(111)
    ax1.plot([j[0] for j in x0_m],
             [j[1] for j in x0_m],
             color='k',label='$R$',lw=3)
    ax1.plot([1e3,1e4],[1,1],
             color='r',ls='--',
             label='$P_c$',lw=3)
    ax2 = ax1.twiny()
    ax2.plot([j[0] for j in pc_m],
             [j[1] for j in pc_m],
             color='r',label='$P_c$',
             ls='--',lw=3)
    ax1.set_xlabel(r'$R$ [km]')
    ax1.set_ylabel(r'$M$ [M$_\odot$]')
    ax2.set_xlabel(r'$P_c$ [GeV/fm$^3$]',labelpad=15)
    ax2.set_xscale('log')
    ax1.set_xlim(5,20)
    
    max_mass = max([j[1] for j in x0_m])
    ax1.plot([5,20],[max_mass,max_mass],ls='dashdot',
             label='$M=0.71$ M$_\odot$')
    ax1.legend(framealpha=0)

    #ax2 = fig.add_subplot(212)
    #ax2.plot([j[0] for j in pc_m],
    #         [j[1] for j in pc_m],color='k')
    ##ax2.scatter([j[0] for j in pc_m],
    ##            [j[1] for j in pc_m],color='k')
    #ax2.set_xlabel(r'$P_c$ [GeV/fm$^3$]')
    #ax2.set_ylabel(r'$M$ [M$_\odot$]')


    plt.show()
    
def MRP2():
    
    X = get_data('xm.dat')
    pc,x0,m = [],[],[]
    
    for i in X:
        pc.append(i[0])
        x0.append(i[1])
        m.append(i[2])
        print(i)

    pc = [j*1.285 for j in pc]  # GeV/fm^3
    m =  [j*4.63  for j in m]   # Solar masses
    x0 = [j*6.84  for j in x0]  # km

    fig = plt.figure()
    ax2  = fig.add_subplot(121)
    ax1 = fig.add_subplot(122)
    ax1.scatter(x0,m,s=1,c='k')
    ax2.scatter(pc,m,s=1,c='k')
    ax1.set_xlabel(r'$R$ [km]')
    ax2.set_ylabel(r'$M$ [M$_\odot$]')
    ax2.set_xlabel(r'$P_c$ [GeV/fm$^3$]')
    ax2.set_xscale('log')

    ax1.plot([0,25],[0.71,0.71],ls='--',c='k')
    ax2.plot([1e-4,1e4],[0.71,0.71],ls='--',c='k',
             label='$M = 0.71$ M$_\odot$')
    ax1.set_xlim(2,23)
    ax2.set_xlim(5e-4,5e2)
    ax2.legend(framealpha=0,loc='lower center')
    
    #ax2 = fig.add_subplot(212)
    #ax2.plot([j[0] for j in pc_m],
    #         [j[1] for j in pc_m],color='k')
    ##ax2.scatter([j[0] for j in pc_m],
    ##            [j[1] for j in pc_m],color='k')
    #ax2.set_xlabel(r'$P_c$ [GeV/fm$^3$]')
    #ax2.set_ylabel(r'$M$ [M$_\odot$]')


    plt.show()

    
MRP2()


