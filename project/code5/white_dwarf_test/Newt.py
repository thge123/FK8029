from PLOTTING import *
from FILES import *
from scipy.integrate import odeint
from numpy import pi,sqrt,log

def ode_solver(x,init):
    
    E = 1e-3
    def TOV(y,x):
        p,m = y
        dydt = [-(E+p)*(m+x**3*p)/(x*(x-2*m)),
                x**2*E]
        return dydt

    sol = odeint(TOV,init,x)
    
    p = []
    m = []
    for i in sol:
        p.append(i[0])
        m.append(i[1])

    return p,m


def main():

    X = get_data('pm.dat')
    x,p,m,E = [],[],[],[]
    x_Newt,p_Newt,m_Newt,E_Newt = [],[],[],[]
    
    for i in X[:-1]:
        x.append(6.84*i[0])
        p.append(1.285*i[1])
        m.append(4.63*i[2])
        E.append(1.285*i[3])

    X = get_data('pm_Newt.dat')
    for i in X[:-1]:
        x_Newt.append(6.84*i[0])
        p_Newt.append(1.285*i[1])
        m_Newt.append(4.63*i[2])
        E_Newt.append(1.285*i[3])

    fig_p = plt.figure()
    ax_p  = fig_p.add_subplot(211)
    ax_p.plot(x,p,lw=3,c='k',label='GR')
    ax_p.plot(x_Newt,p_Newt,lw=3,c='r',label='Newton',ls='--')
    #ax_p.set_xlabel('$r$ [km]')
    ax_p.set_ylabel('$P$ [GeV/fm$^3$]')

    ax_m  = fig_p.add_subplot(212)
    ax_m.plot(x,m,lw=3,c='k',label='GR')
    ax_m.plot(x_Newt,m_Newt,lw=3,c='r',label='Newton',ls='--')
    ax_m.set_xlabel('$r$ [km]')
    ax_m.set_ylabel('$M$ [M$_\odot$]')

    fig_E = plt.figure()
    ax_E  = fig_E.add_subplot(111)
    ax_E.plot(p,E,lw=3,c='k',label='GR')
    ax_E.plot(p_Newt,E_Newt,lw=3,c='r',label='Newton',ls='--')
    ax_E.set_xlabel('$P$ [GeV/fm$^3$]')
    ax_E.set_ylabel('$\epsilon$')

    # Test for const. density
    #p,m = ode_solver(x,[p[0],m[0]])
    #ax_p.plot(x,p,ls='--')
    #ax_m.plot(x,m,ls='--')
    
    ax_p.legend(framealpha=0)
    #ax_m.legend(framealpha=0)
    ax_E.legend(framealpha=0)
    
    print("x = ", x[-1])
    print("M = ", m[-1])

    plt.show()

def EOS_test():
        
    def bissection(f,a,b):
        max_iters = 1000
        iters = 0
        while iters<max_iters:
            c = (a+b)/2
            if abs(f(c))<1e-6:
                return c
            elif f(c)*f(a) < 0:
                b = c
            else:
                a = c
            iters+=1

    def f(x):
        a = 2.1002e-1
        xF = a*(3*pi**2*x)**(1/3)
        yF = sqrt(1+xF**2)
        p  = 2*xF**3*yF/3
        p += -xF*yF
        p += log(xF+yF)
        return p

    x = [j/100 for j in range(100001)]
    plt.plot(x,[f(j) for j in x])
    plt.show()

def MR_test():

    x = [1.7067,1.77,1.8337,
         1.8955,1.9763,2.0925,
         2.297,2.5088,2.6695,
         3.5939,1.6478,1.6755,
         1.623,1.6006,1.5612,
         1.3082,1.2448,1.2018,
         1.1798]
    m = [0.143264,0.139357,0.136506,
         0.13286,0.127827,0.120248,
         0.106555,0.0927495,0.0829332,
         0.0417883,0.145955,0.144729,
         0.14699,0.147873,0.149286,
         0.153224,0.152714,0.151951,
         0.151084]

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.scatter(x,m,facecolor='none',edgecolor='k')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$m$')

    plt.show()

def MR():
    
    X = get_data('xm_Newt.dat')
    pc,x0,m = [],[],[]
    
    for i in X:
        pc.append(i[0])
        x0.append(i[1])
        m.append(i[2])
        print(i)

    pc = [j*1.285 for j in pc]  # GeV/fm^3
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
    
    X = get_data('xm_Newt.dat')
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
    
    X = get_data('xm_Newt.dat')
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
    ax1  = fig.add_subplot(121)
    ax1.scatter(x0,m,s=1,c='k')
    ax2 = fig.add_subplot(122)
    ax2.scatter(pc,m,s=1,c='k')
    ax1.set_xlabel(r'$R$ [km]')
    ax1.set_ylabel(r'$M$ [M$_\odot$]')
    ax2.set_xlabel(r'$P_c$ [GeV/fm$^3$]')
    ax2.set_xscale('log')

    #ax1.plot([0,25],[0.71,0.71],ls='--',c='k')
    #ax2.plot([1e-4,1e4],[0.71,0.71],ls='--',c='k',
    #         label='$M = 0.71$ M$_\odot$')
    #ax1.set_xlim(2,23)
    #ax2.set_xlim(5e-4,5e2)
    ax2.legend(framealpha=0,loc='lower center')
    
    #ax2 = fig.add_subplot(212)
    #ax2.plot([j[0] for j in pc_m],
    #         [j[1] for j in pc_m],color='k')
    ##ax2.scatter([j[0] for j in pc_m],
    ##            [j[1] for j in pc_m],color='k')
    #ax2.set_xlabel(r'$P_c$ [GeV/fm$^3$]')
    #ax2.set_ylabel(r'$M$ [M$_\odot$]')


    plt.show()

    
main()
