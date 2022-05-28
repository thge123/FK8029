from PLOTTING import *
from FILES import *
from scipy.integrate import odeint

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

    def TOV_Newt(y,x):
        p,m = y
        dydt = [-E*m/(x**2),
                x**2*E]
        return dydt

    sol = odeint(TOV_Newt,init,x)
    
    p_Newt = []
    m_Newt = []
    for i in sol:
        p_Newt.append(i[0])
        m_Newt.append(i[1])

    return p,m,p_Newt,m_Newt


def main():

    X = get_data('pm.dat')
    x,p,p_Newt,m,m_Newt = [],[],[],[],[]
    
    for i in X:
        x.append(i[0])
        p.append(i[1])
        p_Newt.append(i[2])
        m.append(i[3])
        m_Newt.append(i[4])

    fig_p = plt.figure()
    ax_p  = fig_p.add_subplot(111)
    ax_p.scatter(x,p,facecolor='none',edgecolor='k')
    ax_p.scatter(x,p_Newt,facecolor='none',edgecolor='r')
    ax_p.set_xlabel('$x$')
    ax_p.set_ylabel('$p$')

    fig_m = plt.figure()
    ax_m  = fig_m.add_subplot(111)
    ax_m.scatter(x,m,facecolor='none',edgecolor='k')
    ax_m.scatter(x,m_Newt,facecolor='none',edgecolor='r')
    ax_m.set_xlabel('$x$')
    ax_m.set_ylabel('$m$')

    # Test for const. density
    #p,m,p_Newt,m_Newt = ode_solver(x,[p[0],m[0]])
    #ax_p.plot(x,p,ls='--')
    #ax_p.plot(x,p_Newt,ls='--')
    #ax_m.plot(x,m,ls='--')
    #ax_m.plot(x,m_Newt,ls='--')

    plt.show()

main()
