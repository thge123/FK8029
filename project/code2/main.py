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
    x,p,m = [],[],[]
    
    for i in X:
        x.append(i[0])
        p.append(i[1])
        m.append(i[2])

    fig_p = plt.figure()
    ax_p  = fig_p.add_subplot(111)
    ax_p.scatter(x,p,facecolor='none',edgecolor='k')
    ax_p.set_xlabel('$x$')
    ax_p.set_ylabel('$p$')

    fig_m = plt.figure()
    ax_m  = fig_m.add_subplot(111)
    ax_m.scatter(x,m,facecolor='none',edgecolor='k')
    ax_m.set_xlabel('$x$')
    ax_m.set_ylabel('$m$')

    # Test for const. density
    p,m = ode_solver(x,[p[0],m[0]])
    ax_p.plot(x,p,ls='--')
    ax_m.plot(x,m,ls='--')

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
        p  = 2*xF**3*sqrt(1+xF)/3
        p += -xF*sqrt(1+xF)
        p += log(xF+sqrt(1+xF))
        return p

    print(bissection(lambda x: f(x)-1,0,100))

    x = [j/100 for j in range(1001)]
    plt.plot(x,[f(j)-1 for j in x])
    plt.show()

    

main()
