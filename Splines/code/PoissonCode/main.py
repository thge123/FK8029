from PLOTTING import *
from BSpline3 import *

def main():

    
    i = int(input("Write sigma [1,2,3,4,5]: ")) - 1 
    N = int(input("Write number of inner points: "))
    X   = array([j/1000 for j in range(1001)])
    exact_sol = [sigma1_exact,sigma2_exact,sigma3_exact,sigma4_exact,sigma5_exact][i]

    if i == 1:
        r = float(input("Write r: "))
    elif i in [2,3,4]:
        r = float(input("Write alpha: "))
    else:
        r = 1

    spl = numsol(N,i,r)
    # Potential times distance
    u   = exact_sol(X,r) + X 
    uu  = spl(X) + X

    # Potentials
    V   = u[1:]/X[1:]
    VV  = uu[1:]/X[1:]

    # Electric fields
    E  = array([0] + [-(V[j+1]-V[j])/(X[1]-X[0]) for j in range(len(V)-1)])
    #EE = (-1 + VV - spl.derivative()(X)[1:])/X[1:]
    EE = array([0] + [-(VV[j+1]-VV[j])/(X[1]-X[0]) for j in range(len(V)-1)])
    #

    fig = plt.figure()
    ax  = fig.add_subplot(211)
    axx  = fig.add_subplot(212)

    #ax.plot(X,u,label='Potential times distance')
    axx.plot(X[1:],V,lw=2,c='b')
    axx.plot(X[:-1],E,c='k',lw=2)
    if i==2:
        axx.text(s='Exact solution',x=0.7,y=75)
    else:
        axx.text(s='Exact solution',x=0,y=0.75)
    axx.grid(linestyle='--')

    xlabels = {0:r'$r/R$',
               1:r'$r/R_2$', 
               2:r'$r/({}a_0)$'.format(int(r)),
               3:r'$r/({}a_0)$'.format(int(r)),
               4:r'$r/({}a_0)$'.format(int(r))}
    axx.set_xlabel(xlabels[i])

    legends = {0:r'$\varphi/(q/4\pi\epsilon_0 R)$',
               1:r'$\varphi/(q/4\pi\epsilon_0 R_2)$',  
               2:r'$\varphi/(q/4\pi\epsilon_0 ({}a_0))$'.format(int(r)),
               3:r'$\varphi/(q/4\pi\epsilon_0 ({}a_0))$'.format(int(r)),
               4:r'$\varphi/(q/4\pi\epsilon_0 ({}a_0))$'.format(int(r)),
               5:r'$E/(q/4\pi\epsilon_0 R^2)$',
               6:r'$E/(q/4\pi\epsilon_0 R_2^2)$',  
               7:r'$E/(q/4\pi\epsilon_0 ({}a_0)^2)$'.format(int(r)),
               8:r'$E/(q/4\pi\epsilon_0 ({}a_0)^2)$'.format(int(r)),
               9:r'$E/(q/4\pi\epsilon_0 ({}a_0)^2)$'.format(int(r))}

    #ax.plot(X,uu,label='Potential times distance')
    ax.plot(X[1:],VV,label=legends[i],lw=2,c='b',ls='--')
    ax.plot(X[:-1],EE,label=legends[i+5],lw=2,ls='dashdot',c='k')
    if i==2:
        ax.text(s='Numerical solution',x=0.7,y=75)
        ax.text(s='1s state',x=0.2,y=125)
    else:
        ax.text(s='Numerical solution',x=0,y=0.75)
    ax.legend(framealpha=0,loc='upper right')
    ax.grid(linestyle='--')
    plt.show()

main()
