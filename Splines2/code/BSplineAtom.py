from numpy import array,zeros,pi,log,argmax
from Numerics import *
from scipy.interpolate import BSpline

Exact = {'u00': lambda r,D: 2*D**1.5*r*exp(-D*r)                  ,
         'u10': lambda r,D: D**1.5*r*(2-r*D)*exp(-D*r/2)/2**1.5   ,
         'u11': lambda r,D: D**2.5*r**2*exp(-D*r/2)/(2*6**0.5)    ,
         'u20': lambda r,D: 2*D**1.5*r*(1-2*r*D/3+2*(r*D)**2/27)
                            *exp(-D*r/3)/(3*3**0.5)               ,
         'u21': lambda r,D: 8*D**2.5*(1-r*D/6)*r**2
                            *exp(-r*D/3)/(27*6**0.5)              ,
         'u22': lambda r,D: 4*D**3.5*r**3*exp(-r*D/3)/(81*30**0.5),
         'u30': lambda r,D: D**1.5*r
                            *(1-3*D*r/4+(D*r)**2/8-(D*r)**3/192)
                            *exp(-r*D/4)/4                        ,
         'u31': lambda r,D: 5*D**2.5*r**2*(1-D*r/4+(D*r)**2/80)
                            *exp(-r*D/4)/(16*15**0.5)             ,
         'u32': lambda r,D: D**3.5*r**3*(1-D*r/12)
                            *exp(-D*r/4)/(64*5**0.5)              ,
         'u33': lambda r,D: D**4.5*r**4*exp(-r*D/4)/(768*35**0.5)
        }

def UniformSphere(x,R0):
    if x<R0:
        return -(3-(x/R0)**2)/(2*R0)
    else:
        return -1/x

def getA(t,K,R0,ell,beta,quadN=20,quadDeg=8):

    ''' The "hamiltonian matrix". The matrix elements are
        given by 
            A_ij = -beta^2*( -int(B'_ik*B'_jk) 
                             -l(l+1)*int(B_ik*B_jk/x^2)
                             +(2/beta)*int(B_ik*B_jk/x) )
    '''

    n = len(t)-K-1
    A = zeros((n-2,n-2))

    for i in range(1,n-1):
        c_i    = zeros(n)
        c_i[i] = 1
        spl_i  = BSpline(t,c_i,K,extrapolate=False)
        dspl_i = spl_i.derivative()
        if i+K+1>n-1:
            m = n-1
        else:
            m = i+K+1
        for j in range(i,m):
            # B_(>i+K) is non-overlapping with B_i

            c_j    = zeros(n)
            c_j[j] = 1
            spl_j  = BSpline(t,c_j,K,extrapolate=False)
            dspl_j = spl_j.derivative()

            f = lambda x: spl_i(x)*spl_j(x)*UniformSphere(x,R0)
            A[i-1,j-1]  = -2*Quadrature(f,t[i],t[i+K+1],N=quadN,deg=quadDeg)/beta
            #A[i-1,j-1]  = -2*Simpson(f,t[i],t[i+K+1])/beta

            f = lambda x: spl_i(x)*spl_j(x)/x**2
            A[i-1,j-1] += -ell*(ell+1)*Quadrature(f,t[j],t[i+K+1],N=quadN,deg=quadDeg)
            #A[i-1,j-1] += -ell*(ell+1)*Simpson(f,t[j],t[i+K+1])


            f = lambda x: dspl_i(x)*dspl_j(x)
            A[i-1,j-1] += -Quadrature(f,t[i],t[i+K+1],N=quadN,deg=quadDeg)
            #A[i-1,j-1] += -Simpson(f,t[i],t[i+K+1])

            xx = array([j/100 for j in range(1001)])
    
            A[i-1,j-1] *= -100*beta**2
            A[j-1,i-1]  = A[i-1,j-1]

    return A


def getB(t,K,quadN=20,quadDeg=8):
    n = len(t)-K-1
    B = zeros((n-2,n-2))

    #fig = plt.figure()
    #ax  = fig.add_subplot(111)
    #xx  = array([j/1000 for j in range(1001)])
    for i in range(1,n-1):
        c_i    = zeros(n)
        c_i[i] = 1
        spl_i  = BSpline(t,c_i,K,extrapolate=False)
        #ax.plot(xx,spl_i(xx))
        if i+K+1>n-1:
            m = n-1
        else:
            m = i+K+1
        for j in range(i,m):
            # B_(>i+K) is non-overlapping with B_i
            c_j    = zeros(n)
            c_j[j] = 1
            spl_j  = BSpline(t,c_j,K,extrapolate=False)
            spl_ij = lambda x: spl_i(x)*spl_j(x)
            B[i-1,j-1] = Quadrature(spl_ij,t[i],t[i+K+1],N=quadN,deg=quadDeg)
            #B[i-1,j-1] = Simpson(spl_ij,t[i],t[i+K+1])
            B[j-1,i-1] = B[i-1,j-1]

    #c_i = zeros(n); c_i[0]=1; c_i[n-1]=1
    #spl = BSpline(t,c_i,K,extrapolate=False)
    #ax.plot(xx,spl(xx))
    #ax.plot
    #ax.scatter(t,[0 for i in t],c='k')
    #plt.show()

    return B

def knots1(X,a):
    return [1-j**a for j in X]

def knots2(X,D):
    return [exp(-D*j) for j in X]

def knots3(X,a):
    return [2**(j**a)-1 for j in X]

def numsol(dct,PlotSplines=False):

    E    = dct['E']
    Z    = dct['Z']
    ell  = dct['ell']
    R0   = dct['R0']
    N    = dct['N']
    mume = dct['mass_ratio']
    D    = dct['D']
    beta = 1/(Z*mume*D)
    quadN = dct['quadN']
    quadDeg = dct['quadDeg']
    
    # Construct inner points
    X = [j/N for j in range(1,N)]
    #X = knots1(X,0.7)
    #X = knots2(X,60*D)
    X = knots3(X,2)
    X.sort()

    # Degree K and knots t.
    K = 3
    t = array((K+1)*[0] + X + (K+1)*[1])

    # Total number of splines
    n = len(t)-K-1
    print("Number of splines: ", n)

    # We wish to solve eigenvalue problem
    #           Ac = lam*Bc.
    # The method uses power iteration.

    A = getA(t,K,R0,ell,beta,quadN=quadN,quadDeg=quadDeg)
    B = getB(t,K,quadN,quadDeg=quadDeg)

    x0 = rand(len(A))
    eig = InvPower(A,B,x0,shift=E)
    print("Energy: ", eig[0])
    if 'eigs' in dct:
        dct['eigs'].append(eig[0])
    else:
        dct['eigs'] = [eig[0]]
    sol = array([0] + [j for j in eig[1]] + [0])
    u_spl = BSpline(t,sol,K,extrapolate=False)

    # normalizing
    u_splSqrd = lambda x: u_spl(x)**2
    A   = Quadrature(u_splSqrd,0,1,N=quadN,deg=quadDeg)
    #A   = Simpson(u_splSqrd,0,1)  # normalizing
    u   = lambda x: u_spl(x)/sqrt(A)

    if PlotSplines:
        ### plot splines
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        xx = array([j/1000 for j in range(1001)])
        ticks = 5 
        xticks      = [int(10*j/ticks)/10 for j in range(ticks+1)]
        xticklabels = [int(j*dct['D']/ticks) for j in range(ticks+1)]
        ax.set_xticks(xticks,xticklabels)
        ax.set_yticks([-2,-1,0,1,2])
        #ax.set_ylabel('$u_{n\ell}$')
        ax.set_ylabel('$u_{30}(r/a_0)$')
        ax.set_xlabel('$r/a_0$')

        for i in range(n):
            c_i    = zeros(n)
            c_i[i] = 1
            spl_i  = BSpline(t,c_i,K,extrapolate=False)
            ax.plot(xx,spl_i(xx),lw=1)

        ax.plot(xx,u(xx),lw=4,c='k',label='Numerical')
        ax.plot(xx,Exact['u20'](xx,dct['D']),lw=4,c='purple',ls='--',
                label = 'Exact')
        ax.scatter(t,[0 for i in t],c='k')
        ax.legend(framealpha=0)
        plt.show()
        ###

    return u

def plot_radials(Atom,exact=False):

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    axs = [ax1,ax2,ax3,ax4]
    lss = ['-','--','dashdot','dotted']
    xx = array([j/1000 for j in range(1001)])

    E   = [-99,-25,-11,-6]
    ell = [0,1,2,3]
    D   = Atom['D']
    D   = [D*(j+1) for j in range(len(E))]
    for i in range(len(E)):
        Atom['E'] = E[i]
        Atom['D'] = D[i]
        ticks       = 5
        xticks      = [int(10*j/ticks)/10 for j in range(ticks+1)]
        xticklabels = [int(j*D[i]/ticks) for j in range(ticks+1)]
        axs[i].set_xticks(xticks,xticklabels)
        axs[i].set_title('$n={}$'.format(i+1),fontsize=18)
        for j in range(i+1):
            Atom['ell'] = ell[j]
            if not exact:
                u = numsol(Atom)
                axs[i].plot(xx,u(xx)**2/D[i],lw=3,label=r'$\ell = {}$'.format(ell[j]),ls=lss[j])
            else:
                name = 'u{}{}'.format(i,j)
                axs[i].plot(xx,Exact[name](xx,Atom['D'])**2/D[i],
                            lw=2,label=r'$\ell = {}$'.format(ell[j]),ls=lss[j])

            name = 'u{}{}'.format(i,j)
            print(name)
            q = argmax(u(xx))
            print(abs(abs(u(xx[q]))-abs(Exact[name](xx[q],Atom['D']))))

        axs[i].legend(framealpha=0)
        axs[i].grid(linestyle='dotted',alpha=0.5)

    axs[2].set_ylabel('$u^2_{n,\ell}$')
    axs[2].set_xlabel('$r/a_0$')
    
    plt.show()
    return 0

def u20_test(Atom):

    xx = array([j/1000 for j in range(1001)])
    u = numsol(Atom)
    u_exact = Exact['u20'](xx,Atom['D'])
    j = argmax(u_exact)
    print('error estimated in: ', xx[j]*Atom['D'], ' a0')

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.plot(xx,u_exact)
    ax.scatter(xx[j],u_exact[j])
    plt.show()

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    
    markers = ['.','v','^','s']
    colors  = ['b','r','purple','k']

    Ns = [j for j in range(2,60,4)]
    quadDegs = [4,5,6]

    for i in range(len(quadDegs)):
        Atom['quadDeg'] = quadDegs[i]
        Y  = []
        for N in Ns:
            Atom['N'] = N
            u = numsol(Atom)
            Y.append(abs(u(xx[j])-u_exact[j]))

        ax.scatter(Ns,Y,facecolor='none',edgecolor=colors[i],
                   s=400,lw=2,marker=markers[i],
                   label='$d = {}$'.format(quadDegs[i]))

    ax.set_xlabel('$N$')
    ax.set_ylabel('$|u_{20}(13a_0)-P_{20}(13a_0)|$')
    ax.set_yscale('log')
    ax.legend(framealpha=0)
    ax.grid(ls='--',alpha=0.5)
    plt.show()

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    
    for i in range(len(quadDegs)):
        Atom['eigs'] = []
        Atom['quadDeg'] = quadDegs[i]
        for N in Ns:
            Atom['N'] = N
            u = numsol(Atom)
        Y = array([abs(j/100+1/9) for j in Atom['eigs']])
        ax.scatter(Ns,Y,
                   facecolor='none',edgecolor=colors[i],
                   marker=markers[i],s=400,lw=2,
                   label='$d={}$'.format(quadDegs[i]))
    
    ax.legend(framealpha=0)
    ax.set_xlabel('$N$')
    ax.set_ylabel("Error $E'$")
    ax.set_yscale('log')
    ax.grid(ls='--',alpha=0.5)
    plt.show()
    

def convergence_test(Atom):

    xx = array([j/1000 for j in range(1001)])
    u = numsol(Atom)
    j = argmax(u(xx))

    Ns = [j*8 for j in range(1,10)]
    Y  = []
    for N in Ns:
        Atom['N'] = N
        u = numsol(Atom)
        Y.append(u(xx[j]))

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.scatter(Ns,Y,facecolor='none',edgecolor='k',s=100)

    
    ax.set_xlabel('$N$')
    ax.set_ylabel('$u(a_0)$')
    plt.show()
    return 0
        

def plot_radial(Atom):
    
    u = numsol(Atom)

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    xx = array([j/1000 for j in range(1001)])
    ticks = 5 
    xticks      = [int(10*j/ticks)/10 for j in range(ticks+1)]
    xticklabels = [int(j*Atom['D']/ticks) for j in range(ticks+1)]
    ax.set_xticks(xticks,xticklabels)
    ax.set_ylabel('$r^2 R_{n,\ell}^2$')
    ax.set_xlabel('$r/a_0$')

    ax.plot(xx,u(xx)**2/Atom['D'],lw=1,c='k')
    plt.show()

    return 0

def spectrum():

    ell = [0,0,0,0,
           1,1,1,
           2,2, 
           3]
    E   = [-66.67,-20.30,-9.66,-5.62,
           -24.67,-11.00,-6.20,
           -11.11,-6.25,
           -6.25]
    ell = array(ell)
    E   = array(E)

    fig = plt.figure()
    ax  = fig.add_subplot(111)

    for n in range(1,5):
        ax.plot([min(ell),max(ell)],2*[1/n**2],ls='--',c='k'
                ,lw=1.5,alpha=0.5)
        if n != 1:
            ax.text(x=2.5,y=1/(n-0.1)**2,s='$1/{}^2$'.format(n))
            

    ax.scatter(ell,-E/100,lw=4,s=500,c='k',marker='_')
    ax.set_xticks([n for n in range(4)])
    ax.set_xlabel(r"$\ell$")
    ax.set_ylabel("$-E'$")
    ax.set_yscale('log')
    ax.grid(alpha=0.7,ls='--',which='both')

    plt.show()

def UniformShell_radials1(Atom):

    dist = 1.5    # bohr radius
    Atom['D'] = 40
    Atom['ell'] = 0
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    xx = array([j/1000 for j in range(1001)])
    rr = Atom['D']*xx
    ax.set_ylabel('$u_{n0}^2$')
    ax.set_xlabel('$r/a_0$')

    Atom['E']  = -60
    Atom['R0'] = dist/Atom['D']
    u = numsol(Atom)
    ax.plot(rr,u(xx)**2/Atom['D'],lw=2,ls='--',c='purple')

    Atom['E']  = -99
    Atom['R0'] = 0
    u = numsol(Atom)
    ax.plot(rr,u(xx)**2/Atom['D'],lw=2,c='k',
            label='$V_C$')

    Atom['E']  = -25
    Atom['R0'] = dist/Atom['D']
    u = numsol(Atom)
    ax.plot(rr,u(xx)**2/Atom['D'],lw=2,c='purple',ls='--',
            label='$V_S$')

    Atom['E']  = -25
    Atom['R0'] = 0
    u = numsol(Atom)
    ax.plot(rr,u(xx)**2/Atom['D'],lw=2,c='k')

    Atom['E']  = -11
    Atom['R0'] = dist/Atom['D']
    u = numsol(Atom)
    ax.plot(rr,u(xx)**2/Atom['D'],lw=2,c='purple',ls='--')

    Atom['E']  = -11
    Atom['R0'] = 0
    u = numsol(Atom)
    ax.plot(rr,u(xx)**2/Atom['D'],lw=2,c='k')

    ax.text(x=2.8,y=0.4,s='$n=1$')
    ax.text(x=8.8,y=0.19,s='$n=2$')
    ax.text(x=22,y=0.12,s='$n=3$')
    ax.set_yticks([0.0,0.5])

    ax.plot([dist,dist],[0,0.5],
            ls='dashdot',lw=3,c='k',
            label='Sphere end')
    ax.set_xscale('log')
    ax.set_yticks([0.0,0.5])
    ax.legend(framealpha=0,loc='upper left')
    

    plt.show()


def UniformShell_radials2(Atom):

    dist = 1.5    # bohr radius
    Atom['D']   = 60
    Atom['ell'] = 1
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    xx = array([j/1000 for j in range(1001)])
    rr = Atom['D']*xx
    ax.set_ylabel('$u_{n1}^2$')
    ax.set_xlabel('$r/a_0$')

    Atom['E']  = -25
    Atom['R0'] = dist/Atom['D']
    u = numsol(Atom)
    ax.plot(rr,u(xx)**2/Atom['D'],lw=2,ls='--',c='purple')

    Atom['E']  = -25
    Atom['R0'] = 0
    u = numsol(Atom)

    ax.plot(rr,u(xx)**2/Atom['D'],lw=2,c='k',
            label='$V_C$')

    Atom['E']  = -11
    Atom['R0'] = dist/Atom['D']
    u = numsol(Atom)
    ax.plot(rr,u(xx)**2/Atom['D'],lw=2,c='purple',ls='--',
            label='$V_S$')

    Atom['E']  = -11
    Atom['R0'] = 0
    u = numsol(Atom)
    ax.plot(rr,u(xx)**2/Atom['D'],lw=2,c='k')

    Atom['E']  = -6
    Atom['R0'] = dist/Atom['D']
    u = numsol(Atom)
    ax.plot(rr,u(xx)**2/Atom['D'],lw=2,c='purple',ls='--')

    Atom['E']  = -6
    Atom['R0'] = 0
    u = numsol(Atom)
    ax.plot(rr,u(xx)**2/Atom['D'],lw=2,c='k')

    ax.text(x=6,y=0.175,s='$n=2$')
    ax.text(x=15,y=0.09,s='$n=3$')
    ax.text(x=29,y=0.05,s='$n=4$')
    ax.set_yticks([0.0,0.2])
    
    ax.plot([dist,dist],[0,0.2],
            ls='dashdot',lw=3,c='k',
            label='Sphere end')
    ax.set_xscale('log')
    ax.legend(framealpha=0,loc='upper left')

    plt.show()

def Uranium_test(Atom):

    D = Atom['D']
    Atom['R0'] = 1.618e-4/Atom['D']
    Atom['R0'] = 0.01

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    xx = array([j/1000 for j in range(1001)])
    
    u = numsol(Atom)
    uu = u(xx)
    rr = D*xx

    ax.plot(xx,u(xx),lw=3,c='k',alpha=0.7)
    ax.set_ylabel('$u_{n,\ell}$')
    ax.set_xlabel(r'$r/a_0$')

    Atom['R0'] = 0
    u = numsol(Atom)
    ax.plot(xx,u(xx),lw=2,c='b',ls='--')

    plt.show()
