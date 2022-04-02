from PLOTTING import *
import FILES
from numpy import array,exp,log
from numpy import abs as ABS

psi = {1: lambda t: exp(-t**2/2),
       2: lambda t: t*exp(-t**2/2),
       3: lambda t: (2*t**2-1)*exp(-t**2/2),
       4: lambda t: (2*t**3-3*t)*exp(-t**2/2)
      }

def plot(ax,params,y,density=False,shift=0,pot_shift=0,Color='k'):
    
    h = params['h']
    N = params['N']

    x = array([j*h for j in range(int(N+1))])
    y0 = max(ABS(y))
    if params['parity'] == 'odd':
        y = array([0] + [j/y0 for j in y] + [0])
        if density:
            y = normalized_density(x,y)
            ax.plot(x,y+shift,c=Color)
            ax.plot(-x,y+shift,c=Color)
        else:
            ax.plot(x,y+shift,c=Color)
            ax.plot(-x,-y+shift,c=Color)
    else:
        y = array([j/y0 for j in y] + [0,0])
        if density:
            y = normalized_density(x,y)
            ax.plot(x,y+shift,c=Color)
            ax.plot(-x,y+shift,c=Color)
        else:
            ax.plot(x,y+shift,c=Color)
            ax.plot(-x,y+shift,c=Color)
    E = params['E']
    print('E = ', E)
    print('err = ', params['err'])
    Pot = [pot(j) for j in x]
    Pot = [pot_shift*j/max(Pot) for j in Pot]
    #ax.plot(x,Pot,ls='--',c=Color)
    #ax.plot(-x,Pot,ls='--',c=Color)

    #ax.plot([x[0],x[-1]],[shift,shift],ls='--',c=Color)
    #ax.plot([-x[-1],x[0]],[shift,shift],ls='--',c=Color)
    

def pot(x):
    alpha = 5
    alpha = 10
    return x
    return x**2+alpha*exp(-x**2)-0.5*(1+log(2*alpha))


def normalized_density(X,y):

    I = 0
    y_out = y**2
    for i in range(len(X)-1):
        y_avg = (y_out[i+1]+y_out[i])/2
        I += y_avg*(X[i+1]-X[i])
    return y_out/(2*I)


def main():

    fig = plt.figure()
    ax  = fig.add_subplot(121)
    ax1  = fig.add_subplot(122)
    M = 5
    for j in [2,7,8]:
        A = 'A/A'+FILES.number2string(j)+'.dat'
        B = 'B/B'+FILES.number2string(j)+'.dat'
        for i in range(1,M):
            state = i    # state 1 = ground state

            X = FILES.get_data(A,fileline=state)
            if (j,i) == (2,4):
                X = -array(X)
            if (j,i) == (2,2):
                X = -array(X)
            if (j,i) == (7,2):
                X = -array(X)
            if (j,i) == (7,3):
                X = -array(X)
            if (j,i) == (8,3):
                X = -array(X)
            E = FILES.get_data(B, fileline=state)
            params = {'E': E[0], 'N': E[1], 'h': E[2], 'err': E[3],}
            y1 = [j*params['h'] for j in range(1,int(params['N']))] 
            y1 = [psi[i](j) for j in y1] 

            if state%2 == 1:
                params['parity'] = 'even'
            else:
                params['parity'] = 'odd'
            
            plot(ax,params,X,density=False,shift=params['E'],pot_shift=2*M)
            if j==2:
                plot(ax1,params,y1,density=False,shift=params['E'],pot_shift=2*M,Color='b')
    ax.set_xlabel(r'$x\sqrt{m\omega/\hbar}$',fontsize=32)
    ax1.set_xlabel(r'$x\sqrt{m\omega/\hbar}$',fontsize=32)
    ax.set_ylabel(r"$2E/\hbar\omega$",fontsize=32)
    ax.plot([1e6,1e6],[1e6,1e6],c='k',label=r'Numerical solution eigenstates')
    ax1.plot([1e6,1e6],[1e6,1e6],c='b',label=r'Exact solution eigenstates')
    ax.set_xlim(-5.2,5.2)
    ax.set_ylim(-0,9)
    ax1.set_xlim(-5.2,5.2)
    ax1.set_ylim(-0,9)
    ax.legend(framealpha=0)
    ax1.legend(framealpha=0)
    #ax.legend(loc='upper right',framealpha=0)


    plt.show()
    

main()


