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
    ax.plot(x,Pot,ls='--',c=Color)
    ax.plot(-x,Pot,ls='--',c=Color)

    ax.plot([x[0],x[-1]],[shift,shift],ls='--',c=Color)
    ax.plot([-x[-1],x[0]],[shift,shift],ls='--',c=Color)
    

def pot(x):
    alpha = 10
    return x**2
    return x**2+alpha*exp(-x**2)-(1+log(alpha))


def normalized_density(X,y):

    I = 0
    y_out = y**2
    for i in range(len(X)-1):
        y_avg = (y_out[i+1]+y_out[i])/2
        I += y_avg*(X[i+1]-X[i])
    return y_out/(2*I)


def main():

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    M = 10
    j = 2 
    A = 'A/A'+FILES.number2string(j)+'.dat'
    B = 'B/B'+FILES.number2string(j)+'.dat'
    for i in range(1,M):
        state = i    # state 1 = ground state

        X = FILES.get_data(A,fileline=state)
        E = FILES.get_data(B, fileline=state)
        params = {'E': E[0], 'N': E[1], 'h': E[2], 'err': E[3]}

        if state%2 == 1:
            params['parity'] = 'even'
        else:
            params['parity'] = 'odd'
        
        plot(ax,params,X,density=True,shift=params['E'],pot_shift=2*M)
    ax.set_xlabel(r'$x\sqrt{m\omega/\hbar}$',fontsize=32)
    ax.set_ylabel(r"$2E/\hbar\omega$",fontsize=32)
    ax.plot([1e6,1e6],[1e6,1e6],c='k',label=r'Probability densities for $\alpha=0$.')
    ax.set_xlim(-8.2,8.2)
    ax.set_xticks([-6,-4,-2,0,2,4,6])
    ax.set_ylim(-0,13)
    ax.set_yticks([1,3,5,7,9,11])
    ax.legend(loc='upper center',framealpha=0)


    plt.show()
    

main()



