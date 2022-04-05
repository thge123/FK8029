from PLOTTING import *
import FILES
from numpy import array,exp,log
from numpy import abs as ABS

psi = {1: lambda t: exp(-t**2/2),
       2: lambda t: t*exp(-t**2/2),
       3: lambda t: (2*t**2-1)*exp(-t**2/2),
       4: lambda t: (2*t**3-3*t)*exp(-t**2/2)
      }

def plot(ax,params,y,density=False,shift=0,pot_shift=0,Color='k',LS='-'):
    
    h = params['h']
    N = params['N']

    x = array([j*h for j in range(int(N+1))])
    y0 = max(ABS(y))
    if params['parity'] == 'odd':
        y = array([0] + [j/y0 for j in y] + [0])
        if density:
            y = normalized_density(x,y)
            ax.plot(x,y+shift,c=Color,lw=2,ls=LS)
            ax.plot(-x,y+shift,c=Color,lw=2,ls=LS)
        else:
            ax.plot(x,y+shift,c=Color,lw=2,ls=LS)
            ax.plot(-x,-y+shift,c=Color,lw=2,ls=LS)
    else:
        y = array([j/y0 for j in y] + [0,0])
        if density:
            y = normalized_density(x,y)
            ax.plot(x,y+shift,c=Color,lw=2,ls=LS)
            ax.plot(-x,y+shift,c=Color,lw=2,ls=LS)
        else:
            ax.plot(x,y+shift,c=Color,lw=2,ls=LS)
            ax.plot(-x,y+shift,c=Color,lw=2,ls=LS)
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
    alpha = 10
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
    ax0 = fig.add_subplot(111)
    ax  = fig.add_subplot(221)
    ax2  = fig.add_subplot(222)
    ax3  = fig.add_subplot(223)
    ax4  = fig.add_subplot(224)
    axs = [ax,ax2,ax3,ax4]
    j = 0 
    Colors = {2:'b',11:'k',12:'r'}
    Ls =     {2:'-',11:'--',12:'dashdot'}
    labels = {2: r'$\alpha=0$',
              11: r'$\alpha=5$',
              12: r'$\alpha=10$'}
    ax0.spines['top'].set_color('none')
    ax0.spines['bottom'].set_color('none')
    ax0.spines['left'].set_color('none')
    ax0.spines['right'].set_color('none')
    ax0.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

    for state in [1,2]:
        A = 'A/A'+FILES.number2string(11)+'.dat'
        B = 'B/B'+FILES.number2string(11)+'.dat'
        X = FILES.get_data(A,fileline=state)
        if state==2:
            X = -array(X)
        E = FILES.get_data(B, fileline=state)
        params = {'E': E[0], 'N': E[1], 'h': E[2], 'err': E[3]}
        if state%2 == 1:
            params['parity'] = 'even'
        else:
            params['parity'] = 'odd'
        plot(axs[state-1],params,X,density=False,shift=0,LS=Ls[11],Color=Colors[11])
    for state in [1,2]:
        A = 'A/A'+FILES.number2string(12)+'.dat'
        B = 'B/B'+FILES.number2string(12)+'.dat'
        X = FILES.get_data(A,fileline=state)
        if state==2:
            X = -array(X)
        E = FILES.get_data(B, fileline=state)
        params = {'E': E[0], 'N': E[1], 'h': E[2], 'err': E[3]}
        if state%2 == 1:
            params['parity'] = 'even'
        else:
            params['parity'] = 'odd'
        plot(axs[state+1],params,X,density=False,shift=0,LS=Ls[12],Color=Colors[12])
    ax.plot([1e6],[1e6],label=r'$\alpha=5$',ls=Ls[11],c=Colors[11])
    ax.plot([1e6],[1e6],label=r'$\alpha=10$',ls=Ls[12],c=Colors[12])
    for i in axs:
        i.set_yticks([-1,0,1])
        i.set_xlim(-8,8)
        i.set_ylim(-1.1,1.1)
    ax.legend(framealpha=0)
    ax.set_title("$E'_1$")
    ax2.set_title("$E'_2$")
    ax3.set_title("$E'_{1}$")
    ax4.set_title("$E'_{2}$")
    ax3.set_xlabel(r'$x\sqrt{m\omega/\hbar}$',fontsize=32)
    ax4.set_xlabel(r'$x\sqrt{m\omega/\hbar}$',fontsize=32)
    ax0.set_ylabel('Eigenstates',fontsize=32)



    plt.show()
    

main()
