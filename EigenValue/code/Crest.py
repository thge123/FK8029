from PLOTTING import *
import FILES
from numpy import array,exp,log
from numpy import abs as ABS
from scipy.stats import linregress

def plot(ax,params,y,Values,density=False,shift=0,pot_shift=0,Color='k'):
    
    h = params['h']
    N = params['N']

    x = array([j*h for j in range(int(N+1))])
    y0 = max(ABS(y))
    if params['parity'] == 'odd':
        y = array([0] + [j/y0 for j in y] + [0])
        if density:
            y = normalized_density(x,y)
            i = Crest_finder(y)
            Avg_dist = avg_dist(x,y)
            Avg_dist_sqrd = avg_dist_sqrd(x,y)
            ax.plot([x[-i],x[-i]],[shift,shift+0.2])
            ax.plot(x,y+shift,c=Color)
            ax.plot(-x,y+shift,c=Color)
        else:
            ax.plot(x,y+shift,c=Color)
            ax.plot(-x,-y+shift,c=Color)
    else:
        y = array([j/y0 for j in y] + [0,0])
        if density:
            y = normalized_density(x,y)
            i = Crest_finder(y)
            Avg_dist = avg_dist(x,y)
            Avg_dist_sqrd = avg_dist_sqrd(x,y)
            ax.plot([x[-i],x[-i]],[shift,shift+0.2])
            ax.plot(x,y+shift,c=Color)
            ax.plot(-x,y+shift,c=Color)
        else:
            ax.plot(x,y+shift,c=Color)
            ax.plot(-x,y+shift,c=Color)
    E = params['E']
    print('E = ', E)
    print('err = ', params['err'])
    Values['Crest'].append(x[-i])
    Values['Avg_dist'].append(Avg_dist)
    Values['Variance'].append(Avg_dist_sqrd - Avg_dist**2)
    #Pot = [pot(j) for j in x]
    #Pot = [pot_shift*j/max(Pot) for j in Pot]
    #ax.plot(x,Pot,ls='--',c=Color)
    #ax.plot(-x,Pot,ls='--',c=Color)

    ax.plot([x[0],x[-1]],[shift,shift],ls='--',c=Color)
    ax.plot([-x[-1],x[0]],[shift,shift],ls='--',c=Color)

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

def avg_dist(X,y_normalized):

    I = 0
    for i in range(len(X)-1):
        x_avg = 0.5*(X[i+1]+X[i])
        y_avg = 0.5*(y_normalized[i+1]+y_normalized[i])
        I += x_avg*y_avg*(X[i+1]-X[i])
    return 2*I

def avg_dist_sqrd(X,y_normalized):

    I = 0
    for i in range(len(X)-1):
        x_avg = 0.5*(X[i+1]+X[i])
        y_avg = 0.5*(y_normalized[i+1]+y_normalized[i])
        I += x_avg**2*y_avg*(X[i+1]-X[i])
    return 2*I


def Crest_finder(X):
    k = 0
    for i in range(1,len(X)):
        if X[-i] >= k:
            k = X[-i]
        else:
            return i-1
def main():

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    M = 100
    j = 2 
    A = 'A/A'+FILES.number2string(j)+'.dat'
    B = 'B/B'+FILES.number2string(j)+'.dat'
    energies = []
    Values = {'Crest'  : [],
             'Avg_dist': [],
             'Variance': []}
    for i in range(2,M):
        state = i    # state 1 = ground state

        X = FILES.get_data(A,fileline=state)
        E = FILES.get_data(B, fileline=state)
        params = {'E': E[0], 'N': E[1], 'h': E[2], 'err': E[3]}
        energies.append(E[0])

        if state%2 == 1:
            params['parity'] = 'even'
        else:
            params['parity'] = 'odd'
        
        plot(ax,params,X,Values,density=True,shift=params['E'],pot_shift=2*M)
    print(Values['Crest'])
    print(Values['Avg_dist'])
    print(Values['Variance'])
    ax.set_xlabel(r'$x\sqrt{m\omega/\hbar}$',fontsize=32)
    ax.set_ylabel(r"$2E/\hbar\omega$",fontsize=32)
    #ax.plot([1e6,1e6],[1e6,1e6],c='k',label=r'Probability densities for $\alpha=10$.')
    #ax.set_xlim(-5.2,5.2)
    #ax.set_xticks([-4,-2,0,2,4])
    #ax.set_ylim(-0,13)
    #ax.set_yticks([1,3,5,7,9,11])
    #ax.legend(loc='upper center',framealpha=0)
    print(energies)
    
    fig2 = plt.figure()
    ax2  = fig2.add_subplot(121)
    ax3  = fig2.add_subplot(122)
    ax2.scatter(energies,Values['Crest'],facecolor='none',edgecolor='k',s=100,label='Wavefront')
    ax2.scatter(energies,Values['Avg_dist'],marker='v',facecolor='none',edgecolor='r',s=100,label=r'Average')
    #ax2.scatter(energies,[j**0.5 for j in Values['Variance']],marker='*',facecolor='none',edgecolor='b',s=100,label=r'Standard deviation')
    ax3.scatter(energies,Values['Crest'],facecolor='none',edgecolor='k',s=100)
    ax3.scatter(energies,Values['Avg_dist'],marker='v',facecolor='none',edgecolor='r',s=100)
    #ax3.scatter(energies,[j**0.5 for j in Values['Variance']],marker='*',facecolor='none',edgecolor='b',s=100)
    #crests_diff = [crests[j+1]-crests[j] for j in range(len(crests)-1)]
    #ax2.scatter(energies[:-1],crests_diff,facecolor='none',edgecolor='k')
    ax2.plot(energies,[j**0.5 for j in energies],label=r'$\sqrt{2E/\hbar\omega}$',c='k')
    ax3.plot(energies,[j**0.5 for j in energies],c='k')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax2.set_xlabel('$2E/\hbar\omega$',fontsize=32)
    ax2.set_ylabel('Distance from origin (r.u.)',fontsize=32)
    ax3.set_xlabel('$2E/\hbar\omega$',fontsize=32)
    ax2.legend(framealpha=0)
    
    plt.show()
    

main()
