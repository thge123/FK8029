from PLOTTING import *
import FILES
from numpy import array

def plot(params,y,density=False):
    
    h = params['h']
    N = params['N']

    x = array([j*h for j in range(int(N+1))])
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    if params['parity'] == 'odd':
        y = array([0] + [j for j in y] + [0])
        if density:
            y = normalized_density(x,y)
            ax.plot(x,y,c='k')
            ax.plot(-x,y,c='k')
        else:
            ax.plot(x,y,c='k')
            ax.plot(-x,-y,c='k')
    else:
        y = array([j for j in y] + [0,0])
        if density:
            y = normalized_density(x,y)
            ax.plot(x,y,c='k')
            ax.plot(-x,y,c='k')
        else:
            ax.plot(x,y,c='k')
            ax.plot(-x,y,c='k')
    #pot = [test_pot(j) for j in x]
    #pot = [j/max(pot) for j in pot]
    #ax.plot(x,pot,ls='--',c='k')
    #ax.plot(-x,pot,ls='--',c='k')


def main():
    state = 2 
    X = FILES.get_data('A/A000.dat',fileline=state)
    E = FILES.get_data('B/B000.dat', fileline=state)
    params = {'E': E[0], 'N': E[1], 'h': E[2], 'err': E[3],}

    if state%2 == 1:
        params['parity'] = 'even'
    else:
        params['parity'] = 'odd'
    
    plot(params,X)
    plt.show()
    

main()














