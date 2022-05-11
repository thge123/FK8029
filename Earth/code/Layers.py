from numpy import meshgrid
from PLOTTING import *

def get_data(filename):
    X = []
    with open(filename) as File:
        for i in File.readlines():
            X.append([float(j) for j in i.split(';')])
    return X
            

def OneLayer():
    filename = "OneLayer.dat"
    
    eps = []
    T0  = []
    T1  = []
    
    e1,e2,T0,T1,T2 = [],[],[],[],[]
    for i in get_data(filename):
        eps.append(i[0])
        T0.append(i[1])
        T1.append(i[2])

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    
    ax.plot(eps,T0,c='k',label='$T_0$',lw=2)
    ax.plot(eps[1:],T1[1:],c='b',ls='--',label='$T_1$',lw=2)

    ax.set_xlabel("$\epsilon_1$")
    ax.set_ylabel("Temperature [K]")

    ax.legend(framealpha=0)
    ax.grid(alpha=0.5,ls='--')

    plt.show()
    
def TwoLayer():
    filename = "TwoLayer.dat"
    
    e1,e2,T0,T1,T2 = [],[],[],[],[]
    for i in get_data(filename):
        e1.append(i[0])
        e2.append(i[1])
        T0.append(i[2])
        T1.append(i[3])
        T2.append(i[4])

    fig = plt.figure()
    ax  = fig.add_subplot(111)

    for e2_val in [0.1,0.2,0.3,0.4,0.5]:
        x = []
        y = []
        for i in range(len(e1)):
            if e2[i] == e2_val:        
                x.append(e1[i])
                y.append(T0[i])
        ax.plot(x,y)
    
    plt.show()
    
    
def main():
    
    OneLayer()

main()


