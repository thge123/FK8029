from PLOTTING import *

def get_data(filename):
    X = []
    with open(filename) as File:
        for i in File.readlines():
            X.append([float(j) for j in i.split(';')])
    return X

def measured_temp():
    T = [15,8.5,2,-4.49,-10.98,-17.47,
         -23.96,-30.45,-36.94,-43.42,
         -49.9,-56.5,-56.5,-51.6,-46.64,
         -22.8,-2.5,-26.13,-53.57,-74.51]
    return [273.15+j for j in T]



def US_STD_ATM():
    
    filename = "US_STD_ATM.dat"
    x = get_data(filename)
    
    z = []
    t = []
    for i in x:
        z.append(i[0])
        t.append(i[1])

    plt.scatter(z,t,facecolor='none',edgecolor='k',s=100)
    plt.scatter(z,measured_temp(),marker='x',s=100)

    plt.show()

def Temp():
    
    filename = "test.dat"
    X = get_data(filename)
    
    rslt = {}
    for i in X:
        if i[0] in rslt:
            rslt[i[0]].append(i[1])
        else:
            rslt[i[0]] = [i[1]]

    plt.scatter([j for j in rslt],
                [rslt[j][-1] for j in rslt],
                facecolor='none',edgecolor='k',s=200,
                label='Multi-layer model')
    plt.scatter([j for j in rslt],
                measured_temp(),marker='x',s=200,
                label='US standard')

    plt.xlabel('Altitude $z$ [m]')
    plt.ylabel('Temperature $T$ [K]')
    filename = "US_STD_ATM.dat"
    x = get_data(filename)
    
    z = []
    t = []
    for i in x:
        z.append(i[0])
        t.append(i[1])

    #plt.scatter(z,t,facecolor='none',edgecolor='b',
    #            marker='v',s=100)
    
    plt.legend(framealpha=0,loc='lower left')
    plt.grid(alpha=0.5,ls='--')
    plt.show()

def Dens():
    
    filename = "test.dat"
    X = get_data(filename)
    
    rslt = {}
    for i in X:
        if len(i) > 2:
            if i[0] in rslt:
                rslt[i[0]].append(i[2])
            else:
                rslt[i[0]] = [i[2]]

    plt.scatter([j for j in rslt],
                [rslt[j][-1] for j in rslt],
                facecolor='none',edgecolor='k',s=200,
                label='Multi-layer model')
    #plt.scatter([j for j in rslt],
    #            measured_temp(),marker='x',s=100)

    filename = "US_STD_ATM.dat"
    x = get_data(filename)
    
    z = []
    t = []
    for i in x:
        z.append(i[0])
        t.append(i[2])

    plt.scatter(z,t,c='b',
                marker='x',s=200,label='US standard')
    plt.grid(alpha=0.5,ls='--')
    plt.xlabel('Altitude $z$ [m]')
    plt.ylabel(r'Density $\rho$ [kg/m$^3$]')
    plt.legend(framealpha=0)
    

    plt.show()

def Iter():
    
    filename = "test.dat"
    X = get_data(filename)
    
    rslt = {}
    for i in X:
        if i[0] in rslt:
            rslt[i[0]].append(i[1])
        else:
            rslt[i[0]] = [i[1]]

    for i in rslt:
        plt.plot([j for j in range(len(rslt[i]))],
                 [j for j in rslt[i]],c='k',lw=0.5)
        plt.scatter([j for j in range(len(rslt[i]))],
                    [j for j in rslt[i]],facecolor='none',
                    edgecolor='k',s=100)

    plt.xlabel('Iteration')
    plt.ylabel('Temperature $T$ [K]')
    filename = "US_STD_ATM.dat"
    x = get_data(filename)
    
    z = []
    t = []
    for i in x:
        z.append(i[0])
        t.append(i[1])

    #plt.scatter(z,t,facecolor='none',edgecolor='b',
    #            marker='v',s=100)
    
    plt.grid(alpha=0.5,ls='--')
    plt.show()

def main():

    Iter()

main()
