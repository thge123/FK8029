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

def test():
    
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
                facecolor='none',edgecolor='k',s=100)
    plt.scatter([j for j in rslt],
                measured_temp(),marker='x',s=100)

    filename = "US_STD_ATM.dat"
    x = get_data(filename)
    
    z = []
    t = []
    for i in x:
        z.append(i[0])
        t.append(i[1])

    plt.scatter(z,t,facecolor='none',edgecolor='b',
                marker='v',s=100)

    plt.show()

def main():

    test()

main()
