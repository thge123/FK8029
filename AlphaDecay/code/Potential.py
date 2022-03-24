from XY import *
from PLOTTING import *

def get_energy():
    return float(input("Write E: "))

def main():
    Potential = XY()
    Potential.get_X('r.dat')
    Potential.get_Y('V.dat')

    r = [Potential.X[0]]
    V = [Potential.Y[0]]
    for i in range(1,len(Potential.Y)):
        r.append(Potential.X[i]-1e-9)
        r.append(Potential.X[i])
        
        V.append(Potential.Y[i-1])
        V.append(Potential.Y[i])

    r.append(Potential.X[-1])
    V.append(V[-1])

    E = get_energy()

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.plot(r,V,c='k',lw=2)
    ax.plot([r[0],r[-1]],[E,E])
    

    plt.show()

main()

