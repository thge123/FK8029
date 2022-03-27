import FILES
from XY import *
from PLOTTING import *
from Density import get_data

def get_energy():
    return float(input("Write E: "))

def main():

    filenumber = int(input("Write filenumber: "))

    r,omega,Coeffs,V = get_data(filenumber)
    Potential = XY()

    r_new = [r[0]]
    V_new = [V[0]]
    for i in range(1,len(V)):
        r_new.append(r[i]-1e-9)
        r_new.append(r[i])
        
        V_new.append(V[i-1])
        V_new.append(V[i])

    r_new.append(r[-1])
    V_new.append(V_new[-1])
    r_new.append(2*r[-1] - r[-2])
    V_new.append(V_new[-1])

    E = get_energy()

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.plot(r_new,V_new,c='k',lw=2)
    ax.plot([r_new[0],r_new[-1]],[E,E])
    

    plt.show()

main()

