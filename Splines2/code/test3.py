from PLOTTING import *
from Numerics import *
from BSplineAtom import *

def main():

    Atom = {'Z'         : 1,     # Atomic number
            'ell'       : 0,     # Angular momentum
            'E'         : -25,   # Guess at eigenvalue
            'N'         : 20,    # Number of inner points
            'D'         : 40,    # Scaling of interval
            'mass_ratio': 1,     # red. mass/mass of elect.
            'R0'        : 0.2}   # Uniform shell cut-off

    u = numsol(Atom)

    xx = array([j/1000 for j in range(1001)])
    plt.plot(xx,u(xx))

    plt.plot(xx,Exact['u10'](xx,Atom['D']))
    plt.show()

main()
