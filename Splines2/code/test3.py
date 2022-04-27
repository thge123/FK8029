from PLOTTING import *
from Numerics import *
from BSplineAtom import *


def main():

    Atom = {'Z'         : 1,     # Atomic number
            'ell'       : 0,     # Angular momentum
            'E'         : -11,   # Guess at eigenvalue
            'N'         : 30,    # Number of inner points
            'D'         : 50,    # Scaling of interval
            'mass_ratio': 1,     # red. mass/mass of elect.
            'R0'        : 0.0,   # Uniform shell cut-off
            'quadN'     : 20,   # Quadrature nodes 
            'quadDeg'   : 8}   # Quadrature nodes 

    #plot_radial(Atom)
    u20_test(Atom)

main()
