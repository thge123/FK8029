from PLOTTING import *
from Numerics import *
from BSplineAtom import *


def main():

    Atom = {'Z'         : 1,     # Atomic number
            'ell'       : 0,     # Angular momentum
            'E'         : -99,   # Guess at eigenvalue
            'N'         : 31,    # Number of inner points
            'D'         : 60,    # Scaling of interval
            'mass_ratio': 1,     # red. mass/mass of elect.
            'R0'        : 0.0,   # Uniform shell cut-off
            'quadN'     : 20,   # Quadrature nodes 
            'quadDeg'   : 5}     # Quadrature degree 

    u20_test(Atom)

main()
