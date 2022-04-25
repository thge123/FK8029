from PLOTTING import *
from Numerics import *
from BSplineAtom import *


def main():

    Atom = {'Z'         : 2,     # Atomic number
            'ell'       : 0,     # Angular momentum
            'E'         : -99,   # Guess at eigenvalue
            'N'         : 20,    # Number of inner points
            'D'         : 50,    # Scaling of interval
            'mass_ratio': 1,     # red. mass/mass of elect.
            'R0'        : 0.0}   # Uniform shell cut-off

    plot_radial(Atom,exact=False)

main()
