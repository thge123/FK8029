from PLOTTING import *
from BSpline2 import *

def main():

    
    i = int(input("Write sigma [1,2,3]: ")) - 1 
    N = int(input("Write number of inner points: "))
    X   = array([j/1000 for j in range(1001)])
    exact_sol = [sigma1_exact,sigma2_exact,sigma3_exact,sigma4_exact,sigma5_exact][i]

    if i == 1:
        r = float(input("Write r: "))
    elif i == 2:
        r = float(input("Write alpha: "))
    else:
        r = 1

    spl = numsol(N,i,r)

main()
