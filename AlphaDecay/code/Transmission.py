from PLOTTING import *
from XY       import *

def main():

    ''' Calculate transmission coefficient. 
        Incoming wave Aexp(iw0x)+Bexp(-iw0x)
        and outgoing wave Fexp(iw1x). Trans-
        mission coefficient is equal to 
        
                T = abs(F/A)^2*w1/w0.
    
        Reflection coefficient:

                R = abs(B/A)^2

        If all right: T+R = 1.
    '''
    
    Coeffs = XY()
    Omegas = XY()
    Coeffs.get_X('X.dat',Complex=True)
    Omegas.get_X('omega.dat')
    
    A = Coeffs.X[0]
    A_sqrd = A[0]**2 + A[1]**2
    print(A_sqrd)
    B = Coeffs.X[1]
    B_sqrd = B[0]**2 + B[1]**2
    F = Coeffs.X[-1]
    F_sqrd = F[0]**2 + F[1]**2

    w0 = Omegas.X[0]
    w1 = Omegas.X[-1]

    T = F_sqrd*w1/(A_sqrd*w0)
    R = B_sqrd/A_sqrd
    print("T = ", T)
    print("R = ", R)
    print("T+R = ", T+R)
    

main()
