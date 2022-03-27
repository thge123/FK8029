from PLOTTING import *
import FILES
from Density import get_data 
from os import listdir

def scatter_TR(ax,filenumber):
    
    r,omega,Coeffs,V = get_data(filenumber)
    A = Coeffs[0]
    A_sqrd = A[0]**2 + A[1]**2
    B = Coeffs[1]
    B_sqrd = B[0]**2 + B[1]**2
    F = Coeffs[-1]
    F_sqrd = F[0]**2 + F[1]**2

    w0 = omega[0]
    w1 = omega[-1]

    T = F_sqrd*w1/(A_sqrd*w0)
    R = B_sqrd/A_sqrd
    
    print("N = ", filenumber, ": ", "T = ",   T)
    print("N = ", filenumber, ": ", "R = ",   R)
    print("N = ", filenumber, ": ", "T+R = ", T+R)
    
    ax.scatter(filenumber,T,edgecolor='k',facecolor='none',s=60)
    #ax.scatter(filenumber,R,marker='v',c='green',s=60)
    #ax.scatter(filenumber,T+R,edgecolor='k',facecolor='none',s=60)
    

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

    filenames = listdir('O')
    filenames.sort()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in filenames:
        scatter_TR(ax,FILES.string2number(i))


    ax.grid(axis='y',ls='--')    
    plt.show()

main()
