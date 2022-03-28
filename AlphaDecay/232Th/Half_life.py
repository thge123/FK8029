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

    #232Th Z = 90
    E = 4.081
    v = 0.69*E**0.5*1e7    #m/s
    A = 232
    R = 1.5*(4**0.3333 + (A-4)**0.3333)    # fm
    f = v/(2*R) * 1e15
    L = f*T
    T_half = 0.6931/L
    
    print("N = ", filenumber, ": ", "T = ",   T_half)
    
    ax.scatter(filenumber,T_half,edgecolor='k',facecolor='none',s=60)
    ax.set_yscale('log')
    
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
    T_half = 4.43e17
    ax.plot([0,100],[T_half,T_half],lw=2,ls='--',c='k')
    #ax.plot([0,500],[0.25e-6,0.25e-6],lw=2,ls='--',c='k')
    ax.text(70,5.43e18,s='$T_{1/2}=4.43\cdot 10^{17}$ s')
    ax.set_xlabel('Number of barriers $N$')
    ax.set_ylabel('$^{232}$Th  $T_{1/2}$ [s]')

    ax.grid(axis='y',ls='--')    
    plt.show()

main()
