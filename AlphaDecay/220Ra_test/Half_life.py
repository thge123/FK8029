from PLOTTING import *
import FILES
from Density import get_data 
from os import listdir

def scatter_TR(ax,Dir,filenumber,Marker,Color):
    
    r,omega,Coeffs,V = get_data(Dir,filenumber)
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

    #220Ra Z = 88
    E = 6.405
    v = 0.69*E**0.*1e7    #m/s
    A = 220
    R = 1.5*(4**0.3333 + (A-4)**0.3333)    # fm
    f = v/(2*R) * 1e15
    L = f*T
    T_half = 0.6931/L
    
    print("N = ", filenumber, ": ", "T = ",   T_half)
    
    ax.scatter(filenumber,T_half,edgecolor=Color,facecolor='none',s=100,marker=Marker)
    
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

    Dirs = ['R0=1.5','R0=1.32','R0=1.2']
    Markers = ['o','v','*',]
    Colors  = ['b','r','purple']
    filenames = listdir(Dirs[0]+'/O')
    filenames.sort()
    

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for j in range(len(Dirs)):
        for i in filenames:
            scatter_TR(ax,Dirs[j],FILES.string2number(i),Markers[j],Colors[j])
    ax.plot([0,100],[56,56],lw=1,ls='--',c='k',label='$T_{1/2}=56$ s')
    #ax.text(70,0.585e-6,s='$T_{1/2}=0.25$ $\mu$s')
    ax.set_xlabel('Number of barriers, $N$')
    ax.set_ylabel('$^{220}$Ra  $T_{1/2}$ [s]')

    ax.set_xticks([1,10,20,30,40,50,60,70,80,90,100])
    ax.set_yscale('log')
    ax.grid(axis='y',ls='--')    

    # Dummy
    ax.scatter([1e20],[1e20],marker=Markers[0],edgecolor=Colors[0],facecolor='none',label='$R_0=1.50$ fm')
    ax.scatter([1e20],[1e20],marker=Markers[1],edgecolor=Colors[1],facecolor='none',label='$R_0=1.32$ fm')
    ax.scatter([1e20],[1e20],marker=Markers[2],edgecolor=Colors[2],facecolor='none',label='$R_0=1.20$ fm')
    ax.legend(framealpha=0)
    
    ax.set_xlim(0,101)
    ax.set_ylim(1e-1,1e13)
    plt.show()

main()
