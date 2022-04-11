from PLOTTING import *
import FILES


def main():
    filenumber = int(input("Write filenumber: "))
    filename = 'A/A'+FILES.number2string(filenumber)+'.dat'
    alpha = []
    AvgE  = []
    VarE  = []
    for j in range(1,10):
        X = FILES.get_data(filename,fileline=j)
        alpha.append(X[0])
        AvgE.append(X[1])
        VarE.append(X[2])
    
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.scatter(alpha,AvgE,edgecolor='k',facecolor='none',s=100)
    ax2.scatter(alpha,VarE,edgecolor='k',facecolor='none',s=100)

    ax2.set_xlabel(r'$\alpha$')
    ax1.set_ylabel(r'$\langle E_L\rangle/\hbar\omega$')
    ax2.set_ylabel(r'Var$(E_L/\hbar\omega)$')
    ax1.grid()
    ax2.grid()
    
    
    plt.show()

main()
