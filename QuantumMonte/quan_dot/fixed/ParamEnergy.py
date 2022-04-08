from PLOTTING import *
import FILES


def main():
    filename = 'A/A001.dat'
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
    ax1.scatter(alpha,AvgE)
    ax2.scatter(alpha,VarE)
    #ax.errorbar(alpha,AvgE,yerr=[j**0.5 for j in VarE],capsize=5,linestyle='none')
    
    plt.show()

main()
