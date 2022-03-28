from PLOTTING import *

def main():

    fig = plt.figure()
    ax = fig.add_subplot(111)

    x = [1,10,20,30,40,50,60,70,80,90,100,110,120]
    y1 = [4.29e-9,1.1e-7,2.03e-7,2.4e-7,2.57e-7,2.67e-7,2.72e-7,2.75e-7,2.78e-7,2.79e-7,2.81e-7,2.82e-7,2.83e-7]    # N = 5 half lifes
    y2 = [8.75e-10,3.12e-9,4.11e-9,4.68e-9,4.98e-9,5.19e-9,5.3e-9,5.39e-9,5.45e-9,5.50e-9,5.53e-9,5.56e-9,5.58e-9]    # N = 10 half lifes
    y3 = [7.2e-10,4e-10,5e-10,6e-10
    # Transmission coeff. y1 = [6.49,3.87,3.71,3.68,3.58,3.76,3.77,3.72,3.75,3.74]    # 50
    # Transmission coeff. y2 = [3.12,5.65,4.23,3.81,4.11,4.30,4.20,4.06,4.16,4.16]    # 100
    #y1 = [j*1e-12 for j in y1]
    #y2 = [j*1e-12 for j in y2]
    
    ax.scatter(x,y1,facecolor='none',edgecolor='b',label='$N=5$',s=100)
    ax.scatter(x,y2,facecolor='none',edgecolor='r',marker='v',label='$N=10$',s=100)
    ax.set_yscale('log')
    ax.legend(framealpha=0)
    #ax.plot(x,y,color='b')
    #ax.scatter(x,y2,facecolor='none',edgecolor='r',marker='v')
    #ax.plot(x,y2,color='r')
    ax.set_xlabel('Number of intervals after barrier, $M$')
    ax.set_ylabel('$^{212}$Po $T_{1/2}$ [s]')
    ax.grid(axis='y',ls='--')
    plt.show()

main()
