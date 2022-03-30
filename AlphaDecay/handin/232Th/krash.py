from PLOTTING import *

def main():

    fig = plt.figure()
    ax = fig.add_subplot(111)

    x = [1,10,20,30,40,50,60,75,100,125]
    y1 = [6.49,3.87,3.71,3.68,3.58,3.76,3.77,3.72,3.75,3.74]
    y2 = [3.12,5.65,4.23,3.81,4.11,4.30,4.20,4.06,4.16,4.16]
    y1 = [j*1e-12 for j in y1]
    y2 = [j*1e-12 for j in y2]
    
    ax.scatter(x,y1,facecolor='none',edgecolor='b')
    ax.plot(x,y1,color='b')
    ax.scatter(x,y2,facecolor='none',edgecolor='r',marker='v')
    ax.plot(x,y2,color='r')
    plt.show()

main()
