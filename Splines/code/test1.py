from PLOTTING import *
from numpy import array,zeros
from numpy.linalg import solve
from scipy.interpolate import BSpline

def B(x, k, i, t):

   if k == 0:
      return 1.0 if t[i] <= x < t[i+1] else 0.0
   if t[i+k] == t[i]:
      c1 = 0.0
   else:
      c1 = (x - t[i])/(t[i+k] - t[i]) * B(x, k-1, i, t)
   if t[i+k+1] == t[i+1]:
      c2 = 0.0
   else:
      c2 = (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * B(x, k-1, i+1, t)

   return c1 + c2


def main1():

    k=5
    n = 14
    N = n-k-3
    t = [0,0,0,0,0,0] + [10*j/(N+1) for j in range(1,N+1)] + [10,10,10,10,10,10]
    print(len(t))
    print(N)
    
    c = zeros(n)
    xx = [j/100 for j in range(1001)]
    for i in range(n):
        c[i] = 1
        spl = BSpline(t,c,k,extrapolate=False)
        yy = [spl(x) for x in xx]
        plt.plot(xx,yy)
        c = zeros(n)
    plt.show()

def main2():

    k = 4
    n = 10
    N = 12
    t = [0,0,0,0,0] + [10*j/(N+1) for j in range(1,N+1)] + [10,10,10,10,10]
    print(N)
    
    xx = [j/100 for j in range(1001)]
    for i in range(n+7):
        yy = [B(x,k,i,t) for x in xx]
        plt.plot(xx,yy)
    plt.show()
    

main1()
        
    



    
