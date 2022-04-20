from PLOTTING import *
from numpy import array,zeros,pi,cos
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

def F(k,N):
    phi = (2*k-1)*pi/(2*N)
    return 5+5*cos(phi)

def main1():

    K = 3
    N = 10
    #X = [10*j/(N+1) for j in range(1,N+1)]
    X = [F(j,N) for j in range(1,N+1)]
    X.sort()
    print(X)
    t = (K+1)*[0] + X + (K+1)*[10]
    n = len(t)-K-1
    print(len(t))
    print(N)
    
    c = zeros(n)
    xx = [j/100 for j in range(1001)]
    for i in range(n):
        c[i] = 1
        spl = BSpline(t,c,K,extrapolate=False)
        yy = [spl(x) for x in xx]
        plt.plot(xx,yy)
        c[i] = 0
    plt.show()

main1()
        
    



    
