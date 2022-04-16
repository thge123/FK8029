from PLOTTING import *
import FILES
from numpy import array,exp,log
from numpy import abs as ABS
from seaborn import kdeplot
from numpy.random import normal

psi = {1: lambda t: exp(-t**2/2),
       2: lambda t: t*exp(-t**2/2),
       3: lambda t: (2*t**2-1)*exp(-t**2/2),
       4: lambda t: (2*t**3-3*t)*exp(-t**2/2)
      }

r1 = []
r2 = []
for i in range(1,15000):
    X = FILES.get_data('A/A001.dat',fileline=i)
    r1.append((X[0],X[1]))
    r2.append((X[2],X[3]))
r1 = array(r1)
r2 = array(r2)

gaussian1 = normal(0,0.5**0.5,size=len(r1))
gaussian2 = normal(0,0.5**0.5,size=len(r2))

fig = plt.figure()
ax1  = fig.add_subplot(111)

x1 = array([j[0] for j in r1])
y1 = array([j[1] for j in r1])
x2 = array([j[0] for j in r2])
y2 = array([j[1] for j in r2])
kdeplot(x=x1,y=y1,fill=True,levels=20,ax=ax1)
kdeplot(x=x2,y=y2,fill=True,levels=20,ax=ax1)
#kdeplot(x=gaussian1,y=gaussian2,fill=False,levels=20,ax=ax1)
ax1.set_aspect('equal')

plt.show()




