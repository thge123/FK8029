from PLOTTING import *
from numpy import log

levels = [(-66.85,0),
          (-24.67,1),
          (-20.33,0),
          (-9.666,0),
          (-2.594,0)]

y = []
for i in levels:
    if i[1] == 0:
        y.append(log(-i[0]))
plt.scatter([log(j+1) for j in range(len(y))],y)
plt.show()

