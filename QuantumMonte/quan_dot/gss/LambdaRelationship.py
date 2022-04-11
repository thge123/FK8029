from PLOTTING import *

fig = plt.figure()
ax  = fig.add_subplot(111)

#x = [-3,-2,-1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,10,20,30]
#y = [-7.97,-2.9,0.27,2,2.122,2.24,2.35,2.45,2.55,2.65,2.74,2.83,2.92,3,3.72,4.33,4.86,5.35,5.81,7.42,10.68,13.43]
x = [-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,20,30]
y = [-7.97,-2.9,0.27,2,3,3.72,4.33,4.86,5.35,5.81,6.24,6.65,7.04,7.42,10.68,13.43]


ax.scatter(x,y,facecolor='none',edgecolor='k',s=100)
ax.set_xticks([-5]+[5*j for j in range(7)])
ax.set_xticks([-5+j for j in range(40)],minor=True)
ax.set_yticks([-10,-5]+[5*j for j in range(4)])
ax.set_yticks([-10+j for j in range(25)],minor=True)
ax.set_xlabel(r'$\lambda$')
ax.set_ylabel(r'$E/\hbar\omega$')
ax.grid(which='both')

for i in range(len(x)):
    print('Lam: ',x[i], '  ', 'E: ', y[i])
plt.show()
