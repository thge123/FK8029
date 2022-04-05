from PLOTTING import *
import FILES

Colors = {11: 'k',
          12: 'b',
          9: 'purple'}
Markers = {11: 'o',
           12: 'v',
           9: 'P'}
Labels  = {11: r'$\alpha = 5$',
           12: r'$\alpha = 10$',
           9: r'$\alpha = 0$'}
fig = plt.figure()
ax  = fig.add_subplot(111)
for j in [9,11,12]:
    B = 'B/B'+FILES.number2string(j)+'.dat'
    y = []
    for i in range(1,24):
        X = FILES.get_data(B,fileline=i)    # h=0.025 N =400
        y.append(X[0])

    x = [j for j in y[:-1]]
    y = [y[j+1]-y[j] for j in range(len(y)-1)]
    print(y)
    ax.scatter(x,y,c=Colors[j],marker=Markers[j],label=Labels[j])
    ax.plot(x,y,c=Colors[j])
ax.legend(framealpha=0)
ax.set_xlabel(r"$E'_j$",fontsize=32)
ax.set_ylabel(r"$E'_{j+1}-E'_j$",fontsize=32)
plt.show()
