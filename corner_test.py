import numpy as np
import corner
import matplotlib.pyplot as plt
import sys

A_d_mean = 6.5352711
A_s_mean = 23.1543445587
b_s_mean = -3.10000

means = [A_s_mean, A_d_mean, b_s_mean]

def getminmax(x):
    mx = np.percentile(x, 97.5)
    mn = np.percentile(x, 2.5)  
    return (mn, mx)
files = sys.argv[1:4]
# burnin = int(sys.argv[-1])
burnin = int(1000)
data = [] 
labels = [] 
for file in files:
    data.append(np.loadtxt(file)[burnin:])
    try:
        labels.append((file.replace("pixel_350_","")).replace(".dat",""))
    except:
        try:
            labels.append(file.replace(".dat",""))
        except:
            labels.append(file)
# ranges=[getminmax(i) for i in data]
ranges = [ (23.0,23.3), (6.53,6.54), (-3.105, -3.095)]
#print(ranges)
samples = np.vstack(data).T
ndim = np.shape(samples)[1]
figure = corner.corner(samples, labels=labels, bins=50, range=ranges)
axes   = np.array(figure.axes).reshape((ndim, ndim))

col = ['r','g','b']

for i in range(ndim):
    ax = axes[i, i]
    ax.axvline(means[i],color=col[i])

for xi in range(ndim):
    for yi in range(xi):
        ax = axes[xi, yi]
        ax.axvline(means[yi],color=col[yi])
        ax.axhline(means[xi],color=col[xi])

plt.show()