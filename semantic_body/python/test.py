import numpy as np
import matplotlib.pyplot as plt

rows = 1400
cols = 26
Exact = np.fromfile('../data/train/roughExact')[2:]
Exact.resize(rows, cols)

y1 = np.amax(Exact, axis=0)
y2 = np.amin(Exact, axis=0)
y3 = (y1-y2)/2
y4 = (y1+y2)/2

plt.xticks(range(0, cols), ["bustline", "waistline", "hipline", "midhipline", "armhole", "head", "collar", "arm", "wrist", "hand", "thigh", "knee", "calf", "ankle", "shoulderWidth", "backWidth", "chestWidth", "breastDist", "chestHeight", "backHeight", "waistLength", "sleeveLength", "frontLength", "backLength", "crotchLength", "outseam"]
            ,rotation=90)
x = np.arange(0, cols)
plt.errorbar(x, y4, yerr=y3, c='r', linewidth=1)
plt.scatter(x, y1, c='b')
plt.scatter(x, y2, c='b')

print(y1)
print(y2)

plt.title('Range of each measurement label')
plt.savefig('../images/chartgraph/measurerange.png')
plt.show()