import numpy as np
import matplotlib.pyplot as plt 


rows = 111
cols = 26

rawExact = np.fromfile('../data/test/roughExact')[2:]
rawExact.resize(rows, cols)
print(rawExact)

rsadExact = np.fromfile('../data/recover/testroughExact')[2:]
rsadExact.resize(rows, cols)
print(rsadExact)

plt.ylabel('error')
plt.xticks(range(0, cols), ["bustline", "waistline", "hipline", "midhipline", "armhole", "head", "collar", "arm", "wrist", "hand", "thigh", "knee", "calf", "ankle", "shoulderWidth", "backWidth", "chestWidth", "breastDist", "chestHeight", "backHeight", "waistLength", "sleeveLength", "frontLength", "backLength", "crotchLength", "outseam"]
            ,rotation=90)

def MSE(a, b, ax):
    return ((a-b)**2).mean(axis=ax)

def Mean(a, b, ax):
    return (abs(a-b)).mean(axis=ax)

def median(a, b, ax):
    return np.median(abs(a-b), axis=ax)

from sklearn import metrics
print("MSE dv:")
print(np.sum(rawExact-rsadExact))

x = np.arange(0, cols)

plt.plot(x, MSE(rawExact, rsadExact, 0), label='dV_pca_rf', color='b')

plt.title('Median')
plt.legend(loc='best')
plt.show()

