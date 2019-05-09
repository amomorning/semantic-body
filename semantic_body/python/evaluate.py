import numpy as np
import matplotlib.pyplot as plt 


rows = 111
cols = 26

rawExact = np.fromfile('../data/test/roughExact')[2:]
rawExact.resize(rows, cols)
print(rawExact)

prExact = np.fromfile('../data/recover/pr_roughexact')[2:]
prExact.resize(rows, cols)

rsprExact = np.fromfile('../data/recover/logRS_roughExact')[2:]
rsprExact.resize(rows, cols)
print(rsprExact)

adExact = np.fromfile('../data/recover/dV_aednn_roughexact')[2:]
adExact.resize(rows, cols)

rsadExact = np.fromfile('../data/recover/logRS_aednn_roughexact')[2:]
rsadExact.resize(rows, cols)

plt.xticks(range(0, cols), ["bustline", "waistline", "hipline", "midhipline", "armhole", "head", "collar", "arm", "wrist", "hand", "thigh", "knee", "calf", "ankle", "shoulderWidth", "backWidth", "chestWidth", "breastDist", "chestHeight", "backHeight", "waistLength", "sleeveLength", "frontLength", "backLength", "crotchLength", "outseam"]
            ,rotation=90)

def MSE(a, b, ax):
    return ((a-b)**2).mean(axis=ax)

def Mean(a, b, ax):
    return (abs(a-b)).mean(axis=ax)

def median(a, b, ax):
    return np.median(abs(a-b), axis=ax)

from sklearn import metrics
print(MSE(prExact[:, 0], rawExact[:, 0], 0))
y_dv = median(prExact, rawExact, 0)
y_rs = median(rsprExact, rawExact, 0)
y_dv_ad = median(adExact, rawExact, 0)
y_rs_ad = Mean(rsadExact, rawExact, 0)
print(y_dv.shape)

print("MSE dv:")
print(metrics.mean_squared_error(prExact, rawExact))

print("MSE RS:")
print(metrics.mean_squared_error(rsprExact, rawExact))
print(Mean(rsprExact, rawExact, (0, 1)))
print(Mean(adExact, rawExact, (0, 1)))
print(Mean(rsadExact, rawExact, (0, 1)))
x = np.arange(0, cols)

plt.plot(x, y_dv, label='dV_pca_rf', color='b')
plt.plot(x, y_rs, label='RS_pca_rf', color='g')
plt.plot(x, y_dv_ad, label='dV_ae_dnn', color='c')
# plt.plot(x, y_rs_ad, label='RS_ae_dnn', color='m')
plt.ylabel('error')
# idx = 4
# print(rsprExact[idx, :])
# plt.plot(x, rawExact[idx, :], label='ground truth', color='r')
# plt.plot(x, prExact[idx, :], label='dV', color='b')
# plt.plot(x, rsprExact[idx, :], label='RS', color='g')
# plt.plot(x, adExact[idx, :], label='DNN', color='c')
plt.title('Median')
plt.legend(loc='best')
plt.savefig('../images/chartgraph/Median.png',bbox_inches='tight')
plt.show()