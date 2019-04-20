import numpy as np
import matplotlib.pyplot as plt 

rows = 111
cols = 26

rawExact = np.fromfile('../data/test/roughExact')[2:]
rawExact.resize(rows, cols)

prExact = np.fromfile('../data/recover/pr_roughexact')[2:]
prExact.resize(rows, cols)

rsprExact = np.fromfile('../data/recover/rs_pr_roughExact')[2:]
rsprExact.resize(rows, cols)

feature = ["bustline", "waistline", "hipline", "midhipline", "armhole", "head", "collar", "arm", "wrist", "hand", "thigh", "knee", "calf", "ankle", "shoulderWidth", "backWidth", "chestWidth", "breastDist", "chestHeight", "backHeight", "waistLength", "sleeveLength", "frontLength", "backLength", "crotchLength", "outseam"]

x = np.arange(0, cols)

idx = 100
print(rsprExact[idx, :])
plt.plot(x, rawExact[idx, :], 'b--')
plt.plot(x, prExact[idx, :], 'r--')
plt.plot(x, rsprExact[idx, :], 'g--')
plt.show()