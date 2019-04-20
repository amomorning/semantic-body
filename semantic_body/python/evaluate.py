import numpy as np
import matplotlib.pyplot as plt 

rows = 111
cols = 26

rawExact = np.fromfile('../data/test/roughExact')[2:]
rawExact.resize(rows, cols)

prExact = np.fromfile('../data/recover/pr_roughexact')[2:]
prExact.resize(rows, cols)

feature = ["bustline", "waistline", "hipline", "midhipline", "armhole", "head", "collar", "arm", "wrist", "hand", "thigh", "knee", "calf", "ankle", "shoulderWidth", "backWidth", "chestWidth", "breastDist", "chestHeight", "backHeight", "waistLength", "sleeveLength", "frontLength", "backLength", "crotchLength", "outseam"]

x = np.arange(0, cols)

plt.plot(x, rawExact[1, :], 'b--')
plt.plot(x, prExact[1, :], 'r--')
plt.show()