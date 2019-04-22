import numpy as np

V = np.fromfile('../data/train/dV')[2:]
V.resize(1400, 37500)

from tensorflow.keras.layers import Dense
from sklearn.model_selection import train_test_split
x_train, x_test, y_train, y_test = train_test_split(V, V, test_size=0.1, random_state=2)
x_train.shape
