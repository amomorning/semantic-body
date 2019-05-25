import pandas as pd

data = pd.read_csv('./data.csv').values
print(data)
fin = data[:, 2]/100.0
print(fin)
from sklearn.externals import joblib
pca = joblib.load('C:\\Users\\amomorning\\source\\repos\\semantic_body\\bodyViz\\Assets\\Exports\\logrs_pca_500c.joblib')

import os
os.environ["CUDA_VISIBLE_DEVICES"]="-1"
from tensorflow.keras.models import load_model
model = load_model('C:\\Users\\amomorning\\source\\repos\\semantic_body\\bodyViz\\Assets\\Exports\\rs_dnn.h5')

import numpy as np

# fin = np.loadtxt('C:\\Users\\amomorning\\source\\repos\\semantic_body\\bodyViz\\Assets\\Resources\\measure.txt')
fin.resize(1, 26)

pred = model.predict(fin)

out = np.dot(pred, pca.components_).astype(np.float64)

out.tofile('C:\\Users\\amomorning\\source\\repos\\semantic_body\\semantic_body\\data\\recover\\testRS')