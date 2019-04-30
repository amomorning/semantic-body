
from sklearn.externals import joblib
import numpy as np

model = joblib.load('../export/logrs_randomforest.joblib')
pca = joblib.load('../export/logrs_pca_450c.joblib')

data = []