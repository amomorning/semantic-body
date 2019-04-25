import numpy as np

RS = np.fromfile('../data/train/logRS')[2:]
RS.resize(1400, 300000)

nc = 450
from sklearn.decomposition import TruncatedSVD
pca = TruncatedSVD(n_components=nc)
newf = pca.fit_transform(RS)

D = np.fromfile('../data/train/roughExact')[2:]
D.resize(1400, 26)

from sklearn.ensemble import RandomForestRegressor

rfr0 = RandomForestRegressor(n_estimators=1000, n_jobs=8, oob_score=True, random_state=10)
#rfr0 = RandomForestRegressor( oob_score=True, random_state=10)
rfr0.fit(D, newf)

from sklearn import metrics
predf = rfr0.predict(D)
pred = np.dot(predf, pca.components_)

print("MSE:")
print(metrics.mean_squared_error(predf, newf))

print("PCA_MSE:")
print(metrics.mean_squared_error(RS, pred))

from sklearn.externals import joblib
joblib.dump(rfr0, '../export/logrs_randomforest.joblib')
joblib.dump(pca, '../export/logrs_pca_450c.joblib')


data = np.fromfile('../data/test/roughExact')[2:]
data.resize(111, 26)

pred = rfr0.predict(data)
out = np.dot(pred, pca.components_)

out.tofile('../data/recover/newlogRS')