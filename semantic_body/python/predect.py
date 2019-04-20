from sklearn.externals import joblib
import numpy as np

if __name__ == "__main__":
    model = joblib.load('../export/randomforest.joblib')
    pca = joblib.load('../export/pca_36c.joblib')


    data = np.fromfile('../data/test/roughExact')[2:]
    print(data.shape)
    data.resize(111, 26)

    pred = model.predict(data)
    out = np.dot(pred, pca.components_)

    print(out.shape)
    out.tofile('../data/recover/dV_pca_rf')
