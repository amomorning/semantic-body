from sklearn.externals import joblib
import numpy as np

if __name__ == "__main__":
    model = joblib.load('model.joblib')
    pca = joblib.load('pca.joblib')

    data = [[1.47614, 1.40532, 1.44258, 1.47699, 0.761712, 0.650765, 0.52607, 0.500625, 0.23843, 0.29884, 0.801796, 0.579924, 0.54716, 0.331767, 0.47362, 0.45998, 0.467646, 0.318187, 0.34539, 0.377015, 0.380852, 0.642929, 0.463352, 0.391945, 0.328395, 1.121300]]

    pred = model.predict(data)
    out = np.dot(pred, pca.components_)

    print(out.shape)
    np.savetxt('new.txt', out, delimiter=' ')