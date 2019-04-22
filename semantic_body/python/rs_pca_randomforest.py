if __name__ == "__main__":
    # load data
    import numpy as np
    RS = np.fromfile('../data/train/RS')[2:]
    RS.resize(1400, 225000)

    D = np.fromfile('../data/train/roughExact')[2:]
    D.resize(1400, 26)

    #PCA 
    from sklearn.decomposition import TruncatedSVD
    pca = TruncatedSVD(n_components=60)
    newf = pca.fit_transform(RS)

    # randomforest find best parameter using gridsearch
    from sklearn.ensemble import RandomForestRegressor

    rfr0 = RandomForestRegressor(n_estimators=700, n_jobs=8, oob_score=True, random_state=10)
    #rfr0 = RandomForestRegressor( oob_score=True, random_state=10)
    rfr0.fit(D, newf)

    # evaluate
    from sklearn import metrics
    predf = rfr0.predict(D)
    pred = np.dot(predf, pca.components_)

    print("MSE:")
    print(metrics.mean_squared_error(predf, newf))

    print("PCA_MSE:")
    print(metrics.mean_squared_error(RS, pred))

    # predict
    data = np.fromfile('../data/test/roughExact')[2:]
    print(data.shape)
    data.resize(111, 26)

    pred = rfr0.predict(data)
    out = np.dot(pred, pca.components_)

    print(out.shape)
    # np.savetxt('../data/newRS.txt', out, delimiter=' ')
    out.tofile('../data/recover/newRS')

    # save model
    from sklearn.externals import joblib
    joblib.dump(rfr0, '../export/rs_randomforest.joblib')
    joblib.dump(pca, '../export/rs_pca_60c.joblib')