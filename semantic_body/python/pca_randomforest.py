import numpy as np

if __name__ == "__main__":
    print("ok")
    # load data
    V = np.fromfile('../data/train/dV')[2:]
    D = np.fromfile('../data/train/roughExact')[2:]

    V.resize(1400, 37500)
    D.resize(1400, 26)


    #PCA 
    from sklearn.decomposition import PCA
    pca = PCA(n_components=36)
    newV = pca.fit_transform(V) #(1400, 36)

    #Random Forest
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import GridSearchCV

    rfr = RandomForestRegressor(n_jobs=8)

    param_grid = {'n_estimators': [1000]}
    model = GridSearchCV(estimator=rfr, param_grid=param_grid, n_jobs=8, cv=10)
    model.fit(D, newV)

    print('Random Forest Regression....')

    print('Best Params:')
    print(model.best_params_)

    print('Best CV Score:')
    print(model.best_score_)

    # Metrics 
    from sklearn import metrics
    pred = np.dot(model.predict(D), pca.components_)

    print("MSE:")
    print(metrics.mean_squared_error(V, pred))


    # Save model
    from sklearn.externals import joblib
    joblib.dump(model, '../export/randomforest.joblib')
    joblib.dump(pca, '../export/pca_36c.joblib')