import numpy as np

if __name__ == "__main__":
    print("ok")
    # load data
    V = np.fromfile('../data/Verts')
    D = np.fromfile('../data/dijkstra')

    vv = V[2:]
    vv.resize(1511, 37500)

    dd = D[2:]
    dd.resize(1511, 26)
    print(dd)


    #PCA 
    from sklearn.decomposition import PCA
    pca = PCA(n_components=36)
    newvv = pca.fit_transform(vv) #(1511, 36)

    #Random Forest
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import GridSearchCV

    rfr = RandomForestRegressor(n_jobs=8)

    param_grid = {'n_estimators': [1000]}
    model = GridSearchCV(estimator=rfr, param_grid=param_grid, n_jobs=8, cv=10)
    model.fit(dd, newvv)

    print('Random Forest Regression....')

    print('Best Params:')
    print(model.best_params_)

    print('Best CV Score:')
    print(model.best_score_)

    # Metrics 
    from sklearn import metrics
    predvv = np.dot(model.predict(dd), pca.components_)

    print("MSE:")
    print(metrics.mean_squared_error(vv, predvv))


    # Save model
    from sklearn.externals import joblib
    joblib.dump(model, 'model.joblib')
    joblib.dump(pca, 'pca.joblib')