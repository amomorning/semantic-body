if __name__ == "__main__":
    # load data
    import numpy as np
    RS = np.fromfile('../data/train/RS')[2:]
    RS.resize(1400, 225000)

    D = np.fromfile('../data/train/roughExact')[2:]
    D.resize(1400, 26)

    #PCA 
    from sklearn.decomposition import PCA
    pca = PCA(n_components=50)
    RS_pca = pca.fit_transform(RS) #(1400, 50)

    # randomforest
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import GridSearchCV

    rfr = RandomForestRegressor(n_jobs=8)

    param_grid = {'n_estimators': range(10, 1000, 10)}
    model = GridSearchCV(estimator=rfr, param_grid=param_grid, n_jobs=8, cv=10)
    model.fit(D, RS_pca)