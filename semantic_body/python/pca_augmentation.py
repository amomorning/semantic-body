import numpy as np
if __name__ == '__main__':
    print('ok')
    F = np.fromfile('../data/feature')

    f = F[2:]
    f.resize(1511, 112500)

    from sklearn.decomposition import PCA
    pca = PCA(n_components=120)
    newf = pca.fit_transform(f)

    print(np.sum(pca.explained_variance_ratio_))
    print(newf.shape)
    pass