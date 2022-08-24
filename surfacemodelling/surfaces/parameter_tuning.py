import os
import sys

import numpy as np
from pykrige.rk import Krige
from sklearn.model_selection import GridSearchCV


def evaluate_kriging(xyz):
    """Prøve og valg beste parameterer for Kriging

    A scikit-learn tool for parameter tuning by cross-validation
    https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.GridSearchCV.html

    Args:
        xyz (list): [x, y, z]

    Returns:
        estimator.best_params_ (dict): Dictionary med beste metode, variogram modelle, nummer
                                       av "lags" og hvis vekt skal være brukt i små "lag-er"
    """

    print("Evaluating Kriging Parameters")
    sys.stdout = open(os.devnull, "w")
    param_dict = {
        "method": ["ordinary", "universal"],
        "variogram_model": ["linear", "power", "gaussian", "spherical"],
        "nlags": [6, 8, 10, 12],
        "weight": [True, False],
    }

    estimator = GridSearchCV(Krige(), param_dict, verbose=False, return_train_score=True)

    XY = []
    z = []
    for pt in xyz:
        XY.append([pt[0], pt[1]])
        z.append(pt[1])

    XY = np.asarray(XY)
    z = np.asarray(z)

    estimator.fit(X=XY, y=z)

    sys.stdout = sys.__stdout__
    if hasattr(estimator, "best_score_"):
        print("best_score R² = {:.3f}".format(estimator.best_score_))
        print("best_params = ", estimator.best_params_)

    return estimator.best_params_
