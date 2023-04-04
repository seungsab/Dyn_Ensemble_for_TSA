from typing import Union

import numpy as np
import pandas as pd

from src.utils.proportion import normalize_and_proportion


class WindowLoss:
    """
    Windowing dynamic ensemble

    Weighting the ensemble according to recent past performance.
    """

    def __init__(self, lambda_=50, n_burn: int = 1):
        """

        param lambda_: Window size. How many recent past observations are used to estimate performance
        param n_burn: Number of instances to burn to prevent leakage. This should be equal to the forecasting horizon.
        """
        self.lambda_ = lambda_
        self.n_burn = n_burn

    def fit(self, Y_hat: pd.DataFrame, y: Union[pd.Series, np.ndarray]):
        pass

    def get_weights(self, Y_hat: pd.DataFrame, y: pd.Series):
        """
        Computing the weights of each model in the ensemble

        param Y_hat: forecasts of each model in a pd.DF
        param y: actual values to compute the error

        return: The weights of models in the ensemble
        """

        se = Y_hat.apply(func=lambda x: (x - y) ** 2, axis=0)
        rolling_mse = se.rolling(window=self.lambda_).mean()

        window_weights = rolling_mse.apply(func=lambda x: normalize_and_proportion(-x), axis=1)

        window_weights = window_weights[:-self.n_burn]

        eq_weights = pd.DataFrame(np.ones_like(Y_hat.head(1)) / Y_hat.shape[1], columns=Y_hat.columns)

        weights = pd.concat([eq_weights, window_weights], axis=0)
        weights.index = Y_hat.index

        return weights
