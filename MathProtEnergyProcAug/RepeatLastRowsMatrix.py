import numpy as np


# Повторение последней строки в матрице
def RepeatLastRowsMatrix(matr, nNeedRows):
    # Число строк в матрице
    nRows = matr.shape[0]

    # Размножаем строки
    if nRows < nNeedRows:
        return np.vstack((matr,
                          np.repeat(matr[-1].reshape(1, -1), nNeedRows - nRows, axis=0)))
    else:
        return matr
