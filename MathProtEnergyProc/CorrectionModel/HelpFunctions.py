import numpy as np
from scipy.linalg import block_diag

#Функция масштабирования
def ScaleFun(minValue, scCf, value):
    #Выводим масштабированный результат
    return minValue + scCf*value

#Создание блочно-диагональной матрицы
def CreateBlockMatrix(subMatrixes
                      ):#Подматрицы, из которых строим блочную матрицу
    #Выводим сформированную блочную матрицу
    return block_diag(*subMatrixes)
    
#Конкатенация массива величин из массива массивов
def ConcatArrayValues(arrOfArrValues):
    #Выводим результат конкатенации
    return np.hstack(arrOfArrValues).tolist()
