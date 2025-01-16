import numpy as np
from scipy.linalg import block_diag

#Положительность необратимых составляющих
def ReluFilter(x):
    ax = np.abs(x)#Модуль аргументов
    
    #Обнуление отрицательных аргументов и вывод результата
    return (x + ax) / 2
def PosLinearFilter(x, lam=1000, xmin=0, scy=1):#Линейный в положительной области фильтр
    #Рассчитываем составляющие
    ax = np.abs(x)#Модуль аргументов
    px = (ax + x) / 2#Обнуление отрицательных аргументов
    
    #Получаем и возвращаем результат
    rez = np.log(1 + np.exp(-lam*ax)) / lam + px
    return xmin + scy*rez
def ExpFilter(x, lam=0.1, xmin=0, scy=1):#Экспоннциальный фильтр
    #Получаем и возвращаем результат
    rez = np.exp(lam*x)
    return xmin + scy*rez

#Функции кинетических матриц из их обратимых и необратимых составляющих
def KineticMatrixFromPosSubMatrix(posSubMatrixes,#Положительные определенные составляющие атрицы
                                  subMatrixesBalance#Податрицы баланса
                                  ):
    #Получаем блочно-диагональную кинетическую матрицу
    blockKineticMatrix = block_diag(*posSubMatrixes)
    blockKineticMatrix = np.array(blockKineticMatrix, dtype=np.double)#Приводим к типу вещественных чисел повышенной точности
    
    #Получаем матрицу обратимых состааляющих
    MatrixBalance = np.hstack(subMatrixesBalance, dtype=np.double)
    
    #Выводим результат
    return np.dot(np.dot(MatrixBalance, blockKineticMatrix), MatrixBalance.transpose())
def KineticMatrixSymAsym(irRevCom,#Необратимые компоненты
                         symRevCom,#Симметричные обратимые компоненты
                         asymRevCom#Антисимметричные обратимые компоненты
                         ):
    #Выводим результат
    return KineticMatrixFromPosSubMatrix(irRevCom, symRevCom) + asymRevCom - asymRevCom.transpose()
