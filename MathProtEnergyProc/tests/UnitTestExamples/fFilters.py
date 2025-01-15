import numpy as np

#Функция условий протекания процессов
def fPosLinearFilter(x, lam):#Тестовый линейный в положительной области фильтр
    #Выводим результат
    return np.log(1 + np.exp(lam*x)) / lam
def fExpFilter(x, lam):#Тестовый экспоненциальный фильтр
    #Выводим результат
    return np.exp(lam*x)
    