import numpy as np

from MathProtEnergyProc.CorrectionModel.HelpFunctions import ScaleFun

#Положительность велиин
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
    return ScaleFun(xmin,scy,rez)
def ExpFilter(x, lam=0.1, xmin=0, scy=1):#Экспоннциальный фильтр
    #Получаем и возвращаем результат
    rez = np.exp(lam*x)
    return ScaleFun(xmin,scy,rez)
