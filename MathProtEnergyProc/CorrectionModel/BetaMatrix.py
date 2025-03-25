import numpy as np

from MathProtEnergyProc.IndexFunctions import GetIndexes, GetPairFromArray

#Доли распределения теплоты по конкретному процессу
def BetaProcess(betaProcessValues#Величины долей распределения теплоты по энергетически степеня свободы
                ):
    #Приводим к массиву
    arBetaProcessValues = np.array(betaProcessValues, dtype=np.double)
    
    #Получаем и выводим доли распределения теплоты
    return arBetaProcessValues / np.sum(arBetaProcessValues)
    
#Компоновка долей распределения теплоты
class BetaMatrix(object):
    #Конструктор класса
    def __init__(self,
    
                 energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                 processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                 
                 betasProcesses#Доли распределения некопенсированных теплот каждого процесса по энергетически степеня свободы
                 ):
        #Приводим описание рассчитанных долей распределения некомпенсированных теплот каждого процесса к двойному массиву
        energyPowersExtBetaDbNames = []#Доли распределения некомпенсированной теплоты энергетических степеней свободы
        for betaProcess in betasProcesses:
            ProcName = betaProcess[1]#Имя процесса
            EnergyPowersNames = betaProcess[0]#Имя энергетической степени свободы
            for energyPowerName in EnergyPowersNames:#Формируем массив
                energyPowersExtBetaDbNames.append((energyPowerName,ProcName))
        
        #Приводим описание пользовательских долей распределения некомпенсированных теплот каждого процесса к двойному массиву
        energyPowersUserBetaDbNames = GetPairFromArray(energyPowersVarBetaNames,
                                                       processCoordinatesVarBetaNames)
        
        #Получаем массив индексов
        self.__BetaIndexes = GetIndexes(energyPowersExtBetaDbNames,
                                        energyPowersUserBetaDbNames)
    
    #Получаем доли распределения некопенсированных теплот энергетических степеней свободы
    def __call__(self,
    
                 arrBetas#Массив массивов долей распределения некопенсированных теплот
                 ):
        #Конкатенуем и выводим массив долей распределения некомпенсированных теплот
        return np.hstack(arrBetas)[self.__BetaIndexes]
    