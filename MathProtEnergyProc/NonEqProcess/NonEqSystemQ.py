import numpy as np
import scipy.sparse as sp

import MathProtEnergyProc.IndexFunctions as indf

#Класс физико-химических систем
class NonEqSystemQ(object):
    __Tbase = 293;#К
    
    #Конструктор класса
    def __init__(self,
                 stateCoordinatesNames,#Имена координат состояния
                 processCoordinatesNames,#Имена координат процессов
                 energyPowersNames,#Имена энергетических степеней свободы
                 reducedTemperaturesEnergyPowersNames,#Имена приведенных температур энергетических степеней свободы
                 energyPowersBetNames,#Имена взаимодействий между энергетическими степенями свободы
                 heatTransfersNames,#Имена потоков переноса теплоты
                 heatTransfersOutputEnergyPowersNames,#Имена энергетических степеней свободы, с которых уходит теплота
                 heatTransfersInputEnergyPowersNames,#Имена энергетических степеней свободы, на которые приходит теплота
                 
                 stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                 heatEnergyPowersStreamsNames,#Имена потоков теплоты на энергетические степени свободы
                 
                 stateFunction,#Функция состояния
                 
                 stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                 processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                 energyPowersVarTemperatureNames,#Имена переменных температур энергетических степеней свободы
                 stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                 energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                 stateCoordinatesVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                 energyPowersVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                 energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                 processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                 reducedTemperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                 energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                 reducedTemperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                 stateCoordinatesVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                 
                 varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                 varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                 varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                 varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                 varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                 varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                 varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                 varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                 
                 stateCoordinatesVarStreamsNames,#Имена переменных внешних потоков
                 heatEnergyPowersVarStreamsNames#Имена переменных внешних потоков теплоты
                 ):
        #Заполняем имена
        self.__StateCoordinatesNames = stateCoordinatesNames.copy()#Имена координат состояния
        self.__ProcessCoordinatesNames = processCoordinatesNames.copy()#Имена координат процессов
        self.__EnergyPowersNames = energyPowersNames.copy()#Имена энергетических степеней свободы
        self.__ReducedTemperaturesEnergyPowersNames = reducedTemperaturesEnergyPowersNames.copy()#Имена приведенных температур энергетических степеней свободы
        self.__EnergyPowersBetNames = energyPowersBetNames.copy()#Имена взаимодействий между энергетическими степенями свободы
        self.__HeatTransfersNames = heatTransfersNames.copy()#Имена потоков переноса теплоты
        
        #Задаем функцию состояния
        self.__StateFunction = stateFunction#Функция состояния
        
        #Формируем матрицу баланса
        nStateCoordinates = len(stateCoordinatesNames)#Число координат состояния
        nProcessCoordinates = len(processCoordinatesNames)#Число координат процессов
        self.__BalanceStateCoordinates = np.zeros((nStateCoordinates,nProcessCoordinates), dtype=np.double)#нулевая матрица баланса
        
        #Массивы индексов переменных элементов матрицы баланса
        self.__StateCoordinatesVarBalanceIndexes = indf.GetIndexes(self.__StateCoordinatesNames,
                                                                   stateCoordinatesVarBalanceNames)#Индексы переменных коэффициентов матрицы баланса по координатам состояния
        self.__ProcessCoordinatesVarBalanceIndexes = indf.GetIndexes(self.__ProcessCoordinatesNames,
                                                                     processCoordinatesVarBalanceNames)#Индексы переменных коэффициентов матрицы баланса по координатам процессов
        
        #Формируем нулевой вектор скоростей изменения координат состояния
        self.__vStateCoordinates = np.zeros((nStateCoordinates,), dtype=np.double)
        
        #Формируем нулевой вектор скоростей изменения координат процессов
        self.__vProcessCoordinates = np.zeros((nProcessCoordinates,), dtype=np.double)
        
        #Формируем индексы внешних потоков
        self.__StateCoordinatesStreamsIndexes = indf.GetIndexes(self.__StateCoordinatesNames,
                                                                stateCoordinatesStreamsNames)#Индексы внешних потоков по координатам состояния
        
        #Формируем нулевой вектор внешних потоков
        nStateCoordinatesStreams = len(stateCoordinatesStreamsNames)#Число внешних потоков
        self.__StateCoordinatesStreams = np.zeros((nStateCoordinatesStreams,), dtype=np.double)
        
        #Формируем индексы внешних потоков теплоты
        self.__HeatEnergyPowersStreamsIndexes = indf.GetIndexes(self.__EnergyPowersNames,
                                                                heatEnergyPowersStreamsNames)#Индексы внешних потоков по координатам состояния
        
        #Формируем нулевой вектор внешних потоков теплоты
        nHeatEnergyPowersStreams = len(heatEnergyPowersStreamsNames)#Число внешних потоков теплоты
        self.__HeatEnergyPowersStreams = np.zeros((nHeatEnergyPowersStreams,), dtype=np.double)
        
        #Формируем индексы переменных внешних потоков
        self.__StateCoordinatesVarStreamsIndexes = indf.Get2Indexes(self.__StateCoordinatesNames,
                                                                    self.__StateCoordinatesStreamsIndexes,
                                                                    stateCoordinatesVarStreamsNames)#Индексы внешних потоков по координатам состояния
        
        #Формируем индексы переменных внешних потоков теплоты
        self.__HeatEnergyPowersVarStreamsIndexes = indf.Get2Indexes(self.__EnergyPowersNames,
                                                                    self.__HeatEnergyPowersStreamsIndexes,
                                                                    heatEnergyPowersVarStreamsNames)#Индексы внешних потоков теплоты энергетическим степеням свободы
        
        #Формируем матрицу переноса теплот
        nEnergyPowers = len(energyPowersNames)#Число энергетических степеней свободы
        nHeatTransfers = len(heatTransfersNames)#Число потоков переноса теплоты
        heatTransfersIndexes = np.arange(nHeatTransfers)#Индексы энергетических степеней свободы
        heatTransfersIndexes = np.hstack((heatTransfersIndexes,heatTransfersIndexes))#Индексы энергетических степеней свободы
        energyPowersIndexesOut = indf.GetIndexes(self.__EnergyPowersNames,
                                                 heatTransfersOutputEnergyPowersNames)#Индексы энергетических степеней свободы
        energyPowersIndexesIn = indf.GetIndexes(self.__EnergyPowersNames,
                                                heatTransfersInputEnergyPowersNames)#Индексы энергетических степеней свободы
        heatOut = -np.ones_like(energyPowersIndexesOut)#Индексы ухода теплоты с энергетических степепенй свободы
        heatIn = np.ones_like(energyPowersIndexesIn)#Индексы ухода теплоты с энергетических степепенй свободы
        self.__HeatTransferMatrix = sp.coo_matrix((np.hstack((heatOut,heatIn)),
                                                   (energyPowersIndexesOut + energyPowersIndexesIn,
                                                    heatTransfersIndexes)), 
                                                  shape=(nEnergyPowers,nHeatTransfers))#Матрица переноса теплот
        
        #Формируем массив температур
        self.__TEnergyPowers = NonEqSystemQ.__Tbase*np.ones((nEnergyPowers,), dtype=np.double)
        
        #Формируем массив скоростей сообщения теплот энергетических степеней свободы
        self.__vHeatPower = np.zeros((nEnergyPowers,), dtype=np.double)
        
        #Формируем массив скоростей приращения внутренних энергий энергетических степеней свободы
        self.__vUPower = np.zeros((nEnergyPowers,), dtype=np.double)
        
        #Формируем массив скоростей приращения приведенных температур
        nReducedTemperaturesEnergyPowers = len(reducedTemperaturesEnergyPowersNames)
        self.__vReducedTemperaturesEnergyPowers = np.zeros((nReducedTemperaturesEnergyPowers,), dtype=np.double)
        
        #Формируем массив скоростей передачи теплот между энергетическими степенями свободы
        self.__vHeatTransfer = np.zeros((nHeatTransfers,), dtype=np.double)
        
        #Формируем индексы переменных температур
        self.__EnergyPowersVarTemperatureIndexes = indf.GetIndexes(self.__EnergyPowersNames,
                                                                   energyPowersVarTemperatureNames)#Индексы переменных температур энергетических степеней свободы
        
        #Формируем матрицу потениалов взаимодействия для энергетических степеней свободы
        self.__PotentialsInter = np.zeros((nEnergyPowers,nStateCoordinates), dtype=np.double)
        
        #Формируем индексы переменных потенциалов взаимодействия
        self.__EnergyPowersVarPotentialsInterIndexes = indf.GetIndexes(self.__EnergyPowersNames,
                                                                       energyPowersVarPotentialsInterNames)#Индексы переменных потенциалов взаимодействия по энергетическим степенәм свободы
        self.__StateCoordinatesVarPotentialsInterIndexes = indf.GetIndexes(self.__StateCoordinatesNames,
                                                                           stateCoordinatesVarPotentialsInterNames)#Индексы переменных потенциалов взаимодействия по координатам состояния
        
        #Формируем матрицу потениалов взаимодействия между энергетическими степенями свободы
        nEnergyPowersBet = len(energyPowersBetNames)#Число взаимодействий между энергетическими степенями свободы
        self.__PotentialsInterBet = np.zeros((nEnergyPowersBet,nStateCoordinates), dtype=np.double)
        
        #Формируем индексы переменных потенциалов взаимодействия между энергетическими степенями свободы
        self.__EnergyPowersVarPotentialsInterBetIndexes = indf.GetIndexes(self.__EnergyPowersBetNames,
                                                                          energyPowersVarPotentialsInterBetNames)#Индексы переменных потенциалов взаимодействия по энергетическим степенәм свободы
        self.__StateCoordinatesVarPotentialsInterBetIndexes = indf.GetIndexes(self.__StateCoordinatesNames,
                                                                              stateCoordinatesVarPotentialsInterBetNames)#Индексы переменных потенциалов взаимодействия по координатам состояния
        
        #Формируем матрицу потениалов взаимодействия для энергетических степеней свободы
        self.__Beta = np.zeros((nProcessCoordinates,nEnergyPowers), dtype=np.double)
        
        #Формируем индексы переменных потенциалов взаимодействия
        self.__EnergyPowersVarBetaIndexes = indf.GetIndexes(self.__EnergyPowersNames,
                                                            energyPowersVarBetaNames)#Индексы переменных потенциалов взаимодействия по энергетическим степенәм свободы
        self.__ProcessCoordinatesVarBetaIndexes = indf.GetIndexes(self.__ProcessCoordinatesNames,
                                                                  processCoordinatesVarBetaNames)#Индексы переменных потенциалов взаимодействия по координатам состояния
                                                                                   
        #Формируем массивы термодинамических сил
        self.__HeatAffinity = np.zeros((nHeatTransfers,), dtype=np.double)#Термодинамические силы по переносу теплоты между энергетическими степенями свободы
        self.__Affinity = np.zeros((nProcessCoordinates,), dtype=np.double)#Термодинамические силы
        
        #Формируем массив скоростей выделениә некомпенсированных теплот
        self.__vHeatProcess = np.zeros((nProcessCoordinates,), dtype=np.double)
        
        #Формируем главный блок кинетической матрицы по координатам процессов
        self.__KineticMatrixPCPC = np.zeros((nProcessCoordinates,nProcessCoordinates), dtype=np.double)
        
        #Формируем индексы главного блока кинетической матрицы по координатам процессов
        self.__ProcessCoordinatesVarKineticPCPCIndexes = indf.GetIndexes(self.__ProcessCoordinatesNames,
                                                                         varKineticPCPCNames)#Индексы сопряженностей между собой координат процессов
        self.__ProcessCoordinatesVarKineticPCPCAffIndexes = indf.GetIndexes(self.__ProcessCoordinatesNames,
                                                                            varKineticPCPCAffNames)#Индексы сопряженностей между собой термодинамических сил
        
        #Формируем перекрестный блок кинетической матрицы по координатам процессов и перенесенным теплотам
        self.__KineticMatrixPCHeat = np.zeros((nProcessCoordinates,nHeatTransfers), dtype=np.double)
        
        #Формируем индексы перекрестного блока кинетической матрицы по координатам процессов и перенесенным теплотам
        self.__ProcessCoordinatesVarKineticPCHeatIndexes = indf.GetIndexes(self.__ProcessCoordinatesNames,
                                                                           varKineticPCHeatNames)#Индексы сопряженностей координат процессов с теплопереносами
        self.__ProcessCoordinatesVarKineticPCHeatAffIndexes = indf.GetIndexes(self.__HeatTransfersNames,
                                                                              varKineticPCHeatAffNames)#Индексы сопряженностей термодинамических сил с теплопереносами
        
        
        #Формируем перекрестный блок кинетической матрицы по перенесенным теплотам и координатам процессов
        self.__KineticMatrixHeatPC = np.zeros((nHeatTransfers,nProcessCoordinates), dtype=np.double)
        
        #Формируем индексы перекрестного блока кинетической матрицы по координатам процессов и перенесенным теплотам
        self.__ProcessCoordinatesVarKineticHeatPCIndexes = indf.GetIndexes(self.__HeatTransfersNames,
                                                                           varKineticHeatPCNames)#Индексы сопряженностей теплопереносов с координатами процессов
        self.__ProcessCoordinatesVarKineticHeatPCAffIndexes = indf.GetIndexes(self.__ProcessCoordinatesNames,
                                                                              varKineticHeatPCAffNames)#Индексы сопряженностей теплопереносов с термодинамическими силами
        
        #Формируем главный блок кинетической матрицы по перенесенным теплотам
        self.__KineticMatrixHeatHeat = np.zeros((nHeatTransfers,nHeatTransfers), dtype=np.double)
        
        #Формируем индексы главного блока кинетической матрицы по перенесенным теплотам
        self.__ProcessCoordinatesVarKineticHeatHeatIndexes = indf.GetIndexes(self.__HeatTransfersNames,
                                                                             varKineticHeatHeatNames)#Индексы сопряженностей между собой перенесенных теплот
        self.__ProcessCoordinatesVarKineticHeatHeatAffIndexes = indf.GetIndexes(self.__HeatTransfersNames,
                                                                                varKineticHeatHeatAffNames)#Индексы сопряженностей между собой термодинамических сил по переносу теплот
        
        #Формируем матрицу обратных теплоемкостей
        self.__InvHeatCapacityMatrix = np.zeros((nReducedTemperaturesEnergyPowers,nEnergyPowers), dtype=np.double)
        
        #Формируем индексы переменных коэффициентов матрицы обратных теплоемкостей
        self.__ReducedTemperaturesEnergyPowersVarInvHeatCapacityIndexes = indf.GetIndexes(self.__ReducedTemperaturesEnergyPowersNames,
                                                                                          reducedTemperaturesEnergyPowersVarInvHeatCapacityNames)
        self.__EnergyPowersVarInvHeatCapacity = indf.GetIndexes(self.__EnergyPowersNames,
                                                                energyPowersVarInvHeatCapacityNames)
        
        #Формируем матрицу приведенных тепловых эффектов
        self.__HeatEffectMatrix = np.zeros((nReducedTemperaturesEnergyPowers,nStateCoordinates), dtype=np.double)
        
        #Формируем индексы переменных коэффициентов матрицы приведенных тепловых эффектов
        self.__ReducedTemperaturesEnergyPowersVarHeatEffectIndexes = indf.GetIndexes(self.__ReducedTemperaturesEnergyPowersNames,
                                                                                     reducedTemperaturesEnergyPowersVarHeatEffectNames)#Индексы переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        self.__StateCoordinatesVarHeatEffectIndexes = indf.GetIndexes(self.__StateCoordinatesNames,
                                                                      stateCoordinatesVarHeatEffectNames)#Индексы переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
    
    #Получание опорной температуры
    def GetTbase():
        return NonEqSystemQ.__Tbase
    
    #Задание постоянных элементов в матрицы баланса
    def SetBalanceStateCoordinatesConstElement(self,
                                               stateCoordinateName,#Имя координаты состояния
                                               processCoordinateName,#Имя координаты процесса
                                               elementValue#Значение элемента
                                               ):
        #Получаем индексы
        stateCoordinateNameInd = indf.GetIndex(self.__StateCoordinatesNames,
                                               stateCoordinateName)#Индекс координаты состояния
        processCoordinateNameInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                                 processCoordinateName)#Индекс координаты процесса
        
        #Задаем элемент
        self.__BalanceStateCoordinates[stateCoordinateNameInd, processCoordinateNameInd] = elementValue
    
    #Получение элементов матрицы баланса
    def GetBalanceStateCoordinatesElement(self,
                                          stateCoordinateName,#Имя координаты состояния
                                          processCoordinateName#Имя координаты процесса
                                          ):
        #Получаем индексы
        stateCoordinateInd = indf.GetIndex(self.__StateCoordinatesNames,
                                           stateCoordinateName)#Индекс координаты состояния
        processCoordinateInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                             processCoordinateName)#Индекс координаты процесса
        
        #Задаем элемент
        return self.__BalanceStateCoordinates[stateCoordinateInd, processCoordinateInd]
        
    #Получение матрицы баланса
    def GetBalanceMatrix(self):
        #Выводим матрицу баланса
        return self.__BalanceStateCoordinates.copy()
            
    #Задание постоянных внешних потоков
    def SetStateCoordinatesStreamsConstElement(self,
                                               stateCoordinateStreamName,#Имя координаты состояния
                                               elementValue#Значение элемента
                                               ):
        #Получаем индексы
        stateStreamInd = indf.Get2Index(self.__StateCoordinatesNames,
                                        self.__StateCoordinatesStreamsIndexes,
                                        stateCoordinateStreamName)#Индекс внешнего потока
        
        #Задаем элемент
        self.__StateCoordinatesStreams[stateStreamInd] = elementValue
            
    #Получение постоянных внешних потоков
    def GetStateCoordinatesStreamsElement(self,
                                          stateCoordinateStreamName#Имя координаты состояния
                                          ):
        #Получаем индексы
        stateStreamInd = indf.Get2Index(self.__StateCoordinatesNames,
                                        self.__StateCoordinatesStreamsIndexes,
                                        stateCoordinateStreamName)#Индекс координаты процесса
        
        #Выводим элемент
        return self.__StateCoordinatesStreams[stateStreamInd]
            
    #Получение массива внешних потоков
    def GetStateCoordinatesStreams(self):
        #Выводим внешние потоки
        return self.__StateCoordinatesStreams.copy()
                    
    #Получение массива внешних потоков
    def GetHeatTransferMatrix(self):
        #Выводим внешние потоки
        return self.__HeatTransferMatrix.toarray().copy()
            
    #Задание постоянных внешних потоков
    def SetHeatEnergyPowersStreamsConstElement(self,
                                               heatEnergyPowersStreamName,#Имя энергетической степени свободы
                                               elementValue#Значение элемента
                                               ):
        #Получаем индексы
        heatStreamInd = indf.Get2Index(self.__EnergyPowersNames,
                                       self.__HeatEnergyPowersStreamsIndexes,
                                       heatEnergyPowersStreamName)#Индекс внешнего потока теплоты
        
        #Задаем элемент
        self.__HeatEnergyPowersStreams[heatStreamInd] = elementValue
            
    #Получение постоянных внешних потоков
    def GetHeatEnergyPowersStreamsElement(self,
                                          heatEnergyPowersStreamName,#Имя энергетической степени свободы
                                          ):
        #Получаем индексы
        heatStreamInd = indf.Get2Index(self.__EnergyPowersNames,
                                       self.__HeatEnergyPowersStreamsIndexes,
                                       heatEnergyPowersStreamName)#Индекс внешнего потока теплоты
        
        #Выводим элемент
        return self.__HeatEnergyPowersStreams[heatStreamInd]
            
    #Получение массива внешних потоков
    def GetHeatEnergyPowersStreams(self):
        #Выводим внешние потоки
        return self.__HeatEnergyPowersStreams.copy()
              
    #Получение массива имен координат состояния
    def GetStateCoordinatesNames(self):
        #Выводим имена координат состояния
        return self.__StateCoordinatesNames.copy()
           
    #Получение массива скоростей изменения координат состояния
    def GetVStateCoordinates(self):
        #Выводим скорости изменения координат состояния
        return self.__vStateCoordinates.copy()
            
    #Получение массива имен координат процессов
    def GetProcessCoordinatesNames(self):
        #Выводим имена координат процессов
        return self.__ProcessCoordinatesNames.copy()
            
    #Получение массива скоростей изменения координат процессов
    def GetVProcessCoordinates(self):
        #Выводим скорости изменения координат процессов
        return self.__vProcessCoordinates.copy()
                         
    #Получение массива термодинамических сил
    def GetAffinity(self):
        #Выводим массив термодинамических сил
        return self.__Affinity.copy()
             
    #Получение массива имен энергетических степеней свободы
    def GetEnergyPowerNames(self):
        #Выводим имена энергетических степеней свободы
        return self.__EnergyPowersNames.copy()
             
    #Получение массива имен энергетических степеней свободы
    def GetReducedTemperaturesEnergyPowersNames(self):
        #Выводим имена энергетических степеней свободы
        return self.__ReducedTemperaturesEnergyPowersNames.copy()
    
    #Получение массива имен взаимодействий между энергетическими степенями свободы
    def GetEnergyPowerBetNames(self):
        #Выводим имена взаимодействий между энергетическими степенями свободы
        return self.__EnergyPowersBetNames.copy()
             
    #Получение массива скоростей сообщения теплоты энергетическим степеням свободы
    def GetVHeatPower(self):
        #Выводим скорости сообщения теплоты энергетическим степеням свободы
        return self.__vHeatPower.copy()
                
    #Получение массива скоростей выделения некомпенсированных теплот
    def GetVHeatProcess(self):
        #Выводим скорости выделения некомпенсированных теплот
        return self.__vHeatProcess.copy()
             
    #Получение массива имен процессов переноса теплоты
    def GetHeatTransfersNames(self):
        #Выводим имена процессов переноса теплоты
        return self.__HeatTransfersNames.copy()
             
    #Получение массива скоростей процессов переноса теплоты
    def GetVHeatTransfers(self):
        #Выводим скорости процессов переноса теплоты
        return self.__vHeatTransfer.copy()
                         
    #Получение массива термодинамических сил переноса теплоты
    def GetHeatAffinity(self):
        #Выводим массив термодинамических сил переноса теплоты
        return self.__HeatAffinity.copy()
               
    #Задание постоянных температур энергетических степеней свободы
    def SetTEnergyPowersConstElement(self,
                                     energyPowerName,#Имя энергетической степени свободы
                                     elementValue#Значение элемента
                                     ):
        #Получаем индексы
        energyPowerInd = indf.GetIndex(self.__EnergyPowersNames,
                                       energyPowerName)#Индекс энергетической степени свободы
        
        #Задаем элемент
        self.__TEnergyPowers[energyPowerInd] = elementValue
            
    #Получение температур энергетических степеней свободы
    def GetTEnergyPowersElement(self,
                                energyPowerName#Имя энергетической степени свободы
                                ):
        #Получаем индексы
        energyPowerInd = indf.GetIndex(self.__EnergyPowersNames,
                                       energyPowerName)#Индекс энергетической степени свободы
        
        #Задаем элемент
        return self.__TEnergyPowers[energyPowerInd]
            
    #Получение массива температур энергетических степеней свободы
    def GetTEnergyPowers(self):
        #Задаем элемент
        return self.__TEnergyPowers.copy()
             
    #Получение массива скоростей изменения внутренних энергий энергетических степеней свободы
    def GetVUPower(self):
        #Выводим скорости изменения внутренних энергий энергетических степеней свободы
        return self.__vUPower.copy()
             
    #Получение массива скоростей изменения приведенных температур энергетических степеней свободы
    def GetVReducedTemperaturesEnergyPowers(self):
        #Выводим скорости изменения приведенных температур энергетических степеней свободы
        return self.__vReducedTemperaturesEnergyPowers.copy()
            
    #Задание постоянных потенцталов взаимодействия энергетических степеней свободы
    def SetPotentialsInterConstElement(self,
                                       energyPowerName,#Имя энергетической степени свободы
                                       stateCoordinateName,#Имя координаты состояния
                                       elementValue#Значение элемента
                                       ):
        #Получаем индексы
        energyPowerInd = indf.GetIndex(self.__EnergyPowersNames,
                                       energyPowerName)#Индекс энергетической степени свободы
        stateCoordinateInd = indf.GetIndex(self.__StateCoordinatesNames,
                                           stateCoordinateName)#Индекс координаты состояния
        
        #Задаем элемент
        self.__PotentialsInter[energyPowerInd,stateCoordinateInd] = elementValue
            
    #Получение потенциалов взаимодействия энергетических степеней свободы
    def GetPotentialsInterElement(self,
                                  energyPowerName,#Имя энергетической степени свободы
                                  stateCoordinateName#Имя координаты состояния
                                  ):
        #Получаем индексыиндексы
        energyPowerInd = indf.GetIndex(self.__EnergyPowersNames,
                                       energyPowerName)#Индекс энергетической степени свободы
        stateCoordinateInd = indf.GetIndex(self.__StateCoordinatesNames,
                                           stateCoordinateName)#Индекс координаты состояния
        
        #Задаем элемент
        return self.__PotentialsInter[energyPowerInd,stateCoordinateInd]
            
    #Получение массива потенциалов взаимодействия энергетических степеней свободы
    def GetPotentialsInter(self):
        #Задаем элемент
        return self.__PotentialsInter.copy()
            
    #Задание постоянных потенцталов взаимодействия между энергетическими степенями свободы
    def SetPotentialsInterBetConstElement(self,
                                          energyPowerBetName,#Имя энергетической степени свободы
                                          stateCoordinateName,#Имя координаты состояния
                                          elementValue#Значение элемента
                                          ):
        #Получаем индексы
        energyPowerBetInd = indf.GetIndex(self.__EnergyPowersBetNames,
                                          energyPowerBetName)#Индекс энергетической степени свободы
        stateCoordinateInd = indf.GetIndex(self.__StateCoordinatesNames,
                                           stateCoordinateName)#Индекс координаты состояния
        
        #Задаем элемент
        self.__PotentialsInterBet[energyPowerBetInd,stateCoordinateInd] = elementValue
            
    #Получение температур энергетических степеней свободы
    def GetPotentialsInterBetElement(self,
                                     energyPowerBetName,#Имя энергетической степени свободы
                                     stateCoordinateName#Имя координаты состояния
                                     ):
        #Получаем индексыиндексы
        energyPowerBetInd = indf.GetIndex(self.__EnergyPowersBetNames,
                                          energyPowerBetName)#Индекс энергетической степени свободы
        stateCoordinateInd = indf.GetIndex(self.__StateCoordinatesNames,
                                           stateCoordinateName)#Индекс координаты состояния
        
        #Задаем элемент
        return self.__PotentialsInterBet[energyPowerBetInd,stateCoordinateInd]
            
    #Получение массива температур энергетических степеней свободы
    def GetPotentialsInterBet(self):
        #Задаем элемент
        return self.__PotentialsInterBet.copy()
            
    #Задание постоянных потенцталов взаимодействия энергетических степеней свободы
    def SetBetaConstElement(self,
                            energyPowerName,#Имя энергетической степени свободы
                            processCoordinateName,#Имя координаты процесса
                            elementValue#Значение элемента
                            ):
        #Получаем индексы
        energyPowerInd = indf.GetIndex(self.__EnergyPowersNames,
                                       energyPowerName)#Индекс энергетической степени свободы
        processCoordinateInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                             processCoordinateName)#Индекс координаты состояния
        
        #Задаем элемент
        self.__Beta[processCoordinateInd,energyPowerInd] = elementValue
            
    #Получение температур энергетических степеней свободы
    def GetBetaElement(self,
                       energyPowerName,#Имя энергетической степени свободы
                       processCoordinateName,#Имя координаты процесса
                       ):
        #Получаем индексыиндексы
        energyPowerInd = indf.GetIndex(self.__EnergyPowersNames,
                                       energyPowerName)#Индекс энергетической степени свободы
        processCoordinateInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                             processCoordinateName)#Индекс координаты состояния
        
        #Задаем элемент
        return self.__Beta[processCoordinateInd,energyPowerInd]
            
    #Получение массива температур энергетических степеней свободы
    def GetBeta(self):
        #Задаем элемент
        return self.__Beta.copy()
            
    #Задание постоянных коэффициентов главного блока кинетической матрицы по координатам процессов
    def SetKineticMatrixPCPCConstElement(self,
                                         processCoordinateName,#Имя координаты процесса
                                         processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                         elementValue#Значение элемента
                                         ):
        #Получаем индексы
        processCoordinateInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                             processCoordinateName)#Индекс координаты процесса
        processCoordinateAffInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                                processCoordinateAffName)#Индекс термодинамической силы, сопряженной координате процесса
        
        #Задаем элемент
        self.__KineticMatrixPCPC[processCoordinateInd,processCoordinateAffInd] = elementValue
            
    #Получение постоянных коэффициентов главного блока кинетической матрицы по координатам процессов
    def GetKineticMatrixPCPCElement(self,
                                    processCoordinateName,#Имя координаты процесса
                                    processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                    ):
        #Получаем индексы
        processCoordinateInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                             processCoordinateName)#Индекс координаты процесса
        processCoordinateAffInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                                processCoordinateAffName)#Индекс термодинамической силы, сопряженной координате процесса
        
        #Выводим элемент
        return self.__KineticMatrixPCPC[processCoordinateInd,processCoordinateAffInd]
            
    #Получение блока постоянных коэффициентов главного блока кинетической матрицы по координатам процессов
    def GetKineticMatrixPCPC(self):
        #Выводим элемент
        return self.__KineticMatrixPCPC.copy()
            
    #Задание постоянных коэффициентов перекрестного блока кинетической матрицы по координатам процессов и перенесенным теплотам
    def SetKineticMatrixPCHeatConstElement(self,
                                           processCoordinateName,#Имя координаты процесса
                                           processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                           elementValue#Значение элемента
                                           ):
        #Получаем индексы
        processCoordinateInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                             processCoordinateName)#Индекс координаты процесса
        processCoordinateAffInd = indf.GetIndex(self.__HeatTransfersNames,
                                                processCoordinateAffName)#Индекс термодинамической силы, сопряженной координате процесса
        
        #Задаем элемент
        self.__KineticMatrixPCHeat[processCoordinateInd,processCoordinateAffInd] = elementValue
            
    #Получение коэффициентов перекрестного блока кинетической матрицы по координатам процессов и перенесенным теплотам
    def GetKineticMatrixPCHeatElement(self,
                                      processCoordinateName,#Имя координаты процесса
                                      processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                      ):
        #Получаем индексы
        processCoordinateInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                             processCoordinateName)#Индекс координаты процесса
        processCoordinateAffInd = indf.GetIndex(self.__HeatTransfersNames,
                                                processCoordinateAffName)#Индекс термодинамической силы, сопряженной координате процесса
        
        #Выводим элемент
        return self.__KineticMatrixPCHeat[processCoordinateInd,processCoordinateAffInd]
            
    #Получение блока коэффициентов перекрестного блока кинетической матрицы по координатам процессов и перенесенным теплотам
    def GetKineticMatrixPCHeat(self):
        #Выводим элемент
        return self.__KineticMatrixPCHeat.copy()
            
    #Задание постоянных коэффициентов перекрестного блока кинетической матрицы по перенесенным теплотам и координатам процессов
    def SetKineticMatrixHeatPCConstElement(self,
                                           processCoordinateName,#Имя координаты процесса
                                           processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                           elementValue#Значение элемента
                                           ):
        #Получаем индексы
        processCoordinateInd = indf.GetIndex(self.__HeatTransfersNames,
                                             processCoordinateName)#Индекс координаты процесса
        processCoordinateAffInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                                processCoordinateAffName)#Индекс термодинамической силы, сопряженной координате процесса
        
        #Задаем элемент
        self.__KineticMatrixHeatPC[processCoordinateInd,processCoordinateAffInd] = elementValue
            
    #Получение коэффициентов перекрестного блока кинетической матрицы по перенесенным теплотам и координатам процессов
    def GetKineticMatrixHeatPCElement(self,
                                      processCoordinateName,#Имя координаты процесса
                                      processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                      ):
        #Получаем индексы
        processCoordinateInd = indf.GetIndex(self.__HeatTransfersNames,
                                             processCoordinateName)#Индекс координаты процесса
        processCoordinateAffInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                                processCoordinateAffName)#Индекс термодинамической силы, сопряженной координате процесса
        
        #Выводим элемент
        return self.__KineticMatrixHeatPC[processCoordinateInd,processCoordinateAffInd]
            
    #Получение блока коэффициентов перекрестного блока кинетической матрицы по перенесенным теплотам и координатам процессов
    def GetKineticMatrixHeatPC(self):
        #Выводим элемент
        return self.__KineticMatrixHeatPC.copy()
            
    #Задание постоянных коэффициентов главного блока кинетической матрицы по перенесенным теплотам
    def SetKineticMatrixHeatHeatConstElement(self,
                                             processCoordinateName,#Имя координаты процесса
                                             processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                             elementValue#Значение элемента
                                             ):
        #Получаем индексы
        processCoordinateInd = indf.GetIndex(self.__HeatTransfersNames,
                                             processCoordinateName)#Индекс координаты процесса
        processCoordinateAffInd = indf.GetIndex(self.__HeatTransfersNames,
                                                processCoordinateAffName)#Индекс термодинамической силы, сопряженной координате процесса
        
        #Задаем элемент
        self.__KineticMatrixHeatHeat[processCoordinateInd,processCoordinateAffInd] = elementValue
            
    #Получение постоянных коэффициентов главного блока кинетической матрицы по перенесенным теплотам
    def GetKineticMatrixHeatHeatElement(self,
                                        processCoordinateName,#Имя координаты процесса
                                        processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                        ):
        #Получаем индексы
        processCoordinateInd = indf.GetIndex(self.__HeatTransfersNames,
                                             rocessCoordinateName)#Индекс координаты процесса
        processCoordinateAffInd = indf.GetIndex(self.__HeatTransfersNames,
                                                processCoordinateAffName)#Индекс термодинамической силы, сопряженной координате процесса
        
        #Выводим элемент
        return self.__KineticMatrixHeatHeat[processCoordinateInd,processCoordinateAffInd]
            
    #Получение блока постоянных коэффициентов главного блока кинетической матрицы по перенесенным теплотам
    def GetKineticMatrixHeatHeat(self):
        #Выводим элемент
        return self.__KineticMatrixHeatHeat.copy()
            
    #Задание постоянных коэффциентов матрицы обратных теплоемкостей
    def SetInvHeatCapacityMatrixConstElement(self,
                                             reducedTemperatureEnergyPowerName,#Имя приведенной температуры
                                             energyPowerName,#Имя энергетической степени свободы
                                             elementValue#Значение элемента
                                             ):
        #Получаем индексы
        reducedTemperatureEnergyPowerInd = indf.GetIndex(self.__ReducedTemperaturesEnergyPowersNames,
                                                         reducedTemperatureEnergyPowerName)#Индекс энергетической степени свободы
        energyPowerInd = indf.GetIndex(self.__EnergyPowersNames,
                                       energyPowerName)#Индекс координаты состояния
        
        #Задаем элемент
        self.__InvHeatCapacityMatrix[reducedTemperatureEnergyPowerInd,energyPowerInd] = elementValue
            
    #Получение коэффциентов матрицы обратных теплоемкостей
    def GetInvHeatCapacityMatrixElement(self,
                                        reducedTemperatureEnergyPowerName,#Имя приведенной температуры
                                        energyPowerName#Имя энергетической степени свободы
                                        ):
        #Получаем индексы
        reducedTemperatureEnergyPowerInd = indf.GetIndex(self.__ReducedTemperaturesEnergyPowersNames,
                                                         reducedTemperatureEnergyPowerName)#Индекс энергетической степени свободы
        energyPowerInd = indf.GetIndex(self.__EnergyPowersNames,
                                       energyPowerName)#Индекс координаты состояния
        
        #Задаем элемент
        return self.__InvHeatCapacityMatrix[reducedTemperatureEnergyPowerInd,energyPowerInd]
            
    #Получение массива коэффциентов матрицы обратных теплоемкостей
    def GetInvHeatCapacityMatrix(self):
        #Задаем элемент
        return self.__InvHeatCapacityMatrix.copy()
            
    #Задание постоянных коэффциентов матрицы тепловых эффектов
    def SetHeatEffectMatrixConstElement(self,
                                        reducedTemperatureEnergyPowerName,#Имя приведенной температуры
                                        stateCoordinateName,#Имя координаты состояния
                                        elementValue#Значение элемента
                                        ):
        #Получаем индексы
        reducedTemperatureEnergyPowerInd = indf.GetIndex(self.__ReducedTemperaturesEnergyPowersNames,
                                                         reducedTemperatureEnergyPowerName)#Индекс энергетической степени свободы
        stateCoordinateInd = indf.GetIndex(self.__StateCoordinatesNames,
                                           stateCoordinateName)#Индекс координаты состояния
        
        #Задаем элемент
        self.__HeatEffectMatrix[reducedTemperatureEnergyPowerInd,stateCoordinateInd] = elementValue
            
    #Получение коэффциентов матрицы тепловых эффектов
    def GetHeatEffectMatrixElement(self,
                                   reducedTemperatureEnergyPowerName,#Имя приведенной температуры
                                   stateCoordinateName#Имя координаты состояния
                                   ):
        #Получаем индексы
        reducedTemperatureEnergyPowerInd = indf.GetIndex(self.__ReducedTemperaturesEnergyPowersNames,
                                                         reducedTemperatureEnergyPowerName)#Индекс энергетической степени свободы
        stateCoordinateInd = indf.GetIndex(self.__StateCoordinatesNames,
                                           stateCoordinateName)#Индекс координаты состояния
        
        #Задаем элемент
        return self.__HeatEffectMatrix[reducedTemperatureEnergyPowerInd,stateCoordinateInd]
            
    #Получение массива коэффциентов матрицы тепловых эффектов
    def GetHeatEffectMatrix(self):
        #Задаем элемент
        return self.__HeatEffectMatrix.copy()
            
    #Расчет величин
    def CountSystem(self,
                    stateCoordinates,#Координаты состояния
                    reducedTemp,#Приведенные температуры энергетических степеней свободы
                    systemParameters#Параметры системы
                    ):
        #Рассчитываем свойства веществ и процессов
        (balanceMatrix,
         stateCoordinatesStreams,
         heatEnergyPowersStreams,
         energyPowerTemperatures,
         potentialInter,
         potentialInterBet,
         beta,kineticMatrixPCPC,
         kineticMatrixPCHeat,
         kineticMatrixHeatPC,
         kineticMatrixHeatHeat,
         invHeatCapacityMatrixCf,
         heatEffectMatrixCf) = self.__StateFunction(stateCoordinates,
                                                    reducedTemp,
                                                    systemParameters)
         
        #Задаем переменные коэффициенты матрицы баланса
        self.__BalanceStateCoordinates[self.__StateCoordinatesVarBalanceIndexes,self.__ProcessCoordinatesVarBalanceIndexes] = balanceMatrix
         
        #Задаем переменные внешние потоки
        self.__StateCoordinatesStreams[self.__StateCoordinatesVarStreamsIndexes] = stateCoordinatesStreams
         
        #Задаем переменные внешние потоки теплоты
        self.__HeatEnergyPowersStreams[self.__HeatEnergyPowersVarStreamsIndexes] = heatEnergyPowersStreams
         
        #Задаем переменные температуры
        self.__TEnergyPowers[self.__EnergyPowersVarTemperatureIndexes] = energyPowerTemperatures
         
        #Задаем переменные потенциалы взаимодействия
        self.__PotentialsInter[self.__EnergyPowersVarPotentialsInterIndexes,self.__StateCoordinatesVarPotentialsInterIndexes] = potentialInter
                 
        #Задаем переменные потенциалы взаимодействия между энергетическими степенями свободы
        self.__PotentialsInterBet[self.__EnergyPowersVarPotentialsInterBetIndexes,self.__StateCoordinatesVarPotentialsInterBetIndexes] = potentialInterBet
         
        #Задаем переменные потенциалы взаимодействия
        self.__Beta[self.__ProcessCoordinatesVarBetaIndexes,self.__EnergyPowersVarBetaIndexes] = beta
        
        #Задаем главный блок кинетической матрицы по координатам процессов
        self.__KineticMatrixPCPC[self.__ProcessCoordinatesVarKineticPCPCIndexes,
                                 self.__ProcessCoordinatesVarKineticPCPCAffIndexes] = kineticMatrixPCPC
        
        #Задаем перекрестный блок кинетической матрицы по координатам процессов и перенесенным теплотам
        self.__KineticMatrixPCHeat[self.__ProcessCoordinatesVarKineticPCHeatIndexes,
                                   self.__ProcessCoordinatesVarKineticPCHeatAffIndexes] = kineticMatrixPCHeat
        
        #Задаем перекрестный блок кинетической матрицы по перенесенным теплотам и координатам процессов
        self.__KineticMatrixHeatPC[self.__ProcessCoordinatesVarKineticHeatPCIndexes,
                                   self.__ProcessCoordinatesVarKineticHeatPCAffIndexes] = kineticMatrixHeatPC
        
        #Задаем главный блок кинетической матрицы по перенесенным теплотам
        self.__KineticMatrixHeatHeat[self.__ProcessCoordinatesVarKineticHeatHeatIndexes,
                                     self.__ProcessCoordinatesVarKineticHeatHeatAffIndexes] = kineticMatrixHeatHeat
        
        #Задаем матрицу обратных теплоемкостей
        self.__InvHeatCapacityMatrix[self.__ReducedTemperaturesEnergyPowersVarInvHeatCapacityIndexes,
                                     self.__EnergyPowersVarInvHeatCapacity] = invHeatCapacityMatrixCf
        
        #Задаем матрицу приведенных тепловых эффектов
        self.__HeatEffectMatrix[self.__ReducedTemperaturesEnergyPowersVarHeatEffectIndexes,
                                self.__StateCoordinatesVarHeatEffectIndexes] = heatEffectMatrixCf
        
        #Определяем термодинамические силы
        AllPotentialsInter = np.sum(self.__PotentialsInter, axis=0) + np.sum(self.__PotentialsInterBet, axis=0)#Суммарные потенциалы взаимодействия
        invTEnergyPowers = NonEqSystemQ.__Tbase / self.__TEnergyPowers#Приведенные температуры
        self.__vHeatProcess = np.dot(AllPotentialsInter,
                                     self.__BalanceStateCoordinates)#Некомпенсированная теплота
        self.__Affinity = self.__vHeatProcess * np.dot(invTEnergyPowers,self.__Beta.transpose())#Термодинамические силы
        self.__HeatAffinity = np.dot(invTEnergyPowers,self.__HeatTransferMatrix.toarray())#Термодинамические силы переноса теплоты
        
        #Определяем скорости протекания физико-химических процессов
        self.__vProcessCoordinates = np.dot(self.__Affinity,self.__KineticMatrixPCPC.transpose()) + \
                                     np.dot(self.__HeatAffinity,self.__KineticMatrixPCHeat.transpose())
        
        #Определяем скорости переноса теплот между энергетическими степенями свободы
        self.__vHeatTransfer = np.dot(self.__Affinity,self.__KineticMatrixHeatPC.transpose()) + \
                               np.dot(self.__HeatAffinity,self.__KineticMatrixHeatHeat.transpose())
        
        #Определяем мощность выделения некомпенсированной теплоты
        self.__vHeatProcess *= self.__vProcessCoordinates
        
        #Определяем скорости изменения координат состояния
        self.__vStateCoordinates = np.dot(self.__vProcessCoordinates,
                                          self.__BalanceStateCoordinates.transpose())#Учитываем внутренние процессы
        self.__vStateCoordinates[self.__StateCoordinatesStreamsIndexes] += self.__StateCoordinatesStreams#Учитываем внешние потоки
        
        #Определяем скорость сообщения теплоты энергетической степени свободы
        self.__vHeatPower = np.dot(self.__vHeatTransfer,
                                   self.__HeatTransferMatrix.toarray().transpose()) + \
                            np.dot(self.__vHeatProcess,self.__Beta)#Учитываем внутренние процессы переноса теплоты
        self.__vHeatPower[self.__HeatEnergyPowersStreamsIndexes] += self.__HeatEnergyPowersStreams#Учитываем внешние потоки теплоты
        
        #Определяем скорость приращения внутренней энергии
        self.__vUPower = self.__vHeatPower - np.dot(self.__vStateCoordinates,
                                                    self.__PotentialsInter.transpose())
        
        #Определяем скорости приращения приведенных температур
        self.__vReducedTemperaturesEnergyPowers = np.dot(self.__vUPower,
                                                         self.__InvHeatCapacityMatrix.transpose()) + \
                                                  np.dot(self.__vStateCoordinates,
                                                         self.__HeatEffectMatrix.transpose())

