import numpy as np

import MathProtEnergyProc.IndexFunctions as indf

#Класс физико-химических систем
class NonEqSystem(object):
    #Конструктор класса
    def __init__(self,
                 stateCoordinatesNames,#Имена координат состояния
                 processCoordinatesNames,#Имена координат процессов
                 
                 stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                 
                 stateFunction,#Функция состояния
                 
                 stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                 processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                 stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                 
                 varKineticNames,#Имена сопряженностей между собой координат процессов
                 varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                 
                 stateCoordinatesVarStreamsNames#Имена переменных внешних потоков
                 ):
        #Заполняем имена
        self.__StateCoordinatesNames = stateCoordinatesNames.copy()#Имена координат состояния
        self.__ProcessCoordinatesNames = processCoordinatesNames.copy()#Имена координат процессов
        
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
        
        #Формируем индексы переменных внешних потоков
        self.__StateCoordinatesVarStreamsIndexes = indf.Get2Indexes(self.__StateCoordinatesNames,
                                                                    self.__StateCoordinatesStreamsIndexes,
                                                                    stateCoordinatesVarStreamsNames)#Индексы внешних потоков по координатам состояния
        
        #Формируем матрицу потениалов взаимодействия
        self.__PotentialsInter = np.zeros((nStateCoordinates,), dtype=np.double)
        
        #Формируем индексы переменных потенциалов взаимодействия
        self.__StateCoordinatesVarPotentialsInterIndexes = indf.GetIndexes(self.__StateCoordinatesNames,
                                                                           stateCoordinatesVarPotentialsInterNames)#Индексы переменных потенциалов взаимодействия по координатам состояния
                                                                                   
        #Формируем массивы термодинамических сил
        self.__Affinity = np.zeros((nProcessCoordinates,), dtype=np.double)#Термодинамические силы
        
        #Формируем кинетическую матрицу
        self.__KineticMatrix = np.zeros((nProcessCoordinates,nProcessCoordinates), dtype=np.double)
        
        #Формируем индексы кинетической матрицы
        self.__ProcessCoordinatesVarKineticIndexes = indf.GetIndexes(self.__ProcessCoordinatesNames,
                                                                     varKineticNames)#Индексы сопряженностей между собой координат процессов
        self.__ProcessCoordinatesVarKineticAffIndexes = indf.GetIndexes(self.__ProcessCoordinatesNames,
                                                                        varKineticAffNames)#Индексы сопряженностей между собой термодинамических сил
                 
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
            
    #Задание постоянных потенцталов взаимодействия
    def SetPotentialsInterConstElement(self,
                                       energyPowerName,#Имя энергетической степени свободы
                                       stateCoordinateName,#Имя координаты состояния
                                       elementValue#Значение элемента
                                       ):
        #Получаем индексы
        stateCoordinateInd = indf.GetIndex(self.__StateCoordinatesNames,
                                           stateCoordinateName)#Индекс координаты состояния
        
        #Задаем элемент
        self.__PotentialsInter[stateCoordinateInd] = elementValue
            
    #Получение потенциалов взаимодействия
    def GetPotentialsInterElement(self,
                                  energyPowerName,#Имя энергетической степени свободы
                                  stateCoordinateName#Имя координаты состояния
                                  ):
        #Получаем индексыиндексы
        stateCoordinateInd = indf.GetIndex(self.__StateCoordinatesNames,
                                           stateCoordinateName)#Индекс координаты состояния
        
        #Задаем элемент
        return self.__PotentialsInter[stateCoordinateInd]
            
    #Получение массива потенциалов взаимодействия энергетических степеней свободы
    def GetPotentialsInter(self):
        #Задаем элемент
        return self.__PotentialsInter.copy()
            
    #Задание постоянных коэффициентов кинетической матрицы
    def SetKineticMatrixConstElement(self,
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
        self.__KineticMatrix[processCoordinateInd,processCoordinateAffInd] = elementValue
            
    #Получение постоянных коэффициентов главного блока кинетической матрицы по координатам процессов
    def GetKineticMatrixElement(self,
                                processCoordinateName,#Имя координаты процесса
                                processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                ):
        #Получаем индексы
        processCoordinateInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                             processCoordinateName)#Индекс координаты процесса
        processCoordinateAffInd = indf.GetIndex(self.__ProcessCoordinatesNames,
                                                processCoordinateAffName)#Индекс термодинамической силы, сопряженной координате процесса
        
        #Выводим элемент
        return self.__KineticMatrix[processCoordinateInd,processCoordinateAffInd]
            
    #Получение блока постоянных коэффициентов главного блока кинетической матрицы по координатам процессов
    def GetKineticMatrix(self):
        #Выводим элемент
        return self.__KineticMatrix.copy()
            
    #Расчет величин
    def CountSystem(self,
                    stateCoordinates,#Координаты состояния
                    systemParameters#Параметры системы
                    ):
        #Рассчитываем свойства веществ и процессов
        (balanceMatrix,
         stateCoordinatesStreams,
         potentialInter,
         kineticMatrix) = self.__StateFunction(stateCoordinates,
                                               systemParameters)
         
        #Задаем переменные коэффициенты матрицы баланса
        self.__BalanceStateCoordinates[self.__StateCoordinatesVarBalanceIndexes,self.__ProcessCoordinatesVarBalanceIndexes] = balanceMatrix
         
        #Задаем переменные внешние потоки
        self.__StateCoordinatesStreams[self.__StateCoordinatesVarStreamsIndexes] = stateCoordinatesStreams
         
        #Задаем переменные потенциалы взаимодействия
        self.__PotentialsInter[self.__StateCoordinatesVarPotentialsInterIndexes] = potentialInter
        
        #Задаем главный блок кинетической матрицы по координатам процессов
        self.__KineticMatrix[self.__ProcessCoordinatesVarKineticIndexes,
                             self.__ProcessCoordinatesVarKineticAffIndexes] = kineticMatrix
        
        #Определяем термодинамические силы
        self.__Affinity = np.dot(self.__PotentialsInter,self.__BalanceStateCoordinates)#Некомпенсированная теплота
        
        #Определяем скорости протекания физико-химических процессов
        self.__vProcessCoordinates = np.dot(self.__Affinity,self.__KineticMatrix.transpose())
        
        #Определяем скорости изменения координат состояния
        self.__vStateCoordinates = np.dot(self.__vProcessCoordinates,
                                          self.__BalanceStateCoordinates.transpose())#Учитываем внутренние процессы
        self.__vStateCoordinates[self.__StateCoordinatesStreamsIndexes] += self.__StateCoordinatesStreams#Учитываем внешние потоки
        