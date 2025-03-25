import numpy as np

from MathProtEnergyProc.IndexFunctions import GetIndexes

#Выходной функтор потенциалов взамодействия
class IntPotentialsBase:
    #Инициализация класса
    def __init__(self,
                 
                 stateCoordinatesNames,#Имена координат состояния
                 energyPowersNames,#Имена энергетических степеней свободы
                 
                 stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                 energyPowersVarPotentialsInterNames#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                 ):
        #Формируем индексы переменных потенциалов взаимодействия
        self.__EnergyPowersVarPotentialsInterIndexes = GetIndexes(energyPowersNames,
                                                                  energyPowersVarPotentialsInterNames)#Индексы переменных потенциалов взаимодействия по энергетическим степенәм свободы
        self.__StateCoordinatesVarPotentialsInterIndexes = GetIndexes(stateCoordinatesNames,
                                                                      stateCoordinatesVarPotentialsInterNames)#Индексы переменных потенциалов взаимодействия по координатам состояния

    #Вывод индексов по координатам состояния
    def _GetStateCoordinatesVarPotentialsInterIndexes(self):
        #Выводим индексы
        return self.__StateCoordinatesVarPotentialsInterIndexes.copy()

    #Вывод индексов по энергетическим степеням свободы
    def _GetEnergyPowersVarPotentialsInterIndexes(self):
        #Выводим индексы
        return self.__EnergyPowersVarPotentialsInterIndexes.copy()

    #Тело функции
    def __call__(self,
    
                 intPotentialsMatrix#Матрица потениалов взаимодействия
                 ):
        #Выводим потенциалы взаимодействия в соответствие с индексами
        return intPotentialsMatrix[self.__EnergyPowersVarPotentialsInterIndexes,
                                   self.__StateCoordinatesVarPotentialsInterIndexes].copy()
class IntPotentialsBaseFun(IntPotentialsBase):
    #Инициализация класса
    def __init__(self,
                 
                 stateCoordinatesNames,#Имена координат состояния
                 energyPowersNames,#Имена энергетических степеней свободы
                 
                 stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                 energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                 
                 intPotentialFunTemp,#Функции потенцалов взаимодействия в зависимости от температуры
                 intPotentialFunInvTemp,#Функции потенцалов взаимодействия в зависимости от обратной температуры
                 
                 isInvTemperatute=False#По температурам или по обратным температурам считаем
                 ):
        #Вызываем конструктор суперкласса
        super().__init__(stateCoordinatesNames,#Имена координат состояния
                         energyPowersNames,#Имена энергетических степеней свободы
                         
                         stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                         energyPowersVarPotentialsInterNames#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                         )
        
        #Задаем функцию расчета потенциалов взаимодействия
        if isInvTemperatute:
            self.__IntPotentialFun = intPotentialFunInvTemp
        else:
            self.__IntPotentialFun = intPotentialFunTemp

    #Функции потенциалов взаимодействия
    def GetIntPotentialFun(self):
        #Выводим потенциалы взаимодействия в соответствие с индексами
        return self.__IntPotentialFun
