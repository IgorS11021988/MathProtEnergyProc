import numpy as np

import MathProtEnergyProc.IndexFunctions as indf

from MathProtEnergyProc.NonEqProcess.NonEqSystemBase import NonEqSystemBase

#Класс физико-химических систем
class NonEqSystem:
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
        #Задаем структуру системы
        self.__baseSystemStructure = NonEqSystemBase(stateCoordinatesNames,#Имена координат состояния
                                                     processCoordinatesNames,#Имена координат процессов
                                                     
                                                     stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                                     
                                                     stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                                     processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                                     stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                                     
                                                     varKineticNames,#Имена сопряженностей между собой координат процессов
                                                     varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                                     
                                                     stateCoordinatesVarStreamsNames#Имена переменных внешних потоков
                                                     )
        
        #Задаем функцию состояния
        self.__StateFunction = stateFunction#Функция состояния
                 
    #Задание постоянных элементов в матрицы баланса
    def SetBalanceStateCoordinatesConstElement(self,
                                               stateCoordinateName,#Имя координаты состояния
                                               processCoordinateName,#Имя координаты процесса
                                               elementValue#Значение элемента
                                               ):
        #Задаем элемент
        self.__baseSystemStructure.SetBalanceStateCoordinatesConstElement(stateCoordinateName,#Имя координаты состояния
                                                                          processCoordinateName,#Имя координаты процесса
                                                                          elementValue#Значение элемента
                                                                          )
    
    #Получение элементов матрицы баланса
    def GetBalanceStateCoordinatesElement(self,
                                          stateCoordinateName,#Имя координаты состояния
                                          processCoordinateName#Имя координаты процесса
                                          ):
        #Получаем элемент
        return self.__baseSystemStructure.GetBalanceStateCoordinatesElement(stateCoordinateName,#Имя координаты состояния
                                                                            processCoordinateName#Имя координаты процесса
                                                                            )
        
    #Получение матрицы баланса
    def GetBalanceMatrix(self):
        #Выводим матрицу баланса
        return self.__baseSystemStructure.GetBalanceMatrix()
            
    #Задание постоянных внешних потоков
    def SetStateCoordinatesStreamsConstElement(self,
                                               stateCoordinateStreamName,#Имя координаты состояния
                                               elementValue#Значение элемента
                                               ):
        #Задаем элемент
        self.__baseSystemStructure.SetStateCoordinatesStreamsConstElement(stateCoordinateStreamName,#Имя координаты состояния
                                                                          elementValue#Значение элемента
                                                                          )
            
    #Получение постоянных внешних потоков
    def GetStateCoordinatesStreamsElement(self,
                                          stateCoordinateStreamName#Имя координаты состояния
                                          ):
        #Выводим элемент
        return self.__baseSystemStructure.GetStateCoordinatesStreamsElement(stateCoordinateStreamName#Имя координаты состояния
                                                                            )
            
    #Получение массива внешних потоков
    def GetStateCoordinatesStreams(self):
        #Выводим внешние потоки
        return self.__baseSystemStructure.GetStateCoordinatesStreams()
              
    #Получение массива имен координат состояния
    def GetStateCoordinatesNames(self):
        #Выводим имена координат состояния
        return self.__baseSystemStructure.GetStateCoordinatesNames()
           
    #Получение массива скоростей изменения координат состояния
    def GetVStateCoordinates(self):
        #Выводим скорости изменения координат состояния
        return self.__baseSystemStructure.GetVStateCoordinates()
            
    #Получение массива имен координат процессов
    def GetProcessCoordinatesNames(self):
        #Выводим имена координат процессов
        return self.__baseSystemStructure.GetProcessCoordinatesNames()
            
    #Получение массива скоростей изменения координат процессов
    def GetVProcessCoordinates(self):
        #Выводим скорости изменения координат процессов
        return self.__baseSystemStructure.GetVProcessCoordinates()
                         
    #Получение массива термодинамических сил
    def GetAffinity(self):
        #Выводим массив термодинамических сил
        return self.__baseSystemStructure.GetAffinity()
            
    #Задание постоянных потенцталов взаимодействия
    def SetPotentialsInterConstElement(self,
                                       energyPowerName,#Имя энергетической степени свободы
                                       stateCoordinateName,#Имя координаты состояния
                                       elementValue#Значение элемента
                                       ):
        #Задаем элемент
        self.__baseSystemStructure.SetPotentialsInterConstElement(energyPowerName,#Имя энергетической степени свободы
                                                                  stateCoordinateName,#Имя координаты состояния
                                                                  elementValue#Значение элемента
                                                                  )
            
    #Получение потенциалов взаимодействия
    def GetPotentialsInterElement(self,
                                  energyPowerName,#Имя энергетической степени свободы
                                  stateCoordinateName#Имя координаты состояния
                                  ):
        #Выводим элемент
        return self.__baseSystemStructure.GetPotentialsInterElement(energyPowerName,#Имя энергетической степени свободы
                                                                    stateCoordinateName#Имя координаты состояния
                                                                    )
            
    #Получение массива потенциалов взаимодействия энергетических степеней свободы
    def GetPotentialsInter(self):
        #Выводим элемент
        return self.__baseSystemStructure.GetPotentialsInter()
            
    #Задание постоянных коэффициентов кинетической матрицы
    def SetKineticMatrixConstElement(self,
                                     processCoordinateName,#Имя координаты процесса
                                     processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                     elementValue#Значение элемента
                                     ):
        #Задаем элемент
        self.__baseSystemStructure.SetKineticMatrixConstElement(processCoordinateName,#Имя координаты процесса
                                                                processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                                                elementValue#Значение элемента
                                                                )
            
    #Получение постоянных коэффициентов главного блока кинетической матрицы по координатам процессов
    def GetKineticMatrixElement(self,
                                processCoordinateName,#Имя координаты процесса
                                processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                ):
        #Выводим элемент
        return self.__baseSystemStructure.GetKineticMatrixElement(processCoordinateName,#Имя координаты процесса
                                                                  processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                                                  )
            
    #Получение блока постоянных коэффициентов главного блока кинетической матрицы по координатам процессов
    def GetKineticMatrix(self):
        #Выводим элемент
        return self.__baseSystemStructure.GetKineticMatrix()
            
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
         
        #Рассчитываем прочие характеристики системы
        self.__baseSystemStructure.CountSystem(balanceMatrix,#Матрица баланса
                                               stateCoordinatesStreams,#Внешние потоки по координата состояния
                                               potentialInter,#Потенциалы взаимодействия
                                               kineticMatrix#Кинетическая матрица
                                               )
        