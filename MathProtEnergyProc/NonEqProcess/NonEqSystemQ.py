import numpy as np
import scipy.sparse as sp

import MathProtEnergyProc.IndexFunctions as indf

from MathProtEnergyProc.NonEqProcess.NonEqSystemQBase import NonEqSystemQBase

#Класс физико-химических систем
class NonEqSystemQ(object):
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
        #Задаем структуру системы
        self.__baseSystemStructure = NonEqSystemQBase(stateCoordinatesNames,#Имена координат состояния
                                                      processCoordinatesNames,#Имена координат процессов
                                                      energyPowersNames,#Имена энергетических степеней свободы
                                                      reducedTemperaturesEnergyPowersNames,#Имена приведенных температур энергетических степеней свободы
                                                      energyPowersBetNames,#Имена взаимодействий между энергетическими степенями свободы
                                                      heatTransfersNames,#Имена потоков переноса теплоты
                                                      heatTransfersOutputEnergyPowersNames,#Имена энергетических степеней свободы, с которых уходит теплота
                                                      heatTransfersInputEnergyPowersNames,#Имена энергетических степеней свободы, на которые приходит теплота
                                                      
                                                      stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                                      heatEnergyPowersStreamsNames,#Имена потоков теплоты на энергетические степени свободы
                                                      
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
        #Выводим элемент
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
                    
    #Получение массива внешних потоков
    def GetHeatTransferMatrix(self):
        #Выводим внешние потоки
        return self.__baseSystemStructure.GetHeatTransferMatrix()
            
    #Задание постоянных внешних потоков
    def SetHeatEnergyPowersStreamsConstElement(self,
                                               heatEnergyPowersStreamName,#Имя энергетической степени свободы
                                               elementValue#Значение элемента
                                               ):
        #Задаем элемент
        self.__baseSystemStructure.SetHeatEnergyPowersStreamsConstElement(heatEnergyPowersStreamName,#Имя энергетической степени свободы
                                                                          elementValue#Значение элемента
                                                                          )
            
    #Получение постоянных внешних потоков
    def GetHeatEnergyPowersStreamsElement(self,
                                          heatEnergyPowersStreamName,#Имя энергетической степени свободы
                                          ):
        #Выводим элемент
        return self.__baseSystemStructure.GetHeatEnergyPowersStreamsElement(heatEnergyPowersStreamName,#Имя энергетической степени свободы
                                                                            )
            
    #Получение массива внешних потоков
    def GetHeatEnergyPowersStreams(self):
        #Выводим внешние потоки
        return self.__baseSystemStructure.GetHeatEnergyPowersStreams()
              
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
             
    #Получение массива имен энергетических степеней свободы
    def GetEnergyPowerNames(self):
        #Выводим имена энергетических степеней свободы
        return self.__baseSystemStructure.GetEnergyPowerNames()
             
    #Получение массива имен энергетических степеней свободы
    def GetReducedTemperaturesEnergyPowersNames(self):
        #Выводим имена энергетических степеней свободы
        return self.__baseSystemStructure.GetReducedTemperaturesEnergyPowersNames()
    
    #Получение массива имен взаимодействий между энергетическими степенями свободы
    def GetEnergyPowerBetNames(self):
        #Выводим имена взаимодействий между энергетическими степенями свободы
        return self.__baseSystemStructure.GetEnergyPowerBetNames()
             
    #Получение массива скоростей сообщения теплоты энергетическим степеням свободы
    def GetVHeatPower(self):
        #Выводим скорости сообщения теплоты энергетическим степеням свободы
        return self.__baseSystemStructure.GetVHeatPower()
                
    #Получение массива скоростей выделения некомпенсированных теплот
    def GetVHeatProcess(self):
        #Выводим скорости выделения некомпенсированных теплот
        return self.__baseSystemStructure.GetVHeatProcess()
             
    #Получение массива имен процессов переноса теплоты
    def GetHeatTransfersNames(self):
        #Выводим имена процессов переноса теплоты
        return self.__baseSystemStructure.GetHeatTransfersNames()
             
    #Получение массива скоростей процессов переноса теплоты
    def GetVHeatTransfers(self):
        #Выводим скорости процессов переноса теплоты
        return self.__baseSystemStructure.GetVHeatTransfers()
                         
    #Получение массива термодинамических сил переноса теплоты
    def GetHeatAffinity(self):
        #Выводим массив термодинамических сил переноса теплоты
        return self.__baseSystemStructure.GetHeatAffinity()
               
    #Задание постоянных температур энергетических степеней свободы
    def SetTEnergyPowersConstElement(self,
                                     energyPowerName,#Имя энергетической степени свободы
                                     elementValue#Значение элемента
                                     ):
        #Задаем элемент
        self.__baseSystemStructure.SetTEnergyPowersConstElement(energyPowerName,#Имя энергетической степени свободы
                                                                elementValue#Значение элемента
                                                                )
            
    #Получение температур энергетических степеней свободы
    def GetTEnergyPowersElement(self,
                                energyPowerName#Имя энергетической степени свободы
                                ):
        #Выводим элемент
        return self.__baseSystemStructure.GetTEnergyPowersElement(energyPowerName#Имя энергетической степени свободы
                                                                  )
            
    #Получение массива температур энергетических степеней свободы
    def GetTEnergyPowers(self):
        #Задаем элемент
        return self.__baseSystemStructure.GetTEnergyPowers()
             
    #Получение массива скоростей изменения внутренних энергий энергетических степеней свободы
    def GetVUPower(self):
        #Выводим скорости изменения внутренних энергий энергетических степеней свободы
        return self.__baseSystemStructure.GetVUPower()
             
    #Получение массива скоростей изменения приведенных температур энергетических степеней свободы
    def GetVReducedTemperaturesEnergyPowers(self):
        #Выводим скорости изменения приведенных температур энергетических степеней свободы
        return self.__baseSystemStructure.GetVReducedTemperaturesEnergyPowers()
            
    #Задание постоянных потенцталов взаимодействия энергетических степеней свободы
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
            
    #Получение потенциалов взаимодействия энергетических степеней свободы
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
        #Задаем элемент
        return self.__baseSystemStructure.GetPotentialsInter()
            
    #Задание постоянных потенцталов взаимодействия между энергетическими степенями свободы
    def SetPotentialsInterBetConstElement(self,
                                          energyPowerBetName,#Имя энергетической степени свободы
                                          stateCoordinateName,#Имя координаты состояния
                                          elementValue#Значение элемента
                                          ):
        #Задаем элемент
        self.__baseSystemStructure.SetPotentialsInterBetConstElement(energyPowerBetName,#Имя энергетической степени свободы
                                                                     stateCoordinateName,#Имя координаты состояния
                                                                     elementValue#Значение элемента
                                                                     )
            
    #Получение температур энергетических степеней свободы
    def GetPotentialsInterBetElement(self,
                                     energyPowerBetName,#Имя энергетической степени свободы
                                     stateCoordinateName#Имя координаты состояния
                                     ):
        #Выводим элемент
        return self.__baseSystemStructure.GetPotentialsInterBetElement(energyPowerBetName,#Имя энергетической степени свободы
                                                                       stateCoordinateName#Имя координаты состояния
                                                                       )
            
    #Получение массива температур энергетических степеней свободы
    def GetPotentialsInterBet(self):
        #Задаем элемент
        return self.__baseSystemStructure.GetPotentialsInterBet()
            
    #Задание постоянных потенцталов взаимодействия энергетических степеней свободы
    def SetBetaConstElement(self,
                            energyPowerName,#Имя энергетической степени свободы
                            processCoordinateName,#Имя координаты процесса
                            elementValue#Значение элемента
                            ):
        #Задаем элемент
        self.__baseSystemStructure.SetBetaConstElement(energyPowerName,#Имя энергетической степени свободы
                                                       processCoordinateName,#Имя координаты процесса
                                                       elementValue#Значение элемента
                                                       )
            
    #Получение температур энергетических степеней свободы
    def GetBetaElement(self,
                       energyPowerName,#Имя энергетической степени свободы
                       processCoordinateName,#Имя координаты процесса
                       ):
        #Выводим элемент
        return self.__baseSystemStructure.GetBetaElement(energyPowerName,#Имя энергетической степени свободы
                                                         processCoordinateName,#Имя координаты процесса
                                                         )
            
    #Получение массива температур энергетических степеней свободы
    def GetBeta(self):
        #Задаем элемент
        return self.__baseSystemStructure.GetBeta()
            
    #Задание постоянных коэффициентов главного блока кинетической матрицы по координатам процессов
    def SetKineticMatrixPCPCConstElement(self,
                                         processCoordinateName,#Имя координаты процесса
                                         processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                         elementValue#Значение элемента
                                         ):
        #Задаем элемент
        self.__baseSystemStructure.SetKineticMatrixPCPCConstElement(processCoordinateName,#Имя координаты процесса
                                                                    processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                                                    elementValue#Значение элемента
                                                                    )
            
    #Получение постоянных коэффициентов главного блока кинетической матрицы по координатам процессов
    def GetKineticMatrixPCPCElement(self,
                                    processCoordinateName,#Имя координаты процесса
                                    processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                    ):
        #Выводим элемент
        return self.__baseSystemStructure.GetKineticMatrixPCPCElement(processCoordinateName,#Имя координаты процесса
                                                                      processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                                                      )
            
    #Получение блока постоянных коэффициентов главного блока кинетической матрицы по координатам процессов
    def GetKineticMatrixPCPC(self):
        #Выводим матрицу
        return self.__baseSystemStructure.GetKineticMatrixPCPC()
            
    #Задание постоянных коэффициентов перекрестного блока кинетической матрицы по координатам процессов и перенесенным теплотам
    def SetKineticMatrixPCHeatConstElement(self,
                                           processCoordinateName,#Имя координаты процесса
                                           processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                           elementValue#Значение элемента
                                           ):
        #Задаем элемент
        self.__baseSystemStructure.SetKineticMatrixPCHeatConstElement(processCoordinateName,#Имя координаты процесса
                                                                      processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                                                      elementValue#Значение элемента
                                                                      )
            
    #Получение коэффициентов перекрестного блока кинетической матрицы по координатам процессов и перенесенным теплотам
    def GetKineticMatrixPCHeatElement(self,
                                      processCoordinateName,#Имя координаты процесса
                                      processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                      ):
        #Выводим элемент
        return self.__baseSystemStructure.GetKineticMatrixPCHeatElement(processCoordinateName,#Имя координаты процесса
                                                                        processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                                                        )
            
    #Получение блока коэффициентов перекрестного блока кинетической матрицы по координатам процессов и перенесенным теплотам
    def GetKineticMatrixPCHeat(self):
        #Выводим матрицу
        return self.__baseSystemStructure.GetKineticMatrixPCHeat()
            
    #Задание постоянных коэффициентов перекрестного блока кинетической матрицы по перенесенным теплотам и координатам процессов
    def SetKineticMatrixHeatPCConstElement(self,
                                           processCoordinateName,#Имя координаты процесса
                                           processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                           elementValue#Значение элемента
                                           ):
        #Задаем элемент
        self.__baseSystemStructure.SetKineticMatrixHeatPCConstElement(processCoordinateName,#Имя координаты процесса
                                                                      processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                                                      elementValue#Значение элемента
                                                                      )
            
    #Получение коэффициентов перекрестного блока кинетической матрицы по перенесенным теплотам и координатам процессов
    def GetKineticMatrixHeatPCElement(self,
                                      processCoordinateName,#Имя координаты процесса
                                      processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                      ):
        #Выводим элемент
        return self.__baseSystemStructure.GetKineticMatrixHeatPCElement(processCoordinateName,#Имя координаты процесса
                                                                        processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                                                        )
            
    #Получение блока коэффициентов перекрестного блока кинетической матрицы по перенесенным теплотам и координатам процессов
    def GetKineticMatrixHeatPC(self):
        #Выводим элемент
        return self.__baseSystemStructure.GetKineticMatrixHeatPC()
            
    #Задание постоянных коэффициентов главного блока кинетической матрицы по перенесенным теплотам
    def SetKineticMatrixHeatHeatConstElement(self,
                                             processCoordinateName,#Имя координаты процесса
                                             processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                             elementValue#Значение элемента
                                             ):
        #Задаем элемент
        self.__baseSystemStructure.SetKineticMatrixHeatHeatConstElement(processCoordinateName,#Имя координаты процесса
                                                                        processCoordinateAffName,#Имя термодинамической силы, сопряженной координате процесса
                                                                        elementValue#Значение элемента
                                                                        )
            
    #Получение постоянных коэффициентов главного блока кинетической матрицы по перенесенным теплотам
    def GetKineticMatrixHeatHeatElement(self,
                                        processCoordinateName,#Имя координаты процесса
                                        processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                        ):
        #Выводим элемент
        return self.__baseSystemStructure.GetKineticMatrixHeatHeatElement(processCoordinateName,#Имя координаты процесса
                                                                          processCoordinateAffName#Имя термодинамической силы, сопряженной координате процесса
                                                                          )
            
    #Получение блока постоянных коэффициентов главного блока кинетической матрицы по перенесенным теплотам
    def GetKineticMatrixHeatHeat(self):
        #Выводим элемент
        return self.__baseSystemStructure.GetKineticMatrixHeatHeat()
            
    #Задание постоянных коэффциентов матрицы обратных теплоемкостей
    def SetInvHeatCapacityMatrixConstElement(self,
                                             reducedTemperatureEnergyPowerName,#Имя приведенной температуры
                                             energyPowerName,#Имя энергетической степени свободы
                                             elementValue#Значение элемента
                                             ):
        #Задаем элемент
        self.__baseSystemStructure.SetInvHeatCapacityMatrixConstElement(reducedTemperatureEnergyPowerName,#Имя приведенной температуры
                                                                        energyPowerName,#Имя энергетической степени свободы
                                                                        elementValue#Значение элемента
                                                                        )
            
    #Получение коэффциентов матрицы обратных теплоемкостей
    def GetInvHeatCapacityMatrixElement(self,
                                        reducedTemperatureEnergyPowerName,#Имя приведенной температуры
                                        energyPowerName#Имя энергетической степени свободы
                                        ):
        #Выводим элемент
        return self.__baseSystemStructure.GetInvHeatCapacityMatrixElement(reducedTemperatureEnergyPowerName,#Имя приведенной температуры
                                                                          energyPowerName#Имя энергетической степени свободы
                                                                          )
            
    #Получение массива коэффциентов матрицы обратных теплоемкостей
    def GetInvHeatCapacityMatrix(self):
        #Задаем элемент
        return self.__baseSystemStructure.GetInvHeatCapacityMatrix()
            
    #Задание постоянных коэффциентов матрицы тепловых эффектов
    def SetHeatEffectMatrixConstElement(self,
                                        reducedTemperatureEnergyPowerName,#Имя приведенной температуры
                                        stateCoordinateName,#Имя координаты состояния
                                        elementValue#Значение элемента
                                        ):
        #Задаем элемент
        self.__baseSystemStructure.SetHeatEffectMatrixConstElement(reducedTemperatureEnergyPowerName,#Имя приведенной температуры
                                                                   stateCoordinateName,#Имя координаты состояния
                                                                   elementValue#Значение элемента
                                                                   )
            
    #Получение коэффциентов матрицы тепловых эффектов
    def GetHeatEffectMatrixElement(self,
                                   reducedTemperatureEnergyPowerName,#Имя приведенной температуры
                                   stateCoordinateName#Имя координаты состояния
                                   ):
        #Выводим элемент
        return self.__baseSystemStructure.GetHeatEffectMatrixElement(reducedTemperatureEnergyPowerName,#Имя приведенной температуры
                                                                     stateCoordinateName#Имя координаты состояния
                                                                     )
            
    #Получение массива коэффциентов матрицы тепловых эффектов
    def GetHeatEffectMatrix(self):
        #Задаем элемент
        return self.__baseSystemStructure.GetHeatEffectMatrix()
            
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
         
        #Рассчитываем прочие характеристики системы
        self.__baseSystemStructure.CountSystem(balanceMatrix,#Матрица баланса (для параметров состояния кроме внутренних энергий энергетических степеней свободы)
                                               stateCoordinatesStreams,#Внешние потоки по координатам состояния (кроме внутренних энергий)
                                               heatEnergyPowersStreams,#Внешние потоки теплоты по энергетическим степеням свободы
                                               energyPowerTemperatures,#Температуры энергетических степеней свободы
                                               potentialInter,#Потенциалы взаимодействия энергтических степеней свободы
                                               potentialInterBet,#Потенциалы взаимодействия, обусловленные взаимодействием между энергтическими степенями свободы
                                               beta,#Доли распределения некомпенсированных теплот ежду энергетическими степенями свободы
                                               kineticMatrixPCPC,#Главный блок кинетической матрицы по координатам процессов (кроме перенесенных теплот)
                                               kineticMatrixPCHeat,#Перекрестный блок кинетической матрицы между координатами процессов (кроме перенесенных теплот) и перенесенными теплотами
                                               kineticMatrixHeatPC,#Перекрестный блок кинетической матрицы между перенесенными теплотами и координатами процессов (кроме перенесенных теплот)
                                               kineticMatrixHeatHeat,#Главный блок кинетической матрицы по перенесенным теплотам
                                               invHeatCapacityMatrixCf,#Приведенные теплоемкости энергетических степеней свободы по координатам состояния (кроме внутренних энергий энергетических степеней свободы)
                                               heatEffectMatrixCf#Приведенные тепловые эффекты энергетических степеней свободы по координатам состояния (кроме внутренних энергий энергетических степеней свободы)
                                               )

