import numpy as np

from MathProtEnergyProc.IndexFunctions import GetIndexes, GetIndex
from MathProtEnergyProc.HeatPowerValues.Base.Temperatures import InvHeatCapacityHeatEffectsFun, HeatEffectsFun, HeatCapacitiesFun

#Выходной функтор тепловых величин (теплоемкотей и тепловых эффектов)
class HeatValues:
    #Инициализация класса
    def __init__(self,
                 
                 stateCoordinatesNames,#Имена координат состояния
                 energyPowersNames,#Имена энергетических степеней свободы
                 
                 temperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                 energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                 temperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                 stateCoordinatesVarHeatEffectNames#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                 ):
        #Сохраняем имена
        self.__StateCoordinatesNames = stateCoordinatesNames.copy()#Имена координат состояния
        self.__EnergyPowersNames = energyPowersNames.copy()#Имена энергетических степеней свободы
        
        #Формируем индексы переменных коэффициентов матрицы обратных теплоемкостей
        self.__TemperaturesEnergyPowersVarInvHeatCapacityIndexes = GetIndexes(energyPowersNames,
                                                                              temperaturesEnergyPowersVarInvHeatCapacityNames)
        self.__EnergyPowersVarInvHeatCapacityIndexes = GetIndexes(energyPowersNames,
                                                                  energyPowersVarInvHeatCapacityNames)
        
        #Формируем индексы переменных коэффициентов матрицы приведенных тепловых эффектов
        self.__TemperaturesEnergyPowersVarHeatEffectIndexes = GetIndexes(energyPowersNames,
                                                                         temperaturesEnergyPowersVarHeatEffectNames)#Индексы переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        self.__StateCoordinatesVarHeatEffectIndexes = GetIndexes(stateCoordinatesNames,
                                                                 stateCoordinatesVarHeatEffectNames)#Индексы переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
    
        #Задаем пустую матрицу теплоемкостей
        self.__HeatCapacities = None
        
        #Задаем пустую матрицу тепловых эффектов
        self.__HeatEffects = None
            
    #Получение коэффциентов матрицы обратных теплоемкостей
    def GetHeatCapacityMatrixElement(self,
    
                                     energyPowerName,#Имя энергетической степени свободы
                                     temperatureEnergyPowerName#Имя температуры энергетической степени свободы
                                     ):
        #Получаем индексы
        temperatureEnergyPowerInd = GetIndex(self.__EnergyPowersNames,
                                             temperatureEnergyPowerName)#Индекс энергетической степени свободы
        energyPowerInd = GetIndex(self.__EnergyPowersNames,
                                  energyPowerName)#Индекс координаты состояния
        
        #Выводим элемент
        return self.__HeatCapacities[temperatureEnergyPowerInd,
                                     energyPowerInd]
            
    #Получение массива коэффциентов матрицы обратных теплоемкостей
    def GetHeatCapacityMatrix(self):
        #Выводим элемент
        return self.__HeatCapacities.copy()
            
    #Получение коэффциентов матрицы тепловых эффектов
    def GetHeatEffectMatrixElement(self,
    
                                   temperaturesEnergyPowerName,#Имя температуры энергетической степени свободы
                                   stateCoordinateName#Имя координаты состояния
                                   ):
        #Получаем индексы
        temperaturesEnergyPowerInd = GetIndex(self.__EnergyPowersNames,
                                              temperaturesEnergyPowerName)#Индекс энергетической степени свободы
        stateCoordinateInd = GetIndex(self.__StateCoordinatesNames,
                                      stateCoordinateName)#Индекс координаты состояния
        
        #Выводим элемент
        return self.__HeatEffects[temperaturesEnergyPowerInd,
                                  stateCoordinateInd]
            
    #Получение массива коэффциентов матрицы тепловых эффектов
    def GetHeatEffectMatrix(self):
        #Выводим элемент
        return self.__HeatEffects.copy()
            
    #Тело функции
    def __call__(self,
    
                 jacSTemperatures,#Якобиан приведенной энтропии по температурам
                 hesSTemperatures,#Матрица Гесса приведенной энтропии по температурам
                 hesSTemperaturesStateCoordinates,#Матрица Гесса приведенной энтропии по температурам и координатам состояния
                 temperatures#Температуры
                 ):
        #Определяем теплоемкости
        self.__HeatCapacities = HeatCapacitiesFun(jacSTemperatures,
                                                  hesSTemperatures,
                                                  temperatures)
        
        #Определяем тепловые эффекты
        self.__HeatEffects = HeatEffectsFun(hesSTemperaturesStateCoordinates,
                                            temperatures)
        
        #Определяем приведенные теплоемкости и приведенные тепловые эффекты
        (invHeatCapacities,
         redHeatEffects) = InvHeatCapacityHeatEffectsFun(self.__HeatCapacities,
                                                         self.__HeatEffects)
        
        #Выводим потенциалы взаимодействия в соответствие с индексами
        return (invHeatCapacities[self.__TemperaturesEnergyPowersVarInvHeatCapacityIndexes,
                                  self.__EnergyPowersVarInvHeatCapacityIndexes].copy(),
                redHeatEffects[self.__TemperaturesEnergyPowersVarHeatEffectIndexes,
                               self.__StateCoordinatesVarHeatEffectIndexes].copy())
