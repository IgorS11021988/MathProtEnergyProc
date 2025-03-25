import numpy as np

from MathProtEnergyProc.HeatPowerValues.Base.IntPotentialsT import IntPotentialFunTOne
from MathProtEnergyProc.HeatPowerValues.Base.IntPotentialsSdU import IntPotentialFunSdUOne

from MathProtEnergyProc.HeatPowerValues.Output.IntPotentialsBase import IntPotentialsBaseFun

#Выходной функтор потенциалов взамодействия
class IntPotentialsOne(IntPotentialsBaseFun):
    #Инициализация класса
    def __init__(self,
                 
                 stateCoordinatesNames,#Имена координат состояния
                 energyPowersNames,#Имена энергетических степеней свободы
                 
                 stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                 energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                 
                 isInvTemperatute=False#По температурам или по обратным температурам считаем
                 ):
        #Вызываем конструктор суперкласса
        super().__init__(stateCoordinatesNames,#Имена координат состояния
                         energyPowersNames,#Имена энергетических степеней свободы
                         
                         stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                         energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                         
                         IntPotentialFunTOne,#Функции потенцалов взаимодействия в зависимости от температуры
                         IntPotentialFunSdUOne,#Функции потенцалов взаимодействия в зависимости от обратной температуры
                         
                         isInvTemperatute=isInvTemperatute#По температурам или по обратным температурам считаем
                         )

    #Тело функции
    def __call__(self,
    
                 jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                 tempratures#Темперетауры по искомым потенциалам взаимодействия
                 ):
        #Выводим потенциалы взаимодействия в соответствие с индексами
        return self.GetIntPotentialFun()(jacSStateCoordinates[self._GetStateCoordinatesVarPotentialsInterIndexes()],
                                         tempratures[self._GetEnergyPowersVarPotentialsInterIndexes()])
