import numpy as np
from scipy import linalg as lg

#Функция температур
def TemperaturesFun(jacSUPowerEnergies#Матрица Якоби производной энтропии по внутренним энергиям энергетических степеней свободы
                    ):
    #Определяем и выводим приведенные температуры
    return 1 / jacSUPowerEnergies

#Функуия приведенных теплоемкостей и приведенных тепловых эффектов
def InvHeatCapacityHeatEffectsFun(heatCapacities,#Теплоемкости
                                  heatEffects#Тепловые эффекты
                                  ):
    #Определем обратные теплоемкости
    invHeatCapacities = lg.inv(heatCapacities)
    
    #Определяем приведенные тепловые эффекты
    redHeatEffects = -lg.solve(heatCapacities, heatEffects)
    
    #Выводим результат
    return (invHeatCapacities,
            redHeatEffects)
def InvHeatCapacityHeatEffectsOneFun(heatCapacities,#Теплоемкости
                                     heatEffects#Тепловые эффекты
                                     ):
    #Определем обратные теплоемкости
    invHeatCapacities = 1 / heatCapacities
    
    #Определяем приведенные тепловые эффекты
    redHeatEffects = -heatEffects / heatCapacities.reshape(-1,1)
    
    #Выводим результат
    return (invHeatCapacities,
            redHeatEffects)

#Функция тепловых эффектов
def HeatEffectsFun(hesSTemperaturesStateCoordinates,#Матрица Гесса приведенной энтропии по температурам и координатам состояния
                   temperatures#Температуры
                   ):
    #Определяем и выводим тепловые эффекты
    return hesSTemperaturesStateCoordinates*np.power(temperatures, 2).reshape(-1,1)

#Функця теплоемкостей
def HeatCapacitiesFun(jacSTemperatures,#Якобиан приведенной энтропии по температурам
                      hesSTemperatures,#Матрица Гесса приведенной энтропии по температурам
                      temperatures#Температуры
                      ):
    #Определяем и выводим теплоемкости
    return np.diag(2*jacSTemperatures*temperatures) + HeatEffectsFun(hesSTemperatures,
                                                                     temperatures)
def HeatCapacitiesOneFun(jacSTemperatures,#Якобиан приведенной энтропии по температурам
                         hesDiagSTemperatures,#Матрица Гесса (диагональная) приведенной энтропии по температурам
                         temperatures#Температуры
                         ):
    #Определяем и выводим теплоемкости
    return 2*jacSTemperatures*temperatures + \
           hesDiagSTemperatures*np.power(temperatures, 2)
