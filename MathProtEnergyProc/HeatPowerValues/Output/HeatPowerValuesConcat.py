import numpy as np

from MathProtEnergyProc.IndexFunctions import GetIndexes, GetPairIndexes

#Конкатенация теплоэнергетических величин
def HeatPowerValuesConcat(energyPowersTemperatureArr,#Температуры энергетических степеней свободы
                          potentialsInterArr,#Потенциалы взаимодействия энергетических степеней свободы
                          potentialsInterBetArr,#Потенциалы взаимодействия между энергетическими степенями свободы
                          invHeatCapacityArr,#Обратные теплоемкости энергетических степеней свободы 
                          redHeatEffectArr#Приведенные тепловые эффекты энергетических степеней свободы 
                          ):
    #Конкатенуем и выводим величины
    return (np.hstack(energyPowersTemperatureArr),#Температуры энергетических степеней свободы
            np.hstack(potentialsInterArr),#Потенциалы взаимодействия энергетических степеней свободы
            np.hstack(potentialsInterBetArr),#Потенциалы взаимодействия между энергетическими степенями свободы
            np.hstack(invHeatCapacityArr),#Обратные теплоемкости энергетических степеней свободы 
            np.hstack(redHeatEffectArr)#Приведенные тепловые эффекты энергетических степеней свободы 
            )

#Выходной функтор тепловых величин (теплоемкотей и тепловых эффектов)
class HeatPowerValuesConcatSelect:
    #Инициализация класса
    def __init__(self,
                 
                 energyPowersTemperatureNamesArr,#Имена температур энергетических степеней свободы
                 energyPowersPotentialsInterNamesArr,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                 stateCoordinatesPotentialsInterNamesArr,#Имена переменных потенциалов взаимодействия по координатам состояния
                 energyPowersPotentialsInterBetNamesArr,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                 stateCoordinatesPotentialsInterBetNamesArr,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                 redTemperaturesEnergyPowersInvHeatCapacityNamesArr,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                 energyPowersInvHeatCapacityNamesArr,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                 redTemperaturesEnergyPowersRedHeatEffectNamesArr,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                 stateCoordinatesRedHeatEffectNamesArr,#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                  
                 energyPowersTemperatureNamesConcat,#Имена температур энергетических степеней свободы
                 energyPowersPotentialsInterNamesConcat,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                 stateCoordinatesPotentialsInterNamesConcat,#Имена переменных потенциалов взаимодействия по координатам состояния
                 energyPowersPotentialsInterBetNamesConcat,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                 stateCoordinatesPotentialsInterBetNamesConcat,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                 redTemperaturesEnergyPowersInvHeatCapacityNamesConcat,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                 energyPowersInvHeatCapacityNamesConcat,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                 redTemperaturesEnergyPowersRedHeatEffectNamesConcat,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                 stateCoordinatesRedHeatEffectNamesConcat#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                 ):
        #Получаем индексы
        self.__EnergyPowersTemperatureIndexes = GetIndexes(np.hstack(energyPowersTemperatureNamesArr).tolist(),
                                                           energyPowersTemperatureNamesConcat)#Индексы температур энергетических степеней свободы
        self.__PotentialsInterIndexes = GetPairIndexes(np.hstack(energyPowersPotentialsInterNamesArr).tolist(),
                                                       np.hstack(stateCoordinatesPotentialsInterNamesArr).tolist(),
                                                       energyPowersPotentialsInterNamesConcat,
                                                       stateCoordinatesPotentialsInterNamesConcat)#Индексы потенциалов взваимодействия энергетических степеней свободы
        self.__PotentialsInterBetIndexes = GetPairIndexes(np.hstack(energyPowersPotentialsInterBetNamesArr).tolist(),
                                                          np.hstack(stateCoordinatesPotentialsInterBetNamesArr).tolist(),
                                                          energyPowersPotentialsInterBetNamesConcat,
                                                          stateCoordinatesPotentialsInterBetNamesConcat)#Индексы потенциалов взваимодействия между энергетическими степенями свободы
        self.__InvHeatCapacityIndexes = GetPairIndexes(np.hstack(redTemperaturesEnergyPowersInvHeatCapacityNamesArr).tolist(),
                                                       np.hstack(energyPowersInvHeatCapacityNamesArr).tolist(),
                                                       redTemperaturesEnergyPowersInvHeatCapacityNamesConcat,
                                                       energyPowersInvHeatCapacityNamesConcat)#Индексы обратных теплоемкостей энергетических степеней свободы 
        self.__RedHeatEffectIndexes = GetPairIndexes(np.hstack(redTemperaturesEnergyPowersRedHeatEffectNamesArr).tolist(),
                                                     np.hstack(stateCoordinatesRedHeatEffectNamesArr).tolist(),
                                                     redTemperaturesEnergyPowersRedHeatEffectNamesConcat,
                                                     stateCoordinatesRedHeatEffectNamesConcat)#Индексы приведенных тепловых эффектов энергетических степеней свободы 
            
    #Тело функции
    def __call__(self,
    
                 energyPowersTemperatureArr,#Температуры энергетических степеней свободы
                 potentialsInterArr,#Потенциалы взаимодействия энергетических степеней свободы
                 potentialsInterBetArr,#Потенциалы взаимодействия между энергетическими степенями свободы
                 invHeatCapacityArr,#Обратные теплоемкости энергетических степеней свободы 
                 redHeatEffectArr#Приведенные тепловые эффекты энергетических степеней свободы 
                 ):
        #Конкатенуем теплоэнергетические величины
        (energyPowersTemperature,#Температуры энергетических степеней свободы
         potentialsInter,#Потенциалы взаимодействия энергетических степеней свободы
         potentialsInterBet,#Потенциалы взаимодействия между энергетическими степенями свободы
         invHeatCapacity,#Обратные теплоемкости энергетических степеней свободы 
         redHeatEffect#Приведенные тепловые эффекты энергетических степеней свободы
         ) = HeatPowerValuesConcat(energyPowersTemperatureArr,#Температуры энергетических степеней свободы
                                   potentialsInterArr,#Потенциалы взаимодействия энергетических степеней свободы
                                   potentialsInterBetArr,#Потенциалы взаимодействия между энергетическими степенями свободы
                                   invHeatCapacityArr,#Обратные теплоемкости энергетических степеней свободы 
                                   redHeatEffectArr#Приведенные тепловые эффекты энергетических степеней свободы 
                                   )
        
        #Выводим теплоэнергетические величины
        return (energyPowersTemperature[self.__EnergyPowersTemperatureIndexes].copy(),#Температуры энергетических степеней свободы
                potentialsInter[self.__PotentialsInterIndexes].copy(),#Потенциалы взаимодействия энергетических степеней свободы
                potentialsInterBet[self.__PotentialsInterBetIndexes].copy(),#Потенциалы взаимодействия между энергетическими степенями свободы
                invHeatCapacity[self.__InvHeatCapacityIndexes].copy(),#Обратные теплоемкости энергетических степеней свободы 
                redHeatEffect[self.__RedHeatEffectIndexes].copy()#Приведенные тепловые эффекты энергетических степеней свободы 
                )
