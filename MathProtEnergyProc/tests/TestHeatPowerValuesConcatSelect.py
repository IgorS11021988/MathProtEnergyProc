import numpy as np

from MathProtEnergyProc.HeatPowerValues import HeatPowerValuesConcatSelect

import unittest

#Модульные тесты
class TestHeatPowerValuesConcatSelect(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testHeatPowerValuesConcatSelect1(self):
        #Исходные данные
        energyPowersTemperatureNamesArr = [["e1","e2","e5"],"e3",["e10","e15"],"e8"]#Имена температур энергетических степеней свободы
        energyPowersPotentialsInterNamesArr     = [["e1","e1"],"e2"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesPotentialsInterNamesArr = [["x1","x2"],"x2"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersPotentialsInterBetNamesArr     = ["e1",["e2","e2"],"e7"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        stateCoordinatesPotentialsInterBetNamesArr = ["x1",["x1","x2"],"x5"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        redTemperaturesEnergyPowersInvHeatCapacityNamesArr = [["e1","e2"],["e2","e2","e7"]]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersInvHeatCapacityNamesArr                = [["e1","e2"],["e1","e5","e5"]]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        redTemperaturesEnergyPowersRedHeatEffectNamesArr = ["e10","e15","e12","e17"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        stateCoordinatesRedHeatEffectNamesArr            = [ "x1", "x2", "x4", "x8"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        energyPowersTemperatureNamesConcat = ["e2","e10","e5","e1","e3","e15","e8"]#Имена температур энергетических степеней свободы
        energyPowersPotentialsInterNamesConcat     = ["e1","e1","e2"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesPotentialsInterNamesConcat = ["x2","x1","x2"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersPotentialsInterBetNamesConcat     = ["e2","e1","e2","e7"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        stateCoordinatesPotentialsInterBetNamesConcat = ["x2","x1","x1","x5"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        redTemperaturesEnergyPowersInvHeatCapacityNamesConcat = ["e2","e1","e2","e2","e7"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersInvHeatCapacityNamesConcat                = ["e1","e1","e2","e5","e5"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        redTemperaturesEnergyPowersRedHeatEffectNamesConcat = ["e12","e10","e15","e17"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        stateCoordinatesRedHeatEffectNamesConcat            = [ "x4", "x1", "x2", "x8"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        energyPowersTemperatureArr = [[123.45, 651.3, 312.3], 453.9, np.array([345.15, 159.183], dtype=np.double), 33.9]#Температуры энергетических степеней свободы
        potentialsInterArr = [[2.19, -5.31], 21.3]#Потенциалы взаимодействия энергетических степеней свободы
        potentialsInterBetArr = [39.63, np.array([-3.33, 0.45], dtype=np.double), -4.5]#Потенциалы взаимодействия между энергетическими степенями свободы
        invHeatCapacityArr = [[75.21, 727.5],[4.59, 6.51, 3.3]]#Обратные теплоемкости энергетических степеней свободы 
        redHeatEffectArr = np.array([3.33, 6.57, 8.1, -6.39], dtype=np.double)#Приведенные тепловые эффекты энергетических степеней свободы 
        
        #Эталонные значения
        energyPowersTemperatureEt = np.array([651.3, 345.15, 312.3, 123.45, 453.9, 159.183, 33.9], dtype=np.double)#Температуры энергетических степеней свободы
        potentialsInterEt = np.array([-5.31, 2.19, 21.3], dtype=np.double)#Потенциалы взаимодействия энергетических степеней свободы
        potentialsInterBetEt = np.array([0.45, 39.63, -3.33, -4.5], dtype=np.double),#Потенциалы взаимодействия между энергетическими степенями свободы
        invHeatCapacityEt = np.array([4.59, 75.21, 727.5,6.51, 3.3], dtype=np.double)#Обратные теплоемкости энергетических степеней свободы 
        redHeatEffectEt = np.array([8.1, 3.33, 6.57, -6.39], dtype=np.double)#Приведенные тепловые эффекты энергетических степеней свободы 
        
        #Создаем конкатенатор
        heatPowerValuesConcatSelect = HeatPowerValuesConcatSelect(energyPowersTemperatureNamesArr,#Имена температур энергетических степеней свободы
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
                                                                  )
        
        #Вызываем функцию
        (energyPowersTemperature,#Температуры энергетических степеней свободы
         potentialsInter,#Потенциалы взаимодействия энергетических степеней свободы
         potentialsInterBet,#Потенциалы взаимодействия между энергетическими степенями свободы
         invHeatCapacity,#Обратные теплоемкости энергетических степеней свободы 
         redHeatEffect#Приведенные тепловые эффекты энергетических степеней свободы
         ) = heatPowerValuesConcatSelect(energyPowersTemperatureArr,#Температуры энергетических степеней свободы
                                         potentialsInterArr,#Потенциалы взаимодействия энергетических степеней свободы
                                         potentialsInterBetArr,#Потенциалы взаимодействия между энергетическими степенями свободы
                                         invHeatCapacityArr,#Обратные теплоемкости энергетических степеней свободы 
                                         redHeatEffectArr#Приведенные тепловые эффекты энергетических степеней свободы 
                                         )
        
        #Проверяем значения
        err = np.max(np.abs(energyPowersTemperature - energyPowersTemperatureEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(potentialsInter - potentialsInterEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(potentialsInterBet - potentialsInterBetEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(invHeatCapacity - invHeatCapacityEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(redHeatEffect - redHeatEffectEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testHeatPowerValuesConcatSelect2(self):
        #Исходные данные
        energyPowersTemperatureNamesArr = [["e1","e2","e5"],"e3","e4",["e10","e15"],"e8","e11"]#Имена температур энергетических степеней свободы
        energyPowersPotentialsInterNamesArr     = [["e1","e1"],"e2",["e3","e4"],"e11","e8"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesPotentialsInterNamesArr = [["x1","x2"],"x2",["x4","x3"], "x2","x2"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersPotentialsInterBetNamesArr     = ["e1",["e2","e2","e8"],"e7",["e2","e8"]]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        stateCoordinatesPotentialsInterBetNamesArr = ["x1",["x1","x2","x4"],"x5",["x7","x7"]]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        redTemperaturesEnergyPowersInvHeatCapacityNamesArr = ["e1","e2","e2","e2","e12","e17"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersInvHeatCapacityNamesArr                = ["e1","e2","e1","e5", "e1", "e2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        redTemperaturesEnergyPowersRedHeatEffectNamesArr = ["e10",["e15","e12"],"e17","e17"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        stateCoordinatesRedHeatEffectNamesArr            = [ "x1",[ "x2", "x4"], "x8", "x1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        energyPowersTemperatureNamesConcat = ["e10","e3","e1","e2","e5","e4","e15","e8","e11"]#Имена температур энергетических степеней свободы
        energyPowersPotentialsInterNamesConcat     = ["e3","e4","e1","e1","e2","e11","e8"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesPotentialsInterNamesConcat = ["x4","x3","x1","x2","x2", "x2","x2"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersPotentialsInterBetNamesConcat     = ["e2","e8","e1","e2","e2","e7","e8"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        stateCoordinatesPotentialsInterBetNamesConcat = ["x7","x4","x1","x1","x2","x5","x7"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        redTemperaturesEnergyPowersInvHeatCapacityNamesConcat = ["e1","e2","e2","e2","e12","e17"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersInvHeatCapacityNamesConcat                = ["e1","e2","e1","e5", "e1", "e2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        redTemperaturesEnergyPowersRedHeatEffectNamesConcat = ["e17","e12","e10","e15","e17"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        stateCoordinatesRedHeatEffectNamesConcat            = [ "x8", "x4", "x1", "x2", "x1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        energyPowersTemperatureArr = [[423.45, 351.9, 612.3], 63.9, 753.9, np.array([345.15, 159.183], dtype=np.double), 33.9, 333.60]#Температуры энергетических степеней свободы
        potentialsInterArr = [[4.563, 63.81], 97.5, [2.19, -5.31], 81.63, 56.1]#Потенциалы взаимодействия энергетических степеней свободы
        potentialsInterBetArr = [39.63, np.array([-3.33, 0.45, -93.69], dtype=np.double), -4.5, [34.5, 87.3]]#Потенциалы взаимодействия между энергетическими степенями свободы
        invHeatCapacityArr = [72.21, 757.5,7.59, 6.51, 3.3, 76.5]#Обратные теплоемкости энергетических степеней свободы 
        redHeatEffectArr = [3.33, [6.57, 8.1], -6.39, -1.23]#Приведенные тепловые эффекты энергетических степеней свободы 
        
        #Эталонные значения
        energyPowersTemperatureEt = np.array([345.15, 63.9, 423.45, 351.9, 612.3, 753.9, 159.183, 33.9, 333.60], dtype=np.double)#Температуры энергетических степеней свободы
        potentialsInterEt = np.array([2.19, -5.31, 4.563, 63.81, 97.5, 81.63, 56.1], dtype=np.double)#Потенциалы взаимодействия энергетических степеней свободы
        potentialsInterBetEt = np.array([34.5, -93.69, 39.63, -3.33, 0.45, -4.5, 87.3], dtype=np.double),#Потенциалы взаимодействия между энергетическими степенями свободы
        invHeatCapacityEt = np.array([72.21, 757.5,7.59, 6.51, 3.3, 76.5], dtype=np.double)#Обратные теплоемкости энергетических степеней свободы 
        redHeatEffectEt = np.array([-6.39, 8.1, 3.33, 6.57, -1.23], dtype=np.double)#Приведенные тепловые эффекты энергетических степеней свободы 
        
        #Создаем конкатенатор
        heatPowerValuesConcatSelect = HeatPowerValuesConcatSelect(energyPowersTemperatureNamesArr,#Имена температур энергетических степеней свободы
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
                                                                  )
        
        #Вызываем функцию
        (energyPowersTemperature,#Температуры энергетических степеней свободы
         potentialsInter,#Потенциалы взаимодействия энергетических степеней свободы
         potentialsInterBet,#Потенциалы взаимодействия между энергетическими степенями свободы
         invHeatCapacity,#Обратные теплоемкости энергетических степеней свободы 
         redHeatEffect#Приведенные тепловые эффекты энергетических степеней свободы
         ) = heatPowerValuesConcatSelect(energyPowersTemperatureArr,#Температуры энергетических степеней свободы
                                         potentialsInterArr,#Потенциалы взаимодействия энергетических степеней свободы
                                         potentialsInterBetArr,#Потенциалы взаимодействия между энергетическими степенями свободы
                                         invHeatCapacityArr,#Обратные теплоемкости энергетических степеней свободы 
                                         redHeatEffectArr#Приведенные тепловые эффекты энергетических степеней свободы 
                                         )
        
        #Проверяем значения
        err = np.max(np.abs(energyPowersTemperature - energyPowersTemperatureEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(potentialsInter - potentialsInterEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(potentialsInterBet - potentialsInterBetEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(invHeatCapacity - invHeatCapacityEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(redHeatEffect - redHeatEffectEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testHeatPowerValuesConcatSelect3(self):
        #Исходные данные
        energyPowersTemperatureNamesArr = [["e1","e2"],"e3",["e10"]]#Имена температур энергетических степеней свободы
        energyPowersPotentialsInterNamesArr     = ["e1","e1","e2"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesPotentialsInterNamesArr = ["x1","x2","x2"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersPotentialsInterBetNamesArr     = ["e1",["e2","e8"]]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        stateCoordinatesPotentialsInterBetNamesArr = ["x1",["x7","x7"]]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        redTemperaturesEnergyPowersInvHeatCapacityNamesArr = [["e1","e2"],["e2","e2","e12"],"e17"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersInvHeatCapacityNamesArr                = [["e1","e2"],["e1","e5", "e1"], "e2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        redTemperaturesEnergyPowersRedHeatEffectNamesArr = [["e15","e12"],"e17","e17"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        stateCoordinatesRedHeatEffectNamesArr            = [[ "x2", "x4"], "x8", "x1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        energyPowersTemperatureNamesConcat = ["e3","e1","e2","e10"]#Имена температур энергетических степеней свободы
        energyPowersPotentialsInterNamesConcat     = ["e2","e1","e1"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesPotentialsInterNamesConcat = ["x2","x1","x2"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersPotentialsInterBetNamesConcat     = ["e1","e2","e8"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        stateCoordinatesPotentialsInterBetNamesConcat = ["x1","x7","x7"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        redTemperaturesEnergyPowersInvHeatCapacityNamesConcat = ["e12","e1","e2","e2","e2","e17"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersInvHeatCapacityNamesConcat                = [ "e1","e1","e2","e1","e5", "e2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        redTemperaturesEnergyPowersRedHeatEffectNamesConcat = ["e15","e12","e17","e17"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        stateCoordinatesRedHeatEffectNamesConcat            = [ "x2", "x4", "x8", "x1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        energyPowersTemperatureArr = [[123.45, 312.3], 453.9, np.array([345.15], dtype=np.double)]#Температуры энергетических степеней свободы
        potentialsInterArr = [5.13, -2.31, 81.3]#Потенциалы взаимодействия энергетических степеней свободы
        potentialsInterBetArr = [39.63, np.array([9.63, -6.75], dtype=np.double)]#Потенциалы взаимодействия между энергетическими степенями свободы
        invHeatCapacityArr = [np.array([75.21, 727.5], dtype=np.double),[4.59, 6.51, 3.3], 88.83]#Обратные теплоемкости энергетических степеней свободы 
        redHeatEffectArr = [[3.33, 6.57], 8.1, -6.39]#Приведенные тепловые эффекты энергетических степеней свободы 
        
        #Эталонные значения
        energyPowersTemperatureEt = np.array([453.9, 123.45, 312.3, 345.15], dtype=np.double)#Температуры энергетических степеней свободы
        potentialsInterEt = np.array([81.3, 5.13, -2.31], dtype=np.double)#Потенциалы взаимодействия энергетических степеней свободы
        potentialsInterBetEt = np.array([39.63, 9.63, -6.75], dtype=np.double),#Потенциалы взаимодействия между энергетическими степенями свободы
        invHeatCapacityEt = np.array([3.3, 75.21, 727.5,4.59, 6.51, 88.83], dtype=np.double)#Обратные теплоемкости энергетических степеней свободы 
        redHeatEffectEt = np.array([3.33, 6.57, 8.1, -6.39], dtype=np.double)#Приведенные тепловые эффекты энергетических степеней свободы 
        
        #Создаем конкатенатор
        heatPowerValuesConcatSelect = HeatPowerValuesConcatSelect(energyPowersTemperatureNamesArr,#Имена температур энергетических степеней свободы
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
                                                                  )
        
        #Вызываем функцию
        (energyPowersTemperature,#Температуры энергетических степеней свободы
         potentialsInter,#Потенциалы взаимодействия энергетических степеней свободы
         potentialsInterBet,#Потенциалы взаимодействия между энергетическими степенями свободы
         invHeatCapacity,#Обратные теплоемкости энергетических степеней свободы 
         redHeatEffect#Приведенные тепловые эффекты энергетических степеней свободы
         ) = heatPowerValuesConcatSelect(energyPowersTemperatureArr,#Температуры энергетических степеней свободы
                                         potentialsInterArr,#Потенциалы взаимодействия энергетических степеней свободы
                                         potentialsInterBetArr,#Потенциалы взаимодействия между энергетическими степенями свободы
                                         invHeatCapacityArr,#Обратные теплоемкости энергетических степеней свободы 
                                         redHeatEffectArr#Приведенные тепловые эффекты энергетических степеней свободы 
                                         )
        
        #Проверяем значения
        err = np.max(np.abs(energyPowersTemperature - energyPowersTemperatureEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(potentialsInter - potentialsInterEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(potentialsInterBet - potentialsInterBetEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(invHeatCapacity - invHeatCapacityEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(redHeatEffect - redHeatEffectEt))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
