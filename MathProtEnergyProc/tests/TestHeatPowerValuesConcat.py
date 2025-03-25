import numpy as np

from MathProtEnergyProc.HeatPowerValues import HeatPowerValuesConcat

import unittest

#Модульные тесты
class TestHeatPowerValuesConcat(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testHeatPowerValuesConcat1(self):
        #Исходные данные
        energyPowersTemperatureArr = [[123.45, 651.3, 312.3], 453.9, np.array([345.15, 159.183], dtype=np.double), 33.9]#Температуры энергетических степеней свободы
        potentialsInterArr = [[2.19, -5.31], 21.3]#Потенциалы взаимодействия энергетических степеней свободы
        potentialsInterBetArr = [39.63, np.array([-3.33, 0.45], dtype=np.double), -4.5]#Потенциалы взаимодействия между энергетическими степенями свободы
        invHeatCapacityArr = [[75.21, 727.5],[4.59, 6.51, 3.3]]#Обратные теплоемкости энергетических степеней свободы 
        redHeatEffectArr = np.array([3.33, 6.57, 8.1, -6.39], dtype=np.double)#Приведенные тепловые эффекты энергетических степеней свободы 
        
        #Эталонные значения
        energyPowersTemperatureEt = np.array([123.45, 651.3, 312.3, 453.9, 345.15, 159.183, 33.9], dtype=np.double)#Температуры энергетических степеней свободы
        potentialsInterEt = np.array([2.19, -5.31, 21.3], dtype=np.double)#Потенциалы взаимодействия энергетических степеней свободы
        potentialsInterBetEt = np.array([39.63, -3.33, 0.45, -4.5], dtype=np.double),#Потенциалы взаимодействия между энергетическими степенями свободы
        invHeatCapacityEt = np.array([75.21, 727.5,4.59, 6.51, 3.3], dtype=np.double)#Обратные теплоемкости энергетических степеней свободы 
        redHeatEffectEt = np.array([3.33, 6.57, 8.1, -6.39], dtype=np.double)#Приведенные тепловые эффекты энергетических степеней свободы 
        
        #Вызываем функцию
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
    def testHeatPowerValuesConcat2(self):
        #Исходные данные
        energyPowersTemperatureArr = [[423.45, 351.9, 612.3], 63.9, 753.9, np.array([345.15, 159.183], dtype=np.double), 33.9, 333.60]#Температуры энергетических степеней свободы
        potentialsInterArr = [[4.563, 63.81], 97.5, [2.19, -5.31], 81.63, 56.1]#Потенциалы взаимодействия энергетических степеней свободы
        potentialsInterBetArr = [39.63, np.array([-3.33, 0.45, -93.69], dtype=np.double), -4.5, [34.5, 87.3]]#Потенциалы взаимодействия между энергетическими степенями свободы
        invHeatCapacityArr = [72.21, 757.5,7.59, 6.51, 3.3, 76.5]#Обратные теплоемкости энергетических степеней свободы 
        redHeatEffectArr = [3.33, [6.57, 8.1], -6.39, -1.23]#Приведенные тепловые эффекты энергетических степеней свободы 
        
        #Эталонные значения
        energyPowersTemperatureEt = np.array([423.45, 351.9, 612.3, 63.9, 753.9, 345.15, 159.183, 33.9, 333.60], dtype=np.double)#Температуры энергетических степеней свободы
        potentialsInterEt = np.array([4.563, 63.81, 97.5, 2.19, -5.31, 81.63, 56.1], dtype=np.double)#Потенциалы взаимодействия энергетических степеней свободы
        potentialsInterBetEt = np.array([39.63, -3.33, 0.45, -93.69, -4.5, 34.5, 87.3], dtype=np.double),#Потенциалы взаимодействия между энергетическими степенями свободы
        invHeatCapacityEt = np.array([72.21, 757.5,7.59, 6.51, 3.3, 76.5], dtype=np.double)#Обратные теплоемкости энергетических степеней свободы 
        redHeatEffectEt = np.array([3.33, 6.57, 8.1, -6.39, -1.23], dtype=np.double)#Приведенные тепловые эффекты энергетических степеней свободы 
        
        #Вызываем функцию
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
    def testHeatPowerValuesConcat3(self):
        #Исходные данные
        energyPowersTemperatureArr = [[123.45, 312.3], 453.9, np.array([345.15], dtype=np.double)]#Температуры энергетических степеней свободы
        potentialsInterArr = [5.13, -2.31, 81.3]#Потенциалы взаимодействия энергетических степеней свободы
        potentialsInterBetArr = [39.63, np.array([9.63, -6.75], dtype=np.double)]#Потенциалы взаимодействия между энергетическими степенями свободы
        invHeatCapacityArr = [np.array([75.21, 727.5], dtype=np.double),[4.59, 6.51, 3.3], 88.83]#Обратные теплоемкости энергетических степеней свободы 
        redHeatEffectArr = [[3.33, 6.57], 8.1, -6.39]#Приведенные тепловые эффекты энергетических степеней свободы 
        
        #Эталонные значения
        energyPowersTemperatureEt = np.array([123.45, 312.3, 453.9, 345.15], dtype=np.double)#Температуры энергетических степеней свободы
        potentialsInterEt = np.array([5.13, -2.31, 81.3], dtype=np.double)#Потенциалы взаимодействия энергетических степеней свободы
        potentialsInterBetEt = np.array([39.63, 9.63, -6.75], dtype=np.double),#Потенциалы взаимодействия между энергетическими степенями свободы
        invHeatCapacityEt = np.array([75.21, 727.5,4.59, 6.51, 3.3, 88.83], dtype=np.double)#Обратные теплоемкости энергетических степеней свободы 
        redHeatEffectEt = np.array([3.33, 6.57, 8.1, -6.39], dtype=np.double)#Приведенные тепловые эффекты энергетических степеней свободы 
        
        #Вызываем функцию
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
