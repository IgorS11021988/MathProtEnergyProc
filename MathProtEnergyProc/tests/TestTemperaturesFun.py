import numpy as np

from MathProtEnergyProc.HeatPowerValues.Base import TemperaturesFun

import unittest

#Модульные тесты
class TestTemperaturesFun(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testTemperaturesFun1(self):
        #Исходные данные
        jacSUPowerEnergies = np.array([1 / 315, 1 / 327, 1 / 213], dtype=np.double)#Матрица Якоби энтропии по внутренним энергиям
        
        #Эталонный результат
        temperaturesEt = np.array([315, 327, 213], dtype=np.double)
        
        #Вызываем функцию
        temperatures = TemperaturesFun(jacSUPowerEnergies)
        
        #Проверяем значения
        err = np.max(np.abs(temperatures - temperaturesEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testTemperaturesFun2(self):
        #Исходные данные
        jacSUPowerEnergies = np.array([1 / 322.5, 1 / 311.1, 1 / 180.3, 1 / 150.9, 1 / 369.3], dtype=np.double)#Матрица Якоби энтропии по внутренним энергиям
        
        #Эталонный результат
        temperaturesEt = np.array([322.5, 311.1, 180.3, 150.9, 369.3], dtype=np.double)
        
        #Вызываем функцию
        temperatures = TemperaturesFun(jacSUPowerEnergies)
        
        #Проверяем значения
        err = np.max(np.abs(temperatures - temperaturesEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testTemperaturesFun3(self):
        #Исходные данные
        jacSUPowerEnergies = np.array([1 / 312.3, 1 / 441.9], dtype=np.double)#Матрица Якоби энтропии по внутренним энергиям
        
        #Эталонный результат
        temperaturesEt = np.array([312.3, 441.9], dtype=np.double)
        
        #Вызываем функцию
        temperatures = TemperaturesFun(jacSUPowerEnergies)
        
        #Проверяем значения
        err = np.max(np.abs(temperatures - temperaturesEt))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
