import numpy as np

from MathProtEnergyProc.HeatPowerValues.Base import HeatCapacitiesOneFun

import unittest

#Модульные тесты
class TestHeatCapacitiesOneFun(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testHeatCapacitiesOneFun1(self):
        #Исходные данные
        jacSTemperatures = np.array([1.53, 0.21, 3.03], dtype=np.double)#Якобиан приведенной энтропии по температурам
        hesDiagSTemperatures = np.array([8.49, 6.93, 9.3], dtype=np.double)#Матрица Гесса приведенной энтропии по температурам
        temperatures = np.array([213.3, 369.3, 930.9], dtype=np.double)#Температуры
        
        #Эталонные значения
        heatCapacitiesEt = np.array([2*1.53*213.3 + 8.49*213.3*213.3, 2*0.21*369.3 + 6.93*369.3*369.3, 2*3.03*930.9 + 9.3*930.9*930.9], dtype=np.double)
        
        #Вызываем функцию
        heatCapacities = HeatCapacitiesOneFun(jacSTemperatures,#Якобиан приведенной энтропии по температурам
                                              hesDiagSTemperatures,#Матрица Гесса приведенной энтропии по температурам
                                              temperatures#Температуры
                                              )
        
        #Проверяем значения
        err = np.max(np.abs(heatCapacities - heatCapacitiesEt))
        self.assertAlmostEqual(err, 0.0, 8)
    def testHeatCapacitiesOneFun2(self):
        #Исходные данные
        jacSTemperatures = np.array([0.33, 2.1, 5.37, 4.35], dtype=np.double)#Якобиан приведенной энтропии по температурам
        hesDiagSTemperatures = np.array([8.49, 6.93, 9.3, 5.31], dtype=np.double)#Матрица Гесса приведенной энтропии по температурам
        temperatures = np.array([243.3, 639.3, 960.3, 153.3], dtype=np.double)#Температуры
        
        #Эталонные значения
        heatCapacitiesEt = np.array([2*0.33*243.3 + 8.49*243.3*243.3, 2*2.1*639.3 + 6.93*639.3*639.3, 2*5.37*960.3 + 9.3*960.3*960.3, 2*4.35*153.3 + 5.31*153.3*153.3], dtype=np.double)
        
        #Вызываем функцию
        heatCapacities = HeatCapacitiesOneFun(jacSTemperatures,#Якобиан приведенной энтропии по температурам
                                              hesDiagSTemperatures,#Матрица Гесса приведенной энтропии по температурам
                                              temperatures#Температуры
                                              )
        
        #Проверяем значения
        err = np.max(np.abs(heatCapacities - heatCapacitiesEt))
        self.assertAlmostEqual(err, 0.0, 8)
    def testHeatCapacitiesOneFun3(self):
        #Исходные данные
        jacSTemperatures = np.array([0.21, 3.9], dtype=np.double)#Якобиан приведенной энтропии по температурам
        hesDiagSTemperatures = np.array([8.49, 6.93], dtype=np.double)#Матрица Гесса приведенной энтропии по температурам
        temperatures = np.array([123.3, 159.3], dtype=np.double)#Температуры
        
        #Эталонные значения
        heatCapacitiesEt = np.array([2*0.21*123.3 + 8.49*123.3*123.3, 2*3.9*159.3 + 6.93*159.3*159.3], dtype=np.double)
        
        #Вызываем функцию
        heatCapacities = HeatCapacitiesOneFun(jacSTemperatures,#Якобиан приведенной энтропии по температурам
                                              hesDiagSTemperatures,#Матрица Гесса приведенной энтропии по температурам
                                              temperatures#Температуры
                                              )
        
        #Проверяем значения
        err = np.max(np.abs(heatCapacities - heatCapacitiesEt))
        self.assertAlmostEqual(err, 0.0, 8)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
