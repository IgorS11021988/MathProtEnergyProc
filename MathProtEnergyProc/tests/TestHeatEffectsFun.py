import numpy as np

from MathProtEnergyProc.HeatPowerValues.Base import HeatEffectsFun

import unittest

#Модульные тесты
class TestHeatEffectsFun(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testHeatEffectsFun1(self):
        #Исходные данные
        hesSTemperaturesStateCoordinates = np.array([[8.49, -0.81, 0.027,  4.5,  -3.9],
                                                     [0.21,  6.93, 0.393, 6.51,   0.3],
                                                     [0.03, 0.063,   9.3, 2.31, -1.23]], dtype=np.double)#Матрица Гесса приведенной энтропии по температурам и координатам состояния
        temperatures = np.array([213.3, 369.3, 930.9], dtype=np.double)#Температуры
        
        #Эталонные значения
        heatEffectsEt = np.array([[8.49*213.3*213.3, -0.81*213.3*213.3, 0.027*213.3*213.3,  4.5*213.3*213.3,  -3.9*213.3*213.3],
                                  [0.21*369.3*369.3,  6.93*369.3*369.3, 0.393*369.3*369.3, 6.51*369.3*369.3,   0.3*369.3*369.3],
                                  [0.03*930.9*930.9, 0.063*930.9*930.9,   9.3*930.9*930.9, 2.31*930.9*930.9, -1.23*930.9*930.9]], dtype=np.double)
        
        #Вызываем функцию
        heatEffects = HeatEffectsFun(hesSTemperaturesStateCoordinates,#Матрица Гесса приведенной энтропии по температурам и координатам состояния
                                     temperatures#Температуры
                                     )
        
        #Проверяем значения
        errHeatEffects = np.max(np.abs(heatEffects - heatEffectsEt))
        self.assertAlmostEqual(errHeatEffects, 0.0, 8)
    def testHeatEffectsFun2(self):
        #Исходные данные
        hesSTemperaturesStateCoordinates = np.array([[8.49, -0.81, 0.027,  4.5,  -3.9,  1.5],
                                                     [0.21,  6.93, 0.393, 6.51,   0.3, 9.51],
                                                     [0.03, 0.063,   9.3, 2.31, -1.23, 5.61],
                                                     [1.65, 2.061,  12.3, 5.31,  4.23, 8.61]], dtype=np.double)#Матрица Гесса приведенной энтропии по температурам и координатам состояния
        temperatures = np.array([243.3, 639.3, 960.3, 153.3], dtype=np.double)#Температуры
        
        #Эталонные значения
        heatEffectsEt = np.array([[8.49*243.3*243.3, -0.81*243.3*243.3, 0.027*243.3*243.3,  4.5*243.3*243.3,  -3.9*243.3*243.3,  1.5*243.3*243.3],
                                  [0.21*639.3*639.3,  6.93*639.3*639.3, 0.393*639.3*639.3, 6.51*639.3*639.3,   0.3*639.3*639.3, 9.51*639.3*639.3],
                                  [0.03*960.3*960.3, 0.063*960.3*960.3,   9.3*960.3*960.3, 2.31*960.3*960.3, -1.23*960.3*960.3, 5.61*960.3*960.3],
                                  [1.65*153.3*153.3, 2.061*153.3*153.3,  12.3*153.3*153.3, 5.31*153.3*153.3,  4.23*153.3*153.3, 8.61*153.3*153.3]], dtype=np.double)
        
        #Вызываем функцию
        heatEffects = HeatEffectsFun(hesSTemperaturesStateCoordinates,#Матрица Гесса приведенной энтропии по температурам и координатам состояния
                                     temperatures#Температуры
                                     )
        
        #Проверяем значения
        errHeatEffects = np.max(np.abs(heatEffects - heatEffectsEt))
        self.assertAlmostEqual(errHeatEffects, 0.0, 9)
    def testHeatEffectsFun3(self):
        #Исходные данные
        hesSTemperaturesStateCoordinates = np.array([[8.49, -0.81, 0.027],
                                                     [0.21,  6.93, 0.393]], dtype=np.double)#Матрица Гесса приведенной энтропии по температурам и координатам состояния
        temperatures = np.array([123.3, 159.3], dtype=np.double)#Температуры
        
        #Эталонные значения
        heatEffectsEt = np.array([[8.49*123.3*123.3, -0.81*123.3*123.3, 0.027*123.3*123.3],
                                  [0.21*159.3*159.3,  6.93*159.3*159.3, 0.393*159.3*159.3]], dtype=np.double)
        
        #Вызываем функцию
        heatEffects = HeatEffectsFun(hesSTemperaturesStateCoordinates,#Матрица Гесса приведенной энтропии по температурам и координатам состояния
                                     temperatures#Температуры
                                     )
        
        #Проверяем значения
        errHeatEffects = np.max(np.abs(heatEffects - heatEffectsEt))
        self.assertAlmostEqual(errHeatEffects, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
