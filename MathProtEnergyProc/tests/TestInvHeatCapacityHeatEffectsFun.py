import numpy as np

from MathProtEnergyProc.HeatPowerValues.Base import InvHeatCapacityHeatEffectsFun

import unittest

#Модульные тесты
class TestInvHeatCapacityHeatEffectsFun(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testInvHeatCapacityHeatEffectsFun1(self):
        #Исходные данные
        heatCapacities = np.array([[ 2.13, -0.51,  0.03],
                                   [-0.51,   3.3, 0.063],
                                   [ 0.03, 0.063,   9.3]], dtype=np.double)#Теплоемкости
        heatEffects = np.array([[8.49, -0.81, 0.027,  4.5,  -3.9],
                                [0.21,  6.93, 0.393, 6.51,   0.3],
                                [0.03, 0.063,   9.3, 2.31, -1.23]], dtype=np.double)#Тепловые эффекты
        
        #Вызываем функцию
        (invHeatCapacities,
         redHeatEffects) = InvHeatCapacityHeatEffectsFun(heatCapacities,#Теплоемкости
                                                         heatEffects#Тепловые эффекты
                                                         )
        
        #Проверяем значения
        err = np.max(np.abs(np.dot(heatCapacities, invHeatCapacities) - np.eye(3)))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(np.dot(heatCapacities, redHeatEffects) + heatEffects))
        self.assertAlmostEqual(err, 0.0, 9)
    def testInvHeatCapacityHeatEffectsFun2(self):
        #Исходные данные
        heatCapacities = np.array([[ 5.13, -0.87,  0.03,  1.11],
                                   [-0.81,   6.9, 0.069,  2.55],
                                   [ 0.03, 0.063,   9.3,  3.69],
                                   [ 1.11,  2.25,  3.69, 30.69]], dtype=np.double)#Теплоемкости
        heatEffects = np.array([[8.49, -0.81, 0.027,  4.5,  -3.9,   6.9],
                                [0.21,  6.93, 0.393, 6.51,   0.3,   3.3],
                                [0.03, 0.063,   9.3, 2.31, -1.23, -7.23],
                                [0.93, 0.963,  15.3, 8.31, -4.23,  1.23]], dtype=np.double)#Тепловые эффекты
        
        #Вызываем функцию
        (invHeatCapacities,
         redHeatEffects) = InvHeatCapacityHeatEffectsFun(heatCapacities,#Теплоемкости
                                                         heatEffects#Тепловые эффекты
                                                         )
        
        #Проверяем значения
        err = np.max(np.abs(np.dot(heatCapacities, invHeatCapacities) - np.eye(4)))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(np.dot(heatCapacities, redHeatEffects) + heatEffects))
        self.assertAlmostEqual(err, 0.0, 9)
    def testInvHeatCapacityHeatEffectsFun3(self):
        #Исходные данные
        heatCapacities = np.array([[20.13, 0.51],
                                   [ 0.51, 69.3]], dtype=np.double)#Теплоемкости
        heatEffects = np.array([[2.49,  6.21, 3.027,  7.5,  -1.5],
                                [3.21, -0.93, 0.393, 3.21, 0.321]], dtype=np.double)#Тепловые эффекты
        
        #Вызываем функцию
        (invHeatCapacities,
         redHeatEffects) = InvHeatCapacityHeatEffectsFun(heatCapacities,#Теплоемкости
                                                         heatEffects#Тепловые эффекты
                                                         )
        
        #Проверяем значения
        err = np.max(np.abs(np.dot(heatCapacities, invHeatCapacities) - np.eye(2)))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(np.dot(heatCapacities, redHeatEffects) + heatEffects))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
