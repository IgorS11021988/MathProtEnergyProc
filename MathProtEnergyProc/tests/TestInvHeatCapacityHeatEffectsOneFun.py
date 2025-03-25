import numpy as np

from MathProtEnergyProc.HeatPowerValues.Base import InvHeatCapacityHeatEffectsOneFun

import unittest

#Модульные тесты
class TestInvHeatCapacityHeatEffectsOneFun(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testInvHeatCapacityHeatEffectsOneFun1(self):
        #Исходные данные
        heatCapacities = np.array([2.13, 3.3, 9.3], dtype=np.double)#Теплоемкости
        heatEffects = np.array([[8.49, -0.81, 0.027,  4.5,  -3.9],
                                [0.21,  6.93, 0.393, 6.51,   0.3],
                                [0.03, 0.063,   9.3, 2.31, -1.23]], dtype=np.double)#Тепловые эффекты
        
        #Эталонные значения
        invHeatCapacitiesEt = np.array([1 / 2.13, 1 / 3.3, 1 / 9.3], dtype=np.double)#Обратные теплоемкости
        redHeatEffectsEt = np.array([[-8.49 / 2.13,   0.81 / 2.13, -0.027 / 2.13,  -4.5 / 2.13,  3.9 / 2.13],
                                     [-0.21 /  3.3,  -6.93 /  3.3, -0.393  / 3.3, -6.51  / 3.3, -0.3  / 3.3],
                                     [-0.03 /  9.3, -0.063 /  9.3,   -9.3 /  9.3, -2.31 /  9.3, 1.23 /  9.3]], dtype=np.double)#Приведенные тепловые эффекты
        
        #Вызываем функцию
        (invHeatCapacities,
         redHeatEffects) = InvHeatCapacityHeatEffectsOneFun(heatCapacities,#Теплоемкости
                                                            heatEffects#Тепловые эффекты
                                                            )
        
        #Проверяем значения
        errInvHeatCapacities = np.max(np.abs(invHeatCapacities - invHeatCapacitiesEt))
        self.assertAlmostEqual(errInvHeatCapacities, 0.0, 9)
        errRedHeatEffects = np.max(np.abs(redHeatEffects - redHeatEffectsEt))
        self.assertAlmostEqual(errRedHeatEffects, 0.0, 9)
    def testInvHeatCapacityHeatEffectsOneFun2(self):
        #Исходные данные
        heatCapacities = np.array([5.13, 6.9, 9.3, 30.69], dtype=np.double)#Теплоемкости
        heatEffects = np.array([[8.49, -0.81, 0.027,  4.5,  -3.9,   6.9],
                                [0.21,  6.93, 0.393, 6.51,   0.3,   3.3],
                                [0.03, 0.063,   9.3, 2.31, -1.23, -7.23],
                                [0.93, 0.963,  15.3, 8.31, -4.23,  1.23]], dtype=np.double)#Тепловые эффекты
        
        #Эталонные значения
        invHeatCapacitiesEt = np.array([1 / 5.13, 1 / 6.9, 1 / 9.3, 1 / 30.69], dtype=np.double)#Обратные теплоемкости
        redHeatEffectsEt = np.array([[-8.49 /  5.13,   0.81 /  5.13, -0.027 /  5.13,  -4.5 /  5.13,  3.9 /  5.13,  -6.9 /  5.13],
                                     [-0.21 /   6.9,  -6.93 /   6.9, -0.393 /   6.9, -6.51 /   6.9, -0.3 /   6.9,  -3.3 /   6.9],
                                     [-0.03 /   9.3, -0.063 /   9.3,   -9.3 /   9.3, -2.31 /   9.3, 1.23 /   9.3,  7.23 /   9.3],
                                     [-0.93 / 30.69, -0.963 / 30.69,  -15.3 / 30.69, -8.31 / 30.69, 4.23 / 30.69, -1.23 / 30.69]], dtype=np.double)#Приведенные тепловые эффекты
        
        #Вызываем функцию
        (invHeatCapacities,
         redHeatEffects) = InvHeatCapacityHeatEffectsOneFun(heatCapacities,#Теплоемкости
                                                            heatEffects#Тепловые эффекты
                                                            )
        
        #Проверяем значения
        errInvHeatCapacities = np.max(np.abs(invHeatCapacities - invHeatCapacitiesEt))
        self.assertAlmostEqual(errInvHeatCapacities, 0.0, 9)
        errRedHeatEffects = np.max(np.abs(redHeatEffects - redHeatEffectsEt))
        self.assertAlmostEqual(errRedHeatEffects, 0.0, 9)
    def testInvHeatCapacityHeatEffectsOneFun3(self):
        #Исходные данные
        heatCapacities = np.array([20.13, 69.3], dtype=np.double)#Теплоемкости
        heatEffects = np.array([[2.49,  6.21, 3.027,  7.5,  -1.5],
                                [3.21, -0.93, 0.393, 3.21, 0.321]], dtype=np.double)#Тепловые эффекты
        
        #Эталонные значения
        invHeatCapacitiesEt = np.array([1 / 20.13, 1 / 69.3], dtype=np.double)#Обратные теплоемкости
        redHeatEffectsEt = np.array([[-2.49 / 20.13, -6.21 / 20.13, -3.027 / 20.13,  -7.5 / 20.13,    1.5 / 20.13],
                                     [-3.21 /  69.3,  0.93 /  69.3, -0.393 /  69.3, -3.21 /  69.3, -0.321 /  69.3]], dtype=np.double)#Приведенные тепловые эффекты
        
        #Вызываем функцию
        (invHeatCapacities,
         redHeatEffects) = InvHeatCapacityHeatEffectsOneFun(heatCapacities,#Теплоемкости
                                                            heatEffects#Тепловые эффекты
                                                            )
        
        #Проверяем значения
        errInvHeatCapacities = np.max(np.abs(invHeatCapacities - invHeatCapacitiesEt))
        self.assertAlmostEqual(errInvHeatCapacities, 0.0, 9)
        errRedHeatEffects = np.max(np.abs(redHeatEffects - redHeatEffectsEt))
        self.assertAlmostEqual(errRedHeatEffects, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
