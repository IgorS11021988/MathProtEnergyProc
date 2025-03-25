import numpy as np

from MathProtEnergyProc.HeatPowerValues.Base import IntPotentialFromCompFun

import unittest

#Модульные тесты
class TestIntPotentialFromCompFun(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testIntPotentialFromCompFun1(self):
        #Исходные данные
        intPotentialComp = np.array([-0.3, 1.5], dtype=np.double)#Компоненты потенциалов взаимодействия
        intPotentialCubMatrix = np.array([[[1.11, -0.9, 2.1],
                                           [-2.7, 0.03, 7.5]], 
                                          [[ 0.81, 0.27, -0.3],
                                           [-0.75, 0.33,  0.9]]], dtype=np.double)
        
        #Эталонный результат
        intPotentialEt = np.array([[ 0.882, 0.675, -1.08],
                                   [-0.315, 0.486,  -0.9]], dtype=np.double)
        
        #Вызываем функцию
        intPotential = IntPotentialFromCompFun(intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                                               intPotentialComp#Независимые составляющие потенциалов взаимодействия
                                               )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialFromCompFun2(self):
        #Исходные данные
        intPotentialComp = np.array([0.45, -4.83, 3.69], dtype=np.double)#Компоненты потенциалов взаимодействия
        intPotentialCubMatrix = np.array([[[1.17, 0.69,  5.1, -1.23, -0.81],
                                           [-9.3, 0.75,  8.1,  0.21, -7.29],
                                           [9.33, 3.51, 11.7,  9.57,  7.53],
                                           [-6.3, 6.45, -8.7,  0.51, -4.65]], 
                                          [[-4.17,  3.69,   3.3,  4.23, -0.21],
                                           [-3.63,  2.25,   5.7, -9.21, -7.11],
                                           [ 9.63, -1.53,  15.3,  6.51,  7.17],
                                           [-4.53, -3.51, -11.7, -3.57, -4.35]], 
                                          [[7.17, 9.69,  8.1, -7.23,  3.81],
                                           [-6.3, 6.75, -2.1,  3.21, -4.29],
                                           [6.33, 6.51, 17.7,  6.57, 10.83],
                                           [-9.3, 3.45,  2.7,  9.21, -1.65]]], dtype=np.double)
        
        #Эталонный результат
        intPotentialEt = np.array([[ 47.1249, 18.2439,  16.245, -47.6631, 14.7087],
                                   [ -9.8991, 14.3775, -31.635,  56.4237, 15.2307],
                                   [-18.9567, 32.9913,  -3.321,  -2.8935,  8.7201],
                                   [-15.2721, 32.5863,  62.559,  51.4575, 12.8295]], dtype=np.double)
        
        #Вызываем функцию
        intPotential = IntPotentialFromCompFun(intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                                               intPotentialComp#Независимые составляющие потенциалов взаимодействия
                                               )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialFromCompFun3(self):
        #Исходные данные
        intPotentialComp = np.array([-0.63, 4.23, 3.9], dtype=np.double)#Компоненты потенциалов взаимодействия
        intPotentialCubMatrix = np.array([[[1.71,  -0.9, 2.1, -8.7]], 
                                          [[3.21,  3.57, 0.3, -5.1]], 
                                          [[6.21, -0.27, 9.3, 11.1]]], dtype=np.double)
        
        #Эталонный результат
        intPotentialEt = np.array([[36.72, 14.6151, 36.216, 27.198]], dtype=np.double)
        
        #Вызываем функцию
        intPotential = IntPotentialFromCompFun(intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                                               intPotentialComp#Независимые составляющие потенциалов взаимодействия
                                               )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
