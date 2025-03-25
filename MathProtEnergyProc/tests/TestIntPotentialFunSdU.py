import numpy as np

from MathProtEnergyProc.HeatPowerValues.Base import IntPotentialFunSdU

import unittest

#Модульные тесты
class TestIntPotentialFunSdU(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testIntPotentialFunSdU1(self):
        #Исходные данные
        intPotentialCubMatrix = np.array([[[1.11, -0.9],
                                           [-2.7, 0.03]], 
                                          [[ 0.81, 0.27],
                                           [-0.75, 0.33]]], dtype=np.double)
        jacSUPowerEnergies = np.array([0.021, 0.045], dtype=np.double)#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы
        jacSStateCoordinates = np.array([ 0.882, 0.675], dtype=np.double) * 0.021 + \
                               np.array([-0.315, 0.486], dtype=np.double) * 0.045
        
        #Эталонный результат
        intPotentialCompEt = np.array([-0.3, 1.5], dtype=np.double)#Компоненты потенциалов взаимодействия
        intPotentialEt = np.array([[ 0.882, 0.675],
                                   [-0.315, 0.486]], dtype=np.double)#Потенциалы взаимодействия энергетических степеней свободы
        
        #Вызываем функцию
        (intPotential,
         intPotentialComp) = IntPotentialFunSdU(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                                jacSUPowerEnergies,#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы
                                                intPotentialCubMatrix#Кубическая матрица потенциалов взаимодействия
                                                )
        
        #Проверяем значения
        err = np.max(np.abs(intPotentialComp - intPotentialCompEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialFunSdU2(self):
        #Исходные данные
        intPotentialCubMatrix = np.array([[[1.17, 0.69,  5.1],
                                           [-9.3, 0.75,  8.1],
                                           [9.33, 3.51, 11.7],
                                           [-6.3, 6.45, -8.7]], 
                                          [[-4.17,  3.69,   3.3],
                                           [-3.63,  2.25,   5.7],
                                           [ 9.63, -1.53,  15.3],
                                           [-4.53, -3.51, -11.7]], 
                                          [[7.17, 9.69,  8.1],
                                           [-6.3, 6.75, -2.1],
                                           [6.33, 6.51, 17.7],
                                           [-9.3, 3.45,  2.7]]], dtype=np.double)
        jacSUPowerEnergies = np.array([0.0243, 0.03183, 0.0753, 0.01563], dtype=np.double)#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы
        jacSStateCoordinates = np.array([ 47.1249, 18.2439,  16.245], dtype=np.double) * 0.0243 + \
                               np.array([ -9.8991, 14.3775, -31.635], dtype=np.double) * 0.03183 + \
                               np.array([-18.9567, 32.9913,  -3.321], dtype=np.double) * 0.0753 + \
                               np.array([-15.2721, 32.5863,  62.559], dtype=np.double) * 0.01563
        
        #Эталонный результат
        intPotentialCompEt = np.array([0.45, -4.83, 3.69], dtype=np.double)#Компоненты потенциалов взаимодействия
        intPotentialEt = np.array([[ 47.1249, 18.2439,  16.245],
                                   [ -9.8991, 14.3775, -31.635],
                                   [-18.9567, 32.9913,  -3.321],
                                   [-15.2721, 32.5863,  62.559]], dtype=np.double)
        
        #Вызываем функцию
        (intPotential,
         intPotentialComp) = IntPotentialFunSdU(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                                jacSUPowerEnergies,#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы
                                                intPotentialCubMatrix#Кубическая матрица потенциалов взаимодействия
                                                )
        
        #Проверяем значения
        err = np.max(np.abs(intPotentialComp - intPotentialCompEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialFunSdU3(self):
        #Исходные данные
        intPotentialCubMatrix = np.array([[[1.71,  -0.9, 2.1]], 
                                          [[3.21,  3.57, 0.3]], 
                                          [[6.21, -0.27, 9.3]]], dtype=np.double)
        jacSUPowerEnergies = np.array([0.1659], dtype=np.double)#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы
        jacSStateCoordinates = np.array([36.72, 14.6151, 36.216], dtype=np.double) * 0.1659
        
        #Эталонный результат
        intPotentialCompEt = np.array([-0.63, 4.23, 3.9], dtype=np.double)#Компоненты потенциалов взаимодействия
        intPotentialEt = np.array([[36.72, 14.6151, 36.216]], dtype=np.double)
        
        #Вызываем функцию
        (intPotential,
         intPotentialComp) = IntPotentialFunSdU(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                                jacSUPowerEnergies,#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы
                                                intPotentialCubMatrix#Кубическая матрица потенциалов взаимодействия
                                                )
        
        #Проверяем значения
        err = np.max(np.abs(intPotentialComp - intPotentialCompEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
