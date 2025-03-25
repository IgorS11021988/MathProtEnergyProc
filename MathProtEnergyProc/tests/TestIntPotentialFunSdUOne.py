import numpy as np

from MathProtEnergyProc.HeatPowerValues.Base import IntPotentialFunSdUOne

import unittest

#Модульные тесты
class TestIntPotentialFunSdUOne(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testIntPotentialFunSdUOne1(self):
        #Исходные данные
        jacSStateCoordinates = np.array([2.1, -3.9, 4.5], dtype=np.double)#Якобиан энтропии по координатам состояния
        jacSUPowerEnergiesIntPotentials = [0.0183, 0.0315, 0.0213]#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы, соответствующая искомым потенциалам взаимодействия
        
        #Эталонный результат
        intPotentialEt = np.array([2.1/0.0183, -3.9/0.0315, 4.5/0.0213], dtype=np.double)
        
        #Вызываем функцию
        intPotential = IntPotentialFunSdUOne(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                             jacSUPowerEnergiesIntPotentials#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы, соответствующая искомым потенциалам взаимодействия
                                             )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialFunSdUOne2(self):
        #Исходные данные
        jacSStateCoordinates = np.array([-5.7, 3.69, 7.23, -1.5, 8.123], dtype=np.double)#Якобиан энтропии по координатам состояния
        jacSUPowerEnergiesIntPotentials = [0.0186, 0.3153, 0.027, 0.027, 0.3153]#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы, соответствующая искомым потенциалам взаимодействия
        
        #Эталонный результат
        intPotentialEt = np.array([-5.7/0.0186, 3.69/0.3153, 7.23/0.027, -1.5/0.027, 8.123/0.3153], dtype=np.double)
        
        #Вызываем функцию
        intPotential = IntPotentialFunSdUOne(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                             jacSUPowerEnergiesIntPotentials#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы, соответствующая искомым потенциалам взаимодействия
                                             )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialFunSdUOne3(self):
        #Исходные данные
        jacSStateCoordinates = np.array([2.7, 3.3], dtype=np.double)#Якобиан энтропии по координатам состояния
        jacSUPowerEnergiesIntPotentials = [0.2193, 0.3759]#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы, соответствующая искомым потенциалам взаимодействия
        
        #Эталонный результат
        intPotentialEt = np.array([2.7/0.2193, 3.3/0.3759], dtype=np.double)
        
        #Вызываем функцию
        intPotential = IntPotentialFunSdUOne(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                             jacSUPowerEnergiesIntPotentials#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы, соответствующая искомым потенциалам взаимодействия
                                             )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
