import numpy as np

from MathProtEnergyProc.HeatPowerValues.Base import IntPotentialFunTOne

import unittest

#Модульные тесты
class TestIntPotentialFunTOne(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testIntPotentialFunTOne1(self):
        #Исходные данные
        jacSStateCoordinates = np.array([2.1, -3.9, 4.5], dtype=np.double)#Якобиан энтропии по координатам состояния
        tempraturesIntPotentials = [186.3, 315.3, 210.3]#Темперетауры по искомым потенциалам взаимодействия
        
        #Эталонный результат
        intPotentialEt = np.array([2.1*186.3, -3.9*315.3, 4.5*210.3], dtype=np.double)
        
        #Вызываем функцию
        intPotential = IntPotentialFunTOne(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                           tempraturesIntPotentials#Темперетауры по искомым потенциалам взаимодействия
                                           )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialFunTOne2(self):
        #Исходные данные
        jacSStateCoordinates = np.array([-5.7, 3.69, 7.23, -1.5, 8.123], dtype=np.double)#Якобиан энтропии по координатам состояния
        tempraturesIntPotentials = [186.3, 315.3, 210.3, 210.3, 315.3]#Темперетауры по искомым потенциалам взаимодействия
        
        #Эталонный результат
        intPotentialEt = np.array([-5.7*186.3, 3.69*315.3, 7.23*210.3, -1.5*210.3, 8.123*315.3], dtype=np.double)
        
        #Вызываем функцию
        intPotential = IntPotentialFunTOne(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                           tempraturesIntPotentials#Темперетауры по искомым потенциалам взаимодействия
                                           )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialFunTOne3(self):
        #Исходные данные
        jacSStateCoordinates = np.array([2.7, 3.3], dtype=np.double)#Якобиан энтропии по координатам состояния
        tempraturesIntPotentials = [219.3, 375.9]#Темперетауры по искомым потенциалам взаимодействия
        
        #Эталонный результат
        intPotentialEt = np.array([2.7*219.3, 3.3*375.9], dtype=np.double)
        
        #Вызываем функцию
        intPotential = IntPotentialFunTOne(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                           tempraturesIntPotentials#Темперетауры по искомым потенциалам взаимодействия
                                           )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
