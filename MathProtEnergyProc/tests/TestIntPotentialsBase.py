import numpy as np

from MathProtEnergyProc.HeatPowerValues import IntPotentialsBase

import unittest

#Модульные тесты
class TestIntPotentialsBase(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testIntPotentialsBase1(self):
        #Исходные данные
        stateCoordinatesNames = ["x1", "x2", "x5", "x7", "x10"]#Имена координат состояния
        energyPowersNames = ["e1", "e2", "e4"]#Имена энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1", "x7", "x5", "x5", "x2", "x10", "x7", "x1", "x10", "x5", "x7", "x10", "x2", "x1", "x2"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames     = ["e1", "e1", "e1", "e2", "e4",  "e4", "e4", "e4",  "e1", "e4", "e2",  "e2", "e2", "e2", "e1"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        intPotentialsMatrix = np.array([[ 1.11, -2.1, 3.69,  5.1, 7.23],
                                        [ 4.23, 5.73, 1.47, -2.1,  7.5],
                                        [-1.53,  5.7, 4.17,  3.9, 8.13]], dtype=np.double)#Матрица потенциалов взаимодействия
        
        #Эталонные значения
        intPotentialsEt = np.array([1.11, 5.1, 3.69, 1.47, 5.7, 8.13, 3.9, -1.53, 7.23, 4.17, -2.1, 7.5, 5.73, 4.23, -2.1], dtype=np.double)#Эталонные потенциалы взаимодействия
        
        #Создаем функтор потенциалов взаимодействия
        intPotentialsFunct = IntPotentialsBase(stateCoordinatesNames,#Имена координат состояния
                                               energyPowersNames,#Имена энергетических степеней свободы
                                               
                                               stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                               energyPowersVarPotentialsInterNames#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                               )
        
        #Вызываем функцию
        intPotentials = intPotentialsFunct(intPotentialsMatrix)#Потенциалы взаимодействия
        
        #Проверяем значения
        err = np.max(np.abs(intPotentials - intPotentialsEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialsBase2(self):
        #Исходные данные
        stateCoordinatesNames = ["x1", "x2", "x5", "x3", "x15", "x7", "x17"]#Имена координат состояния
        energyPowersNames = ["e1", "e8", "e2", "e4"]#Имена энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x2", "x15", "x17", "x5", "x1", "x7", "x7", "x3", "x5", "x3", "x15", "x7", "x17", "x15", "x15", "x17", "x17", "x5", "x7", "x1", "x2"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames     = ["e1",  "e8",  "e2", "e1", "e1", "e2", "e4", "e8", "e8", "e1",  "e1", "e1",  "e1",  "e2",  "e4",  "e4",  "e8", "e2", "e8", "e8", "e8"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        intPotentialsMatrix = np.array([[ 1.11, -3.9,  1.41, 3.69,  5.1,  7.23, -11.1],
                                        [ 4.23, 8.13, -7.53, 1.47, -4.5,  10.5,   4.5],
                                        [-1.53,  5.7, -4.83, 4.17,  3.9,  8.13,  33.9],
                                        [ 7.83,  8.7,  1.83, 1.77, -6.9, 11.13, -36.9]], dtype=np.double)#Матрица потенциалов взаимодействия
        
        #Эталонные значения
        intPotentialsEt = np.array([-3.9, -4.5, 33.9, 1.41, 1.11, 8.13, 11.13, 1.47, -7.53, 3.69, 5.1, 7.23, -11.1, 3.9, -6.9, -36.9, 4.5, -4.83, 10.5, 4.23, 8.13], dtype=np.double)#Эталонные потенциалы взаимодействия
        
        #Создаем функтор потенциалов взаимодействия
        intPotentialsFunct = IntPotentialsBase(stateCoordinatesNames,#Имена координат состояния
                                               energyPowersNames,#Имена энергетических степеней свободы
                                               
                                               stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                               energyPowersVarPotentialsInterNames#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                               )
        
        #Вызываем функцию
        intPotentials = intPotentialsFunct(intPotentialsMatrix)#Потенциалы взаимодействия
        
        #Проверяем значения
        err = np.max(np.abs(intPotentials - intPotentialsEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialsBase3(self):
        #Исходные данные
        stateCoordinatesNames = ["x1", "x5", "x3"]#Имена координат состояния
        energyPowersNames = ["e2", "e7"]#Имена энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1", "x5", "x3", "x5", "x1", "x3"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames     = ["e7", "e2", "e7", "e7", "e2", "e2"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        intPotentialsMatrix = np.array([[ 1.11, -2.1, 3.69],
                                        [ 4.23, 5.73, 1.47]], dtype=np.double)#Матрица потенциалов взаимодействия
        
        #Эталонные значения
        intPotentialsEt = np.array([4.23, -2.1, 1.47, 5.73, 1.11, 3.69], dtype=np.double)#Эталонные потенциалы взаимодействия
        
        #Создаем функтор потенциалов взаимодействия
        intPotentialsFunct = IntPotentialsBase(stateCoordinatesNames,#Имена координат состояния
                                               energyPowersNames,#Имена энергетических степеней свободы
                                               
                                               stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                               energyPowersVarPotentialsInterNames#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                               )
        
        #Вызываем функцию
        intPotentials = intPotentialsFunct(intPotentialsMatrix)#Потенциалы взаимодействия
        
        #Проверяем значения
        err = np.max(np.abs(intPotentials - intPotentialsEt))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
