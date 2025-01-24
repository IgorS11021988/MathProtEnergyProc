import numpy as np

from MathProtEnergyProc.CorrectionModel import BetaProcess

import unittest

#Модульные тесты
class TestBetaProcess(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testBetaProcess1(self):
        #Исходные данные
        cBetaProcessValues = np.array([1.11, 2.13, 3.3], dtype=np.double)
        
        #Эталонные значения
        etBetaProcessValues = np.array([1.11, 2.13, 3.3], dtype=np.double) / 6.54
        
        #Вызываем функцию
        betaProcessValues = BetaProcess(cBetaProcessValues)
        
        #Проверяем значения
        err = np.max(np.abs(betaProcessValues - etBetaProcessValues))
        self.assertAlmostEqual(err, 0.0, 9)
    def testBetaProcess2(self):
        #Исходные данные
        cBetaProcessValues = np.array([4.11, 5.49, 6.39, 8.1], dtype=np.double)
        
        #Эталонные значения
        etBetaProcessValues = np.array([4.11, 5.49, 6.39, 8.1], dtype=np.double) / 24.09
        
        #Вызываем функцию
        betaProcessValues = BetaProcess(cBetaProcessValues)
        
        #Проверяем значения
        err = np.max(np.abs(betaProcessValues - etBetaProcessValues))
        self.assertAlmostEqual(err, 0.0, 9)
    def testBetaProcess3(self):
        #Исходные данные
        cBetaProcessValues = np.array([11.23, 20.19], dtype=np.double)
        
        #Эталонные значения
        etBetaProcessValues = np.array([11.23, 20.19], dtype=np.double) / 31.42
        
        #Вызываем функцию
        betaProcessValues = BetaProcess(cBetaProcessValues)
        
        #Проверяем значения
        err = np.max(np.abs(betaProcessValues - etBetaProcessValues))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
