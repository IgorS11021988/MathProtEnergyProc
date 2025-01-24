import numpy as np

from MathProtEnergyProc.CorrectionModel import ReluFilter

import unittest

#Модульные тесты
class TestReluFilter(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testReluFilter1(self):
        #Исходные данные
        x = np.array([1.0, 2.0, -5.5, 7.7, -10.1], dtype=np.double)
        
        #Эталонный результат
        yEt = np.array([1.0, 2.0, 0.0, 7.7, 0.0], dtype=np.double)
        
        #Вызываем функцию
        y = ReluFilter(x)
        
        #Проверяем значения
        err = np.max(np.abs(y - yEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testReluFilter2(self):
        #Исходные данные
        x = np.array([3, -6, 0.15, -7.7, 18.3, 2.1, -2.7], dtype=np.double)
        
        #Эталонный результат
        yEt = np.array([3, 0.0, 0.15, 0.0, 18.3, 2.1, 0.0], dtype=np.double)
        
        #Вызываем функцию
        y = ReluFilter(x)
        
        #Проверяем значения
        err = np.max(np.abs(y - yEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testReluFilter3(self):
        #Исходные данные
        x = np.array([3.0, -2, 5.5, -8.1, 11.1], dtype=np.double)
        
        #Эталонный результат
        yEt = np.array([3.0, 0.0, 5.5, 0.0, 11.1], dtype=np.double)
        
        #Вызываем функцию
        y = ReluFilter(x)
        
        #Проверяем значения
        err = np.max(np.abs(y - yEt))
        self.assertAlmostEqual(err, 0.0, 5)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
