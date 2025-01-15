import numpy as np

from MathProtEnergyProc.CorrecionModel import PosLinearFilter

from MathProtEnergyProc.tests.UnitTestExamples.fFilters import fPosLinearFilter

import unittest

#Модульные тесты
class TestPosLinearFilter(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testPosLinearFilter1(self):
        #Исходные данные
        lam = 100
        x = np.array([1.0, 2.0, -5.5, 7.7, -10.1], dtype=np.double)
        
        #Эталонный результат
        yEt = np.array([1.0, 2.0, 0.0, 7.7, 0.0], dtype=np.double)
        
        #Вызываем функцию
        y = PosLinearFilter(x, lam=lam)
        
        #Проверяем значения
        err = np.max(np.abs(y - yEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testPosLinearFilter2(self):
        #Исходные данные
        lam = 0.03
        x = np.array([3, -6, 0.15, -7.7, 18.3, 2.1, -2.7], dtype=np.double)
        
        #Эталонный результат
        yEt = fPosLinearFilter(x, lam)
        
        #Вызываем функцию
        y = PosLinearFilter(x, lam=lam)
        
        #Проверяем значения
        err = np.max(np.abs(y - yEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testPosLinearFilter3(self):
        #Исходные данные
        lam = 10
        x = np.array([3.0, -2, 5.5, -8.1, 11.1], dtype=np.double)
        xmin = 3.3
        scy = 1.53
        
        #Эталонный результат
        yEt = xmin + scy*np.array([3, 0.0, 5.5, 0.0, 11.1], dtype=np.double)
        
        #Вызываем функцию
        y = PosLinearFilter(x, lam=lam, xmin=xmin, scy=scy)
        
        #Проверяем значения
        err = np.max(np.abs(y - yEt))
        self.assertAlmostEqual(err, 0.0, 5)
    def testPosLinearFilter4(self):
        #Исходные данные
        lam = 0.57
        x = np.array([4.5, -1.83, 5.5, -8.1, 11.1, 21.9, -36.9], dtype=np.double)
        xmin = 7.5
        scy = 4.83
        
        #Эталонный результат
        _yEt = np.hstack([fPosLinearFilter(np.array([4.5, -1.83, 5.5, -8.1, 11.1, 21.9], dtype=np.double), lam),
                                           np.array([0.0], dtype=np.double)])
        yEt = xmin + scy*_yEt
        
        #Вызываем функцию
        y = PosLinearFilter(x, lam=lam, xmin=xmin, scy=scy)
        
        #Проверяем значения
        err = np.max(np.abs(y - yEt))
        self.assertAlmostEqual(err, 0.0, 5)
    def testPosLinearFilter5(self):
        #Исходные данные
        lam = 0.111
        x = np.array([-1.5, 4.83, -0.57, -2.1, 17.1, 54.9, -72.9], dtype=np.double)
        scy = 8.1
        
        #Эталонный результат
        yEt = scy*fPosLinearFilter(np.array([-1.5, 4.83, -0.57, -2.1, 17.1, 54.9, -72.9], dtype=np.double), lam)
        
        #Вызываем функцию
        y = PosLinearFilter(x, lam=lam, scy=scy)
        
        #Проверяем значения
        err = np.max(np.abs(y - yEt))
        self.assertAlmostEqual(err, 0.0, 6)
    def testPosLinearFilter6(self):
        #Исходные данные
        lam = 3.27
        x = np.array([-7.5, 1.23, 4.5, -5.1, 14.1, 27.3, -81.9], dtype=np.double)
        xmin = 7.5
        
        #Эталонный результат
        _yEt = np.hstack([fPosLinearFilter(np.array([-7.5, 1.23, 4.5, -5.1], dtype=np.double), lam),
                                           np.array([14.1, 27.3, 0.0], dtype=np.double)])
        yEt = xmin + _yEt
        
        #Вызываем функцию
        y = PosLinearFilter(x, lam=lam, xmin=xmin)
        
        #Проверяем значения
        err = np.max(np.abs(y - yEt))
        self.assertAlmostEqual(err, 0.0, 6)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
