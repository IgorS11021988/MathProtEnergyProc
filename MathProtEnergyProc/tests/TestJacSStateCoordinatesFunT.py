import numpy as np

from MathProtEnergyProc.HeatPowerValues.Base import JacSStateCoordinatesFunT

import unittest

#Модульные тесты
class TestJacSStateCoordinatesFunT(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testJacSStateCoordinatesFunT1(self):
        #Исходные данные
        tempratures = np.array([213.69, 15.9, 453.39], dtype=np.double)#Темперетауры энергетических степеней свободы
        intPotential = np.array([[ 21.3, 45.9,  81.3, -75.9],
                                 [-15.3, 19.5,  73.5,  68.1],
                                 [ 16.5, 20.7, -15.3, -90.3]], dtype=np.double)#Матрица потенциалов взаимодействия
        
        #Эталонный результат
        jacSStateCoordinatesEt = np.array([ 21.3, 45.9,  81.3, -75.9], dtype=np.double) / 213.69 + \
                                 np.array([-15.3, 19.5,  73.5,  68.1], dtype=np.double) / 15.9 + \
                                 np.array([ 16.5, 20.7, -15.3, -90.3], dtype=np.double) / 453.39
        
        #Вызываем функцию
        jacSStateCoordinates = JacSStateCoordinatesFunT(tempratures,#Темперетауры энергетических степеней свободы
                                                        intPotential,#Матрица потенциалов взаимодействия
                                                        )
        
        #Проверяем значения
        err = np.max(np.abs(jacSStateCoordinates - jacSStateCoordinatesEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testJacSStateCoordinatesFunT2(self):
        #Исходные данные
        tempratures = np.array([243.69, 318.9, 753.39, 156.3], dtype=np.double)#Темперетауры энергетических степеней свободы
        intPotential = np.array([[ 51.3, 15.9, 111.3, -75.9, 141.3, -72.9],
                                 [-45.3, 49.5, 103.5,  98.1, 133.5,  95.1],
                                 [ 76.5, 50.7, -45.3, -60.3, -48.3, -63.3],
                                 [106.5, 80.7, -75.9,  39.3, -78.9,  36.3]], dtype=np.double)#Матрица потенциалов взаимодействия
        
        #Эталонный результат
        jacSStateCoordinatesEt = np.array([ 51.3, 15.9, 111.3, -75.9, 141.3, -72.9], dtype=np.double) / 243.69 + \
                                 np.array([-45.3, 49.5, 103.5,  98.1, 133.5,  95.1], dtype=np.double) / 318.9 + \
                                 np.array([ 76.5, 50.7, -45.3, -60.3, -48.3, -63.3], dtype=np.double) / 753.39 + \
                                 np.array([106.5, 80.7, -75.9,  39.3, -78.9,  36.3], dtype=np.double) / 156.3
        
        #Вызываем функцию
        jacSStateCoordinates = JacSStateCoordinatesFunT(tempratures,#Темперетауры энергетических степеней свободы
                                                        intPotential,#Матрица потенциалов взаимодействия
                                                        )
        
        #Проверяем значения
        err = np.max(np.abs(jacSStateCoordinates - jacSStateCoordinatesEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testJacSStateCoordinatesFunT3(self):
        #Исходные данные
        tempratures = np.array([165.9], dtype=np.double)#Темперетауры энергетических степеней свободы
        intPotential = np.array([[321.3, 345.9, 681.3]], dtype=np.double)#Матрица потенциалов взаимодействия
        
        #Эталонный результат
        jacSStateCoordinatesEt = np.array([321.3, 345.9, 681.3], dtype=np.double) / 165.9
        
        #Вызываем функцию
        jacSStateCoordinates = JacSStateCoordinatesFunT(tempratures,#Темперетауры энергетических степеней свободы
                                                        intPotential,#Матрица потенциалов взаимодействия
                                                        )
        
        #Проверяем значения
        err = np.max(np.abs(jacSStateCoordinates - jacSStateCoordinatesEt))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
