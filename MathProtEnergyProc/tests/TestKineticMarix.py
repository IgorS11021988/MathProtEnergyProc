import numpy as np

from MathProtEnergyProc.CorrectionModel import KineticMatrix

import unittest

#Модульные тесты
class TestKineticMarix(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testKineticMarix1(self):
        #Исходные данные
        kineticSubMatrix1 = np.array([[  9.3, 1.11,  0.3],
                                      [ 1.05, 14.7, 0.15],
                                      [-0.21, 1.17, 17.7]], dtype=np.double)
        kineticSubMatrix2 = np.array([[ 21.3, 1.47],
                                      [-0.99, 36.9]], dtype=np.double)
        kineticSubMatrixes = [kineticSubMatrix1,kineticSubMatrix2]
        arrKineticMatrixProcessCoordinatesNames = [["c1","c2","c5"], ["c7","c10"]]
        varKineticNames    = ["c1","c2","c1","c2","c7","c10", "c7","c10","c5","c1","c2","c5","c5"]
        varKineticAffNames = ["c1","c1","c2","c2","c7","c10","c10", "c7","c5","c5","c5","c1","c2"]
        
        #Эталонные значения
        etKineticMatrix = np.array([9.3,1.05,1.11,14.7,21.3,36.9,1.47,-0.99,17.7,0.3,0.15,-0.21,1.17], dtype=np.double)
        
        #Вызываем функцию
        fKineticMatrix = KineticMatrix(varKineticNames,#Имена сопряженностей между собой координат процессов
                                       varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                       
                                       arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                                       )#Функтор кинетической матрицы
        kineticMatrix = fKineticMatrix(kineticSubMatrixes)
        
        #Проверяем значения
        err = np.max(np.abs(etKineticMatrix - kineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMarix2(self):
        #Исходные данные
        kineticSubMatrix1 = np.array([[ 15.9, 1.11],
                                      [-1.05, 29.7]], dtype=np.double)
        kineticSubMatrix2 = np.array([[ 90.3, 7.11,  0.3],
                                      [ 1.05, 11.7, 9.45],
                                      [-0.21, 1.17, 17.7]], dtype=np.double)
        kineticSubMatrix3 = np.array([[ 21.3, 1.47],
                                      [-0.99, 36.9]], dtype=np.double)
        kineticSubMatrixes = [kineticSubMatrix1,kineticSubMatrix2,kineticSubMatrix3]
        arrKineticMatrixProcessCoordinatesNames = [["c1","c2"], ["c7","c10","c5"], ["a7","a3"]]
        varKineticNames    = ["c1","c2","c1","c2","c7","c10", "c7","c10","c5", "c5","c10","c7","c5","a7","a3","a3","a7"]
        varKineticAffNames = ["c1","c1","c2","c2","c7","c10","c10", "c7","c5","c10", "c5","c5","c7","a7","a3","a7","a3"]
        
        #Эталонные значения
        etKineticMatrix = np.array([15.9,-1.05,1.11,29.7,90.3,11.7,7.11,1.05,17.7,1.17,9.45,0.3,-0.21,21.3,36.9,-0.99,1.47], dtype=np.double)
        
        #Вызываем функцию
        fKineticMatrix = KineticMatrix(varKineticNames,#Имена сопряженностей между собой координат процессов
                                       varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                       
                                       arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                                       )#Функтор кинетической матрицы
        kineticMatrix = fKineticMatrix(kineticSubMatrixes)
        
        #Проверяем значения
        err = np.max(np.abs(etKineticMatrix - kineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMarix3(self):
        #Исходные данные
        kineticSubMatrix1 = np.array([[  9.3, 1.11],
                                      [ 1.05, 14.7]], dtype=np.double)
        kineticSubMatrix2 = np.array([[21.3]], dtype=np.double)
        kineticSubMatrixes = [kineticSubMatrix1,kineticSubMatrix2]
        arrKineticMatrixProcessCoordinatesNames = [["c1","c2"], ["c7"]]
        varKineticNames    = ["c1","c2","c1","c2","c7"]
        varKineticAffNames = ["c1","c1","c2","c2","c7"]
        
        #Эталонные значения
        etKineticMatrix = np.array([9.3,1.05,1.11,14.7,21.3], dtype=np.double)
        
        #Вызываем функцию
        fKineticMatrix = KineticMatrix(varKineticNames,#Имена сопряженностей между собой координат процессов
                                       varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                       
                                       arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                                       )#Функтор кинетической матрицы
        kineticMatrix = fKineticMatrix(kineticSubMatrixes)
        
        #Проверяем значения
        err = np.max(np.abs(etKineticMatrix - kineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMarix4(self):
        #Исходные данные
        kineticSubMatrix1 = np.array([[  9.3, 1.11],
                                      [ 1.05, 14.7]], dtype=np.double)
        kineticSubMatrix2 = np.array([21.3], dtype=np.double)
        kineticSubMatrixes = [kineticSubMatrix1,kineticSubMatrix2]
        arrKineticMatrixProcessCoordinatesNames = [["c1","c2"], ["c7"]]
        varKineticNames    = ["c1","c2","c1","c2","c7"]
        varKineticAffNames = ["c1","c1","c2","c2","c7"]
        
        #Эталонные значения
        etKineticMatrix = np.array([9.3,1.05,1.11,14.7,21.3], dtype=np.double)
        
        #Вызываем функцию
        fKineticMatrix = KineticMatrix(varKineticNames,#Имена сопряженностей между собой координат процессов
                                       varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                       
                                       arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                                       )#Функтор кинетической матрицы
        kineticMatrix = fKineticMatrix(kineticSubMatrixes)
        
        #Проверяем значения
        err = np.max(np.abs(etKineticMatrix - kineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMarix5(self):
        #Исходные данные
        kineticSubMatrix = np.array([[ 90.3, 7.11,  0.3],
                                     [ 1.05, 11.7, 9.45],
                                     [-0.21, 1.17, 17.7]], dtype=np.double)
        kineticSubMatrixes = [kineticSubMatrix]
        arrKineticMatrixProcessCoordinatesNames = ["c7","c10","c5"]
        varKineticNames    = ["c7","c10", "c7","c10","c5", "c5","c10","c7","c5"]
        varKineticAffNames = ["c7","c10","c10", "c7","c5","c10", "c5","c5","c7"]
        
        #Эталонные значения
        etKineticMatrix = np.array([90.3,11.7,7.11,1.05,17.7,1.17,9.45,0.3,-0.21], dtype=np.double)
        
        #Вызываем функцию
        fKineticMatrix = KineticMatrix(varKineticNames,#Имена сопряженностей между собой координат процессов
                                       varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                       
                                       arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                                       )#Функтор кинетической матрицы
        kineticMatrix = fKineticMatrix(kineticSubMatrixes)
        
        #Проверяем значения
        err = np.max(np.abs(etKineticMatrix - kineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMarix6(self):
        #Исходные данные
        kineticSubMatrix = np.array([[  9.3, 1.11],
                                     [ 1.05, 14.7]], dtype=np.double)
        kineticSubMatrixes = [kineticSubMatrix]
        arrKineticMatrixProcessCoordinatesNames = [["c1","c2"]]
        varKineticNames    = ["c1","c2","c1","c2"]
        varKineticAffNames = ["c1","c1","c2","c2"]
        
        #Эталонные значения
        etKineticMatrix = np.array([9.3,1.05,1.11,14.7], dtype=np.double)
        
        #Вызываем функцию
        fKineticMatrix = KineticMatrix(varKineticNames,#Имена сопряженностей между собой координат процессов
                                       varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                       
                                       arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                                       )#Функтор кинетической матрицы
        kineticMatrix = fKineticMatrix(kineticSubMatrixes)
        
        #Проверяем значения
        err = np.max(np.abs(etKineticMatrix - kineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
