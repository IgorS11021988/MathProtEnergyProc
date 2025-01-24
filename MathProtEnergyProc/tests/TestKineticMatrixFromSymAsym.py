import numpy as np

from MathProtEnergyProc.CorrectionModel import KineticMatrixFromSymAsym

import unittest

#Модульные тесты
class TestKineticMatrixFromSymAsym(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testKineticMatrixFromSymAsym1(self):
        #Исходные данные
        posSubMatrixes = [2.1, 3.3, 6.3]#Составляющие положительно определенные матрицы
        subMatrixBalance1 = np.array([[1.1], [-2.3]], dtype=np.double)
        subMatrixBalance2 = np.array([[0.3], [-3.4]], dtype=np.double)
        subMatrixBalance3 = np.array([[-0.9], [2.5]], dtype=np.double)
        asymKineticMatrix = np.array([[0.0, 3.3],
                                      [0.0, 0.0]], dtype=np.double)
        
        #Составляющие атрицы баланса
        subMatrixesBalance = [subMatrixBalance1,subMatrixBalance2,subMatrixBalance3]
        
        #Эталонный результат
        etSubMatrixBalance = np.array([[ 1.1,  0.3, -0.9],
                                       [-2.3, -3.4,  2.5]], dtype=np.double)
        etPosSubMatrix = np.array([[2.1, 0.0, 0.0],
                                   [0.0, 3.3, 0.0],
                                   [0.0, 0.0, 6.3]], dtype = np.double)
        etAsymKineticMatrix = np.array([[ 0.0, 3.3],
                                        [-3.3, 0.0]], dtype=np.double)
        etKineticMatrix = np.dot(np.dot(etSubMatrixBalance,etPosSubMatrix),etSubMatrixBalance.transpose()) + etAsymKineticMatrix
        
        #Вызываем функцию
        KineticMatrix = KineticMatrixFromSymAsym(posSubMatrixes,#Положительные определенные составляющие атрицы
                                                 subMatrixesBalance,#Податрицы баланса
                                                 asymKineticMatrix#Антисимметричная составляющая
                                                 )
        
        #Проверяем значения
        err = np.max(np.abs(KineticMatrix - etKineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMatrixFromSymAsym2(self):
        #Исходные данные
        posSubMatrix1 = np.array([[2.5, -0.3],
                                  [0.3,  6.3]], dtype=np.double)
        posSubMatrix2 = np.array([5.1], dtype=np.double)
        posSubMatrix3 = np.array([[ 2.5, -0.3, -0.45],
                                  [ 0.3,  6.3,   1.5],
                                  [0.45, -1.5,  20.1]], dtype=np.double)
        subMatrixBalance1 = np.array([[1.1, 3.5], [-2.3, -9.1], [0.3, 0.9],  [4.53, 4.7]], dtype=np.double)
        subMatrixBalance2 = np.array([[0.3], [-3.4], [4.3],  [7.53]], dtype=np.double)
        subMatrixBalance3 = np.array([[-0.9, 0.1, 0.3], [2.5, -2.1, 1.7], [8.3, 9.8, 7.9], [-4.23, 0.23, 0.45]], dtype=np.double)
        asymKineticMatrix = np.array([[0.0, 3.6, -4.5,  2.1],
                                      [0.0, 0.0, 3.45, -5.7],
                                      [0.0, 6.9,  0.0,  6.3],
                                      [0.0, 0.0,  0.0,  0.0]], dtype=np.double)
        
        #Составляющие положительно определенные матрицы
        posSubMatrixes = [posSubMatrix1, posSubMatrix2, posSubMatrix3]
        
        #Составляющие атрицы баланса
        subMatrixesBalance = [subMatrixBalance1,subMatrixBalance2,subMatrixBalance3]
        
        #Эталонный результат
        etSubMatrixBalance = np.array([[ 1.1,  3.5,  0.3,  -0.9,   0.1,  0.3],
                                       [-2.3, -9.1, -3.4,   2.5,  -2.1,  1.7],
                                       [ 0.3,  0.9,  4.3,   8.3,   9.8,  7.9],
                                       [4.53,  4.7, 7.53,  -4.23, 0.23, 0.45]], dtype=np.double)
        etPosSubMatrix = np.array([[2.5, -0.3, 0.0,  0.0,  0.0,   0.0],
                                   [0.3,  6.3, 0.0,  0.0,  0.0,   0.0],
                                   [0.0, 0.0,  5.1,  0.0,  0.0,   0.0],
                                   [0.0, 0.0,  0.0,  2.5, -0.3, -0.45],
                                   [0.0, 0.0,  0.0,  0.3,  6.3,   1.5],
                                   [0.0, 0.0,  0.0, 0.45, -1.5,  20.1]], dtype = np.double)
        etAsymKineticMatrix = np.array([[ 0.0,  3.6,  -4.5,  2.1],
                                        [-3.6,  0.0, -3.45, -5.7],
                                        [ 4.5, 3.45,   0.0,  6.3],
                                        [-2.1,  5.7,  -6.3,  0.0]], dtype=np.double)
        etKineticMatrix = np.dot(np.dot(etSubMatrixBalance,etPosSubMatrix),etSubMatrixBalance.transpose()) + etAsymKineticMatrix
        
        #Вызываем функцию
        KineticMatrix = KineticMatrixFromSymAsym(posSubMatrixes,#Положительные определенные составляющие атрицы
                                                 subMatrixesBalance,#Податрицы баланса
                                                 asymKineticMatrix#Антисимметричная составляющая
                                                 )
        
        #Проверяем значения
        err = np.max(np.abs(KineticMatrix - etKineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMatrixFromSymAsym3(self):
        #Исходные данныеданные
        posSubMatrix1 = np.array([[2.5, 0.45],
                                  [0.3,  6.3]], dtype=np.double)
        posSubMatrix2 = np.array([8.1], dtype=np.double)
        posSubMatrix3 = np.array([[2.7, -0.45],
                                  [0.51,  3.9]], dtype=np.double)
        subMatrixBalance1 = np.array([[1.3, 3.1], [-2.4, -9.2], [0.45, 0.8],  [4.5, 4.6]], dtype=np.double)
        subMatrixBalance2 = np.array([[0.21], [-3.33], [4.1],  [7.51]], dtype=np.double)
        subMatrixBalance3 = np.array([[-0.9, -0.2], [2.2, -1.8], [-8.5, 9.7], [4.43, -0.33]], dtype=np.double)
        asymKineticMatrix = np.array([[3.0,  3.3, -4.5,  2.1],
                                      [2.1,  0.0, 3.45, -5.7],
                                      [2.7, -1.5,  2.7,  6.3],
                                      [0.0, -4.5,  2.1,  0.0]], dtype=np.double)
        
        #Составляющие положительно определенные матрицы
        posSubMatrixes = [posSubMatrix1, posSubMatrix2, posSubMatrix3]
        
        #Составляющие атрицы баланса
        subMatrixesBalance = [subMatrixBalance1,subMatrixBalance2,subMatrixBalance3]
        
        #Эталонный результат
        etSubMatrixBalance = np.array([[ 1.3,  3.1,  0.21,  -0.9,  -0.2],
                                       [-2.4, -9.2, -3.33,   2.2,  -1.8],
                                       [0.45,  0.8,   4.1,  -8.5,   9.7],
                                       [ 4.5,  4.6,  7.51,  4.43, -0.33]], dtype=np.double)
        etPosSubMatrix = np.array([[2.5, 0.45, 0.0,  0.0,   0.0],
                                   [0.3,  6.3, 0.0,  0.0,   0.0],
                                   [0.0, 0.0,  8.1,  0.0,   0.0],
                                   [0.0, 0.0,  0.0,  2.7, -0.45],
                                   [0.0, 0.0,  0.0, 0.51,   3.9]], dtype = np.double)
        etAsymKineticMatrix = np.array([[ 0.0,   1.2, -7.2,  2.1],
                                        [-1.2,   0.0, 4.95, -1.2],
                                        [ 7.2, -4.95,  0.0,  4.2],
                                        [-2.1,   1.2, -4.2,  0.0]], dtype=np.double)
        etKineticMatrix = np.dot(np.dot(etSubMatrixBalance,etPosSubMatrix),etSubMatrixBalance.transpose()) + etAsymKineticMatrix
        
        #Вызываем функцию
        KineticMatrix = KineticMatrixFromSymAsym(posSubMatrixes,#Положительные определенные составляющие атрицы
                                                 subMatrixesBalance,#Податрицы баланса
                                                 asymKineticMatrix#Антисимметричная составляющая
                                                 )
        
        #Проверяем значения
        err = np.max(np.abs(KineticMatrix - etKineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
