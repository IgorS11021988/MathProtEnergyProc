import numpy as np

from MathProtEnergyProc.CorrectionModel import SymRevMatrix

import unittest

#Модульные тесты
class TestSymRevMatrix(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testSymRevMatrix1(self):
        #Исходные данные
        kineticSubMatrix = np.array([[1.11,   0.3, 0.51],
                                     [-0.3,   2.1, 0.81],
                                     [0.21, -0.51,  3.3]], dtype=np.double)
        cfFacStream = np.array([[0.33, -4.5,  1.5],
                                [6.93,  7.5, -0.3]], dtype=np.double)
        cfEkvAffinities = np.array([[ 0.0, -2.7],
                                    [ 3.3,  5.1],
                                    [-8.1,  0.3]], dtype=np.double)
        
        #Вызываем функцию
        symRevCom = SymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                                 cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                 kineticSubMatrix#Блок кинетической матрицы
                                 )
        
        #Проверяем значения
        symKineticSubMatrix = np.array([[2.22, 0.0, 0.72],
                                        [-0.0, 4.2,  0.3],
                                        [0.72, 0.3,  6.6]], dtype=np.double)
        crKineticSubMatrixStream = np.dot(cfFacStream, kineticSubMatrix)
        crKineticSubMatrixAffinities = np.dot(kineticSubMatrix, cfEkvAffinities)
        err = np.max(np.abs(np.dot(symRevCom,symKineticSubMatrix) - crKineticSubMatrixStream - crKineticSubMatrixAffinities.transpose()))
        self.assertAlmostEqual(err, 0.0, 9)
    def testSymRevMatrix2(self):
        #Исходные данные
        kineticSubMatrix = np.array([[ 4.11,   5.7,  0.33,   0.0,  -0.9],
                                     [-5.43,   8.1, -0.69, -1.17,   0.3],
                                     [-0.27, -0.57,   6.3,  0.21,  1.23],
                                     [ 0.15,  0.45,  0.33,  9.63,  0.87],
                                     [ 0.09, -0.39, 0.123,  3.03, 11.67]], dtype=np.double)
        cfFacStream = np.array([[0.33, -4.5,  1.5,  3.15,  0.03],
                                [6.93,  7.5, -0.3, -2.13,  -0.9],
                                [3.69,  4.5, -3.9,  5.13, -0.63]], dtype=np.double)
        cfEkvAffinities = np.array([[ 0.0, -2.7,   2.1],
                                    [ 3.3,  5.1,  -3.9],
                                    [-8.1,  0.3,  0.63],
                                    [ 5.1, 0.69,  0.57],
                                    [ 2.1,  0.9, -0.45]], dtype=np.double)
        
        #Вызываем функцию
        symRevCom = SymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                                 cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                 kineticSubMatrix#Блок кинетической матрицы
                                 )
        
        #Проверяем значения
        symKineticSubMatrix = np.array([[ 8.22,  0.27,  0.06,  0.15, -0.81],
                                        [ 0.27,  16.2, -1.26, -0.72, -0.09],
                                        [ 0.06, -1.26,  12.6,  0.54, 1.353],
                                        [ 0.15, -0.72,  0.54, 19.26,   3.9],
                                        [-0.81, -0.09, 1.353,   3.9, 23.34]], dtype=np.double)
        crKineticSubMatrixStream = np.dot(cfFacStream, kineticSubMatrix)
        crKineticSubMatrixAffinities = np.dot(kineticSubMatrix, cfEkvAffinities)
        err = np.max(np.abs(np.dot(symRevCom,symKineticSubMatrix) - crKineticSubMatrixStream - crKineticSubMatrixAffinities.transpose()))
        self.assertAlmostEqual(err, 0.0, 9)
    def testSymRevMatrix3(self):
        #Исходные данные
        kineticSubMatrix = np.array([[ 7.11,   0.3],
                                     [-0.27,   8.1]], dtype=np.double)
        cfFacStream = np.array([[0.63, -1.5]], dtype=np.double)
        cfEkvAffinities = np.array([[0.1],
                                    [3.8]], dtype=np.double)
        
        #Вызываем функцию
        symRevCom = SymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                                 cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                 kineticSubMatrix#Блок кинетической матрицы
                                 )
        
        #Проверяем значения
        symKineticSubMatrix = np.array([[14.22,  0.03],
                                        [ 0.03,  16.2]], dtype=np.double)
        crKineticSubMatrixStream = np.dot(cfFacStream, kineticSubMatrix)
        crKineticSubMatrixAffinities = np.dot(kineticSubMatrix, cfEkvAffinities)
        err = np.max(np.abs(np.dot(symRevCom,symKineticSubMatrix) - crKineticSubMatrixStream - crKineticSubMatrixAffinities.transpose()))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
