import numpy as np

from MathProtEnergyProc.CorrectionModel import KineticMatrixFromFacStreamEkvAff, SymRevMatrix

import unittest

#Модульные тесты
class TestKineticMatrixFromFacStreamEkvAff(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testKineticMatrixFromFacStreamEkvAff1(self):
        #Исходные данные
        irRevCom = np.array([[2.1, 0.75],
                             [0.63, 5.7]], dtype=np.double)
        kineticSubMatrix = np.array([[1.11,   0.3, 0.51],
                                     [-0.3,   2.1, 0.81],
                                     [0.21, -0.51,  3.3]], dtype=np.double)
        cfFacStream = np.array([[0.33, -4.5,  1.5],
                                [6.93,  7.5, -0.3]], dtype=np.double)
        cfEkvAffinities = np.array([[ 0.0, -2.7],
                                    [ 3.3,  5.1],
                                    [-8.1,  0.3]], dtype=np.double)
        
        #Эталонная кинетическая матрица
        symRevCom = SymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                                 cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                 kineticSubMatrix#Блок кинетической матрицы
                                 )
        crKineticSubMatrixStream = np.dot(cfFacStream, kineticSubMatrix)
        crKineticSubMatrixAffinities = np.dot(kineticSubMatrix, cfEkvAffinities)
        etIrRevCom = irRevCom + np.dot(np.dot(symRevCom, kineticSubMatrix),symRevCom.transpose())
        etKineticMatrix = np.block([[                  etIrRevCom, crKineticSubMatrixStream],
                                    [crKineticSubMatrixAffinities,         kineticSubMatrix]])
        
        #Вызываем функцию
        kineticMatrix = KineticMatrixFromFacStreamEkvAff(irRevCom,#Необратимая составляющая кинетической матрицы
                                                         cfFacStream,#Коэффициенты увлечения потоков
                                                         cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                                         kineticSubMatrix#Блок кинетической матрицы
                                                         )
        
        #Проверяем значения
        err = np.max(np.abs(kineticMatrix - etKineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMatrixFromFacStreamEkvAff2(self):
        #Исходные данные
        irRevCom = np.array([[  6.3, 2.31, 15.5],
                             [-3.51,  8.7, -3.9],
                             [-15.5,  3.9, 16.5]], dtype=np.double)
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
        
        #Эталонная кинетическая матрица
        symRevCom = SymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                                 cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                 kineticSubMatrix#Блок кинетической матрицы
                                 )
        crKineticSubMatrixStream = np.dot(cfFacStream, kineticSubMatrix)
        crKineticSubMatrixAffinities = np.dot(kineticSubMatrix, cfEkvAffinities)
        etIrRevCom = irRevCom + np.dot(np.dot(symRevCom, kineticSubMatrix),symRevCom.transpose())
        etKineticMatrix = np.block([[                  etIrRevCom, crKineticSubMatrixStream],
                                    [crKineticSubMatrixAffinities,         kineticSubMatrix]])
        
        #Вызываем функцию
        kineticMatrix = KineticMatrixFromFacStreamEkvAff(irRevCom,#Необратимая составляющая кинетической матрицы
                                                         cfFacStream,#Коэффициенты увлечения потоков
                                                         cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                                         kineticSubMatrix#Блок кинетической матрицы
                                                         )
        
        #Проверяем значения
        err = np.max(np.abs(kineticMatrix - etKineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMatrixFromFacStreamEkvAff3(self):
        #Исходные данные
        irRevCom = np.array([[9.9]], dtype=np.double)
        kineticSubMatrix = np.array([[ 7.11,   0.3],
                                     [-0.27,   8.1]], dtype=np.double)
        cfFacStream = np.array([[0.63, -1.5]], dtype=np.double)
        cfEkvAffinities = np.array([[0.1],
                                    [3.8]], dtype=np.double)
        
        #Эталонная кинетическая матрица
        symRevCom = SymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                                 cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                 kineticSubMatrix#Блок кинетической матрицы
                                 )
        crKineticSubMatrixStream = np.dot(cfFacStream, kineticSubMatrix)
        crKineticSubMatrixAffinities = np.dot(kineticSubMatrix, cfEkvAffinities)
        etIrRevCom = irRevCom + np.dot(np.dot(symRevCom, kineticSubMatrix),symRevCom.transpose())
        etKineticMatrix = np.block([[                  etIrRevCom, crKineticSubMatrixStream],
                                    [crKineticSubMatrixAffinities,         kineticSubMatrix]])
        
        #Вызываем функцию
        kineticMatrix = KineticMatrixFromFacStreamEkvAff(irRevCom,#Необратимая составляющая кинетической матрицы
                                                         cfFacStream,#Коэффициенты увлечения потоков
                                                         cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                                         kineticSubMatrix#Блок кинетической матрицы
                                                         )
        
        #Проверяем значения
        err = np.max(np.abs(kineticMatrix - etKineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMatrixFromFacStreamEkvAff4(self):
        #Исходные данные
        irRevCom = np.array([[9.9]], dtype=np.double)
        kineticSubMatrix = np.array([[10.11,   3.3],
                                     [-3.27,  11.1]], dtype=np.double)
        cfFacStream = np.array([[0.69, -1.23]], dtype=np.double)
        cfEkvAffinities = np.array([[0.0],
                                    [0.0]], dtype=np.double)
        
        #Эталонная кинетическая матрица
        symRevCom = SymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                                 cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                 kineticSubMatrix#Блок кинетической матрицы
                                 )
        crKineticSubMatrixStream = np.dot(cfFacStream, kineticSubMatrix)
        crKineticSubMatrixAffinities = np.dot(kineticSubMatrix, cfEkvAffinities)
        etIrRevCom = irRevCom + np.dot(np.dot(symRevCom, kineticSubMatrix),symRevCom.transpose())
        etKineticMatrix = np.block([[                  etIrRevCom, crKineticSubMatrixStream],
                                    [crKineticSubMatrixAffinities,         kineticSubMatrix]])
        
        #Вызываем функцию
        kineticMatrix = KineticMatrixFromFacStreamEkvAff(irRevCom,#Необратимая составляющая кинетической матрицы
                                                         cfFacStream,#Коэффициенты увлечения потоков
                                                         cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                                         kineticSubMatrix#Блок кинетической матрицы
                                                         )
        
        #Проверяем значения
        err = np.max(np.abs(kineticMatrix - etKineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMatrixFromFacStreamEkvAff5(self):
        #Исходные данные
        irRevCom = np.array([[  6.3, 2.31, 15.5],
                             [-3.51,  8.7, -3.9],
                             [-15.5,  3.9, 16.5]], dtype=np.double)
        kineticSubMatrix = np.array([[ 4.11,   5.7,  0.33,   0.0,  -0.9],
                                     [-5.43,   8.1, -0.69, -1.17,   0.3],
                                     [-0.27, -0.57,  15.3,  0.21,  1.23],
                                     [ 0.15,  0.45,  0.63,  9.63,  0.87],
                                     [ 0.09, -0.39, 0.123,  3.03, 11.67]], dtype=np.double)
        cfFacStream = np.array([[ 0.0,  0.0,  0.0,   0.0,   0.0],
                                [6.63, -7.5, -0.3, -2.13,  -0.9],
                                [3.69,  4.5, -3.9,  5.13, -0.63]], dtype=np.double)
        cfEkvAffinities = np.array([[ 0.0, 0.0,   2.1],
                                    [ 3.3, 0.0,  -3.9],
                                    [-8.1, 0.0,  0.63],
                                    [ 5.1, 0.0,  0.57],
                                    [ 2.1, 0.0, -0.45]], dtype=np.double)
        
        #Эталонная кинетическая матрица
        symRevCom = SymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                                 cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                 kineticSubMatrix#Блок кинетической матрицы
                                 )
        crKineticSubMatrixStream = np.dot(cfFacStream, kineticSubMatrix)
        crKineticSubMatrixAffinities = np.dot(kineticSubMatrix, cfEkvAffinities)
        etIrRevCom = irRevCom + np.dot(np.dot(symRevCom, kineticSubMatrix),symRevCom.transpose())
        etKineticMatrix = np.block([[                  etIrRevCom, crKineticSubMatrixStream],
                                    [crKineticSubMatrixAffinities,         kineticSubMatrix]])
        
        #Вызываем функцию
        kineticMatrix = KineticMatrixFromFacStreamEkvAff(irRevCom,#Необратимая составляющая кинетической матрицы
                                                         cfFacStream,#Коэффициенты увлечения потоков
                                                         cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                                         kineticSubMatrix#Блок кинетической матрицы
                                                         )
        
        #Проверяем значения
        err = np.max(np.abs(kineticMatrix - etKineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMatrixFromFacStreamEkvAff6(self):
        #Исходные данные
        irRevCom = np.array([[8.1, 0.75],
                             [0.63, 5.7]], dtype=np.double)
        kineticSubMatrix = np.array([[4.11,   0.3, 0.51],
                                     [-0.3,   8.1, 0.81],
                                     [0.21, -0.51,  3.3]], dtype=np.double)
        cfFacStream = np.array([[0.0, 0.0, 0.0],
                                [0.0, 0.0, 0.0]], dtype=np.double)
        cfEkvAffinities = np.array([[ 0.0, -2.7],
                                    [ 3.3,  5.1],
                                    [-8.1,  0.3]], dtype=np.double)
        
        #Эталонная кинетическая матрица
        symRevCom = SymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                                 cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                 kineticSubMatrix#Блок кинетической матрицы
                                 )
        crKineticSubMatrixStream = np.dot(cfFacStream, kineticSubMatrix)
        crKineticSubMatrixAffinities = np.dot(kineticSubMatrix, cfEkvAffinities)
        etIrRevCom = irRevCom + np.dot(np.dot(symRevCom, kineticSubMatrix),symRevCom.transpose())
        etKineticMatrix = np.block([[                  etIrRevCom, crKineticSubMatrixStream],
                                    [crKineticSubMatrixAffinities,         kineticSubMatrix]])
        
        #Вызываем функцию
        kineticMatrix = KineticMatrixFromFacStreamEkvAff(irRevCom,#Необратимая составляющая кинетической матрицы
                                                         cfFacStream,#Коэффициенты увлечения потоков
                                                         cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                                         kineticSubMatrix#Блок кинетической матрицы
                                                         )
        
        #Проверяем значения
        err = np.max(np.abs(kineticMatrix - etKineticMatrix))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
