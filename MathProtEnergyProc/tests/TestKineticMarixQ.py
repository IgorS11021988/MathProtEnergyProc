import numpy as np

from MathProtEnergyProc.CorrectionModel import KineticMatrix, KineticMatrixQ

import unittest

#Модульные тесты
class TestKineticMarixQ(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testKineticMarixQ1(self):
        #Исходные данные
        kineticSubMatrix1 = np.array([[  9.3, 1.11,  0.3],
                                      [ 1.05, 14.7, 0.15],
                                      [-0.21, 1.17, 17.7]], dtype=np.double)
        kineticSubMatrix2 = np.array([[ 21.3, 1.47],
                                      [-0.99, 36.9]], dtype=np.double)
        kineticSubMatrix3 = np.array([[ 27.3, 1.11,  0.9],
                                      [ 1.05, 26.7, -9.15],
                                      [-0.21, 7.17, 38.7]], dtype=np.double)
        kineticSubMatrixes = [kineticSubMatrix1,kineticSubMatrix3,kineticSubMatrix2]
        arrKineticMatrixProcessCoordinatesNames = [["c1","c2","h1"], ["a1","h2","h3"], ["c7","c10"]]
        varKineticPCPCNames    = ["c1","c2","c1","c2","c7","c10", "c7","c10","a1"]
        varKineticPCPCAffNames = ["c1","c1","c2","c2","c7","c10","c10", "c7","a1"]
        varKineticPCHeatNames    = ["c1","c2","a1","a1"]
        varKineticPCHeatAffNames = ["h1","h1","h2","h3"]
        varKineticHeatPCNames    = ["h1","h1","h2","h3"]
        varKineticHeatPCAffNames = ["c1","c2","a1","a1"]
        varKineticHeatHeatNames    = ["h1","h2","h3","h2","h3"]
        varKineticHeatHeatAffNames = ["h1","h2","h3","h3","h2"]
        
        #Эталонные значения
        etKineticMatrixQPCPC = np.array([9.3,1.05,1.11,14.7,21.3,36.9,1.47,-0.99,27.3], dtype=np.double)
        etKineticMatrixQPCHeat = np.array([0.3,0.15,1.11,0.9], dtype=np.double)
        etKineticMatrixQHeatPC = np.array([-0.21,1.17,1.05,-0.21], dtype=np.double)
        etKineticMatrixQHeatHeat = np.array([17.7,26.7,38.7,-9.15,7.17], dtype=np.double)
        
        #Вызываем функцию
        fKineticMatrixQ = KineticMatrixQ(varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                         varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                         varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                         varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                         varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                         varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                         varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                         varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                         
                                         arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                                         )#Функтор кинетической матрицы
        (kineticMatrixQPCPC,
         kineticMatrixQPCHeat,
         kineticMatrixQHeatPC,
         kineticMatrixQHeatHeat
         ) = fKineticMatrixQ(kineticSubMatrixes)
        
        #Проверяем значения
        err = np.max(np.abs(etKineticMatrixQPCPC - kineticMatrixQPCPC))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(etKineticMatrixQPCHeat - kineticMatrixQPCHeat))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(etKineticMatrixQHeatPC - kineticMatrixQHeatPC))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(etKineticMatrixQHeatHeat - kineticMatrixQHeatHeat))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMarixQ2(self):
        #Исходные данные
        kineticSubMatrix1 = np.array([[ 15.9, 1.11],
                                      [-1.05, 29.7]], dtype=np.double)
        kineticSubMatrix2 = np.array([[ 90.3, 7.11,  0.3],
                                      [ 1.05, 11.7, 9.45],
                                      [-0.21, 1.17, 17.7]], dtype=np.double)
        kineticSubMatrix3 = np.array([[ 21.3, 1.47],
                                      [-0.99, 36.9]], dtype=np.double)
        kineticSubMatrixes = [kineticSubMatrix1,kineticSubMatrix2,kineticSubMatrix3]
        arrKineticMatrixProcessCoordinatesNames = [["h1","h2"], ["c7","c10","h3"], ["a7","a3"]]
        varKineticPCPCNames    = ["c7","c10", "c7","c10","a7","a3","a3","a7"]
        varKineticPCPCAffNames = ["c7","c10","c10", "c7","a7","a3","a7","a3"]
        varKineticPCHeatNames    = ["c10","c7"]
        varKineticPCHeatAffNames = [ "h3","h3"]
        varKineticHeatPCNames    = [ "h3","h3"]
        varKineticHeatPCAffNames = ["c10","c7"]
        varKineticHeatHeatNames    = ["h1","h2","h1","h2","h3"]
        varKineticHeatHeatAffNames = ["h1","h1","h2","h2","h3"]
        
        #Эталонные значения
        etKineticMatrixQPCPC = np.array([90.3,11.7,7.11,1.05,21.3,36.9,-0.99,1.47], dtype=np.double)
        etKineticMatrixQPCHeat = np.array([9.45,0.3], dtype=np.double)
        etKineticMatrixQHeatPC = np.array([1.17,-0.21], dtype=np.double)
        etKineticMatrixQHeatHeat = np.array([15.9,-1.05,1.11,29.7,17.7], dtype=np.double)
        
        #Вызываем функцию
        fKineticMatrixQ = KineticMatrixQ(varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                         varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                         varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                         varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                         varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                         varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                         varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                         varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                         
                                         arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                                         )#Функтор кинетической матрицы
        (kineticMatrixQPCPC,
         kineticMatrixQPCHeat,
         kineticMatrixQHeatPC,
         kineticMatrixQHeatHeat
         ) = fKineticMatrixQ(kineticSubMatrixes)
        
        #Проверяем значения
        err = np.max(np.abs(etKineticMatrixQPCPC - kineticMatrixQPCPC))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(etKineticMatrixQPCHeat - kineticMatrixQPCHeat))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(etKineticMatrixQHeatPC - kineticMatrixQHeatPC))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(etKineticMatrixQHeatHeat - kineticMatrixQHeatHeat))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMarixQ3(self):
        #Исходные данные
        kineticSubMatrix1 = np.array([[  9.3, 1.11],
                                      [ 1.05, 14.7]], dtype=np.double)
        kineticSubMatrix2 = np.array([[21.3]], dtype=np.double)
        kineticSubMatrixes = [kineticSubMatrix1,kineticSubMatrix2]
        arrKineticMatrixProcessCoordinatesNames = [["c1","c2"], ["h1"]]
        varKineticPCPCNames    = ["c1","c2","c1","c2"]
        varKineticPCPCAffNames = ["c1","c1","c2","c2"]
        varKineticPCHeatNames    = []
        varKineticPCHeatAffNames = []
        varKineticHeatPCNames    = []
        varKineticHeatPCAffNames = []
        varKineticHeatHeatNames    = ["h1"]
        varKineticHeatHeatAffNames = ["h1"]
        
        #Эталонные значения
        etKineticMatrixQPCPC = np.array([9.3,1.05,1.11,14.7], dtype=np.double)
        etKineticMatrixQHeatHeat = np.array([21.3], dtype=np.double)
        
        #Вызываем функцию
        fKineticMatrixQ = KineticMatrixQ(varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                         varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                         varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                         varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                         varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                         varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                         varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                         varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                         
                                         arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                                         )#Функтор кинетической матрицы
        (kineticMatrixQPCPC,
         kineticMatrixQPCHeat,
         kineticMatrixQHeatPC,
         kineticMatrixQHeatHeat
         ) = fKineticMatrixQ(kineticSubMatrixes)
        
        #Проверяем значения
        err = np.max(np.abs(etKineticMatrixQPCPC - kineticMatrixQPCPC))
        self.assertAlmostEqual(err, 0.0, 9)
        self.assertEqual(kineticMatrixQPCHeat.shape,(0,))
        self.assertEqual(kineticMatrixQHeatPC.shape,(0,))
        err = np.max(np.abs(etKineticMatrixQHeatHeat - kineticMatrixQHeatHeat))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMarixQ4(self):
        #Исходные данные
        kineticSubMatrix1 = np.array([[  9.3, 1.11],
                                      [ 1.05, 14.7]], dtype=np.double)
        kineticSubMatrix2 = np.array([21.3], dtype=np.double)
        kineticSubMatrixes = [kineticSubMatrix1,kineticSubMatrix2]
        arrKineticMatrixProcessCoordinatesNames = [["c1","c2"], ["h1"]]
        varKineticPCPCNames    = ["c1","c2","c1","c2"]
        varKineticPCPCAffNames = ["c1","c1","c2","c2"]
        varKineticPCHeatNames    = []
        varKineticPCHeatAffNames = []
        varKineticHeatPCNames    = []
        varKineticHeatPCAffNames = []
        varKineticHeatHeatNames    = ["h1"]
        varKineticHeatHeatAffNames = ["h1"]
        
        #Эталонные значения
        etKineticMatrixQPCPC = np.array([9.3,1.05,1.11,14.7], dtype=np.double)
        etKineticMatrixQHeatHeat = np.array([21.3], dtype=np.double)
        
        #Вызываем функцию
        fKineticMatrixQ = KineticMatrixQ(varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                         varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                         varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                         varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                         varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                         varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                         varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                         varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                         
                                         arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                                         )#Функтор кинетической матрицы
        (kineticMatrixQPCPC,
         kineticMatrixQPCHeat,
         kineticMatrixQHeatPC,
         kineticMatrixQHeatHeat
         ) = fKineticMatrixQ(kineticSubMatrixes)
        
        #Проверяем значения
        err = np.max(np.abs(etKineticMatrixQPCPC - kineticMatrixQPCPC))
        self.assertAlmostEqual(err, 0.0, 9)
        self.assertEqual(kineticMatrixQPCHeat.shape,(0,))
        self.assertEqual(kineticMatrixQHeatPC.shape,(0,))
        err = np.max(np.abs(etKineticMatrixQHeatHeat - kineticMatrixQHeatHeat))
        self.assertAlmostEqual(err, 0.0, 9)
    def testKineticMarixQ5(self):
        #Исходные данные
        kineticSubMatrix = np.array([[ 90.3, 7.11,  0.3],
                                     [ 1.05, 11.7, 9.45],
                                     [-0.21, 1.17, 17.7]], dtype=np.double)
        kineticSubMatrixes = [kineticSubMatrix]
        arrKineticMatrixProcessCoordinatesNames = ["c7","c10","c5"]
        varKineticPCPCNames    = ["c7","c10", "c7","c10","c5", "c5","c10","c7","c5"]
        varKineticPCPCAffNames = ["c7","c10","c10", "c7","c5","c10", "c5","c5","c7"]
        varKineticPCHeatNames    = []
        varKineticPCHeatAffNames = []
        varKineticHeatPCNames    = []
        varKineticHeatPCAffNames = []
        varKineticHeatHeatNames    = []
        varKineticHeatHeatAffNames = []
        
        #Эталонные значения
        etKineticMatrixQPCPC = np.array([90.3,11.7,7.11,1.05,17.7,1.17,9.45,0.3,-0.21], dtype=np.double)
        
        #Вызываем функциюфункцию
        fKineticMatrixQ = KineticMatrixQ(varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                         varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                         varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                         varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                         varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                         varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                         varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                         varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                         
                                         arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                                         )#Функтор кинетической матрицы
        (kineticMatrixQPCPC,
         kineticMatrixQPCHeat,
         kineticMatrixQHeatPC,
         kineticMatrixQHeatHeat
         ) = fKineticMatrixQ(kineticSubMatrixes)
        
        #Проверяем значения
        err = np.max(np.abs(etKineticMatrixQPCPC - kineticMatrixQPCPC))
        self.assertAlmostEqual(err, 0.0, 9)
        self.assertEqual(kineticMatrixQPCHeat.shape,(0,))
        self.assertEqual(kineticMatrixQHeatPC.shape,(0,))
        self.assertEqual(kineticMatrixQHeatHeat.shape,(0,))
    def testKineticMarixQ6(self):
        #Исходные данные
        kineticSubMatrix = np.array([[  9.3, 1.11],
                                     [ 1.05, 14.7]], dtype=np.double)
        kineticSubMatrixes = [kineticSubMatrix]
        arrKineticMatrixProcessCoordinatesNames = [["h1","h2"]]
        varKineticPCPCNames    = []
        varKineticPCPCAffNames = []
        varKineticPCHeatNames    = []
        varKineticPCHeatAffNames = []
        varKineticHeatPCNames    = []
        varKineticHeatPCAffNames = []
        varKineticHeatHeatNames    = ["h1","h2","h1","h2"]
        varKineticHeatHeatAffNames = ["h1","h1","h2","h2"]
        
        #Эталонные значения
        etKineticMatrixQHeatHeat = np.array([9.3,1.05,1.11,14.7], dtype=np.double)
        
        #Вызываем функцию
        fKineticMatrixQ = KineticMatrixQ(varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                         varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                         varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                         varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                         varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                         varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                         varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                         varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                         
                                         arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                                         )#Функтор кинетической матрицы
        (kineticMatrixQPCPC,
         kineticMatrixQPCHeat,
         kineticMatrixQHeatPC,
         kineticMatrixQHeatHeat
         ) = fKineticMatrixQ(kineticSubMatrixes)
        
        #Проверяем значения
        err = np.max(np.abs(etKineticMatrixQHeatHeat - kineticMatrixQHeatHeat))
        self.assertAlmostEqual(err, 0.0, 9)
        self.assertEqual(kineticMatrixQPCHeat.shape,(0,))
        self.assertEqual(kineticMatrixQHeatPC.shape,(0,))
        self.assertEqual(kineticMatrixQPCPC.shape,(0,))
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
