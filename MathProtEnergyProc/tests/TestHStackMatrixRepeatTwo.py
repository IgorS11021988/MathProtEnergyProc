import numpy as np

from MathProtEnergyProc.DatasAugmentation import HStackMatrixRepeatTwo

import unittest

#Модульные тесты
class TestHStackMatrixRepeatTwo(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testHStackMatrixRepeatTwo1(self):
        #Исходные данные
        matrValues1 = [[1.1, 2.2, 5.5, 7.7],
                       [6.3, 2.3, 4.5, 5.5],
                       [2.7, 2.1, 3.3, 9.3],
                       [1.1, 5.4, 2.8, 1.1],
                       [7.8, 6.9, 1.2, 2.3]]
        matrValues2 = [[4.1, 3.2, 6.5],
                       [8.3, 4.3, 8.5],
                       [6.7, 5.1, 6.3]]
        Indexes1 = [0, 2, 3, 5]
        Indexes2 = [1, 4]
        concatIndexesList = [Indexes1,Indexes2]
        
        #Преобразованные матрицы
        matrValuesRezEt = np.array([[1.1, 2.2, 5.5, 7.7, 4.1, 3.2, 6.5],
                                    [1.1, 2.2, 5.5, 7.7, 8.3, 4.3, 8.5],
                                    [1.1, 2.2, 5.5, 7.7, 6.7, 5.1, 6.3],
                                    [6.3, 2.3, 4.5, 5.5, 4.1, 3.2, 6.5],
                                    [6.3, 2.3, 4.5, 5.5, 8.3, 4.3, 8.5],
                                    [6.3, 2.3, 4.5, 5.5, 6.7, 5.1, 6.3],
                                    [2.7, 2.1, 3.3, 9.3, 4.1, 3.2, 6.5],
                                    [2.7, 2.1, 3.3, 9.3, 8.3, 4.3, 8.5],
                                    [2.7, 2.1, 3.3, 9.3, 6.7, 5.1, 6.3],
                                    [1.1, 5.4, 2.8, 1.1, 4.1, 3.2, 6.5],
                                    [1.1, 5.4, 2.8, 1.1, 8.3, 4.3, 8.5],
                                    [1.1, 5.4, 2.8, 1.1, 6.7, 5.1, 6.3],
                                    [7.8, 6.9, 1.2, 2.3, 4.1, 3.2, 6.5],
                                    [7.8, 6.9, 1.2, 2.3, 8.3, 4.3, 8.5],
                                    [7.8, 6.9, 1.2, 2.3, 6.7, 5.1, 6.3]])
        matrValuesRezEtTuple = (np.array([[1.1, 2.2, 5.5, 7.7],
                                          [1.1, 2.2, 5.5, 7.7],
                                          [1.1, 2.2, 5.5, 7.7],
                                          [6.3, 2.3, 4.5, 5.5],
                                          [6.3, 2.3, 4.5, 5.5],
                                          [6.3, 2.3, 4.5, 5.5],
                                          [2.7, 2.1, 3.3, 9.3],
                                          [2.7, 2.1, 3.3, 9.3],
                                          [2.7, 2.1, 3.3, 9.3],
                                          [1.1, 5.4, 2.8, 1.1],
                                          [1.1, 5.4, 2.8, 1.1],
                                          [1.1, 5.4, 2.8, 1.1],
                                          [7.8, 6.9, 1.2, 2.3],
                                          [7.8, 6.9, 1.2, 2.3],
                                          [7.8, 6.9, 1.2, 2.3]]),
                                np.array([[4.1, 3.2, 6.5],
                                          [8.3, 4.3, 8.5],
                                          [6.7, 5.1, 6.3],
                                          [4.1, 3.2, 6.5],
                                          [8.3, 4.3, 8.5],
                                          [6.7, 5.1, 6.3],
                                          [4.1, 3.2, 6.5],
                                          [8.3, 4.3, 8.5],
                                          [6.7, 5.1, 6.3],
                                          [4.1, 3.2, 6.5],
                                          [8.3, 4.3, 8.5],
                                          [6.7, 5.1, 6.3],
                                          [4.1, 3.2, 6.5],
                                          [8.3, 4.3, 8.5],
                                          [6.7, 5.1, 6.3]]))
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2#Величины 2
                                              )
        
        #Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez - matrValuesRezEt))
        self.assertEqual(dMatrValuesRez, 0.0)
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2,#Величины 2
                                              ConcatValues = True#Нужно ли конкатеновать величины
                                              )
        
        #Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez - matrValuesRezEt))
        self.assertEqual(dMatrValuesRez, 0.0)
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2,#Величины 2
                                              ConcatValues = False#Нужно ли конкатеновать величины
                                              )
        
        #Проверяем значения
        dMatrValuesRez0 = np.max(np.abs(matrValuesRez[0] - matrValuesRezEtTuple[0]))
        dMatrValuesRez1 = np.max(np.abs(matrValuesRez[1] - matrValuesRezEtTuple[1]))
        self.assertEqual(dMatrValuesRez0, 0.0)
        self.assertEqual(dMatrValuesRez1, 0.0)
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2,#Величины 2
                                              ConcatValues = True,#Нужно ли конкатеновать величины
                                              concatValuesIndexesList = concatIndexesList#Список индексов выборок
                                              )
        
        #Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez[0] - matrValuesRezEt[:,Indexes1]))
        self.assertEqual(dMatrValuesRez, 0.0)
        dMatrValuesRez = np.max(np.abs(matrValuesRez[1] - matrValuesRezEt[:,Indexes2]))
        self.assertEqual(dMatrValuesRez, 0.0)
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2,#Величины 2
                                              ConcatValues = False,#Нужно ли конкатеновать величины
                                              concatValuesIndexesList = concatIndexesList#Список индексов выборок
                                              )
        
        #Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez[0] - matrValuesRezEt[:,Indexes1]))
        self.assertEqual(dMatrValuesRez, 0.0)
        dMatrValuesRez = np.max(np.abs(matrValuesRez[1] - matrValuesRezEt[:,Indexes2]))
        self.assertEqual(dMatrValuesRez, 0.0)
    def testHStackMatrixRepeatTwo2(self):
        #Исходные данные
        matrValues1 = [[1.1, 1.2, 4.5],
                       [5.3, 1.3, 3.5],
                       [1.7, 1.1, 1.3],
                       [2.1, 3.4, 1.8],
                       [4.8, 2.9, 0.2],
                       [3.7, 2.5, 6.2]]
        matrValues2 = [[3.1, 3.2, 6.5, 4.5, 5.6],
                       [8.3, 2.3, 8.5, 1.1, 7.7]]
        Indexes1 = [0, 3, 6, 2]
        Indexes2 = [1, 5, 4]
        Indexes3 = [1, 7]
        concatIndexesList = [Indexes1,Indexes2,Indexes3]
        
        #Преобразованные матрицы
        matrValuesRezEt = np.array([[1.1, 1.2, 4.5, 3.1, 3.2, 6.5, 4.5, 5.6],
                                    [1.1, 1.2, 4.5, 8.3, 2.3, 8.5, 1.1, 7.7],
                                    [5.3, 1.3, 3.5, 3.1, 3.2, 6.5, 4.5, 5.6],
                                    [5.3, 1.3, 3.5, 8.3, 2.3, 8.5, 1.1, 7.7],
                                    [1.7, 1.1, 1.3, 3.1, 3.2, 6.5, 4.5, 5.6],
                                    [1.7, 1.1, 1.3, 8.3, 2.3, 8.5, 1.1, 7.7],
                                    [2.1, 3.4, 1.8, 3.1, 3.2, 6.5, 4.5, 5.6],
                                    [2.1, 3.4, 1.8, 8.3, 2.3, 8.5, 1.1, 7.7],
                                    [4.8, 2.9, 0.2, 3.1, 3.2, 6.5, 4.5, 5.6],
                                    [4.8, 2.9, 0.2, 8.3, 2.3, 8.5, 1.1, 7.7],
                                    [3.7, 2.5, 6.2, 3.1, 3.2, 6.5, 4.5, 5.6],
                                    [3.7, 2.5, 6.2, 8.3, 2.3, 8.5, 1.1, 7.7]])
        matrValuesRezEtTuple = (np.array([[1.1, 1.2, 4.5],
                                          [1.1, 1.2, 4.5],
                                          [5.3, 1.3, 3.5],
                                          [5.3, 1.3, 3.5],
                                          [1.7, 1.1, 1.3],
                                          [1.7, 1.1, 1.3],
                                          [2.1, 3.4, 1.8],
                                          [2.1, 3.4, 1.8],
                                          [4.8, 2.9, 0.2],
                                          [4.8, 2.9, 0.2],
                                          [3.7, 2.5, 6.2],
                                          [3.7, 2.5, 6.2]]),
                                np.array([[3.1, 3.2, 6.5, 4.5, 5.6],
                                          [8.3, 2.3, 8.5, 1.1, 7.7],
                                          [3.1, 3.2, 6.5, 4.5, 5.6],
                                          [8.3, 2.3, 8.5, 1.1, 7.7],
                                          [3.1, 3.2, 6.5, 4.5, 5.6],
                                          [8.3, 2.3, 8.5, 1.1, 7.7],
                                          [3.1, 3.2, 6.5, 4.5, 5.6],
                                          [8.3, 2.3, 8.5, 1.1, 7.7],
                                          [3.1, 3.2, 6.5, 4.5, 5.6],
                                          [8.3, 2.3, 8.5, 1.1, 7.7],
                                          [3.1, 3.2, 6.5, 4.5, 5.6],
                                          [8.3, 2.3, 8.5, 1.1, 7.7]]))
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2#Величины 2
                                              )
        
        #Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez - matrValuesRezEt))
        self.assertEqual(dMatrValuesRez, 0.0)
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2,#Величины 2
                                              ConcatValues = True#Нужно ли конкатеновать величины
                                              )
        
        #Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez - matrValuesRezEt))
        self.assertEqual(dMatrValuesRez, 0.0)
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2,#Величины 2
                                              ConcatValues = False#Нужно ли конкатеновать величины
                                              )
        
        #Проверяем значения
        dMatrValuesRez0 = np.max(np.abs(matrValuesRez[0] - matrValuesRezEtTuple[0]))
        dMatrValuesRez1 = np.max(np.abs(matrValuesRez[1] - matrValuesRezEtTuple[1]))
        self.assertEqual(dMatrValuesRez0, 0.0)
        self.assertEqual(dMatrValuesRez1, 0.0)
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2,#Величины 2
                                              ConcatValues = True,#Нужно ли конкатеновать величины
                                              concatValuesIndexesList = concatIndexesList#Список индексов выборок
                                              )
        
        #Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez[0] - matrValuesRezEt[:,Indexes1]))
        self.assertEqual(dMatrValuesRez, 0.0)
        dMatrValuesRez = np.max(np.abs(matrValuesRez[1] - matrValuesRezEt[:,Indexes2]))
        self.assertEqual(dMatrValuesRez, 0.0)
        dMatrValuesRez = np.max(np.abs(matrValuesRez[2] - matrValuesRezEt[:,Indexes3]))
        self.assertEqual(dMatrValuesRez, 0.0)
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2,#Величины 2
                                              ConcatValues = False,#Нужно ли конкатеновать величины
                                              concatValuesIndexesList = concatIndexesList#Список индексов выборок
                                              )
        
        #Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez[0] - matrValuesRezEt[:,Indexes1]))
        self.assertEqual(dMatrValuesRez, 0.0)
        dMatrValuesRez = np.max(np.abs(matrValuesRez[1] - matrValuesRezEt[:,Indexes2]))
        self.assertEqual(dMatrValuesRez, 0.0)
        dMatrValuesRez = np.max(np.abs(matrValuesRez[2] - matrValuesRezEt[:,Indexes3]))
        self.assertEqual(dMatrValuesRez, 0.0)
    def testHStackMatrixRepeatTwo3(self):
        #Исходные данные
        matrValues1 = [[3.1, 2.2, 5.5, 7.7, 10.1],
                       [5.1, 2.3, 4.5, 4.5, 11.3],
                       [3.5, 2.1, 3.3, 6.3, 13.5]]
        matrValues2 = [[4.7, 3.2],
                       [9.3, 4.3],
                       [7.7, 5.1],
                       [4.5, 3.4]]
        Indexes1 = [0, 3, 6, 2]
        concatIndexesList = [Indexes1]
        
        #Преобразованные матрицы
        matrValuesRezEt = np.array([[3.1, 2.2, 5.5, 7.7, 10.1, 4.7, 3.2],
                                    [3.1, 2.2, 5.5, 7.7, 10.1, 9.3, 4.3],
                                    [3.1, 2.2, 5.5, 7.7, 10.1, 7.7, 5.1],
                                    [3.1, 2.2, 5.5, 7.7, 10.1, 4.5, 3.4],
                                    [5.1, 2.3, 4.5, 4.5, 11.3, 4.7, 3.2],
                                    [5.1, 2.3, 4.5, 4.5, 11.3, 9.3, 4.3],
                                    [5.1, 2.3, 4.5, 4.5, 11.3, 7.7, 5.1],
                                    [5.1, 2.3, 4.5, 4.5, 11.3, 4.5, 3.4],
                                    [3.5, 2.1, 3.3, 6.3, 13.5, 4.7, 3.2],
                                    [3.5, 2.1, 3.3, 6.3, 13.5, 9.3, 4.3],
                                    [3.5, 2.1, 3.3, 6.3, 13.5, 7.7, 5.1],
                                    [3.5, 2.1, 3.3, 6.3, 13.5, 4.5, 3.4]])
        matrValuesRezEtTuple = (np.array([[3.1, 2.2, 5.5, 7.7, 10.1],
                                          [3.1, 2.2, 5.5, 7.7, 10.1],
                                          [3.1, 2.2, 5.5, 7.7, 10.1],
                                          [3.1, 2.2, 5.5, 7.7, 10.1],
                                          [5.1, 2.3, 4.5, 4.5, 11.3],
                                          [5.1, 2.3, 4.5, 4.5, 11.3],
                                          [5.1, 2.3, 4.5, 4.5, 11.3],
                                          [5.1, 2.3, 4.5, 4.5, 11.3],
                                          [3.5, 2.1, 3.3, 6.3, 13.5],
                                          [3.5, 2.1, 3.3, 6.3, 13.5],
                                          [3.5, 2.1, 3.3, 6.3, 13.5],
                                          [3.5, 2.1, 3.3, 6.3, 13.5]]),
                                np.array([[4.7, 3.2],
                                          [9.3, 4.3],
                                          [7.7, 5.1],
                                          [4.5, 3.4],
                                          [4.7, 3.2],
                                          [9.3, 4.3],
                                          [7.7, 5.1],
                                          [4.5, 3.4],
                                          [4.7, 3.2],
                                          [9.3, 4.3],
                                          [7.7, 5.1],
                                          [4.5, 3.4]]))
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2#Величины 2
                                              )
        
        #Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez - matrValuesRezEt))
        self.assertEqual(dMatrValuesRez, 0.0)
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2,#Величины 2
                                              ConcatValues = True#Нужно ли конкатеновать величины
                                              )
        
        #Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez - matrValuesRezEt))
        self.assertEqual(dMatrValuesRez, 0.0)
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2,#Величины 2
                                              ConcatValues = False#Нужно ли конкатеновать величины
                                              )
        
        #Проверяем значения
        dMatrValuesRez0 = np.max(np.abs(matrValuesRez[0] - matrValuesRezEtTuple[0]))
        dMatrValuesRez1 = np.max(np.abs(matrValuesRez[1] - matrValuesRezEtTuple[1]))
        self.assertEqual(dMatrValuesRez0, 0.0)
        self.assertEqual(dMatrValuesRez1, 0.0)
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2,#Величины 2
                                              ConcatValues = True,#Нужно ли конкатеновать величины
                                              concatValuesIndexesList = concatIndexesList#Список индексов выборок
                                              )
        
        #Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez[0] - matrValuesRezEt[:,Indexes1]))
        self.assertEqual(dMatrValuesRez, 0.0)
        
        #Генерируем случайные данные
        matrValuesRez = HStackMatrixRepeatTwo(matrValues1,#Величины 1
                                              matrValues2,#Величины 2
                                              ConcatValues = False,#Нужно ли конкатеновать величины
                                              concatValuesIndexesList = concatIndexesList#Список индексов выборок
                                              )
        
        #Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez[0] - matrValuesRezEt[:,Indexes1]))
        self.assertEqual(dMatrValuesRez, 0.0)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
