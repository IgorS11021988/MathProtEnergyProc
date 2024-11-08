import numpy as np

from MathProtEnergyProc.IndexedNames import IndexedNamesFromIndexes

import unittest

#Модульные тесты
class TestIndexedNamesFromIndexes(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testIndexedNamesFromIndexes1(self):
        #Исходные данные
        Indexes = np.array([[0, 1,  2],
                            [0, 1,  5],
                            [0, 1,  8],
                            [0, 1, 11],
                            [0, 3,  2],
                            [0, 3,  5],
                            [0, 3,  8],
                            [0, 3, 11],
                            [0, 5,  2],
                            [0, 5,  5],
                            [0, 5,  8],
                            [0, 5, 11],
                            [1, 1,  2],
                            [1, 1,  5],
                            [1, 1,  8],
                            [1, 1, 11],
                            [1, 3,  2],
                            [1, 3,  5],
                            [1, 3,  8],
                            [1, 3, 11],
                            [1, 5,  2],
                            [1, 5,  5],
                            [1, 5,  8],
                            [1, 5, 11]])
        beginName = "v"#Начало имени
        
        #Индексированные имена
        IndexedNamesEt = ["v0_1_2",
                          "v0_1_5", 
                          "v0_1_8",
                          "v0_1_11",
                          "v0_3_2",
                          "v0_3_5",
                          "v0_3_8",
                          "v0_3_11",
                          "v0_5_2",
                          "v0_5_5",
                          "v0_5_8",
                          "v0_5_11",
                          "v1_1_2",
                          "v1_1_5",
                          "v1_1_8",
                          "v1_1_11",
                          "v1_3_2",
                          "v1_3_5",
                          "v1_3_8",
                          "v1_3_11",
                          "v1_5_2",
                          "v1_5_5",
                          "v1_5_8",
                          "v1_5_11"]
        
        #Генерируем индексированные имена
        IndexedNames = IndexedNamesFromIndexes(Indexes,#Индексы
                                               beginName#Начало имени
                                               )
        
        #Проверяем значения
        self.assertTrue(np.all(IndexedNames == IndexedNamesEt))
    def testIndexedNamesFromIndexes2(self):
        #Исходные данные
        Indexes = [[0, 1, 1,  3],
                   [0, 1, 1,  8],
                   [0, 1, 1, 13],
                   [0, 1, 1, 18],
                   [0, 1, 1, 23],
                   [0, 1, 2,  3],
                   [0, 1, 2,  8],
                   [0, 1, 2, 13],
                   [0, 1, 2, 18],
                   [0, 1, 2, 23],
                   [0, 1, 3,  3],
                   [0, 1, 3,  8],
                   [0, 1, 3, 13],
                   [0, 1, 3, 18],
                   [0, 1, 3, 23],
                   [0, 1, 4,  3],
                   [0, 1, 4,  8],
                   [0, 1, 4, 13],
                   [0, 1, 4, 18],
                   [0, 1, 4, 23],
                   [0, 1, 5,  3],
                   [0, 1, 5,  8],
                   [0, 1, 5, 13],
                   [0, 1, 5, 18],
                   [0, 1, 5, 23],
                   [0, 1, 6,  3],
                   [0, 1, 6,  8],
                   [0, 1, 6, 13],
                   [0, 1, 6, 18],
                   [0, 1, 6, 23],
                   [0, 4, 1,  3],
                   [0, 4, 1,  8],
                   [0, 4, 1, 13],
                   [0, 4, 1, 18],
                   [0, 4, 1, 23],
                   [0, 4, 2,  3],
                   [0, 4, 2,  8],
                   [0, 4, 2, 13],
                   [0, 4, 2, 18],
                   [0, 4, 2, 23],
                   [0, 4, 3,  3],
                   [0, 4, 3,  8],
                   [0, 4, 3, 13],
                   [0, 4, 3, 18],
                   [0, 4, 3, 23],
                   [0, 4, 4,  3],
                   [0, 4, 4,  8],
                   [0, 4, 4, 13],
                   [0, 4, 4, 18],
                   [0, 4, 4, 23],
                   [0, 4, 5,  3],
                   [0, 4, 5,  8],
                   [0, 4, 5, 13],
                   [0, 4, 5, 18],
                   [0, 4, 5, 23],
                   [0, 4, 6,  3],
                   [0, 4, 6,  8],
                   [0, 4, 6, 13],
                   [0, 4, 6, 18],
                   [0, 4, 6, 23],
                   [0, 7, 1,  3],
                   [0, 7, 1,  8],
                   [0, 7, 1, 13],
                   [0, 7, 1, 18],
                   [0, 7, 1, 23],
                   [0, 7, 2,  3],
                   [0, 7, 2,  8],
                   [0, 7, 2, 13],
                   [0, 7, 2, 18],
                   [0, 7, 2, 23],
                   [0, 7, 3,  3],
                   [0, 7, 3,  8],
                   [0, 7, 3, 13],
                   [0, 7, 3, 18],
                   [0, 7, 3, 23],
                   [0, 7, 4,  3],
                   [0, 7, 4,  8],
                   [0, 7, 4, 13],
                   [0, 7, 4, 18],
                   [0, 7, 4, 23],
                   [0, 7, 5,  3],
                   [0, 7, 5,  8],
                   [0, 7, 5, 13],
                   [0, 7, 5, 18],
                   [0, 7, 5, 23],
                   [0, 7, 6,  3],
                   [0, 7, 6,  8],
                   [0, 7, 6, 13],
                   [0, 7, 6, 18],
                   [0, 7, 6, 23],
                   [2, 1, 1,  3],
                   [2, 1, 1,  8],
                   [2, 1, 1, 13],
                   [2, 1, 1, 18],
                   [2, 1, 1, 23],
                   [2, 1, 2,  3],
                   [2, 1, 2,  8],
                   [2, 1, 2, 13],
                   [2, 1, 2, 18],
                   [2, 1, 2, 23],
                   [2, 1, 3,  3],
                   [2, 1, 3,  8],
                   [2, 1, 3, 13],
                   [2, 1, 3, 18],
                   [2, 1, 3, 23],
                   [2, 1, 4,  3],
                   [2, 1, 4,  8],
                   [2, 1, 4, 13],
                   [2, 1, 4, 18],
                   [2, 1, 4, 23],
                   [2, 1, 5,  3],
                   [2, 1, 5,  8],
                   [2, 1, 5, 13],
                   [2, 1, 5, 18],
                   [2, 1, 5, 23],
                   [2, 1, 6,  3],
                   [2, 1, 6,  8],
                   [2, 1, 6, 13],
                   [2, 1, 6, 18],
                   [2, 1, 6, 23],
                   [2, 4, 1,  3],
                   [2, 4, 1,  8],
                   [2, 4, 1, 13],
                   [2, 4, 1, 18],
                   [2, 4, 1, 23],
                   [2, 4, 2,  3],
                   [2, 4, 2,  8],
                   [2, 4, 2, 13],
                   [2, 4, 2, 18],
                   [2, 4, 2, 23],
                   [2, 4, 3,  3],
                   [2, 4, 3,  8],
                   [2, 4, 3, 13],
                   [2, 4, 3, 18],
                   [2, 4, 3, 23],
                   [2, 4, 4,  3],
                   [2, 4, 4,  8],
                   [2, 4, 4, 13],
                   [2, 4, 4, 18],
                   [2, 4, 4, 23],
                   [2, 4, 5,  3],
                   [2, 4, 5,  8],
                   [2, 4, 5, 13],
                   [2, 4, 5, 18],
                   [2, 4, 5, 23],
                   [2, 4, 6,  3],
                   [2, 4, 6,  8],
                   [2, 4, 6, 13],
                   [2, 4, 6, 18],
                   [2, 4, 6, 23],
                   [2, 7, 1,  3],
                   [2, 7, 1,  8],
                   [2, 7, 1, 13],
                   [2, 7, 1, 18],
                   [2, 7, 1, 23],
                   [2, 7, 2,  3],
                   [2, 7, 2,  8],
                   [2, 7, 2, 13],
                   [2, 7, 2, 18],
                   [2, 7, 2, 23],
                   [2, 7, 3,  3],
                   [2, 7, 3,  8],
                   [2, 7, 3, 13],
                   [2, 7, 3, 18],
                   [2, 7, 3, 23],
                   [2, 7, 4,  3],
                   [2, 7, 4,  8],
                   [2, 7, 4, 13],
                   [2, 7, 4, 18],
                   [2, 7, 4, 23],
                   [2, 7, 5,  3],
                   [2, 7, 5,  8],
                   [2, 7, 5, 13],
                   [2, 7, 5, 18],
                   [2, 7, 5, 23],
                   [2, 7, 6,  3],
                   [2, 7, 6,  8],
                   [2, 7, 6, 13],
                   [2, 7, 6, 18],
                   [2, 7, 6, 23],
                   [4, 1, 1,  3],
                   [4, 1, 1,  8],
                   [4, 1, 1, 13],
                   [4, 1, 1, 18],
                   [4, 1, 1, 23],
                   [4, 1, 2,  3],
                   [4, 1, 2,  8],
                   [4, 1, 2, 13],
                   [4, 1, 2, 18],
                   [4, 1, 2, 23],
                   [4, 1, 3,  3],
                   [4, 1, 3,  8],
                   [4, 1, 3, 13],
                   [4, 1, 3, 18],
                   [4, 1, 3, 23],
                   [4, 1, 4,  3],
                   [4, 1, 4,  8],
                   [4, 1, 4, 13],
                   [4, 1, 4, 18],
                   [4, 1, 4, 23],
                   [4, 1, 5,  3],
                   [4, 1, 5,  8],
                   [4, 1, 5, 13],
                   [4, 1, 5, 18],
                   [4, 1, 5, 23],
                   [4, 1, 6,  3],
                   [4, 1, 6,  8],
                   [4, 1, 6, 13],
                   [4, 1, 6, 18],
                   [4, 1, 6, 23],
                   [4, 4, 1,  3],
                   [4, 4, 1,  8],
                   [4, 4, 1, 13],
                   [4, 4, 1, 18],
                   [4, 4, 1, 23],
                   [4, 4, 2,  3],
                   [4, 4, 2,  8],
                   [4, 4, 2, 13],
                   [4, 4, 2, 18],
                   [4, 4, 2, 23],
                   [4, 4, 3,  3],
                   [4, 4, 3,  8],
                   [4, 4, 3, 13],
                   [4, 4, 3, 18],
                   [4, 4, 3, 23],
                   [4, 4, 4,  3],
                   [4, 4, 4,  8],
                   [4, 4, 4, 13],
                   [4, 4, 4, 18],
                   [4, 4, 4, 23],
                   [4, 4, 5,  3],
                   [4, 4, 5,  8],
                   [4, 4, 5, 13],
                   [4, 4, 5, 18],
                   [4, 4, 5, 23],
                   [4, 4, 6,  3],
                   [4, 4, 6,  8],
                   [4, 4, 6, 13],
                   [4, 4, 6, 18],
                   [4, 4, 6, 23],
                   [4, 7, 1,  3],
                   [4, 7, 1,  8],
                   [4, 7, 1, 13],
                   [4, 7, 1, 18],
                   [4, 7, 1, 23],
                   [4, 7, 2,  3],
                   [4, 7, 2,  8],
                   [4, 7, 2, 13],
                   [4, 7, 2, 18],
                   [4, 7, 2, 23],
                   [4, 7, 3,  3],
                   [4, 7, 3,  8],
                   [4, 7, 3, 13],
                   [4, 7, 3, 18],
                   [4, 7, 3, 23],
                   [4, 7, 4,  3],
                   [4, 7, 4,  8],
                   [4, 7, 4, 13],
                   [4, 7, 4, 18],
                   [4, 7, 4, 23],
                   [4, 7, 5,  3],
                   [4, 7, 5,  8],
                   [4, 7, 5, 13],
                   [4, 7, 5, 18],
                   [4, 7, 5, 23],
                   [4, 7, 6,  3],
                   [4, 7, 6,  8],
                   [4, 7, 6, 13],
                   [4, 7, 6, 18],
                   [4, 7, 6, 23]]
        beginName = "b"#Начало имени
        endName = "e"#Конец имени
        
        #Индексированные имена
        IndexedNamesEt = ["b0_1_1_3e",
                          "b0_1_1_8e",
                          "b0_1_1_13e",
                          "b0_1_1_18e",
                          "b0_1_1_23e",
                          "b0_1_2_3e",
                          "b0_1_2_8e",
                          "b0_1_2_13e",
                          "b0_1_2_18e",
                          "b0_1_2_23e",
                          "b0_1_3_3e",
                          "b0_1_3_8e",
                          "b0_1_3_13e",
                          "b0_1_3_18e",
                          "b0_1_3_23e",
                          "b0_1_4_3e",
                          "b0_1_4_8e",
                          "b0_1_4_13e",
                          "b0_1_4_18e",
                          "b0_1_4_23e",
                          "b0_1_5_3e",
                          "b0_1_5_8e",
                          "b0_1_5_13e",
                          "b0_1_5_18e",
                          "b0_1_5_23e",
                          "b0_1_6_3e",
                          "b0_1_6_8e",
                          "b0_1_6_13e",
                          "b0_1_6_18e",
                          "b0_1_6_23e",
                          "b0_4_1_3e",
                          "b0_4_1_8e",
                          "b0_4_1_13e",
                          "b0_4_1_18e",
                          "b0_4_1_23e",
                          "b0_4_2_3e",
                          "b0_4_2_8e",
                          "b0_4_2_13e",
                          "b0_4_2_18e",
                          "b0_4_2_23e",
                          "b0_4_3_3e",
                          "b0_4_3_8e",
                          "b0_4_3_13e",
                          "b0_4_3_18e",
                          "b0_4_3_23e",
                          "b0_4_4_3e",
                          "b0_4_4_8e",
                          "b0_4_4_13e",
                          "b0_4_4_18e",
                          "b0_4_4_23e",
                          "b0_4_5_3e",
                          "b0_4_5_8e",
                          "b0_4_5_13e",
                          "b0_4_5_18e",
                          "b0_4_5_23e",
                          "b0_4_6_3e",
                          "b0_4_6_8e",
                          "b0_4_6_13e",
                          "b0_4_6_18e",
                          "b0_4_6_23e",
                          "b0_7_1_3e",
                          "b0_7_1_8e",
                          "b0_7_1_13e",
                          "b0_7_1_18e",
                          "b0_7_1_23e",
                          "b0_7_2_3e",
                          "b0_7_2_8e",
                          "b0_7_2_13e",
                          "b0_7_2_18e",
                          "b0_7_2_23e",
                          "b0_7_3_3e",
                          "b0_7_3_8e",
                          "b0_7_3_13e",
                          "b0_7_3_18e",
                          "b0_7_3_23e",
                          "b0_7_4_3e",
                          "b0_7_4_8e",
                          "b0_7_4_13e",
                          "b0_7_4_18e",
                          "b0_7_4_23e",
                          "b0_7_5_3e",
                          "b0_7_5_8e",
                          "b0_7_5_13e",
                          "b0_7_5_18e",
                          "b0_7_5_23e",
                          "b0_7_6_3e",
                          "b0_7_6_8e",
                          "b0_7_6_13e",
                          "b0_7_6_18e",
                          "b0_7_6_23e",
                          "b2_1_1_3e",
                          "b2_1_1_8e",
                          "b2_1_1_13e",
                          "b2_1_1_18e",
                          "b2_1_1_23e",
                          "b2_1_2_3e",
                          "b2_1_2_8e",
                          "b2_1_2_13e",
                          "b2_1_2_18e",
                          "b2_1_2_23e",
                          "b2_1_3_3e",
                          "b2_1_3_8e",
                          "b2_1_3_13e",
                          "b2_1_3_18e",
                          "b2_1_3_23e",
                          "b2_1_4_3e",
                          "b2_1_4_8e",
                          "b2_1_4_13e",
                          "b2_1_4_18e",
                          "b2_1_4_23e",
                          "b2_1_5_3e",
                          "b2_1_5_8e",
                          "b2_1_5_13e",
                          "b2_1_5_18e",
                          "b2_1_5_23e",
                          "b2_1_6_3e",
                          "b2_1_6_8e",
                          "b2_1_6_13e",
                          "b2_1_6_18e",
                          "b2_1_6_23e",
                          "b2_4_1_3e",
                          "b2_4_1_8e",
                          "b2_4_1_13e",
                          "b2_4_1_18e",
                          "b2_4_1_23e",
                          "b2_4_2_3e",
                          "b2_4_2_8e",
                          "b2_4_2_13e",
                          "b2_4_2_18e",
                          "b2_4_2_23e",
                          "b2_4_3_3e",
                          "b2_4_3_8e",
                          "b2_4_3_13e",
                          "b2_4_3_18e",
                          "b2_4_3_23e",
                          "b2_4_4_3e",
                          "b2_4_4_8e",
                          "b2_4_4_13e",
                          "b2_4_4_18e",
                          "b2_4_4_23e",
                          "b2_4_5_3e",
                          "b2_4_5_8e",
                          "b2_4_5_13e",
                          "b2_4_5_18e",
                          "b2_4_5_23e",
                          "b2_4_6_3e",
                          "b2_4_6_8e",
                          "b2_4_6_13e",
                          "b2_4_6_18e",
                          "b2_4_6_23e",
                          "b2_7_1_3e",
                          "b2_7_1_8e",
                          "b2_7_1_13e",
                          "b2_7_1_18e",
                          "b2_7_1_23e",
                          "b2_7_2_3e",
                          "b2_7_2_8e",
                          "b2_7_2_13e",
                          "b2_7_2_18e",
                          "b2_7_2_23e",
                          "b2_7_3_3e",
                          "b2_7_3_8e",
                          "b2_7_3_13e",
                          "b2_7_3_18e",
                          "b2_7_3_23e",
                          "b2_7_4_3e",
                          "b2_7_4_8e",
                          "b2_7_4_13e",
                          "b2_7_4_18e",
                          "b2_7_4_23e",
                          "b2_7_5_3e",
                          "b2_7_5_8e",
                          "b2_7_5_13e",
                          "b2_7_5_18e",
                          "b2_7_5_23e",
                          "b2_7_6_3e",
                          "b2_7_6_8e",
                          "b2_7_6_13e",
                          "b2_7_6_18e",
                          "b2_7_6_23e",
                          "b4_1_1_3e",
                          "b4_1_1_8e",
                          "b4_1_1_13e",
                          "b4_1_1_18e",
                          "b4_1_1_23e",
                          "b4_1_2_3e",
                          "b4_1_2_8e",
                          "b4_1_2_13e",
                          "b4_1_2_18e",
                          "b4_1_2_23e",
                          "b4_1_3_3e",
                          "b4_1_3_8e",
                          "b4_1_3_13e",
                          "b4_1_3_18e",
                          "b4_1_3_23e",
                          "b4_1_4_3e",
                          "b4_1_4_8e",
                          "b4_1_4_13e",
                          "b4_1_4_18e",
                          "b4_1_4_23e",
                          "b4_1_5_3e",
                          "b4_1_5_8e",
                          "b4_1_5_13e",
                          "b4_1_5_18e",
                          "b4_1_5_23e",
                          "b4_1_6_3e",
                          "b4_1_6_8e",
                          "b4_1_6_13e",
                          "b4_1_6_18e",
                          "b4_1_6_23e",
                          "b4_4_1_3e",
                          "b4_4_1_8e",
                          "b4_4_1_13e",
                          "b4_4_1_18e",
                          "b4_4_1_23e",
                          "b4_4_2_3e",
                          "b4_4_2_8e",
                          "b4_4_2_13e",
                          "b4_4_2_18e",
                          "b4_4_2_23e",
                          "b4_4_3_3e",
                          "b4_4_3_8e",
                          "b4_4_3_13e",
                          "b4_4_3_18e",
                          "b4_4_3_23e",
                          "b4_4_4_3e",
                          "b4_4_4_8e",
                          "b4_4_4_13e",
                          "b4_4_4_18e",
                          "b4_4_4_23e",
                          "b4_4_5_3e",
                          "b4_4_5_8e",
                          "b4_4_5_13e",
                          "b4_4_5_18e",
                          "b4_4_5_23e",
                          "b4_4_6_3e",
                          "b4_4_6_8e",
                          "b4_4_6_13e",
                          "b4_4_6_18e",
                          "b4_4_6_23e",
                          "b4_7_1_3e",
                          "b4_7_1_8e",
                          "b4_7_1_13e",
                          "b4_7_1_18e",
                          "b4_7_1_23e",
                          "b4_7_2_3e",
                          "b4_7_2_8e",
                          "b4_7_2_13e",
                          "b4_7_2_18e",
                          "b4_7_2_23e",
                          "b4_7_3_3e",
                          "b4_7_3_8e",
                          "b4_7_3_13e",
                          "b4_7_3_18e",
                          "b4_7_3_23e",
                          "b4_7_4_3e",
                          "b4_7_4_8e",
                          "b4_7_4_13e",
                          "b4_7_4_18e",
                          "b4_7_4_23e",
                          "b4_7_5_3e",
                          "b4_7_5_8e",
                          "b4_7_5_13e",
                          "b4_7_5_18e",
                          "b4_7_5_23e",
                          "b4_7_6_3e",
                          "b4_7_6_8e",
                          "b4_7_6_13e",
                          "b4_7_6_18e",
                          "b4_7_6_23e"]
        
        #Генерируем индексированные имена
        IndexedNames = IndexedNamesFromIndexes(Indexes,#Индексы
                                               beginName,#Начало имени
                                               endName = endName#Конец имени
                                               )
        
        #Проверяем значения
        self.assertTrue(np.all(IndexedNames == IndexedNamesEt))
    def testIndexedNamesFromIndexes3(self):
        #Исходные данные
        Indexes = [[2, 1],
                   [2, 5],
                   [2, 9],
                   [4, 1],
                   [4, 5],
                   [4, 9],
                   [6, 1],
                   [6, 5],
                   [6, 9],
                   [8, 1],
                   [8, 5],
                   [8, 9]]
        beginName = "beg"#Начало имени
        endName = "end"#Конец имени
        sepName = "@"#Раздлитель имени
        
        #Индексированные имена
        IndexedNamesEt = ["beg2@1end",
                          "beg2@5end",
                          "beg2@9end",
                          "beg4@1end",
                          "beg4@5end",
                          "beg4@9end",
                          "beg6@1end",
                          "beg6@5end",
                          "beg6@9end",
                          "beg8@1end",
                          "beg8@5end",
                          "beg8@9end"]
        
        #Генерируем индексированные имена
        IndexedNames = IndexedNamesFromIndexes(Indexes,#Индексы
                                               beginName,#Начало имени
                                               endName = endName,#Конец имени
                                               sepName = sepName#Раздлитель имени
                                               )
        
        #Проверяем значения
        self.assertTrue(np.all(IndexedNames == IndexedNamesEt))
    def testIndexedNamesFromIndexes4(self):
        #Исходные данные
        Indexes = [[1],
                   [5],
                   [9],
                   [1],
                   [5],
                   [9],
                   [1],
                   [5],
                   [9],
                   [1],
                   [5],
                   [9]]
        beginName = "beg"#Начало имени
        endName = "end"#Конец имени
        sepName = "@"#Раздлитель имени
        
        #Индексированные имена
        IndexedNamesEt = ["beg1end",
                          "beg5end",
                          "beg9end",
                          "beg1end",
                          "beg5end",
                          "beg9end",
                          "beg1end",
                          "beg5end",
                          "beg9end",
                          "beg1end",
                          "beg5end",
                          "beg9end"]
        
        #Генерируем индексированные имена
        IndexedNames = IndexedNamesFromIndexes(Indexes,#Индексы
                                               beginName,#Начало имени
                                               endName = endName,#Конец имени
                                               sepName = sepName#Раздлитель имени
                                               )
        
        #Проверяем значения
        self.assertTrue(np.all(IndexedNames == IndexedNamesEt))
    def testIndexedNamesFromIndexes5(self):
        #Исходные данные
        Indexes = [2,
                   2,
                   2,
                   4,
                   4,
                   4,
                   6,
                   6,
                   6,
                   8,
                   8,
                   8]
        beginName = "beg"#Начало имени
        endName = "end"#Конец имени
        sepName = "@"#Раздлитель имени
        
        #Индексированные имена
        IndexedNamesEt = ["beg2end",
                          "beg2end",
                          "beg2end",
                          "beg4end",
                          "beg4end",
                          "beg4end",
                          "beg6end",
                          "beg6end",
                          "beg6end",
                          "beg8end",
                          "beg8end",
                          "beg8end"]
        
        #Генерируем индексированные имена
        IndexedNames = IndexedNamesFromIndexes(Indexes,#Индексы
                                               beginName,#Начало имени
                                               endName = endName,#Конец имени
                                               sepName = sepName#Раздлитель имени
                                               )
        
        #Проверяем значения
        self.assertTrue(np.all(IndexedNames == IndexedNamesEt))
    def testIndexedNamesFromIndexes6(self):
        #Исходные данные
        Indexes = [[1],
                   [5],
                   [9],
                   [1],
                   [5],
                   [9],
                   [1],
                   [5],
                   [9],
                   [1],
                   [5],
                   [9]]
        beginName = "beg"#Начало имени
        endName = "e"#Конец имени
        
        #Индексированные имена
        IndexedNamesEt = ["beg1e",
                          "beg5e",
                          "beg9e",
                          "beg1e",
                          "beg5e",
                          "beg9e",
                          "beg1e",
                          "beg5e",
                          "beg9e",
                          "beg1e",
                          "beg5e",
                          "beg9e"]
        
        #Генерируем индексированные имена
        IndexedNames = IndexedNamesFromIndexes(Indexes,#Индексы
                                               beginName,#Начало имени
                                               endName = endName#Конец имени
                                               )
        
        #Проверяем значения
        self.assertTrue(np.all(IndexedNames == IndexedNamesEt))
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()