import numpy as np

from MathProtEnergyProc import NonEqSystemBase

from MathProtEnergyProc.tests.UnitTestExamples.TestNonEq1 import *
from MathProtEnergyProc.tests.UnitTestExamples.TestNonEq2 import *
from MathProtEnergyProc.tests.UnitTestExamples.TestNonEq3 import *

import unittest

#Модульные тесты
class TestNonEqSystemBase(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testNonEqBase1(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81])#Координаты состояния
        systemParameters = [  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                              35,   55,   75,  81,   210,  300,  150,  153,  123, 15.3,
                            45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,  1.5, 33.9,
                            27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,  159, 21.9,
                            43.5, 7.53, 37.5, 1.5, 120.3]#Параметры системы
        
        #Стехиометрические коэффициенты химических реакций
        nu1_1 = 3
        nu1_2 = 2
        nu1_3 = 5
        nu1_4 = 1
        nu1_5 = 2
        nu1_6 = 3
        nu1_7 = 4
        nu1_8 = 1
        nu1_9 = 5
        nu1_10 = 7
        nu1_11 = 2.2
        nu2_1 = 3
        nu2_2 = 5
        nu2_3 = 4
        nu2_4 = 1
        nu2_5 = 2
        nu2_6 = 3
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3
        
        #Эталонный результат
        (etChemPot,etAff,etKineticMatrix,
         etVProcesses,etBalanceMatrix,etVx) = CountSystem1(stateCoordinates,systemParameters,
                                                           nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                           nu1_8,nu1_9,nu1_10,nu1_11,nu2_1,nu2_2,nu2_3,
                                                           nu2_4,nu2_5,nu2_6,Adiff1_1,Adiff1_2,Adiff2_1)
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff1","vDiff2"]#Имена координат процессов
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена переменных потенциалов взаимодействия по координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой термодинамических сил
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []#Имена переменных внешних потоков
        
        #Задаем систему
        nonEqSystem = NonEqSystemBase(stateCoordinatesNames,#Имена координат состояния
                                      processCoordinatesNames,#Имена координат процессов
                                      stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                      
                                      #Задаем рассчитываемые параметры системы
                                      stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                      processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                      stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                        
                                      #Задаем переменные параметры кинетической матрицы
                                      varKineticNames,#Имена сопряженностей между собой координат процессов
                                      varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                        
                                      #Задаем внешние потоки
                                      stateCoordinatesVarStreamsNames#Имена переменных внешних потоков
                                      )
        
        #Задаем постоянные параметры системы
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_1", -nu1_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_1", -nu1_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4", "vChem1_1",  nu1_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_2", -nu1_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_2", -nu1_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_6", "vChem1_2",  nu1_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_2",  nu1_7)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_3", -nu1_8)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_3", -nu1_9)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_3", -nu1_10)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_3",  nu1_11)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vChem2_1", -nu2_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7", "vChem2_1", -nu2_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_8", "vChem2_1",  nu2_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4", "vChem2_2", -nu2_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vChem2_2", -nu2_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7", "vChem2_2",  nu2_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vDiff1", -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4", "vDiff2", -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vDiff1",  1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4", "vDiff2",  1.0)
        nonEqSystem.SetKineticMatrixConstElement("vDiff1", "vDiff1", Adiff1_1)
        nonEqSystem.SetKineticMatrixConstElement("vDiff1", "vDiff2", Adiff1_2)
        nonEqSystem.SetKineticMatrixConstElement("vDiff2", "vDiff1", Adiff2_1)
        
        #Выполняем расчеты
        (balanceMatrix,
         stateCoordinatesStreams,
         potentialInter,
         kineticMatrix) = CountState1(stateCoordinates,
                                      systemParameters)
        nonEqSystem.CountSystem(balanceMatrix,#Матрица баланса
                                stateCoordinatesStreams,#Внешние потоки по координата состояния
                                potentialInter,#Потенциалы взаимодействия
                                kineticMatrix#Кинетическая матрица
                                )
        
        #Проверяем значения
        self.assertEqual(etChemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(etBalanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(etKineticMatrix.tolist(), nonEqSystem.GetKineticMatrix().tolist())#Проверяем главный блок кинетической матрицы
        deltaAff = np.max(np.abs(etAff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaVProcesses = np.max(np.abs(etVProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVXPows = np.max(np.abs(etVx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 9)#Проверяем скорость сообщения теплоты энергетической степени свободы
    def testNonEqBase2(self):
        #Исходные данные
        stateCoordinates = np.array([14, 45, 57, 90, 141, 183, 43, 41, 85, 79])#Координаты состояния
        systemParameters = [  13,   21,   57,  71,   104,  203,   32,   48,   11,   23,
                              37,   53,   75,  85,   240,  310,  183,  159,  133, 16.3,
                            47.9, 8.63,  140, 7.5,  55.4,  415, 31.3, 21.5,  1.6, 35.9,
                            28.3, 4.43, 22.9, 8.5,  52.6, 23.9, 47.3, 35.9,  189, 21.9,
                            43.5, 7.53, 37.5, 1.5, 120.3]#Параметры системы
        
        #Стехиометрические коэффициенты химических реакций
        nu1_1 = 3
        nu1_2 = 1
        nu1_3 = 4
        nu1_4 = 1
        nu1_5 = 1
        nu1_6 = 3
        nu1_7 = 5
        nu1_8 = 1
        nu1_9 = 5.3
        nu1_10 = 6.7
        nu1_11 = 2.3
        nu2_1 = 3
        nu2_2 = 5.2
        nu2_3 = 4.3
        nu2_4 = 1.6
        nu2_5 = 2
        nu2_6 = 3
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        Adiff1_1 = 4.83
        Adiff1_2 = 0.45
        Adiff2_1 = 0.45
        
        #Эталонный результатрезультат
        (etChemPot,etAff,etKineticMatrix,
         etVProcesses,etBalanceMatrix,etVx) = CountSystem1(stateCoordinates,systemParameters,
                                                           nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                           nu1_8,nu1_9,nu1_10,nu1_11,nu2_1,nu2_2,nu2_3,
                                                           nu2_4,nu2_5,nu2_6,Adiff1_1,Adiff1_2,Adiff2_1)
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff1","vDiff2"]#Имена координат процессов
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена переменных потенциалов взаимодействия по координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой термодинамических сил
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []#Имена переменных внешних потоков
        
        #Задаем систему
        nonEqSystem = NonEqSystemBase(stateCoordinatesNames,#Имена координат состояния
                                      processCoordinatesNames,#Имена координат процессов
                                      stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                       
                                      #Задаем рассчитываемые параметры системы
                                      stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                      processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                      stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                        
                                      #Задаем переменные параметры кинетической матрицы
                                      varKineticNames,#Имена сопряженностей между собой координат процессов
                                      varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                        
                                      #Задаем внешние потоки
                                      stateCoordinatesVarStreamsNames#Имена переменных внешних потоков
                                      )
        
        #Задаем постоянные параметры системы
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_1", -nu1_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_1", -nu1_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4", "vChem1_1",  nu1_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_2", -nu1_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_2", -nu1_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_6", "vChem1_2",  nu1_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_2",  nu1_7)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_3", -nu1_8)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_3", -nu1_9)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_3", -nu1_10)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_3",  nu1_11)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vChem2_1", -nu2_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7", "vChem2_1", -nu2_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_8", "vChem2_1",  nu2_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4", "vChem2_2", -nu2_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vChem2_2", -nu2_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7", "vChem2_2",  nu2_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vDiff1", -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4", "vDiff2", -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vDiff1",  1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4", "vDiff2",  1.0)
        nonEqSystem.SetKineticMatrixConstElement("vDiff1", "vDiff1", Adiff1_1)
        nonEqSystem.SetKineticMatrixConstElement("vDiff1", "vDiff2", Adiff1_2)
        nonEqSystem.SetKineticMatrixConstElement("vDiff2", "vDiff1", Adiff2_1)
        
        #Выполняем расчеты
        (balanceMatrix,
         stateCoordinatesStreams,
         potentialInter,
         kineticMatrix) = CountState1(stateCoordinates,
                                      systemParameters)
        nonEqSystem.CountSystem(balanceMatrix,#Матрица баланса
                                stateCoordinatesStreams,#Внешние потоки по координата состояния
                                potentialInter,#Потенциалы взаимодействия
                                kineticMatrix#Кинетическая матрица
                                )
        
        #Проверяем значения
        self.assertEqual(etChemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(etBalanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(etKineticMatrix.tolist(), nonEqSystem.GetKineticMatrix().tolist())#Проверяем главный блок кинетической матрицы
        deltaAff = np.max(np.abs(etAff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaVProcesses = np.max(np.abs(etVProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVXPows = np.max(np.abs(etVx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 8)#Проверяем скорость сообщения теплоты энергетической степени свободы
    def testNonEqBase3(self):
        #Исходные данные
        stateCoordinates = np.array([34, 14, 57, 81, 441, 18, 39, 43, 82, 73])#Координаты состояния
        systemParameters = [  13,   21,   57,  71,   104,  203,   32,   48,   14,   23,
                              37,   53,   75,  85,   240,  310,  183,  159,  133, 16.3,
                            47.9, 7.63,  140, 7.5,  35.4,  415, 31.3, 21.5,  1.6, 32.9,
                            28.3, 4.73, 22.9, 6.5,  52.6, 21.9, 47.3, 35.9,  179, 23.9,
                            43.5, 7.53, 37.5, 1.5, 120.3]#Параметры системы
        
        #Стехиометрические коэффициенты химических реакций
        nu1_1 = 2
        nu1_2 = 3
        nu1_3 = 7
        nu1_4 = 1
        nu1_5 = 3
        nu1_6 = 8
        nu1_7 = 1
        nu1_8 = 2
        nu1_9 = 5.3
        nu1_10 = 4.7
        nu1_11 = 2.3
        nu2_1 = 3
        nu2_2 = 5.2
        nu2_3 = 4.3
        nu2_4 = 1.3
        nu2_5 = 2
        nu2_6 = 1
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        Adiff1_1 = 4.83
        Adiff1_2 = 0.15
        Adiff2_1 = 0.45
        
        #Эталонный результат
        (etChemPot,etAff,etKineticMatrix,
         etVProcesses,etBalanceMatrix,etVx) = CountSystem1(stateCoordinates,systemParameters,
                                                           nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                           nu1_8,nu1_9,nu1_10,nu1_11,nu2_1,nu2_2,nu2_3,
                                                           nu2_4,nu2_5,nu2_6,Adiff1_1,Adiff1_2,Adiff2_1)
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff1","vDiff2"]#Имена координат процессов
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена переменных потенциалов взаимодействия по координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой термодинамических сил
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []#Имена переменных внешних потоков
        
        #Задаем систему
        nonEqSystem = NonEqSystemBase(stateCoordinatesNames,#Имена координат состояния
                                      processCoordinatesNames,#Имена координат процессов
                                      stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                       
                                      #Задаем рассчитываемые параметры системы
                                      stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                      processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                      stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                        
                                      #Задаем переменные параметры кинетической матрицы
                                      varKineticNames,#Имена сопряженностей между собой координат процессов
                                      varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                        
                                      #Задаем внешние потоки
                                      stateCoordinatesVarStreamsNames#Имена переменных внешних потоков
                                      )
        
        #Задаем постоянные параметры системы
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_1", -nu1_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_1", -nu1_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4", "vChem1_1",  nu1_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_2", -nu1_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_2", -nu1_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_6", "vChem1_2",  nu1_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_2",  nu1_7)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_3", -nu1_8)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_3", -nu1_9)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_3", -nu1_10)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_3",  nu1_11)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vChem2_1", -nu2_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7", "vChem2_1", -nu2_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_8", "vChem2_1",  nu2_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4", "vChem2_2", -nu2_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vChem2_2", -nu2_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7", "vChem2_2",  nu2_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vDiff1", -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4", "vDiff2", -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vDiff1",  1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4", "vDiff2",  1.0)
        nonEqSystem.SetKineticMatrixConstElement("vDiff1", "vDiff1", Adiff1_1)
        nonEqSystem.SetKineticMatrixConstElement("vDiff1", "vDiff2", Adiff1_2)
        nonEqSystem.SetKineticMatrixConstElement("vDiff2", "vDiff1", Adiff2_1)
        
        #Выполняем расчеты
        (balanceMatrix,
         stateCoordinatesStreams,
         potentialInter,
         kineticMatrix) = CountState1(stateCoordinates,
                                      systemParameters)
        nonEqSystem.CountSystem(balanceMatrix,#Матрица баланса
                                stateCoordinatesStreams,#Внешние потоки по координата состояния
                                potentialInter,#Потенциалы взаимодействия
                                kineticMatrix#Кинетическая матрица
                                )
        
        #Проверяем значения
        self.assertEqual(etChemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(etBalanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(etKineticMatrix.tolist(), nonEqSystem.GetKineticMatrix().tolist())#Проверяем главный блок кинетической матрицы
        deltaAff = np.max(np.abs(etAff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaVProcesses = np.max(np.abs(etVProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVXPows = np.max(np.abs(etVx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 9)#Проверяем скорость сообщения теплоты энергетической степени свободы
    def testNonEqBase4(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81, 80, 64])#Координаты состояния
        systemParameters = [  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                              35,   55,   75,  81,   210,  300,  150,  153,  123, 15.3,
                            45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,  1.5, 33.9,
                            27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,  159, 21.9,
                            43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3, 0.15, 74.1]#Параметры системы
        
        #Стехиометрические коэффициенты химических реакций
        nu1_1 = 3
        nu1_2 = 2
        nu1_3 = 5
        nu1_4 = 1
        nu1_5 = 2
        nu1_6 = 3
        nu1_7 = 4
        nu1_8 = 1
        nu1_9 = 5
        nu1_10 = 7
        nu1_11 = 2.2
        nu1_12 = 2.1
        nu1_13 = 3.5
        nu1_14 = 4.1
        nu1_15 = 5.2
        nu2_1 = 3
        nu2_2 = 5
        nu2_3 = 4
        nu2_4 = 1
        nu2_5 = 2
        nu2_6 = 3
        nu2_7 = 5
        nu2_8 = 1
        nu2_9 = 6
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        AChem1_4_4 = 3.3
        AChem2_3_3 = 6.3
        AChem2_1_2 = 0.3
        AChem2_2_1 = 0.15
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3
        Adiff3_3 = 4.83
        
        #Внешние потоки вещества
        xExt1_6 = 30.9
        
        #Потенциал взаимодействия
        mu2_9 = 13.5

        #Эталонный результат
        (etChemPot,etAff,etKineticMatrix,etVProcesses,
         etBalanceMatrix,etVx,etStreams) = CountSystem2(stateCoordinates,systemParameters,
                                                        nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                        nu1_8,nu1_9,nu1_10,nu1_11,nu1_12,nu1_13,nu1_14,
                                                        nu1_15,nu2_1,nu2_2,nu2_3,nu2_4,nu2_5,nu2_6,
                                                        nu2_7,nu2_8,nu2_9,Adiff1_1,Adiff1_2,Adiff2_1,
                                                        AChem1_4_4,AChem2_3_3,AChem2_1_2,AChem2_2_1,
                                                        Adiff3_3,xExt1_6,mu2_9)
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8", "x2_9"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem1_4","vChem2_1","vChem2_2","vChem2_3","vDiff1","vDiff2","vDiff3"]#Имена координат процессов
        stateCoordinatesStreamsNames = ["x1_2","x1_6","x2_8"]#Имена координат состояния, изменяемых в результате внешних потоков
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена переменных потенциалов взаимодействия по координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой термодинамических сил
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = ["x1_2","x2_8"]#Имена переменных внешних потоков
        
        #Задаем систему
        nonEqSystem = NonEqSystemBase(stateCoordinatesNames,#Имена координат состояния
                                      processCoordinatesNames,#Имена координат процессов
                                      stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                       
                                      #Задаем рассчитываемые параметры системы
                                      stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                      processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                      stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                        
                                      #Задаем переменные параметры кинетической матрицы
                                      varKineticNames,#Имена сопряженностей между собой координат процессов
                                      varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                        
                                      #Задаем внешние потоки
                                      stateCoordinatesVarStreamsNames#Имена переменных внешних потоков
                                      )
        
        #Задаем постоянные параметры системы
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_1",  -nu1_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_1",  -nu1_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4", "vChem1_1",   nu1_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_2",  -nu1_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_2",  -nu1_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_6", "vChem1_2",   nu1_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_2",   nu1_7)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_3",  -nu1_8)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_3",  -nu1_9)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_3", -nu1_10)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_3",  nu1_11)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_7", "vChem1_4", -nu1_12)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_4", -nu1_13)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_4",  nu1_14)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_6", "vChem1_4",  nu1_15)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vChem2_1",  -nu2_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7", "vChem2_1",  -nu2_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_8", "vChem2_1",   nu2_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4", "vChem2_2",  -nu2_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vChem2_2",  -nu2_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7", "vChem2_2",   nu2_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_8", "vChem2_3",  -nu2_7)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_9", "vChem2_3",   nu2_8)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4", "vChem2_3",   nu2_9)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1",   "vDiff1",    -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4",   "vDiff2",    -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_7",   "vDiff3",    -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1",   "vDiff1",     1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4",   "vDiff2",     1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7",   "vDiff3",     1.0)
        nonEqSystem.SetKineticMatrixConstElement("vChem1_4", "vChem1_4", AChem1_4_4)
        nonEqSystem.SetKineticMatrixConstElement("vChem2_3", "vChem2_3", AChem2_3_3)
        nonEqSystem.SetKineticMatrixConstElement("vChem2_1", "vChem2_2", AChem2_1_2)
        nonEqSystem.SetKineticMatrixConstElement("vChem2_2", "vChem2_1", AChem2_2_1)
        nonEqSystem.SetKineticMatrixConstElement("vDiff1", "vDiff1", Adiff1_1)
        nonEqSystem.SetKineticMatrixConstElement("vDiff1", "vDiff2", Adiff1_2)
        nonEqSystem.SetKineticMatrixConstElement("vDiff2", "vDiff1", Adiff2_1)
        nonEqSystem.SetKineticMatrixConstElement("vDiff3", "vDiff3", Adiff3_3)
        nonEqSystem.SetStateCoordinatesStreamsConstElement("x1_6",xExt1_6)
        nonEqSystem.SetPotentialsInterConstElement("EnPow2", "x2_9", -mu2_9)
        
        #Выполняем расчеты
        (balanceMatrix,
         stateCoordinatesStreams,
         potentialInter,
         kineticMatrix) = CountState2(stateCoordinates,
                                      systemParameters)
        nonEqSystem.CountSystem(balanceMatrix,#Матрица баланса
                                stateCoordinatesStreams,#Внешние потоки по координата состояния
                                potentialInter,#Потенциалы взаимодействия
                                kineticMatrix#Кинетическая матрица
                                )
        
        #Проверяем значения
        self.assertEqual(etChemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(etBalanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(etKineticMatrix.tolist(), nonEqSystem.GetKineticMatrix().tolist())#Проверяем главный блок кинетической матрицы
        self.assertEqual(etStreams.tolist(), nonEqSystem.GetStateCoordinatesStreams().tolist())#Проверяем внешние тепловые потоки
        deltaAff = np.max(np.abs(etAff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaVProcesses = np.max(np.abs(etVProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVXPows = np.max(np.abs(etVx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 9)#Проверяем скорость сообщения теплоты энергетической степени свободы
    def testNonEqBase5(self):
        #Исходные данные
        stateCoordinates = np.array([17, 19, 49, 75, 121, 233, 36, 42, 94, 86, 61, 68])#Координаты состояния
        systemParameters = [  10,   20,   50,  70,   100,  200,    1,   55,   10,   25,
                              35,   55,   75,  81,   210,  300,  160,  173,  123, 15.3,
                            43.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,  1.5, 33.9,
                            27.3, 4.43, 21.9, 7.5,  52.6, 21.9, 45.3, 33.9,  139, 21.9,
                            41.5, 7.53, 36.5, 1.5, 120.3, 74.5, 15.3, 21.3, 0.15, 74.1]#Параметры системы
        
        #Стехиометрические коэффициенты химических реакций
        nu1_1 = 1
        nu1_2 = 4
        nu1_3 = 7
        nu1_4 = 2
        nu1_5 = 5
        nu1_6 = 6
        nu1_7 = 8
        nu1_8 = 2
        nu1_9 = 1
        nu1_10 = 8
        nu1_11 = 1.2
        nu1_12 = 3.1
        nu1_13 = 2.5
        nu1_14 = 1.1
        nu1_15 = 8.2
        nu2_1 = 10
        nu2_2 = 4
        nu2_3 = 6
        nu2_4 = 2
        nu2_5 = 1
        nu2_6 = 4
        nu2_7 = 5
        nu2_8 = 6
        nu2_9 = 2
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        AChem1_4_4 = 2.3
        AChem2_3_3 = 4.3
        AChem2_1_2 = 1.3
        AChem2_2_1 = 1.15
        Adiff1_1 = 2.83
        Adiff1_2 = 1.3
        Adiff2_1 = 1.3
        Adiff3_3 = 3.83
        
        #Внешние потоки вещества
        xExt1_6 = 31.9
        
        #Потенциал взаимодействия
        mu2_9 = 14.5
        
        #Эталонный результат
        (etChemPot,etAff,etKineticMatrix,etVProcesses,
         etBalanceMatrix,etVx,etStreams) = CountSystem2(stateCoordinates,systemParameters,
                                                        nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                        nu1_8,nu1_9,nu1_10,nu1_11,nu1_12,nu1_13,nu1_14,
                                                        nu1_15,nu2_1,nu2_2,nu2_3,nu2_4,nu2_5,nu2_6,
                                                        nu2_7,nu2_8,nu2_9,Adiff1_1,Adiff1_2,Adiff2_1,
                                                        AChem1_4_4,AChem2_3_3,AChem2_1_2,AChem2_2_1,
                                                        Adiff3_3,xExt1_6,mu2_9)
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8", "x2_9"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem1_4","vChem2_1","vChem2_2","vChem2_3","vDiff1","vDiff2","vDiff3"]#Имена координат процессов
        stateCoordinatesStreamsNames = ["x1_2","x1_6","x2_8"]#Имена координат состояния, изменяемых в результате внешних потоков
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена переменных потенциалов взаимодействия по координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой термодинамических сил
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = ["x1_2","x2_8"]#Имена переменных внешних потоков
        
        #Задаем систему
        nonEqSystem = NonEqSystemBase(stateCoordinatesNames,#Имена координат состояния
                                      processCoordinatesNames,#Имена координат процессов
                                      stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                       
                                      #Задаем рассчитываемые параметры системы
                                      stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                      processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                      stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                        
                                      #Задаем переменные параметры кинетической матрицы
                                      varKineticNames,#Имена сопряженностей между собой координат процессов
                                      varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                        
                                      #Задаем внешние потоки
                                      stateCoordinatesVarStreamsNames#Имена переменных внешних потоков
                                      )
        
        #Задаем постоянные параметры системы
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_1",  -nu1_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_1",  -nu1_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4", "vChem1_1",   nu1_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_2",  -nu1_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_2",  -nu1_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_6", "vChem1_2",   nu1_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_2",   nu1_7)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_3",  -nu1_8)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_3",  -nu1_9)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_3", -nu1_10)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_3",  nu1_11)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_7", "vChem1_4", -nu1_12)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_4", -nu1_13)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_4",  nu1_14)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_6", "vChem1_4",  nu1_15)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vChem2_1",  -nu2_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7", "vChem2_1",  -nu2_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_8", "vChem2_1",   nu2_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4", "vChem2_2",  -nu2_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vChem2_2",  -nu2_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7", "vChem2_2",   nu2_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_8", "vChem2_3",  -nu2_7)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_9", "vChem2_3",   nu2_8)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4", "vChem2_3",   nu2_9)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1",   "vDiff1",    -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4",   "vDiff2",    -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_7",   "vDiff3",    -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1",   "vDiff1",     1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4",   "vDiff2",     1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7",   "vDiff3",     1.0)
        nonEqSystem.SetKineticMatrixConstElement("vChem1_4", "vChem1_4", AChem1_4_4)
        nonEqSystem.SetKineticMatrixConstElement("vChem2_3", "vChem2_3", AChem2_3_3)
        nonEqSystem.SetKineticMatrixConstElement("vChem2_1", "vChem2_2", AChem2_1_2)
        nonEqSystem.SetKineticMatrixConstElement("vChem2_2", "vChem2_1", AChem2_2_1)
        nonEqSystem.SetKineticMatrixConstElement("vDiff1", "vDiff1", Adiff1_1)
        nonEqSystem.SetKineticMatrixConstElement("vDiff1", "vDiff2", Adiff1_2)
        nonEqSystem.SetKineticMatrixConstElement("vDiff2", "vDiff1", Adiff2_1)
        nonEqSystem.SetKineticMatrixConstElement("vDiff3", "vDiff3", Adiff3_3)
        nonEqSystem.SetStateCoordinatesStreamsConstElement("x1_6",xExt1_6)
        nonEqSystem.SetPotentialsInterConstElement("EnPow2", "x2_9", -mu2_9)
        
        #Выполняем расчеты
        (balanceMatrix,
         stateCoordinatesStreams,
         potentialInter,
         kineticMatrix) = CountState2(stateCoordinates,
                                      systemParameters)
        nonEqSystem.CountSystem(balanceMatrix,#Матрица баланса
                                stateCoordinatesStreams,#Внешние потоки по координата состояния
                                potentialInter,#Потенциалы взаимодействия
                                kineticMatrix#Кинетическая матрица
                                )
        
        #Проверяем значения
        self.assertEqual(etChemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(etBalanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(etKineticMatrix.tolist(), nonEqSystem.GetKineticMatrix().tolist())#Проверяем главный блок кинетической матрицы
        self.assertEqual(etStreams.tolist(), nonEqSystem.GetStateCoordinatesStreams().tolist())#Проверяем внешние тепловые потоки
        deltaAff = np.max(np.abs(etAff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaVProcesses = np.max(np.abs(etVProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVXPows = np.max(np.abs(etVx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 8)#Проверяем скорость сообщения теплоты энергетической степени свободы
    def testNonEqBase6(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 56, 75, 111, 223, 39, 45, 93, 81, 80, 64])#Координаты состояния
        systemParameters = [  10,   20,   50,  70,   100,  200,   30,   45,   11,   65,
                              35,   55,   75,  81,   280,  300,  150,  153,  123, 15.3,
                            45.9, 5.53,   10, 6.5,  45.4,  315, 21.3,  4.5,  1.5, 33.9,
                            27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,  159, 21.9,
                            43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3, 0.15, 74.1]#Параметры системы
        
        #Стехиометрические коэффициенты химических реакций
        nu1_1 = 3
        nu1_2 = 2
        nu1_3 = 5
        nu1_4 = 1
        nu1_5 = 1
        nu1_6 = 3
        nu1_7 = 4
        nu1_8 = 1
        nu1_9 = 4
        nu1_10 = 7
        nu1_11 = 2.2
        nu1_12 = 2.1
        nu1_13 = 3.8
        nu1_14 = 4.1
        nu1_15 = 4.7
        nu2_1 = 3
        nu2_2 = 5
        nu2_3 = 1
        nu2_4 = 2
        nu2_5 = 6
        nu2_6 = 3
        nu2_7 = 5
        nu2_8 = 7
        nu2_9 = 6

        #Базовые коэффициенты главной кинетической матрицы процессов]
        AChem1_4_4 = 3.3
        AChem2_3_3 = 6.3
        AChem2_1_2 = 0.3
        AChem2_2_1 = 0.15
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3
        Adiff3_3 = 4.83

        #Внешние потоки вещества
        xExt1_6 = 31.9
        
        #Потенциал взаимодействия
        mu2_9 = 13.5
        
        #Эталонный результат
        (etChemPot,etAff,etKineticMatrix,etVProcesses,
         etBalanceMatrix,etVx,etStreams) = CountSystem2(stateCoordinates,systemParameters,
                                                        nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                        nu1_8,nu1_9,nu1_10,nu1_11,nu1_12,nu1_13,nu1_14,
                                                        nu1_15,nu2_1,nu2_2,nu2_3,nu2_4,nu2_5,nu2_6,
                                                        nu2_7,nu2_8,nu2_9,Adiff1_1,Adiff1_2,Adiff2_1,
                                                        AChem1_4_4,AChem2_3_3,AChem2_1_2,AChem2_2_1,
                                                        Adiff3_3,xExt1_6,mu2_9)
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8", "x2_9"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem1_4","vChem2_1","vChem2_2","vChem2_3","vDiff1","vDiff2","vDiff3"]#Имена координат процессов
        stateCoordinatesStreamsNames = ["x1_2","x1_6","x2_8"]#Имена координат состояния, изменяемых в результате внешних потоков
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена переменных потенциалов взаимодействия по координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой термодинамических сил
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = ["x1_2","x2_8"]#Имена переменных внешних потоков
        
        #Задаем систему
        nonEqSystem = NonEqSystemBase(stateCoordinatesNames,#Имена координат состояния
                                      processCoordinatesNames,#Имена координат процессов
                                      stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                        
                                      #Задаем рассчитываемые параметры системы
                                      stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                      processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                      stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                        
                                      #Задаем переменные параметры кинетической матрицы
                                      varKineticNames,#Имена сопряженностей между собой координат процессов
                                      varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                        
                                      #Задаем внешние потоки
                                      stateCoordinatesVarStreamsNames#Имена переменных внешних потоков
                                      )
        
        #Задаем постоянные параметры системы
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_1",  -nu1_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_1",  -nu1_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4", "vChem1_1",   nu1_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_2",  -nu1_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_2",  -nu1_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_6", "vChem1_2",   nu1_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_2",   nu1_7)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_3",  -nu1_8)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_3",  -nu1_9)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_3", -nu1_10)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_3",  nu1_11)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_7", "vChem1_4", -nu1_12)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_4", -nu1_13)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_4",  nu1_14)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_6", "vChem1_4",  nu1_15)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vChem2_1",  -nu2_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7", "vChem2_1",  -nu2_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_8", "vChem2_1",   nu2_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4", "vChem2_2",  -nu2_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1", "vChem2_2",  -nu2_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7", "vChem2_2",   nu2_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_8", "vChem2_3",  -nu2_7)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_9", "vChem2_3",   nu2_8)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4", "vChem2_3",   nu2_9)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1",   "vDiff1",    -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4",   "vDiff2",    -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_7",   "vDiff3",    -1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_1",   "vDiff1",     1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_4",   "vDiff2",     1.0)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x2_7",   "vDiff3",     1.0)
        nonEqSystem.SetKineticMatrixConstElement("vChem1_4", "vChem1_4", AChem1_4_4)
        nonEqSystem.SetKineticMatrixConstElement("vChem2_3", "vChem2_3", AChem2_3_3)
        nonEqSystem.SetKineticMatrixConstElement("vChem2_1", "vChem2_2", AChem2_1_2)
        nonEqSystem.SetKineticMatrixConstElement("vChem2_2", "vChem2_1", AChem2_2_1)
        nonEqSystem.SetKineticMatrixConstElement("vDiff1", "vDiff1", Adiff1_1)
        nonEqSystem.SetKineticMatrixConstElement("vDiff1", "vDiff2", Adiff1_2)
        nonEqSystem.SetKineticMatrixConstElement("vDiff2", "vDiff1", Adiff2_1)
        nonEqSystem.SetKineticMatrixConstElement("vDiff3", "vDiff3", Adiff3_3)
        nonEqSystem.SetStateCoordinatesStreamsConstElement("x1_6",xExt1_6)
        nonEqSystem.SetPotentialsInterConstElement("EnPow2", "x2_9", -mu2_9)
        
        #Выполняем расчеты
        (balanceMatrix,
         stateCoordinatesStreams,
         potentialInter,
         kineticMatrix) = CountState2(stateCoordinates,
                                      systemParameters)
        nonEqSystem.CountSystem(balanceMatrix,#Матрица баланса
                                stateCoordinatesStreams,#Внешние потоки по координата состояния
                                potentialInter,#Потенциалы взаимодействия
                                kineticMatrix#Кинетическая матрица
                                )
        
        #Проверяем значения
        self.assertEqual(etChemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(etBalanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(etKineticMatrix.tolist(), nonEqSystem.GetKineticMatrix().tolist())#Проверяем главный блок кинетической матрицы
        self.assertEqual(etStreams.tolist(), nonEqSystem.GetStateCoordinatesStreams().tolist())#Проверяем внешние тепловые потоки
        deltaAff = np.max(np.abs(etAff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaVProcesses = np.max(np.abs(etVProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVXPows = np.max(np.abs(etVx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 8)#Проверяем скорость сообщения теплоты энергетической степени свободы
    def testNonEqBase7(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213])#Координаты состояния
        systemParameters = [  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                              35,   55,   75,  81,   210,  300,  150,  153,  123, 15.3, 
                            45.9, 5.63,   10, 6.5,  45.4,  315,  161]#Параметры системы
        
        #Стехиометрические коэффициенты химических реакций
        nu1_1 = 3
        nu1_2 = 2
        nu1_3 = 5
        nu1_4 = 1
        nu1_5 = 2
        nu1_6 = 3
        nu1_7 = 4
        nu1_8 = 1
        nu1_9 = 5
        nu1_10 = 7
        nu1_11 = 2.2
        
        #Эталонный результат
        (etChemPot,etAff,etKineticMatrix,etVProcesses,
         etBalanceMatrix,etVx) = CountSystem3(stateCoordinates,systemParameters,
                                              nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,
                                              nu1_7,nu1_8,nu1_9,nu1_10,nu1_11)
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3"]#Имена координат процессов
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена переменных потенциалов взаимодействия по координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3"]#Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3"]#Имена сопряженностей между собой термодинамических сил
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []#Имена переменных внешних потоков
        
        #Задаем систему
        nonEqSystem = NonEqSystemBase(stateCoordinatesNames,#Имена координат состояния
                                      processCoordinatesNames,#Имена координат процессов
                                      stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                        
                                      #Задаем рассчитываемые параметры системы
                                      stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                      processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                      stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                        
                                      #Задаем переменные параметры кинетической матрицы
                                      varKineticNames,#Имена сопряженностей между собой координат процессов
                                      varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                        
                                      #Задаем внешние потоки
                                      stateCoordinatesVarStreamsNames#Имена переменных внешних потоков
                                      )
        
        #Задаем постоянные параметры системы
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_1", -nu1_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_1", -nu1_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4", "vChem1_1",  nu1_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_2", -nu1_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_2", -nu1_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_6", "vChem1_2",  nu1_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_2",  nu1_7)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_3", -nu1_8)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_3", -nu1_9)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_3", -nu1_10)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_3",  nu1_11)
        
        #Выполняем расчеты
        (balanceMatrix,
         stateCoordinatesStreams,
         potentialInter,
         kineticMatrix) = CountState3(stateCoordinates,
                                      systemParameters)
        nonEqSystem.CountSystem(balanceMatrix,#Матрица баланса
                                stateCoordinatesStreams,#Внешние потоки по координата состояния
                                potentialInter,#Потенциалы взаимодействия
                                kineticMatrix#Кинетическая матрица
                                )
        
        #Проверяем значения
        self.assertEqual(etChemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(etBalanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(etKineticMatrix.tolist(), nonEqSystem.GetKineticMatrix().tolist())#Проверяем главный блок кинетической матрицы
        deltaAff = np.max(np.abs(etAff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaVProcesses = np.max(np.abs(etVProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVXPows = np.max(np.abs(etVx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 9)#Проверяем скорость сообщения теплоты энергетической степени свободы         
    def testNonEqBase8(self):
        #Исходные данные
        stateCoordinates = np.array([13, 26, 52, 73, 131, 223])#Координаты состояния
        systemParameters = [  10,   20,   50,  70,   100,  220,   31,   45,   10,   25,
                              35,   56,   75,  83,   211,  300,  152,  156,  143, 15.4, 
                            45.9, 5.63,   11, 6.5,  44.4,  315, 20.3]#Параметры системы
        
        #Стехиометрические коэффициенты химических реакций
        nu1_1 = 4
        nu1_2 = 1
        nu1_3 = 5
        nu1_4 = 1
        nu1_5 = 4
        nu1_6 = 3
        nu1_7 = 7
        nu1_8 = 1
        nu1_9 = 4
        nu1_10 = 7.5
        nu1_11 = 2.8
        
        #Эталонный результат
        (etChemPot,etAff,etKineticMatrix,etVProcesses,
         etBalanceMatrix,etVx) = CountSystem3(stateCoordinates,systemParameters,
                                              nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,
                                              nu1_7,nu1_8,nu1_9,nu1_10,nu1_11)
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3"]#Имена координат процессов
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена переменных потенциалов взаимодействия по координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3"]#Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3"]#Имена сопряженностей между собой термодинамических сил
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []#Имена переменных внешних потоков
        
        #Задаем систему
        nonEqSystem = NonEqSystemBase(stateCoordinatesNames,#Имена координат состояния
                                      processCoordinatesNames,#Имена координат процессов
                                      stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                        
                                      #Задаем рассчитываемые параметры системы
                                      stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                      processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                      stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                        
                                      #Задаем переменные параметры кинетической матрицы
                                      varKineticNames,#Имена сопряженностей между собой координат процессов
                                      varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                        
                                      #Задаем внешние потоки
                                      stateCoordinatesVarStreamsNames#Имена переменных внешних потоков
                                      )
        
        #Задаем постоянные параметры системы
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_1", -nu1_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_1", -nu1_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4", "vChem1_1",  nu1_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_2", -nu1_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_2", -nu1_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_6", "vChem1_2",  nu1_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_2",  nu1_7)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_3", -nu1_8)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_3", -nu1_9)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_3", -nu1_10)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_3",  nu1_11)
        
        #Выполняем расчеты
        (balanceMatrix,
         stateCoordinatesStreams,
         potentialInter,
         kineticMatrix) = CountState3(stateCoordinates,
                                      systemParameters)
        nonEqSystem.CountSystem(balanceMatrix,#Матрица баланса
                                stateCoordinatesStreams,#Внешние потоки по координата состояния
                                potentialInter,#Потенциалы взаимодействия
                                kineticMatrix#Кинетическая матрица
                                )
        
        #Проверяем значения
        self.assertEqual(etChemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(etBalanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(etKineticMatrix.tolist(), nonEqSystem.GetKineticMatrix().tolist())#Проверяем главный блок кинетической матрицы
        deltaAff = np.max(np.abs(etAff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaVProcesses = np.max(np.abs(etVProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVXPows = np.max(np.abs(etVx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 9)#Проверяем скорость сообщения теплоты энергетической степени свободы             
    def testNonEqBase9(self):
        #Исходные данные
        stateCoordinates = np.array([17, 22, 56, 76, 136, 228])#Координаты состояния
        systemParameters = [  10,   22,   50,  70,   100,  220,   31,   45,   10,   25,
                              35,   58,   75,  84,   211,  300,  132,  156,  146, 15.4, 
                            45.9, 5.63,   11, 6.5,  48.4,  314, 20.3]#Параметры системы
        
        #Стехиометрические коэффициенты химических реакций
        nu1_1 = 4
        nu1_2 = 2
        nu1_3 = 5
        nu1_4 = 1
        nu1_5 = 3
        nu1_6 = 3
        nu1_7 = 7
        nu1_8 = 5
        nu1_9 = 4
        nu1_10 = 6.5
        nu1_11 = 3.8
        
        #Эталонный результат
        (etChemPot,etAff,etKineticMatrix,etVProcesses,
         etBalanceMatrix,etVx) = CountSystem3(stateCoordinates,systemParameters,
                                              nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,
                                              nu1_7,nu1_8,nu1_9,nu1_10,nu1_11)
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3"]#Имена координат процессов
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена переменных потенциалов взаимодействия по координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3"]#Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3"]#Имена сопряженностей между собой термодинамических сил
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []#Имена переменных внешних потоков
        
        #Задаем систему
        nonEqSystem = NonEqSystemBase(stateCoordinatesNames,#Имена координат состояния
                                      processCoordinatesNames,#Имена координат процессов
                                      stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                       
                                      #Задаем рассчитываемые параметры системы
                                      stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                      processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                      stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                        
                                      #Задаем переменные параметры кинетической матрицы
                                      varKineticNames,#Имена сопряженностей между собой координат процессов
                                      varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                                        
                                      #Задаем внешние потоки
                                      stateCoordinatesVarStreamsNames#Имена переменных внешних потоков
                                      )
        
        #Задаем постоянные параметры системы
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_1", -nu1_1)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_1", -nu1_2)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_4", "vChem1_1",  nu1_3)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_2", -nu1_4)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_2", -nu1_5)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_6", "vChem1_2",  nu1_6)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_2",  nu1_7)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_1", "vChem1_3", -nu1_8)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_3", "vChem1_3", -nu1_9)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_5", "vChem1_3", -nu1_10)
        nonEqSystem.SetBalanceStateCoordinatesConstElement("x1_2", "vChem1_3",  nu1_11)
        
        #Выполняем расчеты
        (balanceMatrix,
         stateCoordinatesStreams,
         potentialInter,
         kineticMatrix) = CountState3(stateCoordinates,
                                      systemParameters)
        nonEqSystem.CountSystem(balanceMatrix,#Матрица баланса
                                stateCoordinatesStreams,#Внешние потоки по координата состояния
                                potentialInter,#Потенциалы взаимодействия
                                kineticMatrix#Кинетическая матрица
                                )
        
        #Проверяем значения
        self.assertEqual(etChemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(etBalanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(etKineticMatrix.tolist(), nonEqSystem.GetKineticMatrix().tolist())#Проверяем главный блок кинетической матрицы
        deltaAff = np.max(np.abs(etAff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaVProcesses = np.max(np.abs(etVProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVXPows = np.max(np.abs(etVx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 8)#Проверяем скорость сообщения теплоты энергетической степени свободы          
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
