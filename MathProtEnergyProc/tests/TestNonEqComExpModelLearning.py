import numpy as np
from scipy import integrate as integ
from MathProtEnergyProc.IntegrateDyn import standartIntegrateDyn

from MathProtEnergyProc import NonEqSystem
from MathProtEnergyProc import NonEqSystemDyn
from MathProtEnergyProc import CountDynamics, ModelLearning

from MathProtEnergyProc.tests.UnitTestExamples.TestNonEq1 import *
from MathProtEnergyProc.tests.UnitTestExamples.TestNonEq2 import *
from MathProtEnergyProc.tests.UnitTestExamples.TestNonEq3 import *

from MathProtEnergyProc.tests.UnitTestExamples.fU1 import *
from MathProtEnergyProc.tests.UnitTestExamples.fU2 import *
from MathProtEnergyProc.tests.UnitTestExamples.fU3 import *

from MathProtEnergyProc.tests.UnitTestExamples.FunChar1 import *
from MathProtEnergyProc.tests.UnitTestExamples.FunChar2 import *
from MathProtEnergyProc.tests.UnitTestExamples.FunChar3 import *

import unittest

#Модульные тесты
class TestNonEqSystemComExpModelLearning(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testNonEqDyns1(self):
        #Исходные данные
        stateCoordinates1 = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81])#Координаты состояния
        systemParameters1 = np.array([  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                                        35,   55,   75,  81,   210,  300,  150,  153,  123, 15.3,
                                      45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,  1.5, 33.9,
                                      27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,  159, 21.9,
                                      43.5, 7.53, 37.5, 1.5, 120.3])#Параметры системы
        stateCoordinates2 = np.array([11, 23, 45, 75, 141, 113, 39, 95, 83, 73])#Координаты состояния
        systemParameters2 = np.array([  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                                        35,   54,   75,  81,   210,  300,  153,  153,  123, 15.3,
                                      45.9, 5.43,   10, 6.5,  45.4,  315, 21.3, 4.83,  1.5, 33.9,
                                      27.3, 4.53, 21.9, 7.5,  52.6, 22.9, 45.3, 33.9,  159, 21.9,
                                      43.5, 7.53, 37.5, 1.5, 120.3])#Параметры системы
        stateCoordinates3 = np.array([34, 20, 55, 75, 411, 213, 39, 55, 93, 81])#Координаты состояния
        systemParameters3 = np.array([  10,   20,   51,  70,   100,  200,   30,   45,   10,   25,
                                        35,   55,   85,  81,   210,  300,  150,  153,  123, 15.3,
                                      45.9, 5.63,   10, 6.5,  45.4,  615, 21.3,  4.5,  1.5, 33.9,
                                      21.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,  159, 21.9,
                                      43.5, 7.53, 37.5, 1.5, 120.3])#Параметры системы
        
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
        
        #Метод интегрирования дифференциальных уравнений
        method = "RK45"
        integDynamic = standartIntegrateDyn(method=method)
        
        #Время интегрирования
        Tint = 1e-6
        
        #Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)
        
        #Функция эталонного результата
        def countSystem1(t,stateCoordinates,systemParameters):
            #Вызов функции времени
            _systemParameters = fU1(t,systemParameters)
            
            #Вызов расчета
            (chemPot,Aff,kineticMatrix,
             vProcesses,balanceMatrix,vx) = CountSystem1(stateCoordinates,_systemParameters,
                                                         nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                         nu1_8,nu1_9,nu1_10,nu1_11,nu2_1,nu2_2,nu2_3,
                                                         nu2_4,nu2_5,nu2_6,Adiff1_1,Adiff1_2,Adiff2_1)
            
            #Выводим результат
            return vx
        def countSystem1_1(t,stateCoordinates):
            #Выводим результат
            return countSystem1(t,stateCoordinates,systemParameters1)
        def countSystem1_2(t,stateCoordinates):
            #Выводим результат
            return countSystem1(t,stateCoordinates,systemParameters2)
        def countSystem1_3(t,stateCoordinates):
            #Выводим результат
            return countSystem1(t,stateCoordinates,systemParameters3)
            
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff1","vDiff2"]#Имена координат процессов
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        
        #Задаем функцию состояния системы
        stateFunction = CountState1
        
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
        nonEqSystem = NonEqSystem(stateCoordinatesNames,#Имена координат состояния
                                  processCoordinatesNames,#Имена координат процессов
                                  stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                     
                                  #Задаем функцию состояния системы
                                  stateFunction,
                                    
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
        
        #Создаем динамику системы
        nonEqSystemDyn = NonEqSystemDyn(nonEqSystem,#Система
                 
                                        fU1,#Функция условий протекания процессов
                                         
                                        FunChar1,#Функция внешних параметров
                                         
                                        integDynamic#Метод интегрирования дифференциальных уравнений
                                        )
        
        #Рассчитываем эталонную динамику 1 системы
        stateCoordinatesEt = integ.solve_ivp(countSystem1_1,
                                             (0.0, Tint),
                                             stateCoordinates1,#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Рассчитываем эталонную динамику 1 параметров системы
        (t1,Char1et1,Char2et1) = FunChar1(t,#Моменты времени
                                          stateCoordinatesEt,#Координаты состояния
                                          systemParameters1#Параметры системы  
                                          )
        (char1et1,char2et1) = (np.sum(Char1et1*t1),np.sum(Char2et1*t1))
        
        #Рассчитываем эталонную динамику 2 системы
        stateCoordinatesEt = integ.solve_ivp(countSystem1_2,
                                             (0.0, Tint),
                                             stateCoordinates2,#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Рассчитываем эталонную динамику 2 параметров системы
        (t2,Char1et2,Char2et2) = FunChar1(t,#Моменты времени
                                          stateCoordinatesEt,#Координаты состояния
                                          systemParameters2#Параметры системы  
                                          )
        (char1et2,char2et2) = (np.sum(Char1et2*t2),np.sum(Char2et2*t2))
        
        #Рассчитываем эталонную динамику 3 системы
        stateCoordinatesEt = integ.solve_ivp(countSystem1_3,
                                             (0.0, Tint),
                                             stateCoordinates3,#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Рассчитываем эталонную динамику 3 параметров системы
        (t3,Char1et3,Char2et3) = FunChar1(t,#Моменты времени
                                          stateCoordinatesEt,#Координаты состояния
                                          systemParameters3#Параметры системы  
                                          )
        (char1et3,char2et3) = (np.sum(Char1et3*t3),np.sum(Char2et3*t3))
        
        #Эталонные результаты
        char1et = np.array([char1et1,char1et2,char1et3]).reshape(-1,1)
        char2et = np.array([char2et1,char2et2,char2et3]).reshape(-1,1)
        
        #Функция сжатия
        def CompressFunction(dyn, index):
            return (np.sum(dyn[0]*dyn[1]),np.sum(dyn[0]*dyn[2]))
        
        #Создаем класс вычислительного эксперимента
        cDyn = CountDynamics(nonEqSystemDyn,
                             CompressFunction)
        
        #Функтор модели параметров системы
        class learningSystemPar(object):
            #Функтор параметров
            def ModelSystemParametersCount(self,comExpRez,
                                           cStateCoordinates0,
                                           cSystemParameters):
                self.Par = (comExpRez[0],
                            cStateCoordinates0[0],
                            cSystemParameters[0])
        
        #Создаем класс векторной модели
        funSysPar = learningSystemPar()
        cDynLearn = ModelLearning(funSysPar,cDyn)
        
        #Рассчитываем динамику параметров системы
        stateCoordinates = np.vstack([stateCoordinates1,
                                      stateCoordinates2,
                                      stateCoordinates3])
        systemParameters = np.vstack([systemParameters1,
                                      systemParameters2,
                                      systemParameters3])
        t_evals = [t_eval,
                   t_eval,
                   t_eval]
        Tints = np.array([Tint,
                          Tint])
        cDynLearn.ModelLearning(Tints,
                                stateCoordinates,
                                systemParameters,
                                ct_evals = t_evals)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(funSysPar.Par[0] - char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(funSysPar.Par[1] - stateCoordinates[0]))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(funSysPar.Par[2] - systemParameters[0]))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)
    def testNonEqDyns2(self):
        #Исходные данные
        stateCoordinates1 = np.array([14, 20, 55, 75, 111, 217, 29, 45, 96, 81, 80, 64])#Координаты состояния
        systemParameters1 = np.array([  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                                        35,   55,   75,  71,   210,  300,  160,  153,  123, 15.3,
                                      45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,  1.7, 33.9,
                                      23.3, 4.53, 21.9, 3.5,  53.6, 24.9, 45.3, 36.9,  159, 21.9,
                                      43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3, 0.15, 74.1])#Параметры системы
        stateCoordinates2 = np.array([10, 22, 55, 75, 111, 213, 39, 35, 93, 11, 10, 34])#Координаты состояния
        systemParameters2 = np.array([  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                                        35,   25,   75,  81,   210,  300,  150,  153,  123, 15.3,
                                      45.9, 5.63,   10, 6.5,  41.4,  315, 21.3,  4.5,  1.5, 33.9,
                                      27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,  159, 21.9,
                                      43.5, 7.53, 37.5, 1.5, 121.3, 73.5, 15.3, 21.3, 0.15, 74.1])#Параметры системы
        stateCoordinates3 = np.array([10, 20, 55, 75, 121, 213, 39, 45, 97, 81, 83, 64])#Координаты состояния
        systemParameters3 = np.array([  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                                        35,   55,   74,  81,   210,  300,  150,  153,  123, 15.3,
                                      45.9, 5.63,   10, 6.5,  45.4,  615, 21.3,  4.5,  1.5, 33.9,
                                      27.5, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 31.9,  159, 21.9,
                                      43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3, 0.15, 74.1])#Параметры системы
        stateCoordinates4 = np.array([10, 23, 55, 75, 111, 213, 39, 47, 93, 81, 82, 64])#Координаты состояния
        systemParameters4 = np.array([  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                                        35,   55,   75,  81,   210,  300,  151,  153,  123, 15.3,
                                      45.9, 5.63,   10, 6.5,  41.4,  315, 21.3,  4.5,  1.5, 33.9,
                                      22.3, 5.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,  149, 21.9,
                                      43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 11.3, 21.3, 0.15, 74.1])#Параметры системы
        
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

        #Метод интегрирования дифференциальных уравнений
        method = "RK23"
        integDynamic = standartIntegrateDyn(method=method)
        
        #Время интегрирования
        Tint = 1e-6
        
        #Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)
        
        #Функция эталонного результата
        def countSystem2(t,stateCoordinates,systemParameters):
            #Вызов функции времени
            _systemParameters = fU2(t,systemParameters)
            
            #Вызов расчета
            (chemPot,Aff,kineticMatrix,vProcesses,
             balanceMatrix,vx,Streams) = CountSystem2(stateCoordinates,_systemParameters,
                                                      nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                      nu1_8,nu1_9,nu1_10,nu1_11,nu1_12,nu1_13,nu1_14,
                                                      nu1_15,nu2_1,nu2_2,nu2_3,nu2_4,nu2_5,nu2_6,
                                                      nu2_7,nu2_8,nu2_9,Adiff1_1,Adiff1_2,Adiff2_1,
                                                      AChem1_4_4,AChem2_3_3,AChem2_1_2,AChem2_2_1,
                                                      Adiff3_3,xExt1_6,mu2_9)
            
            #Выводим результат
            return vx
        def countSystem2_1(t,stateCoordinates):
            return countSystem2(t,stateCoordinates,systemParameters1)
        def countSystem2_2(t,stateCoordinates):
            return countSystem2(t,stateCoordinates,systemParameters2)
        def countSystem2_3(t,stateCoordinates):
            return countSystem2(t,stateCoordinates,systemParameters3)
        def countSystem2_4(t,stateCoordinates):
            return countSystem2(t,stateCoordinates,systemParameters4)
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8", "x2_9"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem1_4","vChem2_1","vChem2_2","vChem2_3","vDiff1","vDiff2","vDiff3"]#Имена координат процессов
        stateCoordinatesStreamsNames = ["x1_2","x1_6","x2_8"]#Имена координат состояния, изменяемых в результате внешних потоков
        
        #Задаем функцию состояния системы
        stateFunction = CountState2
        
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
        nonEqSystem = NonEqSystem(stateCoordinatesNames,#Имена координат состояния
                                  processCoordinatesNames,#Имена координат процессов
                                  stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                    
                                  #Задаем функцию состояния системы
                                  stateFunction,
                                    
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
        
        #Создаем динамику системы
        nonEqSystemDyn = NonEqSystemDyn(nonEqSystem,#Система
                 
                                        fU2,#Функция условий протекания процессов
                                         
                                        FunChar2,#Функция внешних параметров
                                         
                                        integDynamic#Метод интегрирования дифференциальных уравнений
                                        )
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem2_1,
                                             (0.0, Tint),
                                             stateCoordinates1,#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Рассчитываем эталонную динамику параметров системы
        (t1,Char1et1,Char2et1,Char3et1) = FunChar2(t,#Моменты времени
                                                   stateCoordinatesEt,#Координаты состояния
                                                   systemParameters1#Параметры системы  
                                                   )
        (char1et1,char2et1,char3et1) = (np.sum(t1*Char1et1*Char3et1),np.sum(t1*Char2et1 + Char3et1),np.sum(t1*Char1et1))
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem2_2,
                                             (0.0, Tint),
                                             stateCoordinates2,#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Рассчитываем эталонную динамику параметров системы
        (t2,Char1et2,Char2et2,Char3et2) = FunChar2(t,#Моменты времени
                                                   stateCoordinatesEt,#Координаты состояния
                                                   systemParameters2#Параметры системы  
                                                   )
        (char1et2,char2et2,char3et2) = (np.sum(t2*Char1et2*Char3et2),np.sum(t2*Char2et2 + Char3et2),np.sum(t2*Char1et2))
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem2_3,
                                             (0.0, Tint),
                                             stateCoordinates3,#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Рассчитываем эталонную динамику параметров системы
        (t3,Char1et3,Char2et3,Char3et3) = FunChar2(t,#Моменты времени
                                                   stateCoordinatesEt,#Координаты состояния
                                                   systemParameters3#Параметры системы  
                                                   )
        (char1et3,char2et3,char3et3) = (np.sum(t3*Char1et3*Char3et3),np.sum(t3*Char2et3 + Char3et3),np.sum(t3*Char1et3))
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem2_4,
                                             (0.0, Tint),
                                             stateCoordinates4,#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Рассчитываем эталонную динамику параметров системы
        (t4,Char1et4,Char2et4,Char3et4) = FunChar2(t,#Моменты времени
                                                   stateCoordinatesEt,#Координаты состояния
                                                   systemParameters4#Параметры системы  
                                                   )
        (char1et4,char2et4,char3et4) = (np.sum(t4*Char1et4*Char3et4),np.sum(t4*Char2et4 + Char3et4),np.sum(t4*Char1et4))
        
        #Эталонные результаты
        char1et = np.array([char1et1,char1et2,char1et3,char1et4]).reshape(-1,1)
        char2et = np.array([char2et1,char2et2,char2et3,char2et4]).reshape(-1,1)
        char3et = np.array([char3et1,char3et2,char3et3,char3et4]).reshape(-1,1)
        
        #Функция сжатия
        def CompressFunction(dyn, index):
            #Выводим результат
            return (np.sum(dyn[0]*dyn[1]*dyn[3]),np.sum(dyn[0]*dyn[2] + dyn[3]),np.sum(dyn[0]*dyn[1]))
        
        #Создаем класс вычислительного эксперимента
        cDyn = CountDynamics(nonEqSystemDyn,
                             CompressFunction)
        
        #Функтор модели параметров системы
        class learningSystemPar(object):
            #Функтор параметров
            def ModelSystemParametersCount(self,comExpRez,
                                           cStateCoordinates0,
                                           cSystemParameters):
                self.Par = (comExpRez[2],
                            cStateCoordinates0[1],
                            cSystemParameters[1])
        
        #Создаем класс векторной модели
        funSysPar = learningSystemPar()
        cDynLearn = ModelLearning(funSysPar,cDyn)
        
        #Рассчитываем динамику параметров системы
        stateCoordinates = np.vstack([stateCoordinates1,
                                      stateCoordinates2,
                                      stateCoordinates3,
                                      stateCoordinates4])
        systemParameters = np.vstack([systemParameters1,
                                      systemParameters2,
                                      systemParameters3,
                                      systemParameters4])
        t_evals = [t_eval,
                   t_eval,
                   t_eval,
                   t_eval]
        Tints = np.array([Tint,
                          Tint])
        cDynLearn.ModelLearning(Tints,
                                stateCoordinates,
                                systemParameters,
                                ct_evals = t_evals)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(funSysPar.Par[0] - char3et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(funSysPar.Par[1] - stateCoordinates[1]))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(funSysPar.Par[2] - systemParameters[1]))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)
    def testNonEqDyns3(self):
        #Исходные данные
        stateCoordinates1 = np.array([10, 22, 54, 75, 117, 211])#Координаты состояния
        systemParameters1 = np.array([  10,   20,   50,  70,   101,  200,   30,   47,   10,   25,
                                        35,   55,   75,  82,   210,  300,  150,  153,  123, 15.3, 
                                      45.9, 5.63,   10, 6.5,  45.4,  315,  181])#Параметры системы
        stateCoordinates2 = np.array([10, 20, 55, 75, 111, 213])#Координаты состояния
        systemParameters2 = np.array([  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                                        35,   55,   75,  81,   210,  300,  150,  153,  123, 15.3, 
                                      45.9, 5.63,   10, 6.5,  45.4,  315,  161])#Параметры системы
        
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
        
        #Метод интегрирования дифференциальных уравнений
        method = "LSODA"
        integDynamic = standartIntegrateDyn(method=method)
        
        #Время интегрирования
        Tint = 1e-6
        
        #Функция эталонного результата
        def countSystem3(t,stateCoordinates,systemParameters):
            #Вызов функции времени
            _systemParameters = fU3(t,systemParameters)
            
            #Вызов расчета
            (chemPot,Aff,kineticMatrix,vProcesses,
             balanceMatrix,vx) = CountSystem3(stateCoordinates,_systemParameters,
                                              nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,
                                              nu1_7,nu1_8,nu1_9,nu1_10,nu1_11)
            
            #Выводим результат
            return vx
        def countSystem3_1(t,stateCoordinates):
            return countSystem3(t,stateCoordinates,systemParameters1)
        def countSystem3_2(t,stateCoordinates):
            return countSystem3(t,stateCoordinates,systemParameters2)
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3"]#Имена координат процессов
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        
        #Задаем функцию состояния системы
        stateFunction = CountState3
        
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
        nonEqSystem = NonEqSystem(stateCoordinatesNames,#Имена координат состояния
                                  processCoordinatesNames,#Имена координат процессов
                                  stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                     
                                  #Задаем функцию состояния системы
                                  stateFunction,
                                    
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
        
        #Создаем динамику системы
        nonEqSystemDyn = NonEqSystemDyn(nonEqSystem,#Система
                 
                                        fU3,#Функция условий протекания процессов
                                         
                                        FunChar3,#Функция внешних параметров
                                         
                                        integDynamic#Метод интегрирования дифференциальных уравнений
                                        )
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem3_1,
                                             (0.0, Tint),
                                             stateCoordinates1,#Начальное состояние
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Рассчитываем эталонную динамику параметров системы
        (t1,Char1et1,Char2et1,Char3et1) = FunChar3(t,#Моменты времени
                                                   stateCoordinatesEt,#Координаты состояния
                                                   systemParameters1#Параметры системы  
                                                   )
        (char1et1,char2et1,char3et1) = (np.sum(t1*Char1et1 + Char3et1),np.sum(t1*Char2et1*Char3et1),np.sum(t1*Char2et1))
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem3_2,
                                             (0.0, Tint),
                                             stateCoordinates2,#Начальное состояние
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Рассчитываем эталонную динамику параметров системы
        (t2,Char1et2,Char2et2,Char3et2) = FunChar3(t,#Моменты времени
                                                   stateCoordinatesEt,#Координаты состояния
                                                   systemParameters2#Параметры системы  
                                                   )
        (char1et2,char2et2,char3et2) = (np.sum(t2*Char1et2 + Char3et2),np.sum(t2*Char2et2*Char3et2),np.sum(t2*Char2et2))
        
        #Эталонные результаты
        char1et = np.array([char1et1,char1et2]).reshape(-1,1)
        char2et = np.array([char2et1,char2et2]).reshape(-1,1)
        char3et = np.array([char3et1,char3et2]).reshape(-1,1)
        
        #Функция сжатия
        def CompressFunction(dyn, index):
            #Выводим результат
            return (np.sum(dyn[0]*dyn[1] + dyn[3]),np.sum(dyn[0]*dyn[2]*dyn[3]),np.sum(dyn[0]*dyn[2]))
        
        #Создаем класс вычислительного эксперимента
        cDyn = CountDynamics(nonEqSystemDyn,
                             CompressFunction)
        
        #Функтор модели параметров системы
        class learningSystemPar(object):
            #Функтор параметров
            def ModelSystemParametersCount(self,comExpRez,
                                           cStateCoordinates0,
                                           cSystemParameters):
                self.Par = (comExpRez[1],
                            cStateCoordinates0[0],
                            cSystemParameters[0])
        
        #Создаем класс векторной модели
        funSysPar = learningSystemPar()
        cDynLearn = ModelLearning(funSysPar,cDyn)
        
        #Рассчитываем динамику параметров системы
        stateCoordinates = np.vstack([stateCoordinates1,
                                      stateCoordinates2])
        systemParameters = np.vstack([systemParameters1,
                                      systemParameters2])
        Tints = np.array([Tint,
                          Tint])
        cDynLearn.ModelLearning(Tints,
                                stateCoordinates,
                                systemParameters)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(funSysPar.Par[0] - char2et))
        self.assertAlmostEqual(deltaChar1, 0.0, 6)
        deltaChar2 = np.max(np.abs(funSysPar.Par[1] - stateCoordinates[0]))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(funSysPar.Par[2] - systemParameters[0]))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
