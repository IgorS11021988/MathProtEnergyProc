import numpy as np
from scipy import integrate as integ
from MathProtEnergyProc.IntegrateDyn import standartIntegrateDyn, stepIntegrateDyn

from MathProtEnergyProc import NonEqSystem
from MathProtEnergyProc import NonEqSystemDyn, Model

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


# Модульные тесты
class TestNonEqSystemDyn(unittest.TestCase):
    def setUp(self):
        # Выполнить настройку тестов (если необходимо)
        pass

    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass

    # Модульные тесты
    def testNonEqDyn1(self):
        # Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81])  # Координаты состояния
        systemParameters = np.array([10, 20, 50, 70, 100, 200, 30, 45, 10, 25, 35, 55, 75, 81, 210, 300, 150, 153, 123, 15.3, 45.9, 5.63, 10, 6.5, 45.4, 315, 21.3, 4.5, 1.5, 33.9, 27.3, 4.53, 21.9, 7.5, 53.6, 21.9, 45.3, 33.9, 159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3])  # Параметры системы

        # Стехиометрические коэффициенты химических реакций
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

        # Базовые коэффициенты главной кинетической матрицы процессов]
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3

        # Метод интегрирования дифференциальных уравнений
        method = "RK45"
        integDynamic = standartIntegrateDyn(method=method)

        # Время интегрирования
        Tint = 1e-6

        # Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)

        # Функция эталонного результата
        def countSystem1(t, stateCoordinates):
            # Вызов функции времени
            _systemParameters = fU1(t, systemParameters)

            # Вызов расчета
            (chemPot, Aff, kineticMatrix,
             vProcesses, balanceMatrix, vx) = CountSystem1(stateCoordinates, _systemParameters,
                                                           nu1_1, nu1_2, nu1_3, nu1_4, nu1_5, nu1_6, nu1_7,
                                                           nu1_8, nu1_9, nu1_10, nu1_11, nu2_1, nu2_2, nu2_3,
                                                           nu2_4, nu2_5, nu2_6, Adiff1_1, Adiff1_2, Adiff2_1)

            # Выводим результат
            return vx

        # Задаем структуру системы
        stateCoordinatesNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]  # Имена координат состояния
        processCoordinatesNames = ["vChem1_1", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff1", "vDiff2"]  # Имена координат процессов
        stateCoordinatesStreamsNames = []  # Имена координат состояния, изменяемых в результате внешних потоков

        # Задаем функцию состояния системы
        stateFunction = CountState1

        # Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]  # Имена переменных потенциалов взаимодействия по координатам состояния

        # Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1", "vChem1_1", "vChem1_2", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff2"]  # Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1", "vChem1_2", "vChem1_1", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff2"]  # Имена сопряженностей между собой термодинамических сил

        # Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []  # Имена переменных внешних потоков

        # Задаем систему
        nonEqSystem = NonEqSystem(stateCoordinatesNames,  # Имена координат состояния
                                  processCoordinatesNames,  # Имена координат процессов
                                  stateCoordinatesStreamsNames,  # Имена координат состояния, изменяемых в результате внешних потоков

                                  # Задаем функцию состояния системы
                                  stateFunction,

                                  # Задаем рассчитываемые параметры системы
                                  stateCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам состояния
                                  processCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам процессов
                                  stateCoordinatesVarPotentialsInterNames,  # Имена переменных потенциалов взаимодействия по координатам состояния

                                  # Задаем переменные параметры кинетической матрицы
                                  varKineticNames,  # Имена сопряженностей между собой координат процессов
                                  varKineticAffNames,  # Имена сопряженностей между собой термодинамических сил

                                  # Задаем внешние потоки
                                  stateCoordinatesVarStreamsNames  # Имена переменных внешних потоков
                                  )

        # Задаем постоянные параметры системы
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

        # Создаем динамику системы
        nonEqSystemDyn = NonEqSystemDyn(nonEqSystem,  # Система

                                        fU1,  # Функция условий протекания процессов

                                        FunChar1,  # Функция внешних параметров

                                        integDynamic  # Метод интегрирования дифференциальных уравнений
                                        )

        # Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystem(0.0 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystem(0.3 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystem(0.9 * Tint,
                                          stateCoordinates,
                                          systemParameters)

        # Выполняем эталонные расчеты
        vxEt1 = countSystem1(0.0 * Tint,
                             stateCoordinates)
        vxEt2 = countSystem1(0.3 * Tint,
                             stateCoordinates)
        vxEt3 = countSystem1(0.9 * Tint,
                             stateCoordinates)

        # Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 9)  # Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 9)  # Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 9)  # Проверяем скорость

        # Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem1,
                                             (0.0, Tint),
                                             stateCoordinates,  # Начальное состояние
                                             t_eval=t_eval,  # Моменты времени
                                             method=method  # Метод численного интегрирования
                                             )  # Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1, 1)  # Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()  # Координаты состояния

        # Рассчитываем эталонную динамику параметров системы
        (t, Char1et, Char2et) = FunChar1(t,  # Моменты времени
                                         stateCoordinatesEt,  # Координаты состояния
                                         systemParameters  # Параметры системы
                                         )

        # Рассчитываем динамику параметров системы
        (t, Char1, Char2) = nonEqSystemDyn.NonEqSystemDynamic(Tint,
                                                              stateCoordinates,
                                                              systemParameters,
                                                              t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)

        # Функтор модели параметров системы
        class functorSystemPar(object):
            # Функтор параметров
            def SystemParametersCount(self, measParameters):
                return (measParameters[0:10],
                        measParameters[10::])

        # Получаем класс модели
        funSystemPar = functorSystemPar()
        model = Model(funSystemPar,
                      nonEqSystemDyn)

        # Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              systemParameters))

        # Запускаем динамику
        (t, Char1, Char2) = model.CountDynamic(Tint, measChar,
                                               t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)

    def testNonEqDyn2(self):
        # Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81])  # Координаты состояния
        systemParameters = np.array([10, 20, 50, 70, 100, 200, 30, 45, 10, 25, 35, 55, 75, 81, 210, 300, 150, 153, 123, 15.3, 45.9, 5.63, 10, 6.5, 45.4, 315, 21.3, 4.5, 1.5, 33.9, 27.3, 4.53, 21.9, 7.5, 53.6, 21.9, 45.3, 33.9, 159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3])  # Параметры системы

        # Стехиометрические коэффициенты химических реакций
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

        # Базовые коэффициенты главной кинетической матрицы процессов]
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3

        # Метод интегрирования дифференциальных уравнений
        method = "RK45"
        atol = 2 * 1e-7  # Абсолютная погрешность
        rtol = 3 * 1e-6  # Относительная погрешность
        integDynamic = standartIntegrateDyn(method=method,
                                            atol=atol,
                                            rtol=rtol)

        # Время интегрирования
        Tint = 1e-6

        # Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)

        # Функция эталонного результата
        def countSystem1(t, stateCoordinates):
            # Вызов функции времени
            _systemParameters = fU1(t, systemParameters)

            # Вызов расчета
            (chemPot, Aff, kineticMatrix,
             vProcesses, balanceMatrix, vx) = CountSystem1(stateCoordinates, _systemParameters,
                                                           nu1_1, nu1_2, nu1_3, nu1_4, nu1_5, nu1_6, nu1_7,
                                                           nu1_8, nu1_9, nu1_10, nu1_11, nu2_1, nu2_2, nu2_3,
                                                           nu2_4, nu2_5, nu2_6, Adiff1_1, Adiff1_2, Adiff2_1)

            # Выводим результат
            return vx

        # Задаем структуру системы
        stateCoordinatesNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]  # Имена координат состояния
        processCoordinatesNames = ["vChem1_1", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff1", "vDiff2"]  # Имена координат процессов
        stateCoordinatesStreamsNames = []  # Имена координат состояния, изменяемых в результате внешних потоков

        # Задаем функцию состояния системы
        stateFunction = CountState1

        # Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]  # Имена переменных потенциалов взаимодействия по координатам состояния

        # Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1", "vChem1_1", "vChem1_2", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff2"]  # Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1", "vChem1_2", "vChem1_1", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff2"]  # Имена сопряженностей между собой термодинамических сил

        # Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []  # Имена переменных внешних потоков

        # Задаем систему
        nonEqSystem = NonEqSystem(stateCoordinatesNames,  # Имена координат состояния
                                  processCoordinatesNames,  # Имена координат процессов
                                  stateCoordinatesStreamsNames,  # Имена координат состояния, изменяемых в результате внешних потоков

                                  # Задаем функцию состояния системы
                                  stateFunction,

                                  # Задаем рассчитываемые параметры системы
                                  stateCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам состояния
                                  processCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам процессов
                                  stateCoordinatesVarPotentialsInterNames,  # Имена переменных потенциалов взаимодействия по координатам состояния

                                  # Задаем переменные параметры кинетической матрицы
                                  varKineticNames,  # Имена сопряженностей между собой координат процессов
                                  varKineticAffNames,  # Имена сопряженностей между собой термодинамических сил

                                  # Задаем внешние потоки
                                  stateCoordinatesVarStreamsNames  # Имена переменных внешних потоков
                                  )

        # Задаем постоянные параметры системы
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

        # Создаем динамику системы
        nonEqSystemDyn = NonEqSystemDyn(nonEqSystem,  # Система

                                        fU1,  # Функция условий протекания процессов

                                        FunChar1,  # Функция внешних параметров

                                        integDynamic  # Метод интегрирования дифференциальных уравнений
                                        )

        # Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystem(0.0 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystem(0.3 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystem(0.9 * Tint,
                                          stateCoordinates,
                                          systemParameters)

        # Выполняем эталонные расчеты
        vxEt1 = countSystem1(0.0 * Tint,
                             stateCoordinates)
        vxEt2 = countSystem1(0.3 * Tint,
                             stateCoordinates)
        vxEt3 = countSystem1(0.9 * Tint,
                             stateCoordinates)

        # Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 9)  # Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 9)  # Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 9)  # Проверяем скорость

        # Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem1,
                                             (0.0, Tint),
                                             stateCoordinates,  # Начальное состояние
                                             t_eval=t_eval,  # Моменты времени
                                             atol=atol,
                                             rtol=rtol,
                                             method=method  # Метод численного интегрирования
                                             )  # Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1, 1)  # Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()  # Координаты состояния

        # Рассчитываем эталонную динамику параметров системы
        (t, Char1et, Char2et) = FunChar1(t,  # Моменты времени
                                         stateCoordinatesEt,  # Координаты состояния
                                         systemParameters  # Параметры системы
                                         )

        # Рассчитываем динамику параметров системы
        (t, Char1, Char2) = nonEqSystemDyn.NonEqSystemDynamic(Tint,
                                                              stateCoordinates,
                                                              systemParameters,
                                                              t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)

        # Функтор модели параметров системы
        class functorSystemPar(object):
            # Функтор параметров
            def SystemParametersCount(self, measParameters):
                return (measParameters[0:10],
                        measParameters[10::])

        # Получаем класс модели
        funSystemPar = functorSystemPar()
        model = Model(funSystemPar,
                      nonEqSystemDyn)

        # Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              systemParameters))

        # Запускаем динамику
        (t, Char1, Char2) = model.CountDynamic(Tint, measChar,
                                               t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)

    def testNonEqDyn3(self):
        # Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81])  # Координаты состояния
        systemParameters = np.array([10, 20, 50, 70, 100, 200, 30, 45, 10, 25, 35, 55, 75, 81, 210, 300, 150, 153, 123, 15.3, 45.9, 5.63, 10, 6.5, 45.4, 315, 21.3, 4.5, 1.5, 33.9, 27.3, 4.53, 21.9, 7.5, 53.6, 21.9, 45.3, 33.9, 159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3])  # Параметры системы

        # Стехиометрические коэффициенты химических реакций
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

        # Базовые коэффициенты главной кинетической матрицы процессов]
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3

        # Метод интегрирования дифференциальных уравнений
        method = "RK45"

        # Аттрибуты интегрирования
        def funIntegAttributes(systemParameters, method):
            # Определяем аттрибуты
            atol = np.abs(systemParameters[3]) * 1e-7  # Абсолютная погрешность
            rtol = np.abs(systemParameters[9]) * 1e-4  # Относительная погрешность
            max_step = np.abs(2 * systemParameters[21] + 3 * systemParameters[15]) * 1e-1  # Максимальный шаг

            # Выводим аттрибуты интегрирования
            return (atol, rtol, max_step)
        integDynamic = stepIntegrateDyn(method=method,
                                        funIntegAttributes=funIntegAttributes)

        # Время интегрирования
        Tint = 1e-6

        # Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)

        # Функция эталонного результата
        def countSystem1(t, stateCoordinates):
            # Вызов функции времени
            _systemParameters = fU1(t, systemParameters)

            # Вызов расчета
            (chemPot, Aff, kineticMatrix,
             vProcesses, balanceMatrix, vx) = CountSystem1(stateCoordinates, _systemParameters,
                                                           nu1_1, nu1_2, nu1_3, nu1_4, nu1_5, nu1_6, nu1_7,
                                                           nu1_8, nu1_9, nu1_10, nu1_11, nu2_1, nu2_2, nu2_3,
                                                           nu2_4, nu2_5, nu2_6, Adiff1_1, Adiff1_2, Adiff2_1)

            # Выводим результат
            return vx

        # Задаем структуру системы
        stateCoordinatesNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]  # Имена координат состояния
        processCoordinatesNames = ["vChem1_1", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff1", "vDiff2"]  # Имена координат процессов
        stateCoordinatesStreamsNames = []  # Имена координат состояния, изменяемых в результате внешних потоков

        # Задаем функцию состояния системы
        stateFunction = CountState1

        # Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]  # Имена переменных потенциалов взаимодействия по координатам состояния

        # Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1", "vChem1_1", "vChem1_2", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff2"]  # Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1", "vChem1_2", "vChem1_1", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff2"]  # Имена сопряженностей между собой термодинамических сил

        # Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []  # Имена переменных внешних потоков

        # Задаем систему
        nonEqSystem = NonEqSystem(stateCoordinatesNames,  # Имена координат состояния
                                  processCoordinatesNames,  # Имена координат процессов
                                  stateCoordinatesStreamsNames,  # Имена координат состояния, изменяемых в результате внешних потоков

                                  # Задаем функцию состояния системы
                                  stateFunction,

                                  # Задаем рассчитываемые параметры системы
                                  stateCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам состояния
                                  processCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам процессов
                                  stateCoordinatesVarPotentialsInterNames,  # Имена переменных потенциалов взаимодействия по координатам состояния

                                  # Задаем переменные параметры кинетической матрицы
                                  varKineticNames,  # Имена сопряженностей между собой координат процессов
                                  varKineticAffNames,  # Имена сопряженностей между собой термодинамических сил

                                  # Задаем внешние потоки
                                  stateCoordinatesVarStreamsNames  # Имена переменных внешних потоков
                                  )

        # Задаем постоянные параметры системы
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

        # Создаем динамику системы
        nonEqSystemDyn = NonEqSystemDyn(nonEqSystem,  # Система

                                        fU1,  # Функция условий протекания процессов

                                        FunChar1,  # Функция внешних параметров

                                        integDynamic  # Метод интегрирования дифференциальных уравнений
                                        )

        # Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystem(0.0 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystem(0.3 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystem(0.9 * Tint,
                                          stateCoordinates,
                                          systemParameters)

        # Выполняем эталонные расчеты
        vxEt1 = countSystem1(0.0 * Tint,
                             stateCoordinates)
        vxEt2 = countSystem1(0.3 * Tint,
                             stateCoordinates)
        vxEt3 = countSystem1(0.9 * Tint,
                             stateCoordinates)

        # Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 9)  # Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 9)  # Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 9)  # Проверяем скорость

        # Рассчитываем эталонную динамику системы
        (atol, rtol, max_step) = funIntegAttributes(systemParameters, method)
        stateCoordinatesEt = integ.solve_ivp(countSystem1,
                                             (0.0, Tint),
                                             stateCoordinates,  # Начальное состояние
                                             t_eval=t_eval,  # Моменты времени
                                             atol=atol,
                                             rtol=rtol,
                                             max_step=max_step,
                                             method=method  # Метод численного интегрирования
                                             )  # Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1, 1)  # Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()  # Координаты состояния

        # Рассчитываем эталонную динамику параметров системы
        (t, Char1et, Char2et) = FunChar1(t,  # Моменты времени
                                         stateCoordinatesEt,  # Координаты состояния
                                         systemParameters  # Параметры системы
                                         )

        # Рассчитываем динамику параметров системы
        (t, Char1, Char2) = nonEqSystemDyn.NonEqSystemDynamic(Tint,
                                                              stateCoordinates,
                                                              systemParameters,
                                                              t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)

        # Функтор модели параметров системы
        class functorSystemPar(object):
            # Функтор параметров
            def SystemParametersCount(self, measParameters):
                return (measParameters[0:10],
                        measParameters[10::])

        # Получаем класс модели
        funSystemPar = functorSystemPar()
        model = Model(funSystemPar,
                      nonEqSystemDyn)

        # Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              systemParameters))

        # Запускаем динамику
        (t, Char1, Char2) = model.CountDynamic(Tint, measChar,
                                               t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)

    def testNonEqDyn4(self):
        # Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81, 80, 64])  # Координаты состояния
        systemParameters = np.array([10, 20, 50, 70, 100, 200, 30, 45, 10, 25, 35, 55, 75, 81, 210, 300, 150, 153, 123, 15.3, 45.9, 5.63, 10, 6.5, 45.4, 315, 21.3, 4.5, 1.5, 33.9, 27.3, 4.53, 21.9, 7.5, 53.6, 21.9, 45.3, 33.9, 159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3, 0.15, 74.1])  # Параметры системы

        # Стехиометрические коэффициенты химических реакций
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

        # Базовые коэффициенты главной кинетической матрицы процессов]
        AChem1_4_4 = 3.3
        AChem2_3_3 = 6.3
        AChem2_1_2 = 0.3
        AChem2_2_1 = 0.15
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3
        Adiff3_3 = 4.83

        # Внешние потоки вещества
        xExt1_6 = 30.9

        # Потенциал взаимодействия
        mu2_9 = 13.5

        # Метод интегрирования дифференциальных уравнений
        method = "RK23"
        integDynamic = standartIntegrateDyn(method=method)

        # Время интегрирования
        Tint = 1e-6

        # Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)

        # Функция эталонного результата
        def countSystem2(t, stateCoordinates):
            # Вызов функции времени
            _systemParameters = fU2(t, systemParameters)

            # Вызов расчета
            (chemPot, Aff, kineticMatrix, vProcesses,
             balanceMatrix, vx, Streams) = CountSystem2(stateCoordinates, _systemParameters,
                                                        nu1_1, nu1_2, nu1_3, nu1_4, nu1_5, nu1_6, nu1_7,
                                                        nu1_8, nu1_9, nu1_10, nu1_11, nu1_12, nu1_13, nu1_14,
                                                        nu1_15, nu2_1, nu2_2, nu2_3, nu2_4, nu2_5, nu2_6,
                                                        nu2_7, nu2_8, nu2_9, Adiff1_1, Adiff1_2, Adiff2_1,
                                                        AChem1_4_4, AChem2_3_3, AChem2_1_2, AChem2_2_1,
                                                        Adiff3_3, xExt1_6, mu2_9)

            # Выводим результат
            return vx

        # Задаем структуру системы
        stateCoordinatesNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8", "x2_9"]  # Имена координат состояния
        processCoordinatesNames = ["vChem1_1", "vChem1_2", "vChem1_3", "vChem1_4", "vChem2_1", "vChem2_2", "vChem2_3", "vDiff1", "vDiff2", "vDiff3"]  # Имена координат процессов
        stateCoordinatesStreamsNames = ["x1_2", "x1_6", "x2_8"]  # Имена координат состояния, изменяемых в результате внешних потоков

        # Задаем функцию состояния системы
        stateFunction = CountState2

        # Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8"]  # Имена переменных потенциалов взаимодействия по координатам состояния

        # Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1", "vChem1_1", "vChem1_2", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff2"]  # Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1", "vChem1_2", "vChem1_1", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff2"]  # Имена сопряженностей между собой термодинамических сил

        # Задаем внешние потоки
        stateCoordinatesVarStreamsNames = ["x1_2", "x2_8"]  # Имена переменных внешних потоков

        # Задаем систему
        nonEqSystem = NonEqSystem(stateCoordinatesNames,  # Имена координат состояния
                                  processCoordinatesNames,  # Имена координат процессов
                                  stateCoordinatesStreamsNames,  # Имена координат состояния, изменяемых в результате внешних потоков

                                  # Задаем функцию состояния системы
                                  stateFunction,

                                  # Задаем рассчитываемые параметры системы
                                  stateCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам состояния
                                  processCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам процессов
                                  stateCoordinatesVarPotentialsInterNames,  # Имена переменных потенциалов взаимодействия по координатам состояния

                                  # Задаем переменные параметры кинетической матрицы
                                  varKineticNames,  # Имена сопряженностей между собой координат процессов
                                  varKineticAffNames,  # Имена сопряженностей между собой термодинамических сил

                                  # Задаем внешние потоки
                                  stateCoordinatesVarStreamsNames  # Имена переменных внешних потоков
                                  )

        # Задаем постоянные параметры системы
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
        nonEqSystem.SetStateCoordinatesStreamsConstElement("x1_6", xExt1_6)
        nonEqSystem.SetPotentialsInterConstElement("EnPow2", "x2_9", -mu2_9)

        # Создаем динамику системы
        nonEqSystemDyn = NonEqSystemDyn(nonEqSystem,  # Система

                                        fU2,  # Функция условий протекания процессов

                                        FunChar2,  # Функция внешних параметров

                                        integDynamic  # Метод интегрирования дифференциальных уравнений
                                        )

        # Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystem(0.0 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystem(0.3 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystem(0.9 * Tint,
                                          stateCoordinates,
                                          systemParameters)

        # Выполняем эталонные расчеты
        vxEt1 = countSystem2(0.0 * Tint,
                             stateCoordinates)
        vxEt2 = countSystem2(0.3 * Tint,
                             stateCoordinates)
        vxEt3 = countSystem2(0.9 * Tint,
                             stateCoordinates)

        # Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 9)  # Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 8)  # Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 8)  # Проверяем скорость

        # Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem2,
                                             (0.0, Tint),
                                             stateCoordinates,  # Начальное состояние
                                             t_eval=t_eval,  # Моменты времени
                                             method=method  # Метод численного интегрирования
                                             )  # Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1, 1)  # Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()  # Координаты состояния

        # Рассчитываем эталонную динамику параметров системы
        (t, Char1et, Char2et, Char3et) = FunChar2(t,  # Моменты времени
                                                  stateCoordinatesEt,  # Координаты состояния
                                                  systemParameters  # Параметры системы
                                                  )

        # Рассчитываем динамику параметров системы
        (t, Char1, Char2, Char3) = nonEqSystemDyn.NonEqSystemDynamic(Tint,
                                                                     stateCoordinates,
                                                                     systemParameters,
                                                                     t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)

        # Функтор модели параметров системы
        class functorSystemPar(object):
            # Функтор параметров
            def SystemParametersCount(self, measParameters):
                return (measParameters[0:12],
                        measParameters[12::])

        # Получаем класс модели
        funSystemPar = functorSystemPar()
        model = Model(funSystemPar,
                      nonEqSystemDyn)

        # Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              systemParameters))

        # Запускаем динамику
        (t, Char1, Char2, Char3) = model.CountDynamic(Tint, measChar,
                                                      t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)

    def testNonEqDyn5(self):
        # Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81, 80, 64])  # Координаты состояния
        systemParameters = np.array([10, 20, 50, 70, 100, 200, 30, 45, 10, 25, 35, 55, 75, 81, 210, 300, 150, 153, 123, 15.3, 45.9, 5.63, 10, 6.5, 45.4, 315, 21.3, 4.5, 1.5, 33.9, 27.3, 4.53, 21.9, 7.5, 53.6, 21.9, 45.3, 33.9, 159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3, 0.15, 74.1])  # Параметры системы

        # Стехиометрические коэффициенты химических реакций
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

        # Базовые коэффициенты главной кинетической матрицы процессов]
        AChem1_4_4 = 3.3
        AChem2_3_3 = 6.3
        AChem2_1_2 = 0.3
        AChem2_2_1 = 0.15
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3
        Adiff3_3 = 4.83

        # Внешние потоки вещества
        xExt1_6 = 30.9

        # Потенциал взаимодействия
        mu2_9 = 13.5

        # Метод интегрирования дифференциальных уравнений
        method = "RK23"
        rtol = 2.7 * 1e-2  # Относительная погрешность
        max_step = 3 * 1e-4  # Максимальный шаг интегрирования
        integDynamic = standartIntegrateDyn(method=method,
                                            rtol=rtol,
                                            max_step=max_step)

        # Время интегрирования
        Tint = 1e-6

        # Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)

        # Функция эталонного результата
        def countSystem2(t, stateCoordinates):
            # Вызов функции времени
            _systemParameters = fU2(t, systemParameters)

            # Вызов расчета
            (chemPot, Aff, kineticMatrix, vProcesses,
             balanceMatrix, vx, Streams) = CountSystem2(stateCoordinates, _systemParameters,
                                                        nu1_1, nu1_2, nu1_3, nu1_4, nu1_5, nu1_6, nu1_7,
                                                        nu1_8, nu1_9, nu1_10, nu1_11, nu1_12, nu1_13, nu1_14,
                                                        nu1_15, nu2_1, nu2_2, nu2_3, nu2_4, nu2_5, nu2_6,
                                                        nu2_7, nu2_8, nu2_9, Adiff1_1, Adiff1_2, Adiff2_1,
                                                        AChem1_4_4, AChem2_3_3, AChem2_1_2, AChem2_2_1,
                                                        Adiff3_3, xExt1_6, mu2_9)

            # Выводим результат
            return vx

        # Задаем структуру системы
        stateCoordinatesNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8", "x2_9"]  # Имена координат состояния
        processCoordinatesNames = ["vChem1_1", "vChem1_2", "vChem1_3", "vChem1_4", "vChem2_1", "vChem2_2", "vChem2_3", "vDiff1", "vDiff2", "vDiff3"]  # Имена координат процессов
        stateCoordinatesStreamsNames = ["x1_2", "x1_6", "x2_8"]  # Имена координат состояния, изменяемых в результате внешних потоков

        # Задаем функцию состояния системы
        stateFunction = CountState2

        # Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8"]  # Имена переменных потенциалов взаимодействия по координатам состояния

        # Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1", "vChem1_1", "vChem1_2", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff2"]  # Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1", "vChem1_2", "vChem1_1", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff2"]  # Имена сопряженностей между собой термодинамических сил

        # Задаем внешние потоки
        stateCoordinatesVarStreamsNames = ["x1_2", "x2_8"]  # Имена переменных внешних потоков

        # Задаем систему
        nonEqSystem = NonEqSystem(stateCoordinatesNames,  # Имена координат состояния
                                  processCoordinatesNames,  # Имена координат процессов
                                  stateCoordinatesStreamsNames,  # Имена координат состояния, изменяемых в результате внешних потоков

                                  # Задаем функцию состояния системы
                                  stateFunction,

                                  # Задаем рассчитываемые параметры системы
                                  stateCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам состояния
                                  processCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам процессов
                                  stateCoordinatesVarPotentialsInterNames,  # Имена переменных потенциалов взаимодействия по координатам состояния

                                  # Задаем переменные параметры кинетической матрицы
                                  varKineticNames,  # Имена сопряженностей между собой координат процессов
                                  varKineticAffNames,  # Имена сопряженностей между собой термодинамических сил

                                  # Задаем внешние потоки
                                  stateCoordinatesVarStreamsNames  # Имена переменных внешних потоков
                                  )

        # Задаем постоянные параметры системы
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
        nonEqSystem.SetStateCoordinatesStreamsConstElement("x1_6", xExt1_6)
        nonEqSystem.SetPotentialsInterConstElement("EnPow2", "x2_9", -mu2_9)

        # Создаем динамику системы
        nonEqSystemDyn = NonEqSystemDyn(nonEqSystem,  # Система

                                        fU2,  # Функция условий протекания процессов

                                        FunChar2,  # Функция внешних параметров

                                        integDynamic  # Метод интегрирования дифференциальных уравнений
                                        )

        # Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystem(0.0 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystem(0.3 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystem(0.9 * Tint,
                                          stateCoordinates,
                                          systemParameters)

        # Выполняем эталонные расчеты
        vxEt1 = countSystem2(0.0 * Tint,
                             stateCoordinates)
        vxEt2 = countSystem2(0.3 * Tint,
                             stateCoordinates)
        vxEt3 = countSystem2(0.9 * Tint,
                             stateCoordinates)

        # Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 9)  # Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 8)  # Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 8)  # Проверяем скорость

        # Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem2,
                                             (0.0, Tint),
                                             stateCoordinates,  # Начальное состояние
                                             t_eval=t_eval,  # Моменты времени
                                             rtol=rtol,
                                             max_step=max_step,
                                             method=method  # Метод численного интегрирования
                                             )  # Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1, 1)  # Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()  # Координаты состояния

        # Рассчитываем эталонную динамику параметров системы
        (t, Char1et, Char2et, Char3et) = FunChar2(t,  # Моменты времени
                                                  stateCoordinatesEt,  # Координаты состояния
                                                  systemParameters  # Параметры системы
                                                  )

        # Рассчитываем динамику параметров системы
        (t, Char1, Char2, Char3) = nonEqSystemDyn.NonEqSystemDynamic(Tint,
                                                                     stateCoordinates,
                                                                     systemParameters,
                                                                     t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)

        # Функтор модели параметров системы
        class functorSystemPar(object):
            # Функтор параметров
            def SystemParametersCount(self, measParameters):
                return (measParameters[0:12],
                        measParameters[12::])

        # Получаем класс модели
        funSystemPar = functorSystemPar()
        model = Model(funSystemPar,
                      nonEqSystemDyn)

        # Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              systemParameters))

        # Запускаем динамику
        (t, Char1, Char2, Char3) = model.CountDynamic(Tint, measChar,
                                                      t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)

    def testNonEqDyn6(self):
        # Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81, 80, 64])  # Координаты состояния
        systemParameters = np.array([10, 20, 50, 70, 100, 200, 30, 45, 10, 25, 35, 55, 75, 81, 210, 300, 150, 153, 123, 15.3, 45.9, 5.63, 10, 6.5, 45.4, 315, 21.3, 4.5, 1.5, 33.9, 27.3, 4.53, 21.9, 7.5, 53.6, 21.9, 45.3, 33.9, 159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3, 0.15, 74.1])  # Параметры системы

        # Стехиометрические коэффициенты химических реакций
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

        # Базовые коэффициенты главной кинетической матрицы процессов]
        AChem1_4_4 = 3.3
        AChem2_3_3 = 6.3
        AChem2_1_2 = 0.3
        AChem2_2_1 = 0.15
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3
        Adiff3_3 = 4.83

        # Внешние потоки вещества
        xExt1_6 = 30.9

        # Потенциал взаимодействия
        mu2_9 = 13.5

        # Метод интегрирования дифференциальных уравнений
        method = "RK23"

        # Функция аттрибутов интегрирования
        def funIntegAttributes(systemParameters, method):
            # Определяем аттрибуты
            atol = np.abs(systemParameters[4]) * 1e-8  # Абсолютная погрешность
            rtol = np.abs(systemParameters[14] + 2.1 * systemParameters[18]) * 1e-5  # Относительная погрешность
            max_step = np.abs(5.1 * systemParameters[33] * systemParameters[12]) * 1e-3  # Максимальный шаг

            # Выводим аттрибуты интегрирования
            return (atol, rtol, max_step)
        integDynamic = stepIntegrateDyn(method=method,
                                        funIntegAttributes=funIntegAttributes)

        # Время интегрирования
        Tint = 1e-6

        # Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)

        # Функция эталонного результата
        def countSystem2(t, stateCoordinates):
            # Вызов функции времени
            _systemParameters = fU2(t, systemParameters)

            # Вызов расчета
            (chemPot, Aff, kineticMatrix, vProcesses,
             balanceMatrix, vx, Streams) = CountSystem2(stateCoordinates, _systemParameters,
                                                        nu1_1, nu1_2, nu1_3, nu1_4, nu1_5, nu1_6, nu1_7,
                                                        nu1_8, nu1_9, nu1_10, nu1_11, nu1_12, nu1_13, nu1_14,
                                                        nu1_15, nu2_1, nu2_2, nu2_3, nu2_4, nu2_5, nu2_6,
                                                        nu2_7, nu2_8, nu2_9, Adiff1_1, Adiff1_2, Adiff2_1,
                                                        AChem1_4_4, AChem2_3_3, AChem2_1_2, AChem2_2_1,
                                                        Adiff3_3, xExt1_6, mu2_9)

            # Выводим результат
            return vx

        # Задаем структуру системы
        stateCoordinatesNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8", "x2_9"]  # Имена координат состояния
        processCoordinatesNames = ["vChem1_1", "vChem1_2", "vChem1_3", "vChem1_4", "vChem2_1", "vChem2_2", "vChem2_3", "vDiff1", "vDiff2", "vDiff3"]  # Имена координат процессов
        stateCoordinatesStreamsNames = ["x1_2", "x1_6", "x2_8"]  # Имена координат состояния, изменяемых в результате внешних потоков

        # Задаем функцию состояния системы
        stateFunction = CountState2

        # Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8"]  # Имена переменных потенциалов взаимодействия по координатам состояния

        # Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1", "vChem1_1", "vChem1_2", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff2"]  # Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1", "vChem1_2", "vChem1_1", "vChem1_2", "vChem1_3", "vChem2_1", "vChem2_2", "vDiff2"]  # Имена сопряженностей между собой термодинамических сил

        # Задаем внешние потоки
        stateCoordinatesVarStreamsNames = ["x1_2", "x2_8"]  # Имена переменных внешних потоков

        # Задаем систему
        nonEqSystem = NonEqSystem(stateCoordinatesNames,  # Имена координат состояния
                                  processCoordinatesNames,  # Имена координат процессов
                                  stateCoordinatesStreamsNames,  # Имена координат состояния, изменяемых в результате внешних потоков

                                  # Задаем функцию состояния системы
                                  stateFunction,

                                  # Задаем рассчитываемые параметры системы
                                  stateCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам состояния
                                  processCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам процессов
                                  stateCoordinatesVarPotentialsInterNames,  # Имена переменных потенциалов взаимодействия по координатам состояния

                                  # Задаем переменные параметры кинетической матрицы
                                  varKineticNames,  # Имена сопряженностей между собой координат процессов
                                  varKineticAffNames,  # Имена сопряженностей между собой термодинамических сил

                                  # Задаем внешние потоки
                                  stateCoordinatesVarStreamsNames  # Имена переменных внешних потоков
                                  )

        # Задаем постоянные параметры системы
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
        nonEqSystem.SetStateCoordinatesStreamsConstElement("x1_6", xExt1_6)
        nonEqSystem.SetPotentialsInterConstElement("EnPow2", "x2_9", -mu2_9)

        # Создаем динамику системы
        nonEqSystemDyn = NonEqSystemDyn(nonEqSystem,  # Система

                                        fU2,  # Функция условий протекания процессов

                                        FunChar2,  # Функция внешних параметров

                                        integDynamic  # Метод интегрирования дифференциальных уравнений
                                        )

        # Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystem(0.0 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystem(0.3 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystem(0.9 * Tint,
                                          stateCoordinates,
                                          systemParameters)

        # Выполняем эталонные расчеты
        vxEt1 = countSystem2(0.0 * Tint,
                             stateCoordinates)
        vxEt2 = countSystem2(0.3 * Tint,
                             stateCoordinates)
        vxEt3 = countSystem2(0.9 * Tint,
                             stateCoordinates)

        # Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 9)  # Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 8)  # Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 8)  # Проверяем скорость

        # Рассчитываем эталонную динамику системы
        (atol, rtol, max_step) = funIntegAttributes(systemParameters, method)
        stateCoordinatesEt = integ.solve_ivp(countSystem2,
                                             (0.0, Tint),
                                             stateCoordinates,  # Начальное состояние
                                             t_eval=t_eval,  # Моменты
                                             atol=atol,
                                             rtol=rtol,
                                             max_step=max_step,
                                             method=method  # Метод численного интегрирования
                                             )  # Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1, 1)  # Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()  # Координаты состояния

        # Рассчитываем эталонную динамику параметров системы
        (t, Char1et, Char2et, Char3et) = FunChar2(t,  # Моменты времени
                                                  stateCoordinatesEt,  # Координаты состояния
                                                  systemParameters  # Параметры системы
                                                  )

        # Рассчитываем динамику параметров системы
        (t, Char1, Char2, Char3) = nonEqSystemDyn.NonEqSystemDynamic(Tint,
                                                                     stateCoordinates,
                                                                     systemParameters,
                                                                     t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)

        # Функтор модели параметров системы
        class functorSystemPar(object):
            # Функтор параметров
            def SystemParametersCount(self, measParameters):
                return (measParameters[0:12],
                        measParameters[12::])

        # Получаем класс модели
        funSystemPar = functorSystemPar()
        model = Model(funSystemPar,
                      nonEqSystemDyn)

        # Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              systemParameters))

        # Запускаем динамику
        (t, Char1, Char2, Char3) = model.CountDynamic(Tint, measChar,
                                                      t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)

    def testNonEqDyn7(self):
        # Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213])  # Координаты состояния
        systemParameters = np.array([10, 20, 50, 70, 100, 200, 30, 45, 10, 25, 35, 55, 75, 81, 210, 300, 150, 153, 123, 15.3, 45.9, 5.63, 10, 6.5, 45.4, 315, 161])  # Параметры системы

        # Стехиометрические коэффициенты химических реакций
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

        # Метод интегрирования дифференциальных уравнений
        method = "LSODA"
        integDynamic = standartIntegrateDyn(method=method)

        # Время интегрирования
        Tint = 1e-6

        # Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)

        # Функция эталонного результата
        def countSystem3(t, stateCoordinates):
            # Вызов функции времени
            _systemParameters = fU3(t, systemParameters)

            # Вызов расчета
            (chemPot, Aff, kineticMatrix, vProcesses,
             balanceMatrix, vx) = CountSystem3(stateCoordinates, _systemParameters,
                                               nu1_1, nu1_2, nu1_3, nu1_4, nu1_5, nu1_6,
                                               nu1_7, nu1_8, nu1_9, nu1_10, nu1_11)

            # Выводим результат
            return vx

        # Задаем структуру системы
        stateCoordinatesNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]  # Имена координат состояния
        processCoordinatesNames = ["vChem1_1", "vChem1_2", "vChem1_3"]  # Имена координат процессов
        stateCoordinatesStreamsNames = []  # Имена координат состояния, изменяемых в результате внешних потоков

        # Задаем функцию состояния системы
        stateFunction = CountState3

        # Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]  # Имена переменных потенциалов взаимодействия по координатам состояния

        # Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1", "vChem1_1", "vChem1_2", "vChem1_2", "vChem1_3"]  # Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1", "vChem1_2", "vChem1_1", "vChem1_2", "vChem1_3"]  # Имена сопряженностей между собой термодинамических сил

        # Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []  # Имена переменных внешних потоков

        # Задаем систему
        nonEqSystem = NonEqSystem(stateCoordinatesNames,  # Имена координат состояния
                                  processCoordinatesNames,  # Имена координат процессов
                                  stateCoordinatesStreamsNames,  # Имена координат состояния, изменяемых в результате внешних потоков

                                  # Задаем функцию состояния системы
                                  stateFunction,

                                  # Задаем рассчитываемые параметры системы
                                  stateCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам состояния
                                  processCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам процессов
                                  stateCoordinatesVarPotentialsInterNames,  # Имена переменных потенциалов взаимодействия по координатам состояния

                                  # Задаем переменные параметры кинетической матрицы
                                  varKineticNames,  # Имена сопряженностей между собой координат процессов
                                  varKineticAffNames,  # Имена сопряженностей между собой термодинамических сил

                                  # Задаем внешние потоки
                                  stateCoordinatesVarStreamsNames  # Имена переменных внешних потоков
                                  )

        # Задаем постоянные параметры системы
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

        # Создаем динамику системы
        nonEqSystemDyn = NonEqSystemDyn(nonEqSystem,  # Система

                                        fU3,  # Функция условий протекания процессов

                                        FunChar3,  # Функция внешних параметров

                                        integDynamic  # Метод интегрирования дифференциальных уравнений
                                        )

        # Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystem(0.0 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystem(0.3 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystem(0.9 * Tint,
                                          stateCoordinates,
                                          systemParameters)

        # Выполняем эталонные расчеты
        vxEt1 = countSystem3(0.0 * Tint,
                             stateCoordinates)
        vxEt2 = countSystem3(0.3 * Tint,
                             stateCoordinates)
        vxEt3 = countSystem3(0.9 * Tint,
                             stateCoordinates)

        # Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 9)  # Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 8)  # Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 7)  # Проверяем скорость

        # Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem3,
                                             (0.0, Tint),
                                             stateCoordinates,  # Начальное состояние
                                             t_eval=t_eval,  # Моменты времени
                                             method=method  # Метод численного интегрирования
                                             )  # Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1, 1)  # Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()  # Координаты состояния

        # Рассчитываем эталонную динамику параметров системы
        (t, Char1et, Char2et, Char3et) = FunChar3(t,  # Моменты времени
                                                  stateCoordinatesEt,  # Координаты состояния
                                                  systemParameters  # Параметры системы
                                                  )

        # Рассчитываем динамику параметров системы
        (t, Char1, Char2, Char3) = nonEqSystemDyn.NonEqSystemDynamic(Tint,
                                                                     stateCoordinates,
                                                                     systemParameters,
                                                                     t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)

        # Функтор модели параметров системы
        class functorSystemPar(object):
            # Функтор параметров
            def SystemParametersCount(self, measParameters):
                return (measParameters[0:6],
                        measParameters[6::])

        # Получаем класс модели
        funSystemPar = functorSystemPar()
        model = Model(funSystemPar,
                      nonEqSystemDyn)

        # Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              systemParameters))

        # Запускаем динамику
        (t, Char1, Char2, Char3) = model.CountDynamic(Tint, measChar,
                                                      t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)

    def testNonEqDyn8(self):
        # Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213])  # Координаты состояния
        systemParameters = np.array([10, 20, 50, 70, 100, 200, 30, 45, 10, 25, 35, 55, 75, 81, 210, 300, 150, 153, 123, 15.3, 45.9, 5.63, 10, 6.5, 45.4, 315, 161])  # Параметры системы

        # Стехиометрические коэффициенты химических реакций
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

        # Метод интегрирования дифференциальных уравнений
        method = "LSODA"
        atol = 6.3 * 1e-8  # Абсолютная погрешность
        rtol = 2.7 * 1e-2  # Относительная погрешность
        max_step = 3 * 1e-4  # Максимальный шаг интегрирования
        integDynamic = standartIntegrateDyn(method=method,
                                            atol=atol,
                                            rtol=rtol,
                                            max_step=max_step)

        # Время интегрирования
        Tint = 1e-6

        # Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)

        # Функция эталонного результата
        def countSystem3(t, stateCoordinates):
            # Вызов функции времени
            _systemParameters = fU3(t, systemParameters)

            # Вызов расчета
            (chemPot, Aff, kineticMatrix, vProcesses,
             balanceMatrix, vx) = CountSystem3(stateCoordinates, _systemParameters,
                                               nu1_1, nu1_2, nu1_3, nu1_4, nu1_5, nu1_6,
                                               nu1_7, nu1_8, nu1_9, nu1_10, nu1_11)

            # Выводим результат
            return vx

        # Задаем структуру системы
        stateCoordinatesNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]  # Имена координат состояния
        processCoordinatesNames = ["vChem1_1", "vChem1_2", "vChem1_3"]  # Имена координат процессов
        stateCoordinatesStreamsNames = []  # Имена координат состояния, изменяемых в результате внешних потоков

        # Задаем функцию состояния системы
        stateFunction = CountState3

        # Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]  # Имена переменных потенциалов взаимодействия по координатам состояния

        # Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1", "vChem1_1", "vChem1_2", "vChem1_2", "vChem1_3"]  # Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1", "vChem1_2", "vChem1_1", "vChem1_2", "vChem1_3"]  # Имена сопряженностей между собой термодинамических сил

        # Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []  # Имена переменных внешних потоков

        # Задаем систему
        nonEqSystem = NonEqSystem(stateCoordinatesNames,  # Имена координат состояния
                                  processCoordinatesNames,  # Имена координат процессов
                                  stateCoordinatesStreamsNames,  # Имена координат состояния, изменяемых в результате внешних потоков

                                  # Задаем функцию состояния системы
                                  stateFunction,

                                  # Задаем рассчитываемые параметры системы
                                  stateCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам состояния
                                  processCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам процессов
                                  stateCoordinatesVarPotentialsInterNames,  # Имена переменных потенциалов взаимодействия по координатам состояния

                                  # Задаем переменные параметры кинетической матрицы
                                  varKineticNames,  # Имена сопряженностей между собой координат процессов
                                  varKineticAffNames,  # Имена сопряженностей между собой термодинамических сил

                                  # Задаем внешние потоки
                                  stateCoordinatesVarStreamsNames  # Имена переменных внешних потоков
                                  )

        # Задаем постоянные параметры системы
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

        # Создаем динамику системы
        nonEqSystemDyn = NonEqSystemDyn(nonEqSystem,  # Система

                                        fU3,  # Функция условий протекания процессов

                                        FunChar3,  # Функция внешних параметров

                                        integDynamic  # Метод интегрирования дифференциальных уравнений
                                        )

        # Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystem(0.0 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystem(0.3 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystem(0.9 * Tint,
                                          stateCoordinates,
                                          systemParameters)

        # Выполняем эталонные расчеты
        vxEt1 = countSystem3(0.0 * Tint,
                             stateCoordinates)
        vxEt2 = countSystem3(0.3 * Tint,
                             stateCoordinates)
        vxEt3 = countSystem3(0.9 * Tint,
                             stateCoordinates)

        # Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 9)  # Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 8)  # Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 7)  # Проверяем скорость

        # Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem3,
                                             (0.0, Tint),
                                             stateCoordinates,  # Начальное состояние
                                             t_eval=t_eval,  # Моменты времени
                                             atol=atol,
                                             rtol=rtol,
                                             max_step=max_step,
                                             method=method  # Метод численного интегрирования
                                             )  # Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1, 1)  # Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()  # Координаты состояния

        # Рассчитываем эталонную динамику параметров системы
        (t, Char1et, Char2et, Char3et) = FunChar3(t,  # Моменты времени
                                                  stateCoordinatesEt,  # Координаты состояния
                                                  systemParameters  # Параметры системы
                                                  )

        # Рассчитываем динамику параметров системы
        (t, Char1, Char2, Char3) = nonEqSystemDyn.NonEqSystemDynamic(Tint,
                                                                     stateCoordinates,
                                                                     systemParameters,
                                                                     t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)

        # Функтор модели параметров системы
        class functorSystemPar(object):
            # Функтор параметров
            def SystemParametersCount(self, measParameters):
                return (measParameters[0:6],
                        measParameters[6::])

        # Получаем класс модели
        funSystemPar = functorSystemPar()
        model = Model(funSystemPar,
                      nonEqSystemDyn)

        # Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              systemParameters))

        # Запускаем динамику
        (t, Char1, Char2, Char3) = model.CountDynamic(Tint, measChar,
                                                      t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)

    def testNonEqDyn9(self):
        # Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213])  # Координаты состояния
        systemParameters = np.array([10, 20, 50, 70, 100, 200, 30, 45, 10, 25, 35, 55, 75, 81, 210, 300, 150, 153, 123, 15.3, 45.9, 5.63, 10, 6.5, 45.4, 315, 161])  # Параметры системы

        # Стехиометрические коэффициенты химических реакций
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

        # Метод интегрирования дифференциальных уравнений
        method = "LSODA"

        # Функция аттрибутов интегрирования
        def funIntegAttributes(systemParameters, method):
            # Определяем аттрибуты
            atol = np.abs(systemParameters[4] + 3.9 * systemParameters[1]) * 1e-9  # Абсолютная погрешность
            rtol = np.abs(systemParameters[14]) * 1e-5  # Относительная погрешность
            max_step = np.abs(8.1 * systemParameters[24] * systemParameters[12]) * 1e-4  # Максимальный шаг

            # Выводим аттрибуты интегрирования
            return (atol, rtol, max_step)
        integDynamic = stepIntegrateDyn(method=method,
                                        funIntegAttributes=funIntegAttributes)

        # Время интегрирования
        Tint = 1e-6

        # Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)

        # Функция эталонного результата
        def countSystem3(t, stateCoordinates):
            # Вызов функции времени
            _systemParameters = fU3(t, systemParameters)

            # Вызов расчета
            (chemPot, Aff, kineticMatrix, vProcesses,
             balanceMatrix, vx) = CountSystem3(stateCoordinates, _systemParameters,
                                               nu1_1, nu1_2, nu1_3, nu1_4, nu1_5, nu1_6,
                                               nu1_7, nu1_8, nu1_9, nu1_10, nu1_11)

            # Выводим результат
            return vx

        # Задаем структуру системы
        stateCoordinatesNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]  # Имена координат состояния
        processCoordinatesNames = ["vChem1_1", "vChem1_2", "vChem1_3"]  # Имена координат процессов
        stateCoordinatesStreamsNames = []  # Имена координат состояния, изменяемых в результате внешних потоков

        # Задаем функцию состояния системы
        stateFunction = CountState3

        # Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам процессов
        stateCoordinatesVarPotentialsInterNames = ["x1_1", "x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]  # Имена переменных потенциалов взаимодействия по координатам состояния

        # Задаем переменные параметры кинетической матрицы
        varKineticNames = ["vChem1_1", "vChem1_1", "vChem1_2", "vChem1_2", "vChem1_3"]  # Имена сопряженностей между собой координат процессов
        varKineticAffNames = ["vChem1_1", "vChem1_2", "vChem1_1", "vChem1_2", "vChem1_3"]  # Имена сопряженностей между собой термодинамических сил

        # Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []  # Имена переменных внешних потоков

        # Задаем систему
        nonEqSystem = NonEqSystem(stateCoordinatesNames,  # Имена координат состояния
                                  processCoordinatesNames,  # Имена координат процессов
                                  stateCoordinatesStreamsNames,  # Имена координат состояния, изменяемых в результате внешних потоков

                                  # Задаем функцию состояния системы
                                  stateFunction,

                                  # Задаем рассчитываемые параметры системы
                                  stateCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам состояния
                                  processCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам процессов
                                  stateCoordinatesVarPotentialsInterNames,  # Имена переменных потенциалов взаимодействия по координатам состояния

                                  # Задаем переменные параметры кинетической матрицы
                                  varKineticNames,  # Имена сопряженностей между собой координат процессов
                                  varKineticAffNames,  # Имена сопряженностей между собой термодинамических сил

                                  # Задаем внешние потоки
                                  stateCoordinatesVarStreamsNames  # Имена переменных внешних потоков
                                  )

        # Задаем постоянные параметры системы
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

        # Создаем динамику системы
        nonEqSystemDyn = NonEqSystemDyn(nonEqSystem,  # Система

                                        fU3,  # Функция условий протекания процессов

                                        FunChar3,  # Функция внешних параметров

                                        integDynamic  # Метод интегрирования дифференциальных уравнений
                                        )

        # Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystem(0.0 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystem(0.3 * Tint,
                                          stateCoordinates,
                                          systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystem(0.9 * Tint,
                                          stateCoordinates,
                                          systemParameters)

        # Выполняем эталонные расчеты
        vxEt1 = countSystem3(0.0 * Tint,
                             stateCoordinates)
        vxEt2 = countSystem3(0.3 * Tint,
                             stateCoordinates)
        vxEt3 = countSystem3(0.9 * Tint,
                             stateCoordinates)

        # Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 9)  # Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 8)  # Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 7)  # Проверяем скорость

        # Рассчитываем эталонную динамику системы
        (atol, rtol, max_step) = funIntegAttributes(systemParameters, method)
        stateCoordinatesEt = integ.solve_ivp(countSystem3,
                                             (0.0, Tint),
                                             stateCoordinates,  # Начальное состояние
                                             t_eval=t_eval,  # Моменты времени
                                             atol=atol,
                                             rtol=rtol,
                                             max_step=max_step,
                                             method=method  # Метод численного интегрирования
                                             )  # Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1, 1)  # Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()  # Координаты состояния

        # Рассчитываем эталонную динамику параметров системы
        (t, Char1et, Char2et, Char3et) = FunChar3(t,  # Моменты времени
                                                  stateCoordinatesEt,  # Координаты состояния
                                                  systemParameters  # Параметры системы
                                                  )

        # Рассчитываем динамику параметров системы
        (t, Char1, Char2, Char3) = nonEqSystemDyn.NonEqSystemDynamic(Tint,
                                                                     stateCoordinates,
                                                                     systemParameters,
                                                                     t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)

        # Функтор модели параметров системы
        class functorSystemPar(object):
            # Функтор параметров
            def SystemParametersCount(self, measParameters):
                return (measParameters[0:6],
                        measParameters[6::])

        # Получаем класс модели
        funSystemPar = functorSystemPar()
        model = Model(funSystemPar,
                      nonEqSystemDyn)

        # Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              systemParameters))

        # Запускаем динамику
        (t, Char1, Char2, Char3) = model.CountDynamic(Tint, measChar,
                                                      t_eval=t_eval)

        # Проверяем рассчитанные динамики выходных характеристик
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 9)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)


# Запустить тестирование
if __name__ == "__main__":
    unittest.main()
