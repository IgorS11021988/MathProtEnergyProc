import numpy as np
from scipy import integrate as integ
from MathProtEnergyProc.IntegrateDyn import standartIntegrateDyn, stepIntegrateDyn

from MathProtEnergyProc import NonEqSystemQ
from MathProtEnergyProc import NonEqSystemQDyn, ModelQ

from MathProtEnergyProc.tests.UnitTestExamples.TestNonEqQ1 import *
from MathProtEnergyProc.tests.UnitTestExamples.TestNonEqQ2 import *
from MathProtEnergyProc.tests.UnitTestExamples.TestNonEqQ3 import *
     
from MathProtEnergyProc.tests.UnitTestExamples.fU1 import *
from MathProtEnergyProc.tests.UnitTestExamples.fU2 import *
from MathProtEnergyProc.tests.UnitTestExamples.fU3 import *
     
from MathProtEnergyProc.tests.UnitTestExamples.FunCharQ1 import *
from MathProtEnergyProc.tests.UnitTestExamples.FunCharQ2 import *
from MathProtEnergyProc.tests.UnitTestExamples.FunCharQ3 import *

import unittest

#Модульные тесты
class TestNonEqSystemDynQ(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testNonEqDyn1(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81])#Координаты состояния
        reducedTemp = np.array([350, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters = np.array([ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                                       10,   25,   35,   55,   75,  81,   210,  300,  150,  153,
                                      123, 15.3, 45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                                      1.5, 33.9, 27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,
                                      159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                                     0.15, 74.1])#Параметры системы
        
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
        
        #Доли распределения некомпенсированных теплот по энергетических степеням свободы
        beta2_1 = 0.25
        beta2_2 = 0.75
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3
            
        #Базовые коэффициенты перекрестной кинетической матрицы процессов относительно переноса теплоты
        ADiffHeat1 = 0.15
        ADiffHeat2 = 0.09
        
        #Базовые коэффициенты перекрестной кинетической матрицы переноса теплоты относительно процессов
        AHeatDiff1 = 0.15
        AHeatDiff2 = 0.09
    
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 2.31
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 3.51
        H2_7 = 0.81
    
        #Метод интегрирования дифференциальных уравнений
        method = "RK45"
        integDynamic = standartIntegrateDyn(method=method)
        
        #Время интегрирования
        Tint = 1e-6
        
        #Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)
        
        #Функция эталонного результата
        def countSystem1(t,allStParameters):
            #Выделение координат состояния и приведенных температур
            _stateCoordinates = allStParameters[0:-2]
            _reducedTemp = allStParameters[-2::]
            
            #Вызов функции времени
            _systemParameters = fU1(t,systemParameters)
            
            #Вызов расчета
            (chemPot,potBet,Temp,Betas,Aff,AffHeat,
             kineticMatrixCPCP,kineticMatrixCPHeat,
             kineticMatrixHeatCP,kineticMatrixHeatHeat,
             vProcesses,vQTransf,vHeatProcesses,
             vQPows,heatTransferMatrix,balanceMatrix,
             vx,vUPows,invC,powH,vT) = CountSystemQ1(_stateCoordinates,_reducedTemp,_systemParameters,
                                                     nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                     nu1_8,nu1_9,nu1_10,nu1_11,nu2_1,nu2_2,nu2_3,
                                                     nu2_4,nu2_5,nu2_6,beta2_1,beta2_2,Adiff1_1,
                                                     Adiff1_2,Adiff2_1,ADiffHeat1,ADiffHeat2,
                                                     AHeatDiff1,AHeatDiff2,invC1,H1_1,H2_7)
            
            #Выводим результат
            return np.hstack((vx,vT))
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff1","vDiff2"]#Имена координат процессов
        energyPowersNames = ["EnPow1","EnPow2"]#Имена энергетических степеней свободы
        reducedTemperaturesEnergyPowersNames = ["T1", "T2"]#Имена приведенных температур энергетических степеней свободы
        energyPowersBetNames = []#Имена взаимодействий между энергетическими степенями свободы
        heatTransfersNames = ["vQTransf"]#Имена потоков переноса теплоты
        heatTransfersOutputEnergyPowersNames = ["EnPow1"]#Имена энергетических степеней свободы, с которых уходит теплота
        heatTransfersInputEnergyPowersNames = ["EnPow2"]#Имена энергетических степеней свободы, на которые приходит теплота
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        heatEnergyPowersStreamsNames = []#Имена потоков теплоты на энергетические степени свободы
        
        #Задаем функцию состояния системы
        stateFunction = CountStateQ1
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        energyPowersVarTemperatureNames = ["EnPow1","EnPow2"]#Имена переменных температур энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames = ["EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow2","EnPow2","EnPow2","EnPow2"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        energyPowersVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        energyPowersVarBetaNames = ["EnPow1","EnPow2"]#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
        processCoordinatesVarBetaNames = ["vDiff1","vDiff1"]#Имена переменных долей распределения некомпенсированной теплоты координат процессов
        reducedTemperaturesEnergyPowersVarInvHeatCapacityNames = ["T2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersVarInvHeatCapacityNames = ["EnPow2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        reducedTemperaturesEnergyPowersVarHeatEffectNames = ["T1", "T1", "T2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        stateCoordinatesVarHeatEffectNames = ["x1_2", "x1_5", "x2_1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticPCPCNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой координат процессов
        varKineticPCPCAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой термодинамических сил
        varKineticPCHeatNames = []#Имена сопряженностей координат процессов с теплопереносами
        varKineticPCHeatAffNames = []#Имена сопряженностей термодинамических сил с теплопереносами
        varKineticHeatPCNames = []#Имена сопряженностей теплопереносов с координатами процессов
        varKineticHeatPCAffNames = []#Имена сопряженностей теплопереносов с термодинамическими силами
        varKineticHeatHeatNames = ["vQTransf"]#Имена сопряженностей между собой перенесенных теплот
        varKineticHeatHeatAffNames = ["vQTransf"]#Имена сопряженностей между собой термодинамических сил по переносу теплот
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []#Имена переменных внешних потоков
        heatEnergyPowersVarStreamsNames = []#Имена переменных внешних потоков теплоты
        
        #Задаем систему
        nonEqSystem = NonEqSystemQ(stateCoordinatesNames,#Имена координат состояния
                                   processCoordinatesNames,#Имена координат процессов
                                   energyPowersNames,#Имена энергетических степеней свободы
                                   reducedTemperaturesEnergyPowersNames,#Имена приведенных температур энергетических степеней свободы
                                   energyPowersBetNames,#Имена взаимодействий между энергетическими степенями свободы
                                   heatTransfersNames,#Имена потоков переноса теплоты
                                   heatTransfersOutputEnergyPowersNames,#Имена энергетических степеней свободы, с которых уходит теплота
                                   heatTransfersInputEnergyPowersNames,#Имена энергетических степеней свободы, на которые приходит теплота
                                   stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                   heatEnergyPowersStreamsNames,#Имена потоков теплоты на энергетические степени свободы
                                     
                                   #Задаем функцию состояния системы
                                   stateFunction,
                                    
                                   #Задаем рассчитываемые параметры системы
                                   stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                   processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                   energyPowersVarTemperatureNames,#Имена переменных температур энергетических степеней свободы
                                   stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                   energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                   stateCoordinatesVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                                   energyPowersVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                                   energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                                   processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                                   reducedTemperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                   reducedTemperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   stateCoordinatesVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                    
                                   #Задаем переменные параметры кинетической матрицы
                                   varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                   varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                   varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                   varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                   varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                   varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                   varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                   varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                    
                                   #Задаем внешние потоки
                                   stateCoordinatesVarStreamsNames,#Имена переменных внешних потоков
                                   heatEnergyPowersVarStreamsNames#Имена переменных внешних потоков теплоты
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
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_3",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vDiff2",beta2_1)
        nonEqSystem.SetBetaConstElement("EnPow2","vDiff2",beta2_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff1", "vDiff1", Adiff1_1)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff1", "vDiff2", Adiff1_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff2", "vDiff1", Adiff2_1)
        nonEqSystem.SetKineticMatrixPCHeatConstElement("vDiff1", "vQTransf", ADiffHeat1)
        nonEqSystem.SetKineticMatrixPCHeatConstElement("vDiff2", "vQTransf", ADiffHeat2)
        nonEqSystem.SetKineticMatrixHeatPCConstElement("vQTransf", "vDiff1", AHeatDiff1)
        nonEqSystem.SetKineticMatrixHeatPCConstElement("vQTransf", "vDiff2", AHeatDiff2)
        nonEqSystem.SetInvHeatCapacityMatrixConstElement("T1", "EnPow1", invC1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T1", "x1_1", H1_1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T2", "x2_7", H2_7)
        
        #Создаем динамику системы
        nonEqSystemDyn = NonEqSystemQDyn(nonEqSystem,#Система
                 
                                         fU1,#Функция условий протекания процессов
                                         
                                         FunCharQ1,#Функция внешних параметров
                                         
                                         integDynamic#Метод интегрирования дифференциальных уравнений
                                         )
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystemQ(0.0*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystemQ(0.3*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystemQ(0.9*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        
        #Выполняем эталонные расчеты
        vxEt1 = countSystem1(0.0*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt2 = countSystem1(0.3*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt3 = countSystem1(0.9*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 6)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 6)#Проверяем скорость
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fSystemDyn(0.0*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx2 = nonEqSystemDyn.fSystemDyn(0.3*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx3 = nonEqSystemDyn.fSystemDyn(0.9*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 6)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 6)#Проверяем скорость
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem1,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates,reducedTemp)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t,Char1et,Char2et) = FunCharQ1(t,#Моменты времени
                                        stateCoordinatesEt,#Координаты состояния
                                        reducedTempEt,#Приведенные температуры
                                        systemParameters#Параметры системы  
                                        )
        
        #Рассчитываем динамику параметров системы
        (t,Char1,Char2) = nonEqSystemDyn.NonEqSystemDynamicQ(Tint,
                                                             stateCoordinates,
                                                             reducedTemp,
                                                             systemParameters,
                                                             t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 6)
        
        #Функтор модели параметров системы
        class functorSystemPar(object):
            #Функтор параметров
            def SystemParametersCount(self,measParameters):
                return (measParameters[0:10],
                        measParameters[10:12],
                        measParameters[12::])
        
        #Получаем класс модели
        funSystemPar = functorSystemPar()
        model = ModelQ(funSystemPar,
                       nonEqSystemDyn)
        
        #Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              reducedTemp,
                              systemParameters))
        
        #Запускаем динамику
        (t,Char1,Char2) = model.CountDynamic(Tint,measChar,
                                             t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 6)
    def testNonEqDyn2(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81])#Координаты состояния
        reducedTemp = np.array([350, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters = np.array([ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                                       10,   25,   35,   55,   75,  81,   210,  300,  150,  153,
                                      123, 15.3, 45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                                      1.5, 33.9, 27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,
                                      159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                                     0.15, 74.1])#Параметры системы
        
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
        
        #Доли распределения некомпенсированных теплот по энергетических степеням свободы
        beta2_1 = 0.25
        beta2_2 = 0.75
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3
            
        #Базовые коэффициенты перекрестной кинетической матрицы процессов относительно переноса теплоты
        ADiffHeat1 = 0.15
        ADiffHeat2 = 0.09
        
        #Базовые коэффициенты перекрестной кинетической матрицы переноса теплоты относительно процессов
        AHeatDiff1 = 0.15
        AHeatDiff2 = 0.09
    
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 2.31
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 3.51
        H2_7 = 0.81
    
        #Метод интегрирования дифференциальных уравнений
        method = "RK45"
        atol = 2*1e-7#Абсолютная погрешность
        rtol = 3*1e-6#Относительная погрешность
        integDynamic = standartIntegrateDyn(method=method,
                                            atol=atol,
                                            rtol=rtol)
        
        #Время интегрирования
        Tint = 1e-6
        
        #Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)
        
        #Функция эталонного результата
        def countSystem1(t,allStParameters):
            #Выделение координат состояния и приведенных температур
            _stateCoordinates = allStParameters[0:-2]
            _reducedTemp = allStParameters[-2::]
            
            #Вызов функции времени
            _systemParameters = fU1(t,systemParameters)
            
            #Вызов расчета
            (chemPot,potBet,Temp,Betas,Aff,AffHeat,
             kineticMatrixCPCP,kineticMatrixCPHeat,
             kineticMatrixHeatCP,kineticMatrixHeatHeat,
             vProcesses,vQTransf,vHeatProcesses,
             vQPows,heatTransferMatrix,balanceMatrix,
             vx,vUPows,invC,powH,vT) = CountSystemQ1(_stateCoordinates,_reducedTemp,_systemParameters,
                                                     nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                     nu1_8,nu1_9,nu1_10,nu1_11,nu2_1,nu2_2,nu2_3,
                                                     nu2_4,nu2_5,nu2_6,beta2_1,beta2_2,Adiff1_1,
                                                     Adiff1_2,Adiff2_1,ADiffHeat1,ADiffHeat2,
                                                     AHeatDiff1,AHeatDiff2,invC1,H1_1,H2_7)
            
            #Выводим результат
            return np.hstack((vx,vT))
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff1","vDiff2"]#Имена координат процессов
        energyPowersNames = ["EnPow1","EnPow2"]#Имена энергетических степеней свободы
        reducedTemperaturesEnergyPowersNames = ["T1", "T2"]#Имена приведенных температур энергетических степеней свободы
        energyPowersBetNames = []#Имена взаимодействий между энергетическими степенями свободы
        heatTransfersNames = ["vQTransf"]#Имена потоков переноса теплоты
        heatTransfersOutputEnergyPowersNames = ["EnPow1"]#Имена энергетических степеней свободы, с которых уходит теплота
        heatTransfersInputEnergyPowersNames = ["EnPow2"]#Имена энергетических степеней свободы, на которые приходит теплота
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        heatEnergyPowersStreamsNames = []#Имена потоков теплоты на энергетические степени свободы
        
        #Задаем функцию состояния системы
        stateFunction = CountStateQ1
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        energyPowersVarTemperatureNames = ["EnPow1","EnPow2"]#Имена переменных температур энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames = ["EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow2","EnPow2","EnPow2","EnPow2"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        energyPowersVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        energyPowersVarBetaNames = ["EnPow1","EnPow2"]#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
        processCoordinatesVarBetaNames = ["vDiff1","vDiff1"]#Имена переменных долей распределения некомпенсированной теплоты координат процессов
        reducedTemperaturesEnergyPowersVarInvHeatCapacityNames = ["T2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersVarInvHeatCapacityNames = ["EnPow2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        reducedTemperaturesEnergyPowersVarHeatEffectNames = ["T1", "T1", "T2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        stateCoordinatesVarHeatEffectNames = ["x1_2", "x1_5", "x2_1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticPCPCNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой координат процессов
        varKineticPCPCAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой термодинамических сил
        varKineticPCHeatNames = []#Имена сопряженностей координат процессов с теплопереносами
        varKineticPCHeatAffNames = []#Имена сопряженностей термодинамических сил с теплопереносами
        varKineticHeatPCNames = []#Имена сопряженностей теплопереносов с координатами процессов
        varKineticHeatPCAffNames = []#Имена сопряженностей теплопереносов с термодинамическими силами
        varKineticHeatHeatNames = ["vQTransf"]#Имена сопряженностей между собой перенесенных теплот
        varKineticHeatHeatAffNames = ["vQTransf"]#Имена сопряженностей между собой термодинамических сил по переносу теплот
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []#Имена переменных внешних потоков
        heatEnergyPowersVarStreamsNames = []#Имена переменных внешних потоков теплоты
        
        #Задаем систему
        nonEqSystem = NonEqSystemQ(stateCoordinatesNames,#Имена координат состояния
                                   processCoordinatesNames,#Имена координат процессов
                                   energyPowersNames,#Имена энергетических степеней свободы
                                   reducedTemperaturesEnergyPowersNames,#Имена приведенных температур энергетических степеней свободы
                                   energyPowersBetNames,#Имена взаимодействий между энергетическими степенями свободы
                                   heatTransfersNames,#Имена потоков переноса теплоты
                                   heatTransfersOutputEnergyPowersNames,#Имена энергетических степеней свободы, с которых уходит теплота
                                   heatTransfersInputEnergyPowersNames,#Имена энергетических степеней свободы, на которые приходит теплота
                                   stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                   heatEnergyPowersStreamsNames,#Имена потоков теплоты на энергетические степени свободы
                                     
                                   #Задаем функцию состояния системы
                                   stateFunction,
                                    
                                   #Задаем рассчитываемые параметры системы
                                   stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                   processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                   energyPowersVarTemperatureNames,#Имена переменных температур энергетических степеней свободы
                                   stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                   energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                   stateCoordinatesVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                                   energyPowersVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                                   energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                                   processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                                   reducedTemperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                   reducedTemperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   stateCoordinatesVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                    
                                   #Задаем переменные параметры кинетической матрицы
                                   varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                   varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                   varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                   varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                   varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                   varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                   varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                   varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                    
                                   #Задаем внешние потоки
                                   stateCoordinatesVarStreamsNames,#Имена переменных внешних потоков
                                   heatEnergyPowersVarStreamsNames#Имена переменных внешних потоков теплоты
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
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_3",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vDiff2",beta2_1)
        nonEqSystem.SetBetaConstElement("EnPow2","vDiff2",beta2_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff1", "vDiff1", Adiff1_1)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff1", "vDiff2", Adiff1_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff2", "vDiff1", Adiff2_1)
        nonEqSystem.SetKineticMatrixPCHeatConstElement("vDiff1", "vQTransf", ADiffHeat1)
        nonEqSystem.SetKineticMatrixPCHeatConstElement("vDiff2", "vQTransf", ADiffHeat2)
        nonEqSystem.SetKineticMatrixHeatPCConstElement("vQTransf", "vDiff1", AHeatDiff1)
        nonEqSystem.SetKineticMatrixHeatPCConstElement("vQTransf", "vDiff2", AHeatDiff2)
        nonEqSystem.SetInvHeatCapacityMatrixConstElement("T1", "EnPow1", invC1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T1", "x1_1", H1_1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T2", "x2_7", H2_7)
        
        #Создаем динамику системы
        nonEqSystemDyn = NonEqSystemQDyn(nonEqSystem,#Система
                 
                                         fU1,#Функция условий протекания процессов
                                         
                                         FunCharQ1,#Функция внешних параметров
                                         
                                         integDynamic#Метод интегрирования дифференциальных уравнений
                                         )
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystemQ(0.0*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystemQ(0.3*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystemQ(0.9*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        
        #Выполняем эталонные расчеты
        vxEt1 = countSystem1(0.0*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt2 = countSystem1(0.3*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt3 = countSystem1(0.9*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 6)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 6)#Проверяем скорость
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fSystemDyn(0.0*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx2 = nonEqSystemDyn.fSystemDyn(0.3*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx3 = nonEqSystemDyn.fSystemDyn(0.9*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 6)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 6)#Проверяем скорость
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem1,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates,reducedTemp)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             atol=atol,
                                             rtol=rtol,
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t,Char1et,Char2et) = FunCharQ1(t,#Моменты времени
                                        stateCoordinatesEt,#Координаты состояния
                                        reducedTempEt,#Приведенные температуры
                                        systemParameters#Параметры системы  
                                        )
        
        #Рассчитываем динамику параметров системы
        (t,Char1,Char2) = nonEqSystemDyn.NonEqSystemDynamicQ(Tint,
                                                             stateCoordinates,
                                                             reducedTemp,
                                                             systemParameters,
                                                             t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 6)
        
        #Функтор модели параметров системы
        class functorSystemPar(object):
            #Функтор параметров
            def SystemParametersCount(self,measParameters):
                return (measParameters[0:10],
                        measParameters[10:12],
                        measParameters[12::])
        
        #Получаем класс модели
        funSystemPar = functorSystemPar()
        model = ModelQ(funSystemPar,
                       nonEqSystemDyn)
        
        #Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              reducedTemp,
                              systemParameters))
        
        #Запускаем динамику
        (t,Char1,Char2) = model.CountDynamic(Tint,measChar,
                                             t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 6)
    def testNonEqDyn3(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81])#Координаты состояния
        reducedTemp = np.array([350, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters = np.array([ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                                       10,   25,   35,   55,   75,  81,   210,  300,  150,  153,
                                      123, 15.3, 45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                                      1.5, 33.9, 27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,
                                      159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                                     0.15, 74.1])#Параметры системы
        
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
        
        #Доли распределения некомпенсированных теплот по энергетических степеням свободы
        beta2_1 = 0.25
        beta2_2 = 0.75
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3
            
        #Базовые коэффициенты перекрестной кинетической матрицы процессов относительно переноса теплоты
        ADiffHeat1 = 0.15
        ADiffHeat2 = 0.09
        
        #Базовые коэффициенты перекрестной кинетической матрицы переноса теплоты относительно процессов
        AHeatDiff1 = 0.15
        AHeatDiff2 = 0.09
    
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 2.31
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 3.51
        H2_7 = 0.81
    
        #Метод интегрирования дифференциальных уравнений
        method = "RK45"
        def funIntegAttributes(systemParameters,method):
            #Определяем аттрибуты
            atol = np.abs(systemParameters[3])*1e-7#Абсолютная погрешность
            rtol = np.abs(systemParameters[9])*1e-4#Относительная погрешность
            max_step = np.abs(2*systemParameters[21] + 3*systemParameters[15])*1e-1#Максимальный шаг
            
            #Выводим аттрибуты интегрирования
            return (atol,rtol,max_step)
        integDynamic = stepIntegrateDyn(method=method,
                                        funIntegAttributes=funIntegAttributes)
        
        #Время интегрирования
        Tint = 1e-6
        
        #Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)
        
        #Функция эталонного результата
        def countSystem1(t,allStParameters):
            #Выделение координат состояния и приведенных температур
            _stateCoordinates = allStParameters[0:-2]
            _reducedTemp = allStParameters[-2::]
            
            #Вызов функции времени
            _systemParameters = fU1(t,systemParameters)
            
            #Вызов расчета
            (chemPot,potBet,Temp,Betas,Aff,AffHeat,
             kineticMatrixCPCP,kineticMatrixCPHeat,
             kineticMatrixHeatCP,kineticMatrixHeatHeat,
             vProcesses,vQTransf,vHeatProcesses,
             vQPows,heatTransferMatrix,balanceMatrix,
             vx,vUPows,invC,powH,vT) = CountSystemQ1(_stateCoordinates,_reducedTemp,_systemParameters,
                                                     nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                     nu1_8,nu1_9,nu1_10,nu1_11,nu2_1,nu2_2,nu2_3,
                                                     nu2_4,nu2_5,nu2_6,beta2_1,beta2_2,Adiff1_1,
                                                     Adiff1_2,Adiff2_1,ADiffHeat1,ADiffHeat2,
                                                     AHeatDiff1,AHeatDiff2,invC1,H1_1,H2_7)
            
            #Выводим результат
            return np.hstack((vx,vT))
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff1","vDiff2"]#Имена координат процессов
        energyPowersNames = ["EnPow1","EnPow2"]#Имена энергетических степеней свободы
        reducedTemperaturesEnergyPowersNames = ["T1", "T2"]#Имена приведенных температур энергетических степеней свободы
        energyPowersBetNames = []#Имена взаимодействий между энергетическими степенями свободы
        heatTransfersNames = ["vQTransf"]#Имена потоков переноса теплоты
        heatTransfersOutputEnergyPowersNames = ["EnPow1"]#Имена энергетических степеней свободы, с которых уходит теплота
        heatTransfersInputEnergyPowersNames = ["EnPow2"]#Имена энергетических степеней свободы, на которые приходит теплота
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        heatEnergyPowersStreamsNames = []#Имена потоков теплоты на энергетические степени свободы
        
        #Задаем функцию состояния системы
        stateFunction = CountStateQ1
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        energyPowersVarTemperatureNames = ["EnPow1","EnPow2"]#Имена переменных температур энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames = ["EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow2","EnPow2","EnPow2","EnPow2"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        energyPowersVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        energyPowersVarBetaNames = ["EnPow1","EnPow2"]#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
        processCoordinatesVarBetaNames = ["vDiff1","vDiff1"]#Имена переменных долей распределения некомпенсированной теплоты координат процессов
        reducedTemperaturesEnergyPowersVarInvHeatCapacityNames = ["T2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersVarInvHeatCapacityNames = ["EnPow2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        reducedTemperaturesEnergyPowersVarHeatEffectNames = ["T1", "T1", "T2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        stateCoordinatesVarHeatEffectNames = ["x1_2", "x1_5", "x2_1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticPCPCNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой координат процессов
        varKineticPCPCAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой термодинамических сил
        varKineticPCHeatNames = []#Имена сопряженностей координат процессов с теплопереносами
        varKineticPCHeatAffNames = []#Имена сопряженностей термодинамических сил с теплопереносами
        varKineticHeatPCNames = []#Имена сопряженностей теплопереносов с координатами процессов
        varKineticHeatPCAffNames = []#Имена сопряженностей теплопереносов с термодинамическими силами
        varKineticHeatHeatNames = ["vQTransf"]#Имена сопряженностей между собой перенесенных теплот
        varKineticHeatHeatAffNames = ["vQTransf"]#Имена сопряженностей между собой термодинамических сил по переносу теплот
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []#Имена переменных внешних потоков
        heatEnergyPowersVarStreamsNames = []#Имена переменных внешних потоков теплоты
        
        #Задаем систему
        nonEqSystem = NonEqSystemQ(stateCoordinatesNames,#Имена координат состояния
                                   processCoordinatesNames,#Имена координат процессов
                                   energyPowersNames,#Имена энергетических степеней свободы
                                   reducedTemperaturesEnergyPowersNames,#Имена приведенных температур энергетических степеней свободы
                                   energyPowersBetNames,#Имена взаимодействий между энергетическими степенями свободы
                                   heatTransfersNames,#Имена потоков переноса теплоты
                                   heatTransfersOutputEnergyPowersNames,#Имена энергетических степеней свободы, с которых уходит теплота
                                   heatTransfersInputEnergyPowersNames,#Имена энергетических степеней свободы, на которые приходит теплота
                                   stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                   heatEnergyPowersStreamsNames,#Имена потоков теплоты на энергетические степени свободы
                                     
                                   #Задаем функцию состояния системы
                                   stateFunction,
                                    
                                   #Задаем рассчитываемые параметры системы
                                   stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                   processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                   energyPowersVarTemperatureNames,#Имена переменных температур энергетических степеней свободы
                                   stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                   energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                   stateCoordinatesVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                                   energyPowersVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                                   energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                                   processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                                   reducedTemperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                   reducedTemperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   stateCoordinatesVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                    
                                   #Задаем переменные параметры кинетической матрицы
                                   varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                   varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                   varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                   varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                   varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                   varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                   varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                   varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                    
                                   #Задаем внешние потоки
                                   stateCoordinatesVarStreamsNames,#Имена переменных внешних потоков
                                   heatEnergyPowersVarStreamsNames#Имена переменных внешних потоков теплоты
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
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_3",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vDiff2",beta2_1)
        nonEqSystem.SetBetaConstElement("EnPow2","vDiff2",beta2_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff1", "vDiff1", Adiff1_1)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff1", "vDiff2", Adiff1_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff2", "vDiff1", Adiff2_1)
        nonEqSystem.SetKineticMatrixPCHeatConstElement("vDiff1", "vQTransf", ADiffHeat1)
        nonEqSystem.SetKineticMatrixPCHeatConstElement("vDiff2", "vQTransf", ADiffHeat2)
        nonEqSystem.SetKineticMatrixHeatPCConstElement("vQTransf", "vDiff1", AHeatDiff1)
        nonEqSystem.SetKineticMatrixHeatPCConstElement("vQTransf", "vDiff2", AHeatDiff2)
        nonEqSystem.SetInvHeatCapacityMatrixConstElement("T1", "EnPow1", invC1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T1", "x1_1", H1_1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T2", "x2_7", H2_7)
        
        #Создаем динамику системы
        nonEqSystemDyn = NonEqSystemQDyn(nonEqSystem,#Система
                 
                                         fU1,#Функция условий протекания процессов
                                         
                                         FunCharQ1,#Функция внешних параметров
                                         
                                         integDynamic#Метод интегрирования дифференциальных уравнений
                                         )
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystemQ(0.0*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystemQ(0.3*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystemQ(0.9*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        
        #Выполняем эталонные расчеты
        vxEt1 = countSystem1(0.0*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt2 = countSystem1(0.3*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt3 = countSystem1(0.9*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 6)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 6)#Проверяем скорость
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fSystemDyn(0.0*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx2 = nonEqSystemDyn.fSystemDyn(0.3*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx3 = nonEqSystemDyn.fSystemDyn(0.9*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 6)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 6)#Проверяем скорость
        
        #Рассчитываем эталонную динамику системы
        (atol,rtol,max_step) = funIntegAttributes(systemParameters,method)
        stateCoordinatesEt = integ.solve_ivp(countSystem1,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates,reducedTemp)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             atol=atol,
                                             rtol=rtol,
                                             max_step=max_step,
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t,Char1et,Char2et) = FunCharQ1(t,#Моменты времени
                                        stateCoordinatesEt,#Координаты состояния
                                        reducedTempEt,#Приведенные температуры
                                        systemParameters#Параметры системы  
                                        )
        
        #Рассчитываем динамику параметров системы
        (t,Char1,Char2) = nonEqSystemDyn.NonEqSystemDynamicQ(Tint,
                                                             stateCoordinates,
                                                             reducedTemp,
                                                             systemParameters,
                                                             t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 6)
        
        #Функтор модели параметров системы
        class functorSystemPar(object):
            #Функтор параметров
            def SystemParametersCount(self,measParameters):
                return (measParameters[0:10],
                        measParameters[10:12],
                        measParameters[12::])
        
        #Получаем класс модели
        funSystemPar = functorSystemPar()
        model = ModelQ(funSystemPar,
                       nonEqSystemDyn)
        
        #Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              reducedTemp,
                              systemParameters))
        
        #Запускаем динамику
        (t,Char1,Char2) = model.CountDynamic(Tint,measChar,
                                             t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 6)
    def testNonEqDyn4(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81, 80, 64])#Координаты состояния
        reducedTemp = np.array([350, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters = systemParameters = np.array([ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                                                          10,   25,   35,   55,   75,  81,   210,  300,  150,  153,
                                                         123, 15.3, 45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                                                         1.5, 33.9, 27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,
                                                         159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                                                        0.15, 74.1, 81.0, 5.13, 31.8, 3.3, 39.63, 54.3, 27.3, 0.15,
                                                        0.15])#Параметры системы
        
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
        
        #Доли распределения некомпенсированных теплот по энергетических степеням свободы
        beta2_1 = 0.25
        beta2_2 = 0.75
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        AChem1_4_4 = 3.3
        AChem2_3_3 = 6.3
        AChem2_1_2 = 0.3
        AChem2_2_1 = 0.15
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3
        Adiff3_3 = 4.83
            
        #Базовые коэффициенты перекрестной кинетической матрицы процессов относительно переноса теплоты
        ADiffHeat2 = 0.09
        ADiffHeat3 = 0.063
        
        #Базовые коэффициенты перекрестной кинетической матрицы переноса теплоты относительно процессов
        AHeatDiff2 = 0.09
        AHeatDiff3 = 0.03
    
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 2.31
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 3.51
        H2_7 = 0.81
    
        #Температура окружающей среды
        Tokr = 200
        
        #Коэффициент теплотдачи в окружающую среду
        AHeatOkr = 105
        
        #Внешние потоки теплоты
        ExtQEnPow1 = 90.0
        
        #Внешние потоки вещества
        xExt1_6 = 30.9
        
        #Потенциал взаимодействия
        mu2_9 = 13.5
        
        #Потенциал взаимодействия между энергетическими степенями свободы
        muBet1_3 = 18
        muBet1_5 = 187
    
        #Метод интегрирования дифференциальных уравнений
        method = "RK23"
        integDynamic = standartIntegrateDyn(method=method)
        
        #Время интегрирования
        Tint = 1e-6
        
        #Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)
        
        #Функция эталонного результата
        def countSystem2(t,allStParameters):
            #Выделение координат состояния и приведенных температур
            _stateCoordinates = allStParameters[0:-2]
            _reducedTemp = allStParameters[-2::]
            
            #Вызов функции времени
            _systemParameters = fU2(t,systemParameters)
            
            #Вызов расчета
            (chemPot,potBet,Temp,Betas,Aff,AffHeat,
             kineticMatrixCPCP,kineticMatrixCPHeat,
             kineticMatrixHeatCP,kineticMatrixHeatHeat,
             vProcesses,vQTransf,vHeatProcesses,
             vQPows,heatTransferMatrix,balanceMatrix,
             vx,vUPows,invC,powH,vT,heatStreams,Streams) = CountSystemQ2(_stateCoordinates,_reducedTemp,_systemParameters,
                                                                         nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                                         nu1_8,nu1_9,nu1_10,nu1_11,nu1_12,nu1_13,nu1_14,
                                                                         nu1_15,nu2_1,nu2_2,nu2_3,nu2_4,nu2_5,nu2_6,
                                                                         nu2_7,nu2_8,nu2_9,beta2_1,beta2_2,Adiff1_1,
                                                                         Adiff1_2,Adiff2_1,ADiffHeat2,AHeatDiff2,invC1,
                                                                         H1_1,H2_7,Tokr,AHeatOkr,ExtQEnPow1,AChem1_4_4,
                                                                         mu2_9,AChem2_3_3,AChem2_1_2,AChem2_2_1,Adiff3_3,
                                                                         ADiffHeat3,AHeatDiff3,xExt1_6,muBet1_3,muBet1_5)
            
            #Выводим результат
            return np.hstack((vx,vT))
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8", "x2_9"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem1_4","vChem2_1","vChem2_2","vChem2_3","vDiff1","vDiff2","vDiff3"]#Имена координат процессов
        energyPowersNames = ["EnPow1","EnPow2","EnOkr"]#Имена энергетических степеней свободы
        reducedTemperaturesEnergyPowersNames = ["T1", "T2"]#Имена приведенных температур энергетических степеней свободы
        energyPowersBetNames = ["EnPowBet1","EnPowBet2","EnPowBet3","EnPowBet4"]#Имена взаимодействий между энергетическими степенями свободы
        heatTransfersNames = ["vQTransf","vQTrOkr"]#Имена потоков переноса теплоты
        heatTransfersOutputEnergyPowersNames = ["EnPow1","EnPow1"]#Имена энергетических степеней свободы, с которых уходит теплота
        heatTransfersInputEnergyPowersNames = ["EnPow2","EnOkr"]#Имена энергетических степеней свободы, на которые приходит теплота
        stateCoordinatesStreamsNames = ["x1_2","x1_6","x2_8"]#Имена координат состояния, изменяемых в результате внешних потоков
        heatEnergyPowersStreamsNames = ["EnPow1","EnOkr"]#Имена потоков теплоты на энергетические степени свободы
        
        #Задаем функцию состояния системы
        stateFunction = CountStateQ2
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        energyPowersVarTemperatureNames = ["EnPow1","EnPow2"]#Имена переменных температур энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames = ["EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow2","EnPow2","EnPow2","EnPow2"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesVarPotentialsInterBetNames = ["x1_2"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        energyPowersVarPotentialsInterBetNames = ["EnPowBet2"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        energyPowersVarBetaNames = ["EnPow1","EnPow2"]#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
        processCoordinatesVarBetaNames = ["vDiff1","vDiff1"]#Имена переменных долей распределения некомпенсированной теплоты координат процессов
        reducedTemperaturesEnergyPowersVarInvHeatCapacityNames = ["T2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersVarInvHeatCapacityNames = ["EnPow2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        reducedTemperaturesEnergyPowersVarHeatEffectNames = ["T1", "T1", "T2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        stateCoordinatesVarHeatEffectNames = ["x1_2", "x1_5", "x2_1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticPCPCNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой координат процессов
        varKineticPCPCAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой термодинамических сил
        varKineticPCHeatNames = ["vDiff1"]#Имена сопряженностей координат процессов с теплопереносами
        varKineticPCHeatAffNames = ["vQTransf"]#Имена сопряженностей термодинамических сил с теплопереносами
        varKineticHeatPCNames = ["vQTransf"]#Имена сопряженностей теплопереносов с координатами процессов
        varKineticHeatPCAffNames = ["vDiff1"]#Имена сопряженностей теплопереносов с термодинамическими силами
        varKineticHeatHeatNames = ["vQTransf"]#Имена сопряженностей между собой перенесенных теплот
        varKineticHeatHeatAffNames = ["vQTransf"]#Имена сопряженностей между собой термодинамических сил по переносу теплот
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = ["x1_2","x2_8"]#Имена переменных внешних потоков
        heatEnergyPowersVarStreamsNames = ["EnOkr"]#Имена переменных внешних потоков теплоты
        
        #Задаем систему
        nonEqSystem = NonEqSystemQ(stateCoordinatesNames,#Имена координат состояния
                                   processCoordinatesNames,#Имена координат процессов
                                   energyPowersNames,#Имена энергетических степеней свободы
                                   reducedTemperaturesEnergyPowersNames,#Имена приведенных температур энергетических степеней свободы
                                   energyPowersBetNames,#Имена взаимодействий между энергетическими степенями свободы
                                   heatTransfersNames,#Имена потоков переноса теплоты
                                   heatTransfersOutputEnergyPowersNames,#Имена энергетических степеней свободы, с которых уходит теплота
                                   heatTransfersInputEnergyPowersNames,#Имена энергетических степеней свободы, на которые приходит теплота
                                   stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                   heatEnergyPowersStreamsNames,#Имена потоков теплоты на энергетические степени свободы
                                    
                                   #Задаем функцию состояния системы
                                   stateFunction,
                                    
                                   #Задаем рассчитываемые параметры системы
                                   stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                   processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                   energyPowersVarTemperatureNames,#Имена переменных температур энергетических степеней свободы
                                   stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                   energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                   stateCoordinatesVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                                   energyPowersVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                                   energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                                   processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                                   reducedTemperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                   reducedTemperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   stateCoordinatesVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                    
                                   #Задаем переменные параметры кинетической матрицы
                                   varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                   varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                   varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                   varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                   varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                   varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                   varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                   varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                    
                                   #Задаем внешние потоки
                                   stateCoordinatesVarStreamsNames,#Имена переменных внешних потоков
                                   heatEnergyPowersVarStreamsNames#Имена переменных внешних потоков теплоты
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
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_3",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_4",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_3",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vDiff2",beta2_1)
        nonEqSystem.SetBetaConstElement("EnPow2","vDiff2",beta2_2)
        nonEqSystem.SetBetaConstElement("EnPow1","vDiff3",beta2_1)
        nonEqSystem.SetBetaConstElement("EnPow2","vDiff3",beta2_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vChem1_4", "vChem1_4", AChem1_4_4)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vChem2_3", "vChem2_3", AChem2_3_3)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vChem2_1", "vChem2_2", AChem2_1_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vChem2_2", "vChem2_1", AChem2_2_1)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff1", "vDiff1", Adiff1_1)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff1", "vDiff2", Adiff1_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff2", "vDiff1", Adiff2_1)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff3", "vDiff3", Adiff3_3)
        nonEqSystem.SetKineticMatrixPCHeatConstElement("vDiff2", "vQTransf", ADiffHeat2)
        nonEqSystem.SetKineticMatrixHeatPCConstElement("vQTransf", "vDiff2", AHeatDiff2)
        nonEqSystem.SetKineticMatrixPCHeatConstElement("vDiff3", "vQTransf", ADiffHeat3)
        nonEqSystem.SetKineticMatrixHeatPCConstElement("vQTransf", "vDiff3", AHeatDiff3)
        nonEqSystem.SetKineticMatrixHeatHeatConstElement("vQTrOkr", "vQTrOkr", AHeatOkr)
        nonEqSystem.SetInvHeatCapacityMatrixConstElement("T1", "EnPow1", invC1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T1", "x1_1", H1_1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T2", "x2_7", H2_7)
        nonEqSystem.SetTEnergyPowersConstElement("EnOkr", Tokr)
        nonEqSystem.SetHeatEnergyPowersStreamsConstElement("EnPow1",ExtQEnPow1)
        nonEqSystem.SetStateCoordinatesStreamsConstElement("x1_6",xExt1_6)
        nonEqSystem.SetPotentialsInterConstElement("EnPow2", "x2_9", -mu2_9)
        nonEqSystem.SetPotentialsInterBetConstElement("EnPowBet1", "x1_5", -muBet1_5)
        nonEqSystem.SetPotentialsInterBetConstElement("EnPowBet3", "x1_3", -muBet1_3)
        
        #Создаем динамику системы
        nonEqSystemDyn = NonEqSystemQDyn(nonEqSystem,#Система
                 
                                         fU2,#Функция условий протекания процессов
                                         
                                         FunCharQ2,#Функция внешних параметров
                                         
                                         integDynamic#Метод интегрирования дифференциальных уравнений
                                         )
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystemQ(0.0*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystemQ(0.3*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystemQ(0.9*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        
        #Выполняем эталонные расчеты
        vxEt1 = countSystem2(0.0*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt2 = countSystem2(0.3*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt3 = countSystem2(0.9*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 5)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 5)#Проверяем скорость
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fSystemDyn(0.0*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx2 = nonEqSystemDyn.fSystemDyn(0.3*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx3 = nonEqSystemDyn.fSystemDyn(0.9*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 5)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 5)#Проверяем скорость
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem2,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates,reducedTemp)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t,Char1et,Char2et,Char3et) = FunCharQ2(t,#Моменты времени
                                                stateCoordinatesEt,#Координаты состояния
                                                reducedTempEt,#Приведенные температуры
                                                systemParameters#Параметры системы  
                                                )
        
        #Рассчитываем динамику параметров системы
        (t,Char1,Char2,Char3) = nonEqSystemDyn.NonEqSystemDynamicQ(Tint,
                                                                   stateCoordinates,
                                                                   reducedTemp,
                                                                   systemParameters,
                                                                   t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 7)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 7)
        
        #Функтор модели параметров системы
        class functorSystemPar(object):
            #Функтор параметров
            def SystemParametersCount(self,measParameters):
                return (measParameters[0:12],
                        measParameters[12:14],
                        measParameters[14::])
        
        #Получаем класс модели
        funSystemPar = functorSystemPar()
        model = ModelQ(funSystemPar,
                       nonEqSystemDyn)
        
        #Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              reducedTemp,
                              systemParameters))
        
        #Запускаем динамику
        (t,Char1,Char2,Char3) = model.CountDynamic(Tint,measChar,
                                                   t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 7)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 7)
    def testNonEqDyn5(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81, 80, 64])#Координаты состояния
        reducedTemp = np.array([350, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters = systemParameters = np.array([ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                                                          10,   25,   35,   55,   75,  81,   210,  300,  150,  153,
                                                         123, 15.3, 45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                                                         1.5, 33.9, 27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,
                                                         159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                                                        0.15, 74.1, 81.0, 5.13, 31.8, 3.3, 39.63, 54.3, 27.3, 0.15,
                                                        0.15])#Параметры системы
        
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
        
        #Доли распределения некомпенсированных теплот по энергетических степеням свободы
        beta2_1 = 0.25
        beta2_2 = 0.75
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        AChem1_4_4 = 3.3
        AChem2_3_3 = 6.3
        AChem2_1_2 = 0.3
        AChem2_2_1 = 0.15
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3
        Adiff3_3 = 4.83
            
        #Базовые коэффициенты перекрестной кинетической матрицы процессов относительно переноса теплоты
        ADiffHeat2 = 0.09
        ADiffHeat3 = 0.063
        
        #Базовые коэффициенты перекрестной кинетической матрицы переноса теплоты относительно процессов
        AHeatDiff2 = 0.09
        AHeatDiff3 = 0.03
    
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 2.31
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 3.51
        H2_7 = 0.81
    
        #Температура окружающей среды
        Tokr = 200
        
        #Коэффициент теплотдачи в окружающую среду
        AHeatOkr = 105
        
        #Внешние потоки теплоты
        ExtQEnPow1 = 90.0
        
        #Внешние потоки вещества
        xExt1_6 = 30.9
        
        #Потенциал взаимодействия
        mu2_9 = 13.5
        
        #Потенциал взаимодействия между энергетическими степенями свободы
        muBet1_3 = 18
        muBet1_5 = 187
    
        #Метод интегрирования дифференциальных уравнений
        method = "RK23"
        rtol = 2.7*1e-2#Относительная погрешность
        max_step = 3*1e-4#Максимальный шаг интегрирования
        integDynamic = standartIntegrateDyn(method=method,
                                            rtol=rtol,
                                            max_step=max_step)
        
        #Время интегрирования
        Tint = 1e-6
        
        #Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)
        
        #Функция эталонного результата
        def countSystem2(t,allStParameters):
            #Выделение координат состояния и приведенных температур
            _stateCoordinates = allStParameters[0:-2]
            _reducedTemp = allStParameters[-2::]
            
            #Вызов функции времени
            _systemParameters = fU2(t,systemParameters)
            
            #Вызов расчета
            (chemPot,potBet,Temp,Betas,Aff,AffHeat,
             kineticMatrixCPCP,kineticMatrixCPHeat,
             kineticMatrixHeatCP,kineticMatrixHeatHeat,
             vProcesses,vQTransf,vHeatProcesses,
             vQPows,heatTransferMatrix,balanceMatrix,
             vx,vUPows,invC,powH,vT,heatStreams,Streams) = CountSystemQ2(_stateCoordinates,_reducedTemp,_systemParameters,
                                                                         nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                                         nu1_8,nu1_9,nu1_10,nu1_11,nu1_12,nu1_13,nu1_14,
                                                                         nu1_15,nu2_1,nu2_2,nu2_3,nu2_4,nu2_5,nu2_6,
                                                                         nu2_7,nu2_8,nu2_9,beta2_1,beta2_2,Adiff1_1,
                                                                         Adiff1_2,Adiff2_1,ADiffHeat2,AHeatDiff2,invC1,
                                                                         H1_1,H2_7,Tokr,AHeatOkr,ExtQEnPow1,AChem1_4_4,
                                                                         mu2_9,AChem2_3_3,AChem2_1_2,AChem2_2_1,Adiff3_3,
                                                                         ADiffHeat3,AHeatDiff3,xExt1_6,muBet1_3,muBet1_5)
            
            #Выводим результат
            return np.hstack((vx,vT))
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8", "x2_9"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem1_4","vChem2_1","vChem2_2","vChem2_3","vDiff1","vDiff2","vDiff3"]#Имена координат процессов
        energyPowersNames = ["EnPow1","EnPow2","EnOkr"]#Имена энергетических степеней свободы
        reducedTemperaturesEnergyPowersNames = ["T1", "T2"]#Имена приведенных температур энергетических степеней свободы
        energyPowersBetNames = ["EnPowBet1","EnPowBet2","EnPowBet3","EnPowBet4"]#Имена взаимодействий между энергетическими степенями свободы
        heatTransfersNames = ["vQTransf","vQTrOkr"]#Имена потоков переноса теплоты
        heatTransfersOutputEnergyPowersNames = ["EnPow1","EnPow1"]#Имена энергетических степеней свободы, с которых уходит теплота
        heatTransfersInputEnergyPowersNames = ["EnPow2","EnOkr"]#Имена энергетических степеней свободы, на которые приходит теплота
        stateCoordinatesStreamsNames = ["x1_2","x1_6","x2_8"]#Имена координат состояния, изменяемых в результате внешних потоков
        heatEnergyPowersStreamsNames = ["EnPow1","EnOkr"]#Имена потоков теплоты на энергетические степени свободы
        
        #Задаем функцию состояния системы
        stateFunction = CountStateQ2
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        energyPowersVarTemperatureNames = ["EnPow1","EnPow2"]#Имена переменных температур энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames = ["EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow2","EnPow2","EnPow2","EnPow2"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesVarPotentialsInterBetNames = ["x1_2"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        energyPowersVarPotentialsInterBetNames = ["EnPowBet2"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        energyPowersVarBetaNames = ["EnPow1","EnPow2"]#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
        processCoordinatesVarBetaNames = ["vDiff1","vDiff1"]#Имена переменных долей распределения некомпенсированной теплоты координат процессов
        reducedTemperaturesEnergyPowersVarInvHeatCapacityNames = ["T2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersVarInvHeatCapacityNames = ["EnPow2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        reducedTemperaturesEnergyPowersVarHeatEffectNames = ["T1", "T1", "T2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        stateCoordinatesVarHeatEffectNames = ["x1_2", "x1_5", "x2_1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticPCPCNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой координат процессов
        varKineticPCPCAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой термодинамических сил
        varKineticPCHeatNames = ["vDiff1"]#Имена сопряженностей координат процессов с теплопереносами
        varKineticPCHeatAffNames = ["vQTransf"]#Имена сопряженностей термодинамических сил с теплопереносами
        varKineticHeatPCNames = ["vQTransf"]#Имена сопряженностей теплопереносов с координатами процессов
        varKineticHeatPCAffNames = ["vDiff1"]#Имена сопряженностей теплопереносов с термодинамическими силами
        varKineticHeatHeatNames = ["vQTransf"]#Имена сопряженностей между собой перенесенных теплот
        varKineticHeatHeatAffNames = ["vQTransf"]#Имена сопряженностей между собой термодинамических сил по переносу теплот
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = ["x1_2","x2_8"]#Имена переменных внешних потоков
        heatEnergyPowersVarStreamsNames = ["EnOkr"]#Имена переменных внешних потоков теплоты
        
        #Задаем систему
        nonEqSystem = NonEqSystemQ(stateCoordinatesNames,#Имена координат состояния
                                   processCoordinatesNames,#Имена координат процессов
                                   energyPowersNames,#Имена энергетических степеней свободы
                                   reducedTemperaturesEnergyPowersNames,#Имена приведенных температур энергетических степеней свободы
                                   energyPowersBetNames,#Имена взаимодействий между энергетическими степенями свободы
                                   heatTransfersNames,#Имена потоков переноса теплоты
                                   heatTransfersOutputEnergyPowersNames,#Имена энергетических степеней свободы, с которых уходит теплота
                                   heatTransfersInputEnergyPowersNames,#Имена энергетических степеней свободы, на которые приходит теплота
                                   stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                   heatEnergyPowersStreamsNames,#Имена потоков теплоты на энергетические степени свободы
                                    
                                   #Задаем функцию состояния системы
                                   stateFunction,
                                    
                                   #Задаем рассчитываемые параметры системы
                                   stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                   processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                   energyPowersVarTemperatureNames,#Имена переменных температур энергетических степеней свободы
                                   stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                   energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                   stateCoordinatesVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                                   energyPowersVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                                   energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                                   processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                                   reducedTemperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                   reducedTemperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   stateCoordinatesVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                    
                                   #Задаем переменные параметры кинетической матрицы
                                   varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                   varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                   varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                   varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                   varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                   varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                   varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                   varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                    
                                   #Задаем внешние потоки
                                   stateCoordinatesVarStreamsNames,#Имена переменных внешних потоков
                                   heatEnergyPowersVarStreamsNames#Имена переменных внешних потоков теплоты
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
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_3",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_4",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_3",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vDiff2",beta2_1)
        nonEqSystem.SetBetaConstElement("EnPow2","vDiff2",beta2_2)
        nonEqSystem.SetBetaConstElement("EnPow1","vDiff3",beta2_1)
        nonEqSystem.SetBetaConstElement("EnPow2","vDiff3",beta2_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vChem1_4", "vChem1_4", AChem1_4_4)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vChem2_3", "vChem2_3", AChem2_3_3)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vChem2_1", "vChem2_2", AChem2_1_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vChem2_2", "vChem2_1", AChem2_2_1)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff1", "vDiff1", Adiff1_1)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff1", "vDiff2", Adiff1_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff2", "vDiff1", Adiff2_1)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff3", "vDiff3", Adiff3_3)
        nonEqSystem.SetKineticMatrixPCHeatConstElement("vDiff2", "vQTransf", ADiffHeat2)
        nonEqSystem.SetKineticMatrixHeatPCConstElement("vQTransf", "vDiff2", AHeatDiff2)
        nonEqSystem.SetKineticMatrixPCHeatConstElement("vDiff3", "vQTransf", ADiffHeat3)
        nonEqSystem.SetKineticMatrixHeatPCConstElement("vQTransf", "vDiff3", AHeatDiff3)
        nonEqSystem.SetKineticMatrixHeatHeatConstElement("vQTrOkr", "vQTrOkr", AHeatOkr)
        nonEqSystem.SetInvHeatCapacityMatrixConstElement("T1", "EnPow1", invC1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T1", "x1_1", H1_1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T2", "x2_7", H2_7)
        nonEqSystem.SetTEnergyPowersConstElement("EnOkr", Tokr)
        nonEqSystem.SetHeatEnergyPowersStreamsConstElement("EnPow1",ExtQEnPow1)
        nonEqSystem.SetStateCoordinatesStreamsConstElement("x1_6",xExt1_6)
        nonEqSystem.SetPotentialsInterConstElement("EnPow2", "x2_9", -mu2_9)
        nonEqSystem.SetPotentialsInterBetConstElement("EnPowBet1", "x1_5", -muBet1_5)
        nonEqSystem.SetPotentialsInterBetConstElement("EnPowBet3", "x1_3", -muBet1_3)
        
        #Создаем динамику системы
        nonEqSystemDyn = NonEqSystemQDyn(nonEqSystem,#Система
                 
                                         fU2,#Функция условий протекания процессов
                                         
                                         FunCharQ2,#Функция внешних параметров
                                         
                                         integDynamic#Метод интегрирования дифференциальных уравнений
                                         )
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystemQ(0.0*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystemQ(0.3*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystemQ(0.9*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        
        #Выполняем эталонные расчеты
        vxEt1 = countSystem2(0.0*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt2 = countSystem2(0.3*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt3 = countSystem2(0.9*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 5)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 5)#Проверяем скорость
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fSystemDyn(0.0*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx2 = nonEqSystemDyn.fSystemDyn(0.3*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx3 = nonEqSystemDyn.fSystemDyn(0.9*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 5)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 5)#Проверяем скорость
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem2,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates,reducedTemp)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             rtol=rtol,
                                             max_step=max_step,
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t,Char1et,Char2et,Char3et) = FunCharQ2(t,#Моменты времени
                                                stateCoordinatesEt,#Координаты состояния
                                                reducedTempEt,#Приведенные температуры
                                                systemParameters#Параметры системы  
                                                )
        
        #Рассчитываем динамику параметров системы
        (t,Char1,Char2,Char3) = nonEqSystemDyn.NonEqSystemDynamicQ(Tint,
                                                                   stateCoordinates,
                                                                   reducedTemp,
                                                                   systemParameters,
                                                                   t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 7)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 7)
        
        #Функтор модели параметров системы
        class functorSystemPar(object):
            #Функтор параметров
            def SystemParametersCount(self,measParameters):
                return (measParameters[0:12],
                        measParameters[12:14],
                        measParameters[14::])
        
        #Получаем класс модели
        funSystemPar = functorSystemPar()
        model = ModelQ(funSystemPar,
                       nonEqSystemDyn)
        
        #Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              reducedTemp,
                              systemParameters))
        
        #Запускаем динамику
        (t,Char1,Char2,Char3) = model.CountDynamic(Tint,measChar,
                                                   t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 7)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 7)
    def testNonEqDyn6(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81, 80, 64])#Координаты состояния
        reducedTemp = np.array([350, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters = systemParameters = np.array([ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                                                          10,   25,   35,   55,   75,  81,   210,  300,  150,  153,
                                                         123, 15.3, 45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                                                         1.5, 33.9, 27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,
                                                         159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                                                        0.15, 74.1, 81.0, 5.13, 31.8, 3.3, 39.63, 54.3, 27.3, 0.15,
                                                        0.15])#Параметры системы
        
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
        
        #Доли распределения некомпенсированных теплот по энергетических степеням свободы
        beta2_1 = 0.25
        beta2_2 = 0.75
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        AChem1_4_4 = 3.3
        AChem2_3_3 = 6.3
        AChem2_1_2 = 0.3
        AChem2_2_1 = 0.15
        Adiff1_1 = 1.83
        Adiff1_2 = 0.3
        Adiff2_1 = 0.3
        Adiff3_3 = 4.83
            
        #Базовые коэффициенты перекрестной кинетической матрицы процессов относительно переноса теплоты
        ADiffHeat2 = 0.09
        ADiffHeat3 = 0.063
        
        #Базовые коэффициенты перекрестной кинетической матрицы переноса теплоты относительно процессов
        AHeatDiff2 = 0.09
        AHeatDiff3 = 0.03
    
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 2.31
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 3.51
        H2_7 = 0.81
    
        #Температура окружающей среды
        Tokr = 200
        
        #Коэффициент теплотдачи в окружающую среду
        AHeatOkr = 105
        
        #Внешние потоки теплоты
        ExtQEnPow1 = 90.0
        
        #Внешние потоки вещества
        xExt1_6 = 30.9
        
        #Потенциал взаимодействия
        mu2_9 = 13.5
        
        #Потенциал взаимодействия между энергетическими степенями свободы
        muBet1_3 = 18
        muBet1_5 = 187
    
        #Метод интегрирования дифференциальных уравнений
        method = "RK23"
        def funIntegAttributes(systemParameters,method):
            #Определяем аттрибуты
            atol = np.abs(systemParameters[4])*1e-8#Абсолютная погрешность
            rtol = np.abs(systemParameters[14] + 2.1*systemParameters[18])*1e-5#Относительная погрешность
            max_step = np.abs(5.1*systemParameters[33] *systemParameters[12])*1e-3#Максимальный шаг
            
            #Выводим аттрибуты интегрирования
            return (atol,rtol,max_step)
        integDynamic = stepIntegrateDyn(method=method,
                                        funIntegAttributes=funIntegAttributes)
        
        #Время интегрирования
        Tint = 1e-6
        
        #Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)
        
        #Функция эталонного результата
        def countSystem2(t,allStParameters):
            #Выделение координат состояния и приведенных температур
            _stateCoordinates = allStParameters[0:-2]
            _reducedTemp = allStParameters[-2::]
            
            #Вызов функции времени
            _systemParameters = fU2(t,systemParameters)
            
            #Вызов расчета
            (chemPot,potBet,Temp,Betas,Aff,AffHeat,
             kineticMatrixCPCP,kineticMatrixCPHeat,
             kineticMatrixHeatCP,kineticMatrixHeatHeat,
             vProcesses,vQTransf,vHeatProcesses,
             vQPows,heatTransferMatrix,balanceMatrix,
             vx,vUPows,invC,powH,vT,heatStreams,Streams) = CountSystemQ2(_stateCoordinates,_reducedTemp,_systemParameters,
                                                                         nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                                         nu1_8,nu1_9,nu1_10,nu1_11,nu1_12,nu1_13,nu1_14,
                                                                         nu1_15,nu2_1,nu2_2,nu2_3,nu2_4,nu2_5,nu2_6,
                                                                         nu2_7,nu2_8,nu2_9,beta2_1,beta2_2,Adiff1_1,
                                                                         Adiff1_2,Adiff2_1,ADiffHeat2,AHeatDiff2,invC1,
                                                                         H1_1,H2_7,Tokr,AHeatOkr,ExtQEnPow1,AChem1_4_4,
                                                                         mu2_9,AChem2_3_3,AChem2_1_2,AChem2_2_1,Adiff3_3,
                                                                         ADiffHeat3,AHeatDiff3,xExt1_6,muBet1_3,muBet1_5)
            
            #Выводим результат
            return np.hstack((vx,vT))
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8", "x2_9"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3","vChem1_4","vChem2_1","vChem2_2","vChem2_3","vDiff1","vDiff2","vDiff3"]#Имена координат процессов
        energyPowersNames = ["EnPow1","EnPow2","EnOkr"]#Имена энергетических степеней свободы
        reducedTemperaturesEnergyPowersNames = ["T1", "T2"]#Имена приведенных температур энергетических степеней свободы
        energyPowersBetNames = ["EnPowBet1","EnPowBet2","EnPowBet3","EnPowBet4"]#Имена взаимодействий между энергетическими степенями свободы
        heatTransfersNames = ["vQTransf","vQTrOkr"]#Имена потоков переноса теплоты
        heatTransfersOutputEnergyPowersNames = ["EnPow1","EnPow1"]#Имена энергетических степеней свободы, с которых уходит теплота
        heatTransfersInputEnergyPowersNames = ["EnPow2","EnOkr"]#Имена энергетических степеней свободы, на которые приходит теплота
        stateCoordinatesStreamsNames = ["x1_2","x1_6","x2_8"]#Имена координат состояния, изменяемых в результате внешних потоков
        heatEnergyPowersStreamsNames = ["EnPow1","EnOkr"]#Имена потоков теплоты на энергетические степени свободы
        
        #Задаем функцию состояния системы
        stateFunction = CountStateQ2
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        energyPowersVarTemperatureNames = ["EnPow1","EnPow2"]#Имена переменных температур энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x2_1", "x2_4", "x2_7", "x2_8"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames = ["EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow2","EnPow2","EnPow2","EnPow2"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesVarPotentialsInterBetNames = ["x1_2"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        energyPowersVarPotentialsInterBetNames = ["EnPowBet2"]#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        energyPowersVarBetaNames = ["EnPow1","EnPow2"]#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
        processCoordinatesVarBetaNames = ["vDiff1","vDiff1"]#Имена переменных долей распределения некомпенсированной теплоты координат процессов
        reducedTemperaturesEnergyPowersVarInvHeatCapacityNames = ["T2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersVarInvHeatCapacityNames = ["EnPow2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        reducedTemperaturesEnergyPowersVarHeatEffectNames = ["T1", "T1", "T2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        stateCoordinatesVarHeatEffectNames = ["x1_2", "x1_5", "x2_1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticPCPCNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой координат процессов
        varKineticPCPCAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3","vChem2_1","vChem2_2","vDiff2"]#Имена сопряженностей между собой термодинамических сил
        varKineticPCHeatNames = ["vDiff1"]#Имена сопряженностей координат процессов с теплопереносами
        varKineticPCHeatAffNames = ["vQTransf"]#Имена сопряженностей термодинамических сил с теплопереносами
        varKineticHeatPCNames = ["vQTransf"]#Имена сопряженностей теплопереносов с координатами процессов
        varKineticHeatPCAffNames = ["vDiff1"]#Имена сопряженностей теплопереносов с термодинамическими силами
        varKineticHeatHeatNames = ["vQTransf"]#Имена сопряженностей между собой перенесенных теплот
        varKineticHeatHeatAffNames = ["vQTransf"]#Имена сопряженностей между собой термодинамических сил по переносу теплот
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = ["x1_2","x2_8"]#Имена переменных внешних потоков
        heatEnergyPowersVarStreamsNames = ["EnOkr"]#Имена переменных внешних потоков теплоты
        
        #Задаем систему
        nonEqSystem = NonEqSystemQ(stateCoordinatesNames,#Имена координат состояния
                                   processCoordinatesNames,#Имена координат процессов
                                   energyPowersNames,#Имена энергетических степеней свободы
                                   reducedTemperaturesEnergyPowersNames,#Имена приведенных температур энергетических степеней свободы
                                   energyPowersBetNames,#Имена взаимодействий между энергетическими степенями свободы
                                   heatTransfersNames,#Имена потоков переноса теплоты
                                   heatTransfersOutputEnergyPowersNames,#Имена энергетических степеней свободы, с которых уходит теплота
                                   heatTransfersInputEnergyPowersNames,#Имена энергетических степеней свободы, на которые приходит теплота
                                   stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                   heatEnergyPowersStreamsNames,#Имена потоков теплоты на энергетические степени свободы
                                    
                                   #Задаем функцию состояния системы
                                   stateFunction,
                                    
                                   #Задаем рассчитываемые параметры системы
                                   stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                   processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                   energyPowersVarTemperatureNames,#Имена переменных температур энергетических степеней свободы
                                   stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                   energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                   stateCoordinatesVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                                   energyPowersVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                                   energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                                   processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                                   reducedTemperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                   reducedTemperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   stateCoordinatesVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                    
                                   #Задаем переменные параметры кинетической матрицы
                                   varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                   varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                   varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                   varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                   varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                   varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                   varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                   varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                    
                                   #Задаем внешние потоки
                                   stateCoordinatesVarStreamsNames,#Имена переменных внешних потоков
                                   heatEnergyPowersVarStreamsNames#Имена переменных внешних потоков теплоты
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
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_3",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_4",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow2","vChem2_3",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vDiff2",beta2_1)
        nonEqSystem.SetBetaConstElement("EnPow2","vDiff2",beta2_2)
        nonEqSystem.SetBetaConstElement("EnPow1","vDiff3",beta2_1)
        nonEqSystem.SetBetaConstElement("EnPow2","vDiff3",beta2_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vChem1_4", "vChem1_4", AChem1_4_4)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vChem2_3", "vChem2_3", AChem2_3_3)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vChem2_1", "vChem2_2", AChem2_1_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vChem2_2", "vChem2_1", AChem2_2_1)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff1", "vDiff1", Adiff1_1)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff1", "vDiff2", Adiff1_2)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff2", "vDiff1", Adiff2_1)
        nonEqSystem.SetKineticMatrixPCPCConstElement("vDiff3", "vDiff3", Adiff3_3)
        nonEqSystem.SetKineticMatrixPCHeatConstElement("vDiff2", "vQTransf", ADiffHeat2)
        nonEqSystem.SetKineticMatrixHeatPCConstElement("vQTransf", "vDiff2", AHeatDiff2)
        nonEqSystem.SetKineticMatrixPCHeatConstElement("vDiff3", "vQTransf", ADiffHeat3)
        nonEqSystem.SetKineticMatrixHeatPCConstElement("vQTransf", "vDiff3", AHeatDiff3)
        nonEqSystem.SetKineticMatrixHeatHeatConstElement("vQTrOkr", "vQTrOkr", AHeatOkr)
        nonEqSystem.SetInvHeatCapacityMatrixConstElement("T1", "EnPow1", invC1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T1", "x1_1", H1_1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T2", "x2_7", H2_7)
        nonEqSystem.SetTEnergyPowersConstElement("EnOkr", Tokr)
        nonEqSystem.SetHeatEnergyPowersStreamsConstElement("EnPow1",ExtQEnPow1)
        nonEqSystem.SetStateCoordinatesStreamsConstElement("x1_6",xExt1_6)
        nonEqSystem.SetPotentialsInterConstElement("EnPow2", "x2_9", -mu2_9)
        nonEqSystem.SetPotentialsInterBetConstElement("EnPowBet1", "x1_5", -muBet1_5)
        nonEqSystem.SetPotentialsInterBetConstElement("EnPowBet3", "x1_3", -muBet1_3)
        
        #Создаем динамику системы
        nonEqSystemDyn = NonEqSystemQDyn(nonEqSystem,#Система
                 
                                         fU2,#Функция условий протекания процессов
                                         
                                         FunCharQ2,#Функция внешних параметров
                                         
                                         integDynamic#Метод интегрирования дифференциальных уравнений
                                         )
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystemQ(0.0*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystemQ(0.3*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystemQ(0.9*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        
        #Выполняем эталонные расчеты
        vxEt1 = countSystem2(0.0*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt2 = countSystem2(0.3*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt3 = countSystem2(0.9*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 5)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 5)#Проверяем скорость
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fSystemDyn(0.0*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx2 = nonEqSystemDyn.fSystemDyn(0.3*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx3 = nonEqSystemDyn.fSystemDyn(0.9*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 5)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 5)#Проверяем скорость
        
        #Рассчитываем эталонную динамику системы
        (atol,rtol,max_step) = funIntegAttributes(systemParameters,method)
        stateCoordinatesEt = integ.solve_ivp(countSystem2,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates,reducedTemp)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени 
                                             atol=atol,
                                             rtol=rtol,
                                             max_step=max_step,
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t,Char1et,Char2et,Char3et) = FunCharQ2(t,#Моменты времени
                                                stateCoordinatesEt,#Координаты состояния
                                                reducedTempEt,#Приведенные температуры
                                                systemParameters#Параметры системы  
                                                )
        
        #Рассчитываем динамику параметров системы
        (t,Char1,Char2,Char3) = nonEqSystemDyn.NonEqSystemDynamicQ(Tint,
                                                                   stateCoordinates,
                                                                   reducedTemp,
                                                                   systemParameters,
                                                                   t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 7)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 7)
        
        #Функтор модели параметров системы
        class functorSystemPar(object):
            #Функтор параметров
            def SystemParametersCount(self,measParameters):
                return (measParameters[0:12],
                        measParameters[12:14],
                        measParameters[14::])
        
        #Получаем класс модели
        funSystemPar = functorSystemPar()
        model = ModelQ(funSystemPar,
                       nonEqSystemDyn)
        
        #Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              reducedTemp,
                              systemParameters))
        
        #Запускаем динамику
        (t,Char1,Char2,Char3) = model.CountDynamic(Tint,measChar,
                                                   t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 7)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 7)
    def testNonEqDyn7(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213])#Координаты состояния
        reducedTemp = np.array([350])#Приведенные температуры энергетических степеней свободы
        systemParameters = systemParameters = np.array([  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                                                          35,   55,   75,  81,   210,  300,  150,  153,  123, 15.3, 
                                                        45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,  1.5])#Параметры системы
        
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
        
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 2.31
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 3.51
    
        #Метод интегрирования дифференциальных уравнений
        method = "LSODA"
        integDynamic = standartIntegrateDyn(method=method)
        
        #Время интегрирования
        Tint = 1e-6
        
        #Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)
        
        #Функция эталонного результата
        def countSystem3(t,allStParameters):
            #Выделение координат состояния и приведенных температур
            _stateCoordinates = allStParameters[0:-1]
            _reducedTemp = allStParameters[-1::]
            
            #Вызов функции времени
            _systemParameters = fU3(t,systemParameters)
            
            #Вызов расчета
            (chemPot,potBet,Temp,Betas,Aff,AffHeat,
             kineticMatrixCPCP,kineticMatrixCPHeat,
             kineticMatrixHeatCP,kineticMatrixHeatHeat,
             vProcesses,vQTransf,vHeatProcesses,
             vQPows,heatTransferMatrix,balanceMatrix,
             vx,vUPows,invC,powH,vT) = CountSystemQ3(_stateCoordinates,_reducedTemp,_systemParameters,
                                                     nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                     nu1_8,nu1_9,nu1_10,nu1_11,invC1,H1_1)
            
            #Выводим результат
            return np.hstack((vx,vT))
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3"]#Имена координат процессов
        energyPowersNames = ["EnPow1"]#Имена энергетических степеней свободы
        reducedTemperaturesEnergyPowersNames = ["T1"]#Имена приведенных температур энергетических степеней свободы
        energyPowersBetNames = []#Имена взаимодействий между энергетическими степенями свободы
        heatTransfersNames = []#Имена потоков переноса теплоты
        heatTransfersOutputEnergyPowersNames = []#Имена энергетических степеней свободы, с которых уходит теплота
        heatTransfersInputEnergyPowersNames = []#Имена энергетических степеней свободы, на которые приходит теплота
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        heatEnergyPowersStreamsNames = []#Имена потоков теплоты на энергетические степени свободы
        
        #Задаем функцию состояния системы
        stateFunction = CountStateQ3
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        energyPowersVarTemperatureNames = ["EnPow1"]#Имена переменных температур энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames = ["EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow1"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        energyPowersVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        energyPowersVarBetaNames = []#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
        processCoordinatesVarBetaNames = []#Имена переменных долей распределения некомпенсированной теплоты координат процессов
        reducedTemperaturesEnergyPowersVarInvHeatCapacityNames = []#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersVarInvHeatCapacityNames = []#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        reducedTemperaturesEnergyPowersVarHeatEffectNames = ["T1", "T1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        stateCoordinatesVarHeatEffectNames = ["x1_2", "x1_5"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticPCPCNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3"]#Имена сопряженностей между собой координат процессов
        varKineticPCPCAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3"]#Имена сопряженностей между собой термодинамических сил
        varKineticPCHeatNames = []#Имена сопряженностей координат процессов с теплопереносами
        varKineticPCHeatAffNames = []#Имена сопряженностей термодинамических сил с теплопереносами
        varKineticHeatPCNames = []#Имена сопряженностей теплопереносов с координатами процессов
        varKineticHeatPCAffNames = []#Имена сопряженностей теплопереносов с термодинамическими силами
        varKineticHeatHeatNames = []#Имена сопряженностей между собой перенесенных теплот
        varKineticHeatHeatAffNames = []#Имена сопряженностей между собой термодинамических сил по переносу теплот
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []#Имена переменных внешних потоков
        heatEnergyPowersVarStreamsNames = []#Имена переменных внешних потоков теплоты
        
        #Задаем систему
        nonEqSystem = NonEqSystemQ(stateCoordinatesNames,#Имена координат состояния
                                   processCoordinatesNames,#Имена координат процессов
                                   energyPowersNames,#Имена энергетических степеней свободы
                                   reducedTemperaturesEnergyPowersNames,#Имена приведенных температур энергетических степеней свободы
                                   energyPowersBetNames,#Имена взаимодействий между энергетическими степенями свободы
                                   heatTransfersNames,#Имена потоков переноса теплоты
                                   heatTransfersOutputEnergyPowersNames,#Имена энергетических степеней свободы, с которых уходит теплота
                                   heatTransfersInputEnergyPowersNames,#Имена энергетических степеней свободы, на которые приходит теплота
                                   stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                   heatEnergyPowersStreamsNames,#Имена потоков теплоты на энергетические степени свободы
                                     
                                   #Задаем функцию состояния системы
                                   stateFunction,
                                    
                                   #Задаем рассчитываемые параметры системы
                                   stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                   processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                   energyPowersVarTemperatureNames,#Имена переменных температур энергетических степеней свободы
                                   stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                   energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                   stateCoordinatesVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                                   energyPowersVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                                   energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                                   processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                                   reducedTemperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                   reducedTemperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   stateCoordinatesVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                    
                                   #Задаем переменные параметры кинетической матрицы
                                   varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                   varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                   varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                   varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                   varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                   varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                   varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                   varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                    
                                   #Задаем внешние потоки
                                   stateCoordinatesVarStreamsNames,#Имена переменных внешних потоков
                                   heatEnergyPowersVarStreamsNames#Имена переменных внешних потоков теплоты
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
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_3",1.0)
        nonEqSystem.SetInvHeatCapacityMatrixConstElement("T1", "EnPow1", invC1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T1", "x1_1", H1_1)
        
        #Создаем динамику системы
        nonEqSystemDyn = NonEqSystemQDyn(nonEqSystem,#Система
                 
                                         fU3,#Функция условий протекания процессов
                                         
                                         FunCharQ3,#Функция внешних параметров
                                         
                                         integDynamic#Метод интегрирования дифференциальных уравнений
                                         )
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystemQ(0.0*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystemQ(0.3*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystemQ(0.9*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        
        #Выполняем эталонные расчеты
        vxEt1 = countSystem3(0.0*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt2 = countSystem3(0.3*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt3 = countSystem3(0.9*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 5)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 5)#Проверяем скорость
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fSystemDyn(0.0*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx2 = nonEqSystemDyn.fSystemDyn(0.3*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx3 = nonEqSystemDyn.fSystemDyn(0.9*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 5)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 5)#Проверяем скорость
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem3,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates,reducedTemp)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-1::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-1]
        
        #Рассчитываем эталонную динамику параметров системы
        (t,Char1et,Char2et,Char3et) = FunCharQ3(t,#Моменты времени
                                                stateCoordinatesEt,#Координаты состояния
                                                reducedTempEt,#Приведенные температуры
                                                systemParameters#Параметры системы  
                                                )
        
        #Рассчитываем динамику параметров системы
        (t,Char1,Char2,Char3) = nonEqSystemDyn.NonEqSystemDynamicQ(Tint,
                                                                   stateCoordinates,
                                                                   reducedTemp,
                                                                   systemParameters,
                                                                   t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 7)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 7)
        
        #Функтор модели параметров системы
        class functorSystemPar(object):
            #Функтор параметров
            def SystemParametersCount(self,measParameters):
                return (measParameters[0:6],
                        measParameters[6:7],
                        measParameters[7::])
        
        #Получаем класс модели
        funSystemPar = functorSystemPar()
        model = ModelQ(funSystemPar,
                       nonEqSystemDyn)
        
        #Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              reducedTemp,
                              systemParameters))
        
        #Запускаем динамику
        (t,Char1,Char2,Char3) = model.CountDynamic(Tint,measChar,
                                                   t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 7)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 7)
    def testNonEqDyn8(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213])#Координаты состояния
        reducedTemp = np.array([350])#Приведенные температуры энергетических степеней свободы
        systemParameters = systemParameters = np.array([  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                                                          35,   55,   75,  81,   210,  300,  150,  153,  123, 15.3, 
                                                        45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,  1.5])#Параметры системы
        
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
        
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 2.31
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 3.51
    
        #Метод интегрирования дифференциальных уравнений
        method = "LSODA"
        atol = 6.3*1e-8#Абсолютная погрешность
        rtol = 2.7*1e-2#Относительная погрешность
        max_step = 3*1e-4#Максимальный шаг интегрирования
        integDynamic = standartIntegrateDyn(method=method,
                                            atol=atol,
                                            rtol=rtol,
                                            max_step=max_step)
        
        #Время интегрирования
        Tint = 1e-6
        
        #Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)
        
        #Функция эталонного результата
        def countSystem3(t,allStParameters):
            #Выделение координат состояния и приведенных температур
            _stateCoordinates = allStParameters[0:-1]
            _reducedTemp = allStParameters[-1::]
            
            #Вызов функции времени
            _systemParameters = fU3(t,systemParameters)
            
            #Вызов расчета
            (chemPot,potBet,Temp,Betas,Aff,AffHeat,
             kineticMatrixCPCP,kineticMatrixCPHeat,
             kineticMatrixHeatCP,kineticMatrixHeatHeat,
             vProcesses,vQTransf,vHeatProcesses,
             vQPows,heatTransferMatrix,balanceMatrix,
             vx,vUPows,invC,powH,vT) = CountSystemQ3(_stateCoordinates,_reducedTemp,_systemParameters,
                                                     nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                     nu1_8,nu1_9,nu1_10,nu1_11,invC1,H1_1)
            
            #Выводим результат
            return np.hstack((vx,vT))
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3"]#Имена координат процессов
        energyPowersNames = ["EnPow1"]#Имена энергетических степеней свободы
        reducedTemperaturesEnergyPowersNames = ["T1"]#Имена приведенных температур энергетических степеней свободы
        energyPowersBetNames = []#Имена взаимодействий между энергетическими степенями свободы
        heatTransfersNames = []#Имена потоков переноса теплоты
        heatTransfersOutputEnergyPowersNames = []#Имена энергетических степеней свободы, с которых уходит теплота
        heatTransfersInputEnergyPowersNames = []#Имена энергетических степеней свободы, на которые приходит теплота
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        heatEnergyPowersStreamsNames = []#Имена потоков теплоты на энергетические степени свободы
        
        #Задаем функцию состояния системы
        stateFunction = CountStateQ3
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        energyPowersVarTemperatureNames = ["EnPow1"]#Имена переменных температур энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames = ["EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow1"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        energyPowersVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        energyPowersVarBetaNames = []#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
        processCoordinatesVarBetaNames = []#Имена переменных долей распределения некомпенсированной теплоты координат процессов
        reducedTemperaturesEnergyPowersVarInvHeatCapacityNames = []#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersVarInvHeatCapacityNames = []#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        reducedTemperaturesEnergyPowersVarHeatEffectNames = ["T1", "T1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        stateCoordinatesVarHeatEffectNames = ["x1_2", "x1_5"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticPCPCNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3"]#Имена сопряженностей между собой координат процессов
        varKineticPCPCAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3"]#Имена сопряженностей между собой термодинамических сил
        varKineticPCHeatNames = []#Имена сопряженностей координат процессов с теплопереносами
        varKineticPCHeatAffNames = []#Имена сопряженностей термодинамических сил с теплопереносами
        varKineticHeatPCNames = []#Имена сопряженностей теплопереносов с координатами процессов
        varKineticHeatPCAffNames = []#Имена сопряженностей теплопереносов с термодинамическими силами
        varKineticHeatHeatNames = []#Имена сопряженностей между собой перенесенных теплот
        varKineticHeatHeatAffNames = []#Имена сопряженностей между собой термодинамических сил по переносу теплот
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []#Имена переменных внешних потоков
        heatEnergyPowersVarStreamsNames = []#Имена переменных внешних потоков теплоты
        
        #Задаем систему
        nonEqSystem = NonEqSystemQ(stateCoordinatesNames,#Имена координат состояния
                                   processCoordinatesNames,#Имена координат процессов
                                   energyPowersNames,#Имена энергетических степеней свободы
                                   reducedTemperaturesEnergyPowersNames,#Имена приведенных температур энергетических степеней свободы
                                   energyPowersBetNames,#Имена взаимодействий между энергетическими степенями свободы
                                   heatTransfersNames,#Имена потоков переноса теплоты
                                   heatTransfersOutputEnergyPowersNames,#Имена энергетических степеней свободы, с которых уходит теплота
                                   heatTransfersInputEnergyPowersNames,#Имена энергетических степеней свободы, на которые приходит теплота
                                   stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                   heatEnergyPowersStreamsNames,#Имена потоков теплоты на энергетические степени свободы
                                     
                                   #Задаем функцию состояния системы
                                   stateFunction,
                                    
                                   #Задаем рассчитываемые параметры системы
                                   stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                   processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                   energyPowersVarTemperatureNames,#Имена переменных температур энергетических степеней свободы
                                   stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                   energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                   stateCoordinatesVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                                   energyPowersVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                                   energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                                   processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                                   reducedTemperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                   reducedTemperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   stateCoordinatesVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                    
                                   #Задаем переменные параметры кинетической матрицы
                                   varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                   varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                   varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                   varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                   varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                   varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                   varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                   varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                    
                                   #Задаем внешние потоки
                                   stateCoordinatesVarStreamsNames,#Имена переменных внешних потоков
                                   heatEnergyPowersVarStreamsNames#Имена переменных внешних потоков теплоты
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
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_3",1.0)
        nonEqSystem.SetInvHeatCapacityMatrixConstElement("T1", "EnPow1", invC1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T1", "x1_1", H1_1)
        
        #Создаем динамику системы
        nonEqSystemDyn = NonEqSystemQDyn(nonEqSystem,#Система
                 
                                         fU3,#Функция условий протекания процессов
                                         
                                         FunCharQ3,#Функция внешних параметров
                                         
                                         integDynamic#Метод интегрирования дифференциальных уравнений
                                         )
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystemQ(0.0*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystemQ(0.3*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystemQ(0.9*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        
        #Выполняем эталонные расчеты
        vxEt1 = countSystem3(0.0*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt2 = countSystem3(0.3*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt3 = countSystem3(0.9*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 5)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 5)#Проверяем скорость
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fSystemDyn(0.0*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx2 = nonEqSystemDyn.fSystemDyn(0.3*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx3 = nonEqSystemDyn.fSystemDyn(0.9*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 5)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 5)#Проверяем скорость
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem3,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates,reducedTemp)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             atol=atol,
                                             rtol=rtol,
                                             max_step=max_step,
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-1::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-1]
        
        #Рассчитываем эталонную динамику параметров системы
        (t,Char1et,Char2et,Char3et) = FunCharQ3(t,#Моменты времени
                                                stateCoordinatesEt,#Координаты состояния
                                                reducedTempEt,#Приведенные температуры
                                                systemParameters#Параметры системы  
                                                )
        
        #Рассчитываем динамику параметров системы
        (t,Char1,Char2,Char3) = nonEqSystemDyn.NonEqSystemDynamicQ(Tint,
                                                                   stateCoordinates,
                                                                   reducedTemp,
                                                                   systemParameters,
                                                                   t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 7)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 7)
        
        #Функтор модели параметров системы
        class functorSystemPar(object):
            #Функтор параметров
            def SystemParametersCount(self,measParameters):
                return (measParameters[0:6],
                        measParameters[6:7],
                        measParameters[7::])
        
        #Получаем класс модели
        funSystemPar = functorSystemPar()
        model = ModelQ(funSystemPar,
                       nonEqSystemDyn)
        
        #Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              reducedTemp,
                              systemParameters))
        
        #Запускаем динамику
        (t,Char1,Char2,Char3) = model.CountDynamic(Tint,measChar,
                                                   t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 7)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 7)
    def testNonEqDyn9(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213])#Координаты состояния
        reducedTemp = np.array([350])#Приведенные температуры энергетических степеней свободы
        systemParameters = systemParameters = np.array([  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                                                          35,   55,   75,  81,   210,  300,  150,  153,  123, 15.3, 
                                                        45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,  1.5])#Параметры системы
        
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
        
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 2.31
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 3.51
    
        #Метод интегрирования дифференциальных уравнений
        method = "LSODA"
        def funIntegAttributes(systemParameters,method):
            #Определяем аттрибуты
            atol = np.abs(systemParameters[4] + 3.9*systemParameters[1])*1e-9#Абсолютная погрешность
            rtol = np.abs(systemParameters[14])*1e-5#Относительная погрешность
            max_step = np.abs(8.1*systemParameters[24] *systemParameters[12])*1e-4#Максимальный шаг
            
            #Выводим аттрибуты интегрирования
            return (atol,rtol,max_step)
        integDynamic = stepIntegrateDyn(method=method,
                                        funIntegAttributes=funIntegAttributes)
        
        #Время интегрирования
        Tint = 1e-6
        
        #Моменты времени
        t_eval = np.linspace(0.0, Tint, 10000)
        
        #Функция эталонного результата
        def countSystem3(t,allStParameters):
            #Выделение координат состояния и приведенных температур
            _stateCoordinates = allStParameters[0:-1]
            _reducedTemp = allStParameters[-1::]
            
            #Вызов функции времени
            _systemParameters = fU3(t,systemParameters)
            
            #Вызов расчета
            (chemPot,potBet,Temp,Betas,Aff,AffHeat,
             kineticMatrixCPCP,kineticMatrixCPHeat,
             kineticMatrixHeatCP,kineticMatrixHeatHeat,
             vProcesses,vQTransf,vHeatProcesses,
             vQPows,heatTransferMatrix,balanceMatrix,
             vx,vUPows,invC,powH,vT) = CountSystemQ3(_stateCoordinates,_reducedTemp,_systemParameters,
                                                     nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                     nu1_8,nu1_9,nu1_10,nu1_11,invC1,H1_1)
            
            #Выводим результат
            return np.hstack((vx,vT))
        
        #Задаем структуру системы
        stateCoordinatesNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена координат состояния
        processCoordinatesNames = ["vChem1_1","vChem1_2","vChem1_3"]#Имена координат процессов
        energyPowersNames = ["EnPow1"]#Имена энергетических степеней свободы
        reducedTemperaturesEnergyPowersNames = ["T1"]#Имена приведенных температур энергетических степеней свободы
        energyPowersBetNames = []#Имена взаимодействий между энергетическими степенями свободы
        heatTransfersNames = []#Имена потоков переноса теплоты
        heatTransfersOutputEnergyPowersNames = []#Имена энергетических степеней свободы, с которых уходит теплота
        heatTransfersInputEnergyPowersNames = []#Имена энергетических степеней свободы, на которые приходит теплота
        stateCoordinatesStreamsNames = []#Имена координат состояния, изменяемых в результате внешних потоков
        heatEnergyPowersStreamsNames = []#Имена потоков теплоты на энергетические степени свободы
        
        #Задаем функцию состояния системы
        stateFunction = CountStateQ3
        
        #Задаем рассчитываемые параметры системы
        stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
        processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
        energyPowersVarTemperatureNames = ["EnPow1"]#Имена переменных температур энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1_1","x1_2", "x1_3", "x1_4", "x1_5", "x1_6"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames = ["EnPow1","EnPow1","EnPow1","EnPow1","EnPow1","EnPow1"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        stateCoordinatesVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
        energyPowersVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
        energyPowersVarBetaNames = []#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
        processCoordinatesVarBetaNames = []#Имена переменных долей распределения некомпенсированной теплоты координат процессов
        reducedTemperaturesEnergyPowersVarInvHeatCapacityNames = []#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersVarInvHeatCapacityNames = []#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        reducedTemperaturesEnergyPowersVarHeatEffectNames = ["T1", "T1"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        stateCoordinatesVarHeatEffectNames = ["x1_2", "x1_5"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        
        #Задаем переменные параметры кинетической матрицы
        varKineticPCPCNames = ["vChem1_1","vChem1_1","vChem1_2","vChem1_2","vChem1_3"]#Имена сопряженностей между собой координат процессов
        varKineticPCPCAffNames = ["vChem1_1","vChem1_2","vChem1_1","vChem1_2","vChem1_3"]#Имена сопряженностей между собой термодинамических сил
        varKineticPCHeatNames = []#Имена сопряженностей координат процессов с теплопереносами
        varKineticPCHeatAffNames = []#Имена сопряженностей термодинамических сил с теплопереносами
        varKineticHeatPCNames = []#Имена сопряженностей теплопереносов с координатами процессов
        varKineticHeatPCAffNames = []#Имена сопряженностей теплопереносов с термодинамическими силами
        varKineticHeatHeatNames = []#Имена сопряженностей между собой перенесенных теплот
        varKineticHeatHeatAffNames = []#Имена сопряженностей между собой термодинамических сил по переносу теплот
        
        #Задаем внешние потоки
        stateCoordinatesVarStreamsNames = []#Имена переменных внешних потоков
        heatEnergyPowersVarStreamsNames = []#Имена переменных внешних потоков теплоты
        
        #Задаем систему
        nonEqSystem = NonEqSystemQ(stateCoordinatesNames,#Имена координат состояния
                                   processCoordinatesNames,#Имена координат процессов
                                   energyPowersNames,#Имена энергетических степеней свободы
                                   reducedTemperaturesEnergyPowersNames,#Имена приведенных температур энергетических степеней свободы
                                   energyPowersBetNames,#Имена взаимодействий между энергетическими степенями свободы
                                   heatTransfersNames,#Имена потоков переноса теплоты
                                   heatTransfersOutputEnergyPowersNames,#Имена энергетических степеней свободы, с которых уходит теплота
                                   heatTransfersInputEnergyPowersNames,#Имена энергетических степеней свободы, на которые приходит теплота
                                   stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                                   heatEnergyPowersStreamsNames,#Имена потоков теплоты на энергетические степени свободы
                                     
                                   #Задаем функцию состояния системы
                                   stateFunction,
                                    
                                   #Задаем рассчитываемые параметры системы
                                   stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                                   processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                                   energyPowersVarTemperatureNames,#Имена переменных температур энергетических степеней свободы
                                   stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                   energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                   stateCoordinatesVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                                   energyPowersVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                                   energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                                   processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                                   reducedTemperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                   reducedTemperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                   stateCoordinatesVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                    
                                   #Задаем переменные параметры кинетической матрицы
                                   varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                   varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                   varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                   varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                   varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                   varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                   varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                   varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                    
                                   #Задаем внешние потоки
                                   stateCoordinatesVarStreamsNames,#Имена переменных внешних потоков
                                   heatEnergyPowersVarStreamsNames#Имена переменных внешних потоков теплоты
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
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_1",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_2",1.0)
        nonEqSystem.SetBetaConstElement("EnPow1","vChem1_3",1.0)
        nonEqSystem.SetInvHeatCapacityMatrixConstElement("T1", "EnPow1", invC1)
        nonEqSystem.SetHeatEffectMatrixConstElement("T1", "x1_1", H1_1)
        
        #Создаем динамику системы
        nonEqSystemDyn = NonEqSystemQDyn(nonEqSystem,#Система
                 
                                         fU3,#Функция условий протекания процессов
                                         
                                         FunCharQ3,#Функция внешних параметров
                                         
                                         integDynamic#Метод интегрирования дифференциальных уравнений
                                         )
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fNonEqSystemQ(0.0*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx2 = nonEqSystemDyn.fNonEqSystemQ(0.3*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        vx3 = nonEqSystemDyn.fNonEqSystemQ(0.9*Tint,
                                           stateCoordinates,
                                           reducedTemp,
                                           systemParameters)
        
        #Выполняем эталонные расчеты
        vxEt1 = countSystem3(0.0*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt2 = countSystem3(0.3*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        vxEt3 = countSystem3(0.9*Tint,
                             np.hstack((stateCoordinates,reducedTemp)))
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 5)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 5)#Проверяем скорость
        
        #Выполняем расчеты
        vx1 = nonEqSystemDyn.fSystemDyn(0.0*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx2 = nonEqSystemDyn.fSystemDyn(0.3*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        vx3 = nonEqSystemDyn.fSystemDyn(0.9*Tint,
                                        np.hstack((stateCoordinates,reducedTemp)),
                                        systemParameters)
        
        #Проверяем значения
        deltaVXPows1 = np.max(np.abs(vx1 - vxEt1))
        self.assertAlmostEqual(deltaVXPows1, 0.0, 5)#Проверяем скорость
        deltaVXPows2 = np.max(np.abs(vx2 - vxEt2))
        self.assertAlmostEqual(deltaVXPows2, 0.0, 5)#Проверяем скорость
        deltaVXPows3 = np.max(np.abs(vx3 - vxEt3))
        self.assertAlmostEqual(deltaVXPows3, 0.0, 5)#Проверяем скорость
        
        #Рассчитываем эталонную динамику системы
        (atol,rtol,max_step) = funIntegAttributes(systemParameters,method)
        stateCoordinatesEt = integ.solve_ivp(countSystem3,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates,reducedTemp)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             atol=atol,
                                             rtol=rtol,
                                             max_step=max_step,
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-1::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-1]
        
        #Рассчитываем эталонную динамику параметров системы
        (t,Char1et,Char2et,Char3et) = FunCharQ3(t,#Моменты времени
                                                stateCoordinatesEt,#Координаты состояния
                                                reducedTempEt,#Приведенные температуры
                                                systemParameters#Параметры системы  
                                                )
        
        #Рассчитываем динамику параметров системы
        (t,Char1,Char2,Char3) = nonEqSystemDyn.NonEqSystemDynamicQ(Tint,
                                                                   stateCoordinates,
                                                                   reducedTemp,
                                                                   systemParameters,
                                                                   t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 7)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 7)
        
        #Функтор модели параметров системы
        class functorSystemPar(object):
            #Функтор параметров
            def SystemParametersCount(self,measParameters):
                return (measParameters[0:6],
                        measParameters[6:7],
                        measParameters[7::])
        
        #Получаем класс модели
        funSystemPar = functorSystemPar()
        model = ModelQ(funSystemPar,
                       nonEqSystemDyn)
        
        #Наблюдаемые параметры
        measChar = np.hstack((stateCoordinates,
                              reducedTemp,
                              systemParameters))
        
        #Запускаем динамику
        (t,Char1,Char2,Char3) = model.CountDynamic(Tint,measChar,
                                                   t_eval = t_eval)
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(Char1 - Char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 9)
        deltaChar2 = np.max(np.abs(Char2 - Char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 7)
        deltaChar3 = np.max(np.abs(Char3 - Char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 7)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
