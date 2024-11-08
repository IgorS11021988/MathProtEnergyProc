import numpy as np
from scipy import integrate as integ
from MathProtEnergyProc.IntegrateDyn import standartIntegrateDyn

from MathProtEnergyProc import NonEqSystemQ
from MathProtEnergyProc import NonEqSystemQDyn
from MathProtEnergyProc import CountDynamicsQ, VectorModelQ

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
class TestNonEqSystemComExpModelQ(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testNonEqDyn1(self):
        #Исходные данные
        stateCoordinates1 = np.array([10, 23, 55, 75, 121, 213, 39, 45, 93, 83])#Координаты состояния
        reducedTemp1 = np.array([350, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters1 = np.array([ 0.3,  0.5,   10,   20,   50,  70,   100,  200,   30,   45,
                                        10,   25,   35,   55,   75,  81,   210,  200,  150,  153,
                                       123, 15.3, 45.9, 5.63,   10, 6.5,  45.4,  615, 21.3,  4.5,
                                       1.5, 33.9, 27.3, 4.53, 21.9, 7.5,  52.6, 21.9, 45.3, 33.9,
                                       159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                                      0.15, 74.1])#Параметры системы
        stateCoordinates2 = np.array([12, 20, 55, 75, 141, 213, 39, 45, 93, 81])#Координаты состояния
        reducedTemp2 = np.array([350, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters2 = np.array([ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                                        10,   25,   35,   35,   75,  81,   210,  310,  150,  153,
                                       123, 15.3, 45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                                       1.5, 33.9, 27.3, 4.53, 21.9, 2.5,  43.6, 21.9, 45.3, 33.9,
                                       159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                                      0.15, 74.1])#Параметры системы
        stateCoordinates3 = np.array([10, 20, 55, 75, 214, 213, 39, 45, 93, 71])#Координаты состояния
        reducedTemp3 = np.array([350, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters3 = np.array([ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                                        10,   25,   35,   55,   65,  81,   210,  150,  150,  153,
                                       123, 15.3, 45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                                       1.5, 33.9, 27.3, 4.53, 21.9, 3.5,  53.6, 21.7, 45.3, 33.9,
                                       159, 21.9, 43.5, 7.53, 37.5, 2.5, 120.3, 73.5, 15.3, 21.3,
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
        def countSystem1(t,allStParameters,systemParameters):
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
        def countSystem1_1(t,allStParameters):
            return countSystem1(t,allStParameters,systemParameters1)
        def countSystem1_2(t,allStParameters):
            return countSystem1(t,allStParameters,systemParameters2)
        def countSystem1_3(t,allStParameters):
            return countSystem1(t,allStParameters,systemParameters3)
        
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
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem1_1,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates1,reducedTemp1)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t1,Char1et1,Char2et1) = FunCharQ1(t,#Моменты времени
                                           stateCoordinatesEt,#Координаты состояния
                                           reducedTempEt,#Приведенные температуры
                                           systemParameters1#Параметры системы  
                                           )
        (char1et1,char2et1) = (np.sum(Char1et1*t1),np.sum(Char2et1*t1))
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem1_2,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates2,reducedTemp2)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t2,Char1et2,Char2et2) = FunCharQ1(t,#Моменты времени
                                           stateCoordinatesEt,#Координаты состояния
                                           reducedTempEt,#Приведенные температуры
                                           systemParameters2#Параметры системы  
                                           )
        (char1et2,char2et2) = (np.sum(Char1et2*t2),np.sum(Char2et2*t2))
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem1_3,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates3,reducedTemp3)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t3,Char1et3,Char2et3) = FunCharQ1(t,#Моменты времени
                                           stateCoordinatesEt,#Координаты состояния
                                           reducedTempEt,#Приведенные температуры
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
        cDyn = CountDynamicsQ(nonEqSystemDyn,
                              CompressFunction)
        
        #Функтор модели параметров системы
        class functorSystemPar(object):
            #Функтор параметров
            def SystemParametersCount(self,measParameters):
                return (measParameters[:,0:10],
                        measParameters[:,10:12],
                        measParameters[:,12::])
        
        #Создаем класс векторной модели
        funSysPar = functorSystemPar()
        cDynVecQ = VectorModelQ(funSysPar,cDyn)
        
        #Рассчитываем динамику параметров системы
        stateCoordinates = np.vstack([stateCoordinates1,
                                      stateCoordinates2,
                                      stateCoordinates3])
        reducedTemp = np.vstack([reducedTemp1,
                                 reducedTemp2,
                                 reducedTemp3])
        systemParameters = np.vstack([systemParameters1,
                                      systemParameters2,
                                      systemParameters3])
        t_evals = [t_eval,
                   t_eval,
                   t_eval]
        Tints = np.array([Tint,
                          Tint])
        dyns = cDynVecQ.ModelingDynamicsQ(Tints,
                                          np.hstack((stateCoordinates,
                                                     reducedTemp,
                                                     systemParameters)),
                                          t_evals = t_evals)
        (char1,char2) = dyns
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(char1 - char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 7)
        self.assertTrue(char1.shape == char1et.shape)
        deltaChar2 = np.max(np.abs(char2 - char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 4)
        self.assertTrue(char2.shape == char2et.shape)
    def testNonEqDyn2(self):
        #Исходные данные
        stateCoordinates1 = np.array([11, 21, 55, 75, 111, 213, 39, 47, 93, 81, 80, 64])#Координаты состояния
        reducedTemp1 = np.array([340, 413])#Приведенные температуры энергетических степеней свободы
        systemParameters1 = systemParameters = np.array([ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                                                           10,   25,   35,   55,   75,  81,   210,  300,  150,  153,
                                                          123, 15.3, 44.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                                                          1.5, 33.9, 27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,
                                                          159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                                                         0.15, 74.1, 81.0, 5.13, 31.8, 3.3, 39.63, 54.3, 27.3, 0.15,
                                                         0.15])#Параметры системы
        stateCoordinates2 = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81, 80, 64])#Координаты состояния
        reducedTemp2 = np.array([350, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters2 = systemParameters = np.array([ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                                                           10,   25,   35,   55,   75,  81,   210,  300,  150,  153,
                                                          123, 15.3, 45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                                                          1.5, 36.9, 27.6, 4.53, 21.5, 7.5,  53.6, 21.9, 45.3, 33.9,
                                                          159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                                                         0.15, 74.1, 81.0, 5.13, 31.8, 3.3, 39.63, 54.3, 27.3, 0.15,
                                                         0.15])#Параметры системы
        stateCoordinates3 = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81, 80, 64])#Координаты состояния
        reducedTemp3 = np.array([390, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters3 = systemParameters = np.array([ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                                                           10,   25,   35,   55,   75,  81,   210,  300,  150,  153,
                                                          123, 15.3, 41.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                                                          1.5, 33.9, 27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,
                                                          159, 21.9, 42.5, 7.53, 31.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                                                         0.15, 74.1, 81.0, 5.13, 31.8, 3.3, 39.63, 54.3, 27.3, 0.15,
                                                         0.15])#Параметры системы
        stateCoordinates4 = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81, 80, 54])#Координаты состояния
        reducedTemp4 = np.array([310, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters4 = systemParameters = np.array([ 0.3,  0.7,   10,   20,   20,  70,   100,  200,   30,   45,
                                                           10,   25,   35,   55,   75,  81,   210,  300,  150,  153,
                                                          123, 11.3, 45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                                                          1.5, 33.9, 27.3, 4.53, 20.9, 7.5,  53.6, 21.9, 45.3, 33.9,
                                                          159, 21.9, 43.5, 7.53, 37.5, 4.5, 120.3, 73.5, 15.3, 21.3,
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
        def countSystem2(t,allStParameters,systemParameters):
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
        def countSystem2_1(t,allStParameters):
            return countSystem2(t,allStParameters,systemParameters1)
        def countSystem2_2(t,allStParameters):
            return countSystem2(t,allStParameters,systemParameters2)
        def countSystem2_3(t,allStParameters):
            return countSystem2(t,allStParameters,systemParameters3)
        def countSystem2_4(t,allStParameters):
            return countSystem2(t,allStParameters,systemParameters4)
        
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
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem2_1,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates1,reducedTemp1)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t1,Char1et1,Char2et1,Char3et1) = FunCharQ2(t,#Моменты времени
                                                    stateCoordinatesEt,#Координаты состояния
                                                    reducedTempEt,#Приведенные температуры
                                                    systemParameters1#Параметры системы  
                                                    )
        (char1et1,char2et1,char3et1) = (np.sum(t1*Char1et1*Char3et1),np.sum(t1*Char2et1 + Char3et1),np.sum(t1*Char1et1))
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem2_2,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates2,reducedTemp2)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t2,Char1et2,Char2et2,Char3et2) = FunCharQ2(t,#Моменты времени
                                                    stateCoordinatesEt,#Координаты состояния
                                                    reducedTempEt,#Приведенные температуры
                                                    systemParameters2#Параметры системы  
                                                    )
        (char1et2,char2et2,char3et2) = (np.sum(t2*Char1et2*Char3et2),np.sum(t2*Char2et2 + Char3et2),np.sum(t2*Char1et2))
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem2_3,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates3,reducedTemp3)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t3,Char1et3,Char2et3,Char3et3) = FunCharQ2(t,#Моменты времени
                                                    stateCoordinatesEt,#Координаты состояния
                                                    reducedTempEt,#Приведенные температуры
                                                    systemParameters3#Параметры системы  
                                                    )
        (char1et3,char2et3,char3et3) = (np.sum(t3*Char1et3*Char3et3),np.sum(t3*Char2et3 + Char3et3),np.sum(t3*Char1et3))
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem2_4,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates4,reducedTemp4)),#Начальное состояние
                                             t_eval=t_eval,#Моменты времени
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-2::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-2]
        
        #Рассчитываем эталонную динамику параметров системы
        (t4,Char1et4,Char2et4,Char3et4) = FunCharQ2(t,#Моменты времени
                                                    stateCoordinatesEt,#Координаты состояния
                                                    reducedTempEt,#Приведенные температуры
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
        cDyn = CountDynamicsQ(nonEqSystemDyn,
                              CompressFunction)
        
        #Функтор модели параметров системы
        class functorSystemPar(object):
            #Функтор параметров
            def SystemParametersCount(self,measParameters):
                return (measParameters[:,0:12],
                        measParameters[:,12:14],
                        measParameters[:,14::])
        
        #Создаем класс векторной модели
        funSysPar = functorSystemPar()
        cDynVecQ = VectorModelQ(funSysPar,cDyn)
        
        #Рассчитываем динамику параметров системы
        stateCoordinates = np.vstack([stateCoordinates1,
                                      stateCoordinates2,
                                      stateCoordinates3,
                                      stateCoordinates4])
        reducedTemp = np.vstack([reducedTemp1,
                                 reducedTemp2,
                                 reducedTemp3,
                                 reducedTemp4])
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
        dyns = cDynVecQ.ModelingDynamicsQ(Tints,
                                          np.hstack((stateCoordinates,
                                                     reducedTemp,
                                                     systemParameters)),
                                          t_evals = t_evals)
        (char1,char2,char3) = dyns
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(char1 - char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 2)
        self.assertTrue(char1.shape == char1et.shape)
        deltaChar2 = np.max(np.abs(char2 - char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, -1)
        self.assertTrue(char2.shape == char2et.shape)
        deltaChar3 = np.max(np.abs(char3 - char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)
        self.assertTrue(char3.shape == char3et.shape)
    def testNonEqDyn3(self):
        #Исходные данные
        stateCoordinates1 = np.array([10, 20, 55, 75, 111, 213])#Координаты состояния
        reducedTemp1 = np.array([350])#Приведенные температуры энергетических степеней свободы
        systemParameters1 = systemParameters = np.array([  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                                                           35,   55,   75,  81,   210,  300,  150,  153,  123, 15.3, 
                                                         45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,  1.5])#Параметры системы
        stateCoordinates2 = np.array([10, 20, 55, 75, 111, 213])#Координаты состояния
        reducedTemp2 = np.array([350])#Приведенные температуры энергетических степеней свободы
        systemParameters2 = systemParameters = np.array([  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
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
        
        #Функция эталонного результата
        def countSystem3(t,allStParameters,systemParameters):
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
        def countSystem3_1(t,allStParameters):
            return countSystem3(t,allStParameters,systemParameters1)
        def countSystem3_2(t,allStParameters):
            return countSystem3(t,allStParameters,systemParameters2)
        
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
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem3_1,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates1,reducedTemp1)),#Начальное состояние
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-1::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-1]
        
        #Рассчитываем эталонную динамику параметров системы
        (t1,Char1et1,Char2et1,Char3et1) = FunCharQ3(t,#Моменты времени
                                                    stateCoordinatesEt,#Координаты состояния
                                                    reducedTempEt,#Приведенные температуры
                                                    systemParameters1#Параметры системы  
                                                    )
        (char1et1,char2et1,char3et1) = (np.sum(t1*Char1et1 + Char3et1),np.sum(t1*Char2et1*Char3et1),np.sum(t1*Char2et1))
        
        #Рассчитываем эталонную динамику системы
        stateCoordinatesEt = integ.solve_ivp(countSystem3_2,
                                             (0.0, Tint),
                                             np.hstack((stateCoordinates2,reducedTemp2)),#Начальное состояние
                                             method = method#Метод численного интегрирования
                                             )#Интегрирование дифференциальных уравнений
        t = stateCoordinatesEt.t.reshape(-1,1)#Моменты времени
        stateCoordinatesEt = stateCoordinatesEt.y.transpose()#Координаты состояния
        
        #Выделяем координаты состояния и приведенную температуру
        reducedTempEt = stateCoordinatesEt[:,-1::]
        stateCoordinatesEt = stateCoordinatesEt[:,0:-1]
        
        #Рассчитываем эталонную динамику параметров системы
        (t2,Char1et2,Char2et2,Char3et2) = FunCharQ3(t,#Моменты времени
                                                    stateCoordinatesEt,#Координаты состояния
                                                    reducedTempEt,#Приведенные температуры
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
        cDyn = CountDynamicsQ(nonEqSystemDyn,
                              CompressFunction)
        
        #Функтор модели параметров системы
        class functorSystemPar(object):
            #Функтор параметров
            def SystemParametersCount(self,measParameters):
                return (measParameters[:,0:6],
                        measParameters[:,6:7],
                        measParameters[:,7::])
        
        #Создаем класс векторной модели
        funSysPar = functorSystemPar()
        cDynVecQ = VectorModelQ(funSysPar,cDyn)
        
        #Рассчитываем динамику параметров системы
        stateCoordinates = np.vstack([stateCoordinates1,
                                      stateCoordinates2])
        reducedTemp = np.vstack([reducedTemp1,
                                 reducedTemp2])
        systemParameters = np.vstack([systemParameters1,
                                      systemParameters2])
        Tints = np.array([Tint,
                          Tint])
        dyns = cDynVecQ.ModelingDynamicsQ(Tints,
                                          np.hstack((stateCoordinates,
                                                     reducedTemp,
                                                     systemParameters)))
        (char1,char2,char3) = dyns
        
        #Проверяем рассчитанные динамики выходных характеристик системы
        deltaChar1 = np.max(np.abs(char1 - char1et))
        self.assertAlmostEqual(deltaChar1, 0.0, 0)
        self.assertTrue(char1.shape == char1et.shape)
        deltaChar2 = np.max(np.abs(char2 - char2et))
        self.assertAlmostEqual(deltaChar2, 0.0, 2)
        self.assertTrue(char2.shape == char2et.shape)
        deltaChar3 = np.max(np.abs(char3 - char3et))
        self.assertAlmostEqual(deltaChar3, 0.0, 9)
        self.assertTrue(char3.shape == char3et.shape)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
