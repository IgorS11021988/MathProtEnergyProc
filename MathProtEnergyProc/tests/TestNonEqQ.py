import numpy as np

from MathProtEnergyProc import NonEqSystemQ

from MathProtEnergyProc.tests.UnitTestExamples.TestNonEqQ1 import *
from MathProtEnergyProc.tests.UnitTestExamples.TestNonEqQ2 import *
from MathProtEnergyProc.tests.UnitTestExamples.TestNonEqQ3 import *

import unittest

#Модульные тесты
class TestNonEqSystemQ(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testNonEq1(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81])#Координаты состояния
        reducedTemp = np.array([350, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters = [ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                              10,   25,   35,   55,   75,  81,   210,  300,  150,  153,
                             123, 15.3, 45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                             1.5, 33.9, 27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,
                             159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                            0.15, 74.1]#Параметры системы
        
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
    
        #Эталонный результат
        (chemPot,potBet,Temp,Betas,Aff,AffHeat,
         kineticMatrixCPCP,kineticMatrixCPHeat,
         kineticMatrixHeatCP,kineticMatrixHeatHeat,
         vProcesses,vQTransf,vHeatProcesses,
         vQPows,heatTransferMatrix,balanceMatrix,
         vx,vUPows,invC,powH,vT) = CountSystemQ1(stateCoordinates,reducedTemp,systemParameters,
                                                 nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                 nu1_8,nu1_9,nu1_10,nu1_11,nu2_1,nu2_2,nu2_3,
                                                 nu2_4,nu2_5,nu2_6,beta2_1,beta2_2,Adiff1_1,
                                                 Adiff1_2,Adiff2_1,ADiffHeat1,ADiffHeat2,
                                                 AHeatDiff1,AHeatDiff2,invC1,H1_1,H2_7)
        
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
        
        #Выполняем расчеты
        nonEqSystem.CountSystem(stateCoordinates,
                                reducedTemp,
                                systemParameters)
        
        #Проверяем значения
        self.assertEqual(chemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(potBet.tolist(), nonEqSystem.GetPotentialsInterBet().tolist())#Проверяем потенциалы взаимодействия между энергетическими степенями свободы
        self.assertEqual(Temp.tolist(), nonEqSystem.GetTEnergyPowers().tolist())#Проверяем температуры энергетических степеней свободы
        self.assertEqual(Betas.tolist(), nonEqSystem.GetBeta().tolist())#Проверяем доли распределения теплоты между энергетическими степенями свободы
        self.assertEqual(balanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(heatTransferMatrix.tolist(), nonEqSystem.GetHeatTransferMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(kineticMatrixCPCP.tolist(), nonEqSystem.GetKineticMatrixPCPC().tolist())#Проверяем главный блок кинетической матрицы
        self.assertEqual(kineticMatrixCPHeat.tolist(), nonEqSystem.GetKineticMatrixPCHeat().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatCP.tolist(), nonEqSystem.GetKineticMatrixHeatPC().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatHeat.tolist(), nonEqSystem.GetKineticMatrixHeatHeat().tolist())#Проверяем главный блок кинетической матрицы по теплообмену
        self.assertEqual(invC.tolist(), nonEqSystem.GetInvHeatCapacityMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(powH.tolist(), nonEqSystem.GetHeatEffectMatrix().tolist())#Проверяем матрицу переноса теплоты
        deltaAff = np.max(np.abs(Aff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaAffHeat = np.max(np.abs(AffHeat - nonEqSystem.GetHeatAffinity()))
        self.assertAlmostEqual(deltaAffHeat, 0.0, 9)#Проверяем термодинамические силы переноса теплотытеплообмену
        deltaVProcesses = np.max(np.abs(vProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVQTransf = np.max(np.abs(vQTransf - nonEqSystem.GetVHeatTransfers()))
        self.assertAlmostEqual(deltaVQTransf, 0.0, 9)#Проверяем скорость переноса теплоты
        deltaVHeatProcesses = np.max(np.abs(vHeatProcesses - nonEqSystem.GetVHeatProcess()))
        self.assertAlmostEqual(deltaVHeatProcesses, 0.0, 6)#Проверяем скорость выделения некомпенсированной теплоты
        deltaVQPows = np.max(np.abs(vQPows - nonEqSystem.GetVHeatPower()))
        self.assertAlmostEqual(deltaVQPows, 0.0, 6)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVXPows = np.max(np.abs(vx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 9)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVUPows = np.max(np.abs(vUPows - nonEqSystem.GetVUPower()))
        self.assertAlmostEqual(deltaVUPows, 0.0, 7)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
        deltaVRTPows = np.max(np.abs(vT - nonEqSystem.GetVReducedTemperaturesEnergyPowers()))
        self.assertAlmostEqual(deltaVRTPows, 0.0, 6)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
    def testNonEq2(self):
        #Исходные данные
        stateCoordinates = np.array([14, 45, 57, 90, 141, 183, 43, 41, 85, 79])#Координаты состояния
        reducedTemp = np.array([320, 723])#Приведенные температуры энергетических степеней свободы
        systemParameters = [ 0.4,  0.6,   13,   21,   57,  71,   104,  203,   32,   48,
                              11,   23,   37,   53,   75,  85,   240,  310,  183,  159,
                             133, 16.3, 47.9, 8.63,  140, 7.5,  55.4,  415, 31.3, 21.5,
                             1.6, 35.9, 28.3, 4.43, 22.9, 8.5,  52.6, 23.9, 47.3, 35.9,
                             189, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                            0.35, 77.1]#Параметры системы
        
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
        
        #Доли распределения некомпенсированных теплот по энергетических степеням свободы
        beta2_1 = 0.35
        beta2_2 = 0.65
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        Adiff1_1 = 4.83
        Adiff1_2 = 0.45
        Adiff2_1 = 0.45
            
        #Базовые коэффициенты перекрестной кинетической матрицы процессов относительно переноса теплоты
        ADiffHeat1 = 0.25
        ADiffHeat2 = 0.11
        
        #Базовые коэффициенты перекрестной кинетической матрицы переноса теплоты относительно процессов
        AHeatDiff1 = 0.35
        AHeatDiff2 = 0.08
    
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 5.31
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 6.51
        H2_7 = 3.81
    
        #Эталонный результат
        (chemPot,potBet,Temp,Betas,Aff,AffHeat,
         kineticMatrixCPCP,kineticMatrixCPHeat,
         kineticMatrixHeatCP,kineticMatrixHeatHeat,
         vProcesses,vQTransf,vHeatProcesses,
         vQPows,heatTransferMatrix,balanceMatrix,
         vx,vUPows,invC,powH,vT) = CountSystemQ1(stateCoordinates,reducedTemp,systemParameters,
                                                 nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                 nu1_8,nu1_9,nu1_10,nu1_11,nu2_1,nu2_2,nu2_3,
                                                 nu2_4,nu2_5,nu2_6,beta2_1,beta2_2,Adiff1_1,
                                                 Adiff1_2,Adiff2_1,ADiffHeat1,ADiffHeat2,
                                                 AHeatDiff1,AHeatDiff2,invC1,H1_1,H2_7)
        
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
        
        #Выполняем расчеты
        nonEqSystem.CountSystem(stateCoordinates,
                                reducedTemp,
                                systemParameters)
        
        #Проверяем значения
        self.assertEqual(chemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(potBet.tolist(), nonEqSystem.GetPotentialsInterBet().tolist())#Проверяем потенциалы взаимодействия между энергетическими степенями свободы
        self.assertEqual(Temp.tolist(), nonEqSystem.GetTEnergyPowers().tolist())#Проверяем температуры энергетических степеней свободы
        self.assertEqual(Betas.tolist(), nonEqSystem.GetBeta().tolist())#Проверяем доли распределения теплоты между энергетическими степенями свободы
        self.assertEqual(balanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(heatTransferMatrix.tolist(), nonEqSystem.GetHeatTransferMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(kineticMatrixCPCP.tolist(), nonEqSystem.GetKineticMatrixPCPC().tolist())#Проверяем главный блок кинетической матрицы
        self.assertEqual(kineticMatrixCPHeat.tolist(), nonEqSystem.GetKineticMatrixPCHeat().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatCP.tolist(), nonEqSystem.GetKineticMatrixHeatPC().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatHeat.tolist(), nonEqSystem.GetKineticMatrixHeatHeat().tolist())#Проверяем главный блок кинетической матрицы по теплообмену
        self.assertEqual(invC.tolist(), nonEqSystem.GetInvHeatCapacityMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(powH.tolist(), nonEqSystem.GetHeatEffectMatrix().tolist())#Проверяем матрицу переноса теплоты
        deltaAff = np.max(np.abs(Aff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaAffHeat = np.max(np.abs(AffHeat - nonEqSystem.GetHeatAffinity()))
        self.assertAlmostEqual(deltaAffHeat, 0.0, 9)#Проверяем термодинамические силы переноса теплотытеплообмену
        deltaVProcesses = np.max(np.abs(vProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVQTransf = np.max(np.abs(vQTransf - nonEqSystem.GetVHeatTransfers()))
        self.assertAlmostEqual(deltaVQTransf, 0.0, 9)#Проверяем скорость переноса теплоты
        deltaVHeatProcesses = np.max(np.abs(vHeatProcesses - nonEqSystem.GetVHeatProcess()))
        self.assertAlmostEqual(deltaVHeatProcesses, 0.0, 6)#Проверяем скорость выделения некомпенсированной теплоты
        deltaVQPows = np.max(np.abs(vQPows - nonEqSystem.GetVHeatPower()))
        self.assertAlmostEqual(deltaVQPows, 0.0, 6)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVXPows = np.max(np.abs(vx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 9)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVUPows = np.max(np.abs(vUPows - nonEqSystem.GetVUPower()))
        self.assertAlmostEqual(deltaVUPows, 0.0, 6)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
        deltaVRTPows = np.max(np.abs(vT - nonEqSystem.GetVReducedTemperaturesEnergyPowers()))
        self.assertAlmostEqual(deltaVRTPows, 0.0, 5)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
    def testNonEq3(self):
        #Исходные данные
        stateCoordinates = np.array([34, 14, 57, 81, 441, 18, 39, 43, 82, 73])#Координаты состояния
        reducedTemp = np.array([320, 723])#Приведенные температуры энергетических степеней свободы
        systemParameters = [ 0.1,  0.9,   13,   21,   57,  71,   104,  203,   32,   48,
                              14,   23,   37,   53,   75,  85,   240,  310,  183,  159,
                             133, 16.3, 47.9, 7.63,  140, 7.5,  35.4,  415, 31.3, 21.5,
                             1.6, 32.9, 28.3, 4.73, 22.9, 6.5,  52.6, 21.9, 47.3, 35.9,
                             179, 23.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 14.7, 21.3,
                            0.35, 77.1]#Параметры системы
        
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
        
        #Доли распределения некомпенсированных теплот по энергетических степеням свободы
        beta2_1 = 0.15
        beta2_2 = 0.85
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        Adiff1_1 = 4.83
        Adiff1_2 = 0.15
        Adiff2_1 = 0.45
            
        #Базовые коэффициенты перекрестной кинетической матрицы процессов относительно переноса теплоты
        ADiffHeat1 = 0.25
        ADiffHeat2 = 0.05
        
        #Базовые коэффициенты перекрестной кинетической матрицы переноса теплоты относительно процессов
        AHeatDiff1 = 0.47
        AHeatDiff2 = 0.03
    
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 4.31
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 5.51
        H2_7 = 2.81
    
        #Эталонный результат
        (chemPot,potBet,Temp,Betas,Aff,AffHeat,
         kineticMatrixCPCP,kineticMatrixCPHeat,
         kineticMatrixHeatCP,kineticMatrixHeatHeat,
         vProcesses,vQTransf,vHeatProcesses,
         vQPows,heatTransferMatrix,balanceMatrix,
         vx,vUPows,invC,powH,vT) = CountSystemQ1(stateCoordinates,reducedTemp,systemParameters,
                                                 nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                 nu1_8,nu1_9,nu1_10,nu1_11,nu2_1,nu2_2,nu2_3,
                                                 nu2_4,nu2_5,nu2_6,beta2_1,beta2_2,Adiff1_1,
                                                 Adiff1_2,Adiff2_1,ADiffHeat1,ADiffHeat2,
                                                 AHeatDiff1,AHeatDiff2,invC1,H1_1,H2_7)
        
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
        
        #Выполняем расчеты
        nonEqSystem.CountSystem(stateCoordinates,
                                reducedTemp,
                                systemParameters)
        
        #Проверяем значения
        self.assertEqual(chemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(potBet.tolist(), nonEqSystem.GetPotentialsInterBet().tolist())#Проверяем потенциалы взаимодействия между энергетическими степенями свободы
        self.assertEqual(Temp.tolist(), nonEqSystem.GetTEnergyPowers().tolist())#Проверяем температуры энергетических степеней свободы
        self.assertEqual(Betas.tolist(), nonEqSystem.GetBeta().tolist())#Проверяем доли распределения теплоты между энергетическими степенями свободы
        self.assertEqual(balanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(heatTransferMatrix.tolist(), nonEqSystem.GetHeatTransferMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(kineticMatrixCPCP.tolist(), nonEqSystem.GetKineticMatrixPCPC().tolist())#Проверяем главный блок кинетической матрицы
        self.assertEqual(kineticMatrixCPHeat.tolist(), nonEqSystem.GetKineticMatrixPCHeat().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatCP.tolist(), nonEqSystem.GetKineticMatrixHeatPC().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatHeat.tolist(), nonEqSystem.GetKineticMatrixHeatHeat().tolist())#Проверяем главный блок кинетической матрицы по теплообмену
        self.assertEqual(invC.tolist(), nonEqSystem.GetInvHeatCapacityMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(powH.tolist(), nonEqSystem.GetHeatEffectMatrix().tolist())#Проверяем матрицу переноса теплоты
        deltaAff = np.max(np.abs(Aff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaAffHeat = np.max(np.abs(AffHeat - nonEqSystem.GetHeatAffinity()))
        self.assertAlmostEqual(deltaAffHeat, 0.0, 9)#Проверяем термодинамические силы переноса теплотытеплообмену
        deltaVProcesses = np.max(np.abs(vProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVQTransf = np.max(np.abs(vQTransf - nonEqSystem.GetVHeatTransfers()))
        self.assertAlmostEqual(deltaVQTransf, 0.0, 9)#Проверяем скорость переноса теплоты
        deltaVHeatProcesses = np.max(np.abs(vHeatProcesses - nonEqSystem.GetVHeatProcess()))
        self.assertAlmostEqual(deltaVHeatProcesses, 0.0, 7)#Проверяем скорость выделения некомпенсированной теплоты
        deltaVQPows = np.max(np.abs(vQPows - nonEqSystem.GetVHeatPower()))
        self.assertAlmostEqual(deltaVQPows, 0.0, 7)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVXPows = np.max(np.abs(vx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 8)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVUPows = np.max(np.abs(vUPows - nonEqSystem.GetVUPower()))
        self.assertAlmostEqual(deltaVUPows, 0.0, 7)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
        deltaVRTPows = np.max(np.abs(vT - nonEqSystem.GetVReducedTemperaturesEnergyPowers()))
        self.assertAlmostEqual(deltaVRTPows, 0.0, 5)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
    def testNonEq4(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213, 39, 45, 93, 81, 80, 64])#Координаты состояния
        reducedTemp = np.array([350, 423])#Приведенные температуры энергетических степеней свободы
        systemParameters = [ 0.3,  0.7,   10,   20,   50,  70,   100,  200,   30,   45,
                              10,   25,   35,   55,   75,  81,   210,  300,  150,  153,
                             123, 15.3, 45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                             1.5, 33.9, 27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,
                             159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                            0.15, 74.1, 81.0, 5.13, 31.8, 3.3, 39.63, 54.3, 27.3, 0.15,
                            0.17]#Параметры системы
        
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
    
        #Эталонный результат
        (chemPot,potBet,Temp,Betas,Aff,AffHeat,
         kineticMatrixCPCP,kineticMatrixCPHeat,
         kineticMatrixHeatCP,kineticMatrixHeatHeat,
         vProcesses,vQTransf,vHeatProcesses,
         vQPows,heatTransferMatrix,balanceMatrix,
         vx,vUPows,invC,powH,vT,heatStreams,Streams) = CountSystemQ2(stateCoordinates,reducedTemp,systemParameters,
                                                                     nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                                     nu1_8,nu1_9,nu1_10,nu1_11,nu1_12,nu1_13,nu1_14,
                                                                     nu1_15,nu2_1,nu2_2,nu2_3,nu2_4,nu2_5,nu2_6,
                                                                     nu2_7,nu2_8,nu2_9,beta2_1,beta2_2,Adiff1_1,
                                                                     Adiff1_2,Adiff2_1,ADiffHeat2,AHeatDiff2,invC1,
                                                                     H1_1,H2_7,Tokr,AHeatOkr,ExtQEnPow1,AChem1_4_4,
                                                                     mu2_9,AChem2_3_3,AChem2_1_2,AChem2_2_1,Adiff3_3,
                                                                     ADiffHeat3,AHeatDiff3,xExt1_6,muBet1_3,muBet1_5)
        
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
        
        #Выполняем расчеты
        nonEqSystem.CountSystem(stateCoordinates,
                                reducedTemp,
                                systemParameters)
        
        #Проверяем значения
        self.assertEqual(chemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(potBet.tolist(), nonEqSystem.GetPotentialsInterBet().tolist())#Проверяем потенциалы взаимодействия между энергетическими степенями свободы
        self.assertEqual(Temp.tolist(), nonEqSystem.GetTEnergyPowers().tolist())#Проверяем температуры энергетических степеней свободы
        self.assertEqual(Betas.tolist(), nonEqSystem.GetBeta().tolist())#Проверяем доли распределения теплоты между энергетическими степенями свободы
        self.assertEqual(balanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(heatTransferMatrix.tolist(), nonEqSystem.GetHeatTransferMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(kineticMatrixCPCP.tolist(), nonEqSystem.GetKineticMatrixPCPC().tolist())#Проверяем главный блок кинетической матрицы
        self.assertEqual(kineticMatrixCPHeat.tolist(), nonEqSystem.GetKineticMatrixPCHeat().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatCP.tolist(), nonEqSystem.GetKineticMatrixHeatPC().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatHeat.tolist(), nonEqSystem.GetKineticMatrixHeatHeat().tolist())#Проверяем главный блок кинетической матрицы по теплообмену
        self.assertEqual(invC.tolist(), nonEqSystem.GetInvHeatCapacityMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(powH.tolist(), nonEqSystem.GetHeatEffectMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(heatStreams.tolist(), nonEqSystem.GetHeatEnergyPowersStreams().tolist())#Проверяем внешние тепловые потоки
        self.assertEqual(Streams.tolist(), nonEqSystem.GetStateCoordinatesStreams().tolist())#Проверяем внешние тепловые потоки
        deltaAff = np.max(np.abs(Aff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaAffHeat = np.max(np.abs(AffHeat - nonEqSystem.GetHeatAffinity()))
        self.assertAlmostEqual(deltaAffHeat, 0.0, 9)#Проверяем термодинамические силы переноса теплотытеплообмену
        deltaVProcesses = np.max(np.abs(vProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVQTransf = np.max(np.abs(vQTransf - nonEqSystem.GetVHeatTransfers()))
        self.assertAlmostEqual(deltaVQTransf, 0.0, 9)#Проверяем скорость переноса теплоты
        deltaVHeatProcesses = np.max(np.abs(vHeatProcesses - nonEqSystem.GetVHeatProcess()))
        self.assertAlmostEqual(deltaVHeatProcesses, 0.0, 6)#Проверяем скорость выделения некомпенсированной теплоты
        deltaVQPows = np.max(np.abs(vQPows - nonEqSystem.GetVHeatPower()))
        self.assertAlmostEqual(deltaVQPows, 0.0, 6)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVXPows = np.max(np.abs(vx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 9)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVUPows = np.max(np.abs(vUPows - nonEqSystem.GetVUPower()))
        self.assertAlmostEqual(deltaVUPows, 0.0, 6)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
        deltaVRTPows = np.max(np.abs(vT - nonEqSystem.GetVReducedTemperaturesEnergyPowers()))
        self.assertAlmostEqual(deltaVRTPows, 0.0, 5)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
    def testNonEq5(self):
        #Исходные данные
        stateCoordinates = np.array([17, 19, 49, 75, 121, 233, 36, 42, 94, 86, 61, 68])#Координаты состояния
        reducedTemp = np.array([347, 613])#Приведенные температуры энергетических степеней свободы
        systemParameters = [ 0.2,  0.8,   10,   20,   50,  70,   100,  200,   1,   55,
                              10,   25,   35,   55,   75,  81,   210,  300,  160,  173,
                             123, 15.3, 43.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,
                             1.5, 33.9, 27.3, 4.43, 21.9, 7.5,  52.6, 21.9, 45.3, 33.9,
                             139, 21.9, 41.5, 7.53, 36.5, 1.5, 120.3, 74.5, 15.3, 21.3,
                            0.15, 74.1, 81.0, 5.13, 34.8, 3.3, 39.63, 54.3, 27.3, 0.25,
                            0.55]#Параметры системы
        
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
        
        #Доли распределения некомпенсированных теплот по энергетических степеням свободы
        beta2_1 = 0.15
        beta2_2 = 0.85
        
        #Базовые коэффициенты главной кинетической матрицы процессов]
        AChem1_4_4 = 2.3
        AChem2_3_3 = 4.3
        AChem2_1_2 = 1.3
        AChem2_2_1 = 1.15
        Adiff1_1 = 2.83
        Adiff1_2 = 1.3
        Adiff2_1 = 1.3
        Adiff3_3 = 3.83
            
        #Базовые коэффициенты перекрестной кинетической матрицы процессов относительно переноса теплоты
        ADiffHeat2 = 0.19
        ADiffHeat3 = 0.163
        
        #Базовые коэффициенты перекрестной кинетической матрицы переноса теплоты относительно процессов
        AHeatDiff2 = 0.19
        AHeatDiff3 = 0.13
    
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 2.41
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 3.61
        H2_7 = 0.71
    
        #Температура окружающей среды
        Tokr = 210
        
        #Коэффициент теплотдачи в окружающую среду
        AHeatOkr = 115
        
        #Внешние потоки теплоты
        ExtQEnPow1 = 93.0
        
        #Внешние потоки вещества
        xExt1_6 = 31.9
        
        #Потенциал взаимодействия
        mu2_9 = 14.5
        
        #Потенциал взаимодействия между энергетическими степенями свободы
        muBet1_3 = 19
        muBet1_5 = 177
    
        #Эталонный результат
        (chemPot,potBet,Temp,Betas,Aff,AffHeat,
         kineticMatrixCPCP,kineticMatrixCPHeat,
         kineticMatrixHeatCP,kineticMatrixHeatHeat,
         vProcesses,vQTransf,vHeatProcesses,
         vQPows,heatTransferMatrix,balanceMatrix,
         vx,vUPows,invC,powH,vT,heatStreams,Streams) = CountSystemQ2(stateCoordinates,reducedTemp,systemParameters,
                                                                     nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                                     nu1_8,nu1_9,nu1_10,nu1_11,nu1_12,nu1_13,nu1_14,
                                                                     nu1_15,nu2_1,nu2_2,nu2_3,nu2_4,nu2_5,nu2_6,
                                                                     nu2_7,nu2_8,nu2_9,beta2_1,beta2_2,Adiff1_1,
                                                                     Adiff1_2,Adiff2_1,ADiffHeat2,AHeatDiff2,invC1,
                                                                     H1_1,H2_7,Tokr,AHeatOkr,ExtQEnPow1,AChem1_4_4,
                                                                     mu2_9,AChem2_3_3,AChem2_1_2,AChem2_2_1,Adiff3_3,
                                                                     ADiffHeat3,AHeatDiff3,xExt1_6,muBet1_3,muBet1_5)
        
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
        
        #Выполняем расчеты
        nonEqSystem.CountSystem(stateCoordinates,
                                reducedTemp,
                                systemParameters)
        
        #Проверяем значения
        self.assertEqual(chemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(potBet.tolist(), nonEqSystem.GetPotentialsInterBet().tolist())#Проверяем потенциалы взаимодействия между энергетическими степенями свободы
        self.assertEqual(Temp.tolist(), nonEqSystem.GetTEnergyPowers().tolist())#Проверяем температуры энергетических степеней свободы
        self.assertEqual(Betas.tolist(), nonEqSystem.GetBeta().tolist())#Проверяем доли распределения теплоты между энергетическими степенями свободы
        self.assertEqual(balanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(heatTransferMatrix.tolist(), nonEqSystem.GetHeatTransferMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(kineticMatrixCPCP.tolist(), nonEqSystem.GetKineticMatrixPCPC().tolist())#Проверяем главный блок кинетической матрицы
        self.assertEqual(kineticMatrixCPHeat.tolist(), nonEqSystem.GetKineticMatrixPCHeat().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatCP.tolist(), nonEqSystem.GetKineticMatrixHeatPC().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatHeat.tolist(), nonEqSystem.GetKineticMatrixHeatHeat().tolist())#Проверяем главный блок кинетической матрицы по теплообмену
        self.assertEqual(invC.tolist(), nonEqSystem.GetInvHeatCapacityMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(powH.tolist(), nonEqSystem.GetHeatEffectMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(heatStreams.tolist(), nonEqSystem.GetHeatEnergyPowersStreams().tolist())#Проверяем внешние тепловые потоки
        self.assertEqual(Streams.tolist(), nonEqSystem.GetStateCoordinatesStreams().tolist())#Проверяем внешние тепловые потоки
        deltaAff = np.max(np.abs(Aff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaAffHeat = np.max(np.abs(AffHeat - nonEqSystem.GetHeatAffinity()))
        self.assertAlmostEqual(deltaAffHeat, 0.0, 9)#Проверяем термодинамические силы переноса теплотытеплообмену
        deltaVProcesses = np.max(np.abs(vProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVQTransf = np.max(np.abs(vQTransf - nonEqSystem.GetVHeatTransfers()))
        self.assertAlmostEqual(deltaVQTransf, 0.0, 9)#Проверяем скорость переноса теплоты
        deltaVHeatProcesses = np.max(np.abs(vHeatProcesses - nonEqSystem.GetVHeatProcess()))
        self.assertAlmostEqual(deltaVHeatProcesses, 0.0, 5)#Проверяем скорость выделения некомпенсированной теплоты
        deltaVQPows = np.max(np.abs(vQPows - nonEqSystem.GetVHeatPower()))
        self.assertAlmostEqual(deltaVQPows, 0.0, 5)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVXPows = np.max(np.abs(vx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 8)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVUPows = np.max(np.abs(vUPows - nonEqSystem.GetVUPower()))
        self.assertAlmostEqual(deltaVUPows, 0.0, 4)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
        deltaVRTPows = np.max(np.abs(vT - nonEqSystem.GetVReducedTemperaturesEnergyPowers()))
        self.assertAlmostEqual(deltaVRTPows, 0.0, 3)#Проверяем скорость изменения внутренней энергии энергетической степенисвободы
    def testNonEq6(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 56, 75, 111, 223, 39, 45, 93, 81, 80, 64])#Координаты состояния
        reducedTemp = np.array([250, 323])#Приведенные температуры энергетических степеней свободы
        systemParameters = [ 0.6,  0.4,   10,   20,   50,  70,   100,  200,   30,   45,
                              11,   65,   35,   55,   75,  81,   280,  300,  150,  153,
                             123, 15.3, 45.9, 5.53,   10, 6.5,  45.4,  315, 21.3,  4.5,
                             1.5, 33.9, 27.3, 4.53, 21.9, 7.5,  53.6, 21.9, 45.3, 33.9,
                             159, 21.9, 43.5, 7.53, 37.5, 1.5, 120.3, 73.5, 15.3, 21.3,
                            0.15, 74.1, 81.0, 5.13, 31.8, 6.9, 39.63, 54.3, 27.3, 0.15,
                            0.15]#Параметры системы
        
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
        
        #Доли распределения некомпенсированных теплот по энергетических степеням свободы
        beta2_1 = 0.65
        beta2_2 = 0.35
        
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
        Tokr = 250
        
        #Коэффициент теплотдачи в окружающую среду
        AHeatOkr = 145
        
        #Внешние потоки теплоты
        ExtQEnPow1 = 90.0
        
        #Внешние потоки вещества
        xExt1_6 = 31.9
        
        #Потенциал взаимодействия
        mu2_9 = 13.5
        
        #Потенциал взаимодействия между энергетическими степенями свободы
        muBet1_3 = 18
        muBet1_5 = 127
    
        #Эталонный результат
        (chemPot,potBet,Temp,Betas,Aff,AffHeat,
         kineticMatrixCPCP,kineticMatrixCPHeat,
         kineticMatrixHeatCP,kineticMatrixHeatHeat,
         vProcesses,vQTransf,vHeatProcesses,
         vQPows,heatTransferMatrix,balanceMatrix,
         vx,vUPows,invC,powH,vT,heatStreams,Streams) = CountSystemQ2(stateCoordinates,reducedTemp,systemParameters,
                                                                     nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                                     nu1_8,nu1_9,nu1_10,nu1_11,nu1_12,nu1_13,nu1_14,
                                                                     nu1_15,nu2_1,nu2_2,nu2_3,nu2_4,nu2_5,nu2_6,
                                                                     nu2_7,nu2_8,nu2_9,beta2_1,beta2_2,Adiff1_1,
                                                                     Adiff1_2,Adiff2_1,ADiffHeat2,AHeatDiff2,invC1,
                                                                     H1_1,H2_7,Tokr,AHeatOkr,ExtQEnPow1,AChem1_4_4,
                                                                     mu2_9,AChem2_3_3,AChem2_1_2,AChem2_2_1,Adiff3_3,
                                                                     ADiffHeat3,AHeatDiff3,xExt1_6,muBet1_3,muBet1_5)
        
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
        
        #Выполняем расчеты
        nonEqSystem.CountSystem(stateCoordinates,
                                reducedTemp,
                                systemParameters)
        
        #Проверяем значения
        self.assertEqual(chemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(potBet.tolist(), nonEqSystem.GetPotentialsInterBet().tolist())#Проверяем потенциалы взаимодействия между энергетическими степенями свободы
        self.assertEqual(Temp.tolist(), nonEqSystem.GetTEnergyPowers().tolist())#Проверяем температуры энергетических степеней свободы
        self.assertEqual(Betas.tolist(), nonEqSystem.GetBeta().tolist())#Проверяем доли распределения теплоты между энергетическими степенями свободы
        self.assertEqual(balanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(heatTransferMatrix.tolist(), nonEqSystem.GetHeatTransferMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(kineticMatrixCPCP.tolist(), nonEqSystem.GetKineticMatrixPCPC().tolist())#Проверяем главный блок кинетической матрицы
        self.assertEqual(kineticMatrixCPHeat.tolist(), nonEqSystem.GetKineticMatrixPCHeat().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatCP.tolist(), nonEqSystem.GetKineticMatrixHeatPC().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatHeat.tolist(), nonEqSystem.GetKineticMatrixHeatHeat().tolist())#Проверяем главный блок кинетической матрицы по теплообмену
        self.assertEqual(invC.tolist(), nonEqSystem.GetInvHeatCapacityMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(powH.tolist(), nonEqSystem.GetHeatEffectMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(heatStreams.tolist(), nonEqSystem.GetHeatEnergyPowersStreams().tolist())#Проверяем внешние тепловые потоки
        self.assertEqual(Streams.tolist(), nonEqSystem.GetStateCoordinatesStreams().tolist())#Проверяем внешние тепловые потоки
        deltaAff = np.max(np.abs(Aff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        deltaAffHeat = np.max(np.abs(AffHeat - nonEqSystem.GetHeatAffinity()))
        self.assertAlmostEqual(deltaAffHeat, 0.0, 9)#Проверяем термодинамические силы переноса теплотытеплообмену
        deltaVProcesses = np.max(np.abs(vProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        deltaVQTransf = np.max(np.abs(vQTransf - nonEqSystem.GetVHeatTransfers()))
        self.assertAlmostEqual(deltaVQTransf, 0.0, 9)#Проверяем скорость переноса теплоты
        deltaVHeatProcesses = np.max(np.abs(vHeatProcesses - nonEqSystem.GetVHeatProcess()))
        self.assertAlmostEqual(deltaVHeatProcesses, 0.0, 6)#Проверяем скорость выделения некомпенсированной теплоты
        deltaVQPows = np.max(np.abs(vQPows - nonEqSystem.GetVHeatPower()))
        self.assertAlmostEqual(deltaVQPows, 0.0, 6)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVXPows = np.max(np.abs(vx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 9)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVUPows = np.max(np.abs(vUPows - nonEqSystem.GetVUPower()))
        self.assertAlmostEqual(deltaVUPows, 0.0, 6)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
        deltaVRTPows = np.max(np.abs(vT - nonEqSystem.GetVReducedTemperaturesEnergyPowers()))
        self.assertAlmostEqual(deltaVRTPows, 0.0, 5)#Проверяем скорость изменения внутренней энергии энергетической степени
    def testNonEq7(self):
        #Исходные данные
        stateCoordinates = np.array([10, 20, 55, 75, 111, 213])#Координаты состояния
        reducedTemp = np.array([350])#Приведенные температуры энергетических степеней свободы
        systemParameters = [  10,   20,   50,  70,   100,  200,   30,   45,   10,   25,
                              35,   55,   75,  81,   210,  300,  150,  153,  123, 15.3, 
                            45.9, 5.63,   10, 6.5,  45.4,  315, 21.3,  4.5,  1.5]#Параметры системы
        
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
    
        #Эталонный результат
        (chemPot,potBet,Temp,Betas,Aff,AffHeat,
         kineticMatrixCPCP,kineticMatrixCPHeat,
         kineticMatrixHeatCP,kineticMatrixHeatHeat,
         vProcesses,vQTransf,vHeatProcesses,
         vQPows,heatTransferMatrix,balanceMatrix,
         vx,vUPows,invC,powH,vT) = CountSystemQ3(stateCoordinates,reducedTemp,systemParameters,
                                                 nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                 nu1_8,nu1_9,nu1_10,nu1_11,invC1,H1_1)
        
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
        
        #Выполняем расчеты
        nonEqSystem.CountSystem(stateCoordinates,
                                reducedTemp,
                                systemParameters)
        
        #Проверяем значения
        self.assertEqual(chemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(potBet.tolist(), nonEqSystem.GetPotentialsInterBet().tolist())#Проверяем потенциалы взаимодействия между энергетическими степенями свободы
        self.assertEqual(Temp.tolist(), nonEqSystem.GetTEnergyPowers().tolist())#Проверяем температуры энергетических степеней свободы
        self.assertEqual(Betas.tolist(), nonEqSystem.GetBeta().tolist())#Проверяем доли распределения теплоты между энергетическими степенями свободы
        self.assertEqual(balanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(heatTransferMatrix.tolist(), nonEqSystem.GetHeatTransferMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(kineticMatrixCPCP.tolist(), nonEqSystem.GetKineticMatrixPCPC().tolist())#Проверяем главный блок кинетической матрицы
        self.assertEqual(kineticMatrixCPHeat.tolist(), nonEqSystem.GetKineticMatrixPCHeat().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatCP.tolist(), nonEqSystem.GetKineticMatrixHeatPC().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatHeat.tolist(), nonEqSystem.GetKineticMatrixHeatHeat().tolist())#Проверяем главный блок кинетической матрицы по теплообмену
        self.assertEqual(invC.tolist(), nonEqSystem.GetInvHeatCapacityMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(powH.tolist(), nonEqSystem.GetHeatEffectMatrix().tolist())#Проверяем матрицу переноса теплоты
        deltaAff = np.max(np.abs(Aff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        self.assertEqual(AffHeat.tolist(), nonEqSystem.GetHeatAffinity().tolist())#Проверяем термодинамические силы переноса теплотытеплообмену
        deltaVProcesses = np.max(np.abs(vProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        self.assertEqual(vQTransf.tolist(), nonEqSystem.GetVHeatTransfers().tolist())#Проверяем скорость переноса теплоты
        deltaVHeatProcesses = np.max(np.abs(vHeatProcesses - nonEqSystem.GetVHeatProcess()))
        self.assertAlmostEqual(deltaVHeatProcesses, 0.0, 7)#Проверяем скорость выделения некомпенсированной теплоты
        deltaVQPows = np.max(np.abs(vQPows - nonEqSystem.GetVHeatPower()))
        self.assertAlmostEqual(deltaVQPows, 0.0, 9)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVXPows = np.max(np.abs(vx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 9)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVUPows = np.max(np.abs(vUPows - nonEqSystem.GetVUPower()))
        self.assertAlmostEqual(deltaVUPows, 0.0, 5)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
        deltaVRTPows = np.max(np.abs(vT - nonEqSystem.GetVReducedTemperaturesEnergyPowers()))
        self.assertAlmostEqual(deltaVRTPows, 0.0, 5)#Проверяем скорость изменения внутренней энергии энергетической степени свободы         
    def testNonEq8(self):
        #Исходные данные
        stateCoordinates = np.array([13, 26, 52, 73, 131, 223])#Координаты состояния
        reducedTemp = np.array([350])#Приведенные температуры энергетических степеней свободы
        systemParameters = [  10,   20,   50,  70,   100,  220,   31,   45,   10,   25,
                              35,   56,   75,  83,   211,  300,  152,  156,  143, 15.4, 
                            45.9, 5.63,   11, 6.5,  44.4,  315, 20.3,  4.5,  1.5]#Параметры системы
        
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
        
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 2.61
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 3.81
    
        #Эталонный результат
        (chemPot,potBet,Temp,Betas,Aff,AffHeat,
         kineticMatrixCPCP,kineticMatrixCPHeat,
         kineticMatrixHeatCP,kineticMatrixHeatHeat,
         vProcesses,vQTransf,vHeatProcesses,
         vQPows,heatTransferMatrix,balanceMatrix,
         vx,vUPows,invC,powH,vT) = CountSystemQ3(stateCoordinates,reducedTemp,systemParameters,
                                                 nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                 nu1_8,nu1_9,nu1_10,nu1_11,invC1,H1_1)
        
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
        
        #Выполняем расчеты
        nonEqSystem.CountSystem(stateCoordinates,
                                reducedTemp,
                                systemParameters)
        
        #Проверяем значения
        self.assertEqual(chemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(potBet.tolist(), nonEqSystem.GetPotentialsInterBet().tolist())#Проверяем потенциалы взаимодействия между энергетическими степенями свободы
        self.assertEqual(Temp.tolist(), nonEqSystem.GetTEnergyPowers().tolist())#Проверяем температуры энергетических степеней свободы
        self.assertEqual(Betas.tolist(), nonEqSystem.GetBeta().tolist())#Проверяем доли распределения теплоты между энергетическими степенями свободы
        self.assertEqual(balanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(heatTransferMatrix.tolist(), nonEqSystem.GetHeatTransferMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(kineticMatrixCPCP.tolist(), nonEqSystem.GetKineticMatrixPCPC().tolist())#Проверяем главный блок кинетической матрицы
        self.assertEqual(kineticMatrixCPHeat.tolist(), nonEqSystem.GetKineticMatrixPCHeat().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatCP.tolist(), nonEqSystem.GetKineticMatrixHeatPC().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatHeat.tolist(), nonEqSystem.GetKineticMatrixHeatHeat().tolist())#Проверяем главный блок кинетической матрицы по теплообмену
        self.assertEqual(invC.tolist(), nonEqSystem.GetInvHeatCapacityMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(powH.tolist(), nonEqSystem.GetHeatEffectMatrix().tolist())#Проверяем матрицу переноса теплоты
        deltaAff = np.max(np.abs(Aff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        self.assertEqual(AffHeat.tolist(), nonEqSystem.GetHeatAffinity().tolist())#Проверяем термодинамические силы переноса теплотытеплообмену
        deltaVProcesses = np.max(np.abs(vProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        self.assertEqual(vQTransf.tolist(), nonEqSystem.GetVHeatTransfers().tolist())#Проверяем скорость переноса теплоты
        deltaVHeatProcesses = np.max(np.abs(vHeatProcesses - nonEqSystem.GetVHeatProcess()))
        self.assertAlmostEqual(deltaVHeatProcesses, 0.0, 6)#Проверяем скорость выделения некомпенсированной теплоты
        deltaVQPows = np.max(np.abs(vQPows - nonEqSystem.GetVHeatPower()))
        self.assertAlmostEqual(deltaVQPows, 0.0, 5)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVXPows = np.max(np.abs(vx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 8)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVUPows = np.max(np.abs(vUPows - nonEqSystem.GetVUPower()))
        self.assertAlmostEqual(deltaVUPows, 0.0, 6)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
        deltaVRTPows = np.max(np.abs(vT - nonEqSystem.GetVReducedTemperaturesEnergyPowers()))
        self.assertAlmostEqual(deltaVRTPows, 0.0, 5)#Проверяем скорость изменения внутренней энергии энергетической степени свободы         
    def testNonEq9(self):
        #Исходные данные
        stateCoordinates = np.array([17, 22, 56, 76, 136, 228])#Координаты состояния
        reducedTemp = np.array([350])#Приведенные температуры энергетических степеней свободы
        systemParameters = [  10,   22,   50,  70,   100,  220,   31,   45,   10,   25,
                              35,   58,   75,  84,   211,  300,  132,  156,  146, 15.4, 
                            45.9, 5.63,   11, 6.5,  48.4,  314, 20.3,  4.1,  1.5]#Параметры системы
        
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
        
        #Обратные теплоемкости энергетических степеней свободы
        invC1 = 1.61
        
        #Приведенные тепловые эффекты энергетических степеней свободы
        H1_1 = 0.81
    
        #Эталонный результат
        (chemPot,potBet,Temp,Betas,Aff,AffHeat,
         kineticMatrixCPCP,kineticMatrixCPHeat,
         kineticMatrixHeatCP,kineticMatrixHeatHeat,
         vProcesses,vQTransf,vHeatProcesses,
         vQPows,heatTransferMatrix,balanceMatrix,
         vx,vUPows,invC,powH,vT) = CountSystemQ3(stateCoordinates,reducedTemp,systemParameters,
                                                 nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                                                 nu1_8,nu1_9,nu1_10,nu1_11,invC1,H1_1)
        
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
        
        #Выполняем расчеты
        nonEqSystem.CountSystem(stateCoordinates,
                                reducedTemp,
                                systemParameters)
        
        #Проверяем значения
        self.assertEqual(chemPot.tolist(), nonEqSystem.GetPotentialsInter().tolist())#Проверяем потенциалы взаимодействия энергетических степеней свободы
        self.assertEqual(potBet.tolist(), nonEqSystem.GetPotentialsInterBet().tolist())#Проверяем потенциалы взаимодействия между энергетическими степенями свободы
        self.assertEqual(Temp.tolist(), nonEqSystem.GetTEnergyPowers().tolist())#Проверяем температуры энергетических степеней свободы
        self.assertEqual(Betas.tolist(), nonEqSystem.GetBeta().tolist())#Проверяем доли распределения теплоты между энергетическими степенями свободы
        self.assertEqual(balanceMatrix.tolist(), nonEqSystem.GetBalanceMatrix().tolist())#Проверяем матрицу баланса
        self.assertEqual(heatTransferMatrix.tolist(), nonEqSystem.GetHeatTransferMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(kineticMatrixCPCP.tolist(), nonEqSystem.GetKineticMatrixPCPC().tolist())#Проверяем главный блок кинетической матрицы
        self.assertEqual(kineticMatrixCPHeat.tolist(), nonEqSystem.GetKineticMatrixPCHeat().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatCP.tolist(), nonEqSystem.GetKineticMatrixHeatPC().tolist())#Проверяем перекрестный блок кинетической матрицы
        self.assertEqual(kineticMatrixHeatHeat.tolist(), nonEqSystem.GetKineticMatrixHeatHeat().tolist())#Проверяем главный блок кинетической матрицы по теплообмену
        self.assertEqual(invC.tolist(), nonEqSystem.GetInvHeatCapacityMatrix().tolist())#Проверяем матрицу переноса теплоты
        self.assertEqual(powH.tolist(), nonEqSystem.GetHeatEffectMatrix().tolist())#Проверяем матрицу переноса теплоты
        deltaAff = np.max(np.abs(Aff - nonEqSystem.GetAffinity()))
        self.assertAlmostEqual(deltaAff, 0.0, 9)#Проверяем термодинамические силы
        self.assertEqual(AffHeat.tolist(), nonEqSystem.GetHeatAffinity().tolist())#Проверяем термодинамические силы переноса теплотытеплообмену
        deltaVProcesses = np.max(np.abs(vProcesses - nonEqSystem.GetVProcessCoordinates()))
        self.assertAlmostEqual(deltaVProcesses, 0.0, 9)#Проверяем скорость процессов
        self.assertEqual(vQTransf.tolist(), nonEqSystem.GetVHeatTransfers().tolist())#Проверяем скорость переноса теплоты
        deltaVHeatProcesses = np.max(np.abs(vHeatProcesses - nonEqSystem.GetVHeatProcess()))
        self.assertAlmostEqual(deltaVHeatProcesses, 0.0, 5)#Проверяем скорость выделения некомпенсированной теплоты
        deltaVQPows = np.max(np.abs(vQPows - nonEqSystem.GetVHeatPower()))
        self.assertAlmostEqual(deltaVQPows, 0.0, 5)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVXPows = np.max(np.abs(vx - nonEqSystem.GetVStateCoordinates()))
        self.assertAlmostEqual(deltaVXPows, 0.0, 8)#Проверяем скорость сообщения теплоты энергетической степени свободы
        deltaVUPows = np.max(np.abs(vUPows - nonEqSystem.GetVUPower()))
        self.assertAlmostEqual(deltaVUPows, 0.0, 5)#Проверяем скорость изменения внутренней энергии энергетической степени свободы
        deltaVRTPows = np.max(np.abs(vT - nonEqSystem.GetVReducedTemperaturesEnergyPowers()))
        self.assertAlmostEqual(deltaVRTPows, 0.0, 5)#Проверяем скорость изменения внутренней энергии энергетической степени свободы         
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
