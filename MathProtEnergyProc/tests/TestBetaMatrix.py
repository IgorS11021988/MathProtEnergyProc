import numpy as np

from MathProtEnergyProc.CorrectionModel import BetaProcess, BetaMatrix

import unittest

#Модульные тесты
class TestBetaMatrix(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testBetaMatrix1(self):
        #Исходные данные
        betasProcesses = [(["Pow1", "Pow2", "Pow5"],
                           "proc1"),
                          (["Pow3", "Pow2", "Pow7", "Pow10"],
                           "proc2"),
                          (["Pow3", "Pow10"],
                           "proc3")
                          ]
        energyPowersVarBetaNames       = ["Pow1",  "Pow2",  "Pow3",  "Pow2",  "Pow7",  "Pow10", "Pow5",  "Pow3",  "Pow10"]
        processCoordinatesVarBetaNames = ["proc1", "proc1", "proc2", "proc2", "proc2", "proc2", "proc1", "proc3", "proc3"]
        Betas1 = np.array([0.15, 0.33, 0.52], dtype=np.double)
        Betas2 = np.array([0.51, 0.33, 0.1, 0.06], dtype=np.double)
        Betas3 = np.array([0.17, 0.83], dtype=np.double)
        arrBetas = [Betas1, Betas2, Betas3]
        
        #Эталонные значения
        etBetaValues = np.array([0.15, 0.33, 0.51, 0.33, 0.1, 0.06, 0.52, 0.17, 0.83], dtype=np.double)
        
        #Вызываем функцию
        fBetaMatrix = BetaMatrix(energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                                 processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                                 
                                 betasProcesses#Доли распределения некопенсированных теплот каждого процесса по энергетически степеня свободы
                                 )
        betaValues = fBetaMatrix(arrBetas)
        
        #Проверяем значения
        err = np.max(np.abs(betaValues - etBetaValues))
        self.assertAlmostEqual(err, 0.0, 9)
    def testBetaMatrix2(self):
        #Исходные данные
        betasProcesses = [(["Pow1", "Pow2"],
                           "proc1"),
                          (["Pow3", "Pow2", "Pow7", "Pow10", "Pow5"],
                           "proc2"),
                          (["Pow3", "Pow10"],
                           "proc3"),
                          (["Pow4", "Pow11"],
                           "proc4"),
                          (["Pow3", "Pow15", "Pow1"],
                           "proc5")
                          ]
        energyPowersVarBetaNames       = ["Pow1",  "Pow3",  "Pow2",  "Pow2",  "Pow7",  "Pow10", "Pow5",  "Pow3",  "Pow10", "Pow4",  "Pow11", "Pow15", "Pow3",  "Pow1"]
        processCoordinatesVarBetaNames = ["proc1", "proc2", "proc1", "proc2", "proc2", "proc2", "proc2", "proc3", "proc3", "proc4", "proc4", "proc5", "proc5", "proc5"]
        Betas1 = np.array([0.45, 0.55], dtype=np.double)
        Betas2 = np.array([0.21, 0.19, 0.3, 0.27, 0.03], dtype=np.double)
        Betas3 = [0.17, 0.83]
        Betas4 = np.array([0.1, 0.9], dtype=np.double)
        Betas5 = np.array([0.27, 0.23, 0.5], dtype=np.double)
        arrBetas = [Betas1, Betas2, Betas3, Betas4, Betas5]
        
        #Эталонные значения
        etBetaValues = np.array([0.45, 0.21, 0.55, 0.19, 0.3, 0.27, 0.03, 0.17, 0.83, 0.1, 0.9, 0.23, 0.27, 0.5], dtype=np.double)
        
        #Вызываем функцию
        fBetaMatrix = BetaMatrix(energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                                 processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                                 
                                 betasProcesses#Доли распределения некопенсированных теплот каждого процесса по энергетически степеня свободы
                                 )
        betaValues = fBetaMatrix(arrBetas)
        
        #Проверяем значения
        err = np.max(np.abs(betaValues - etBetaValues))
        self.assertAlmostEqual(err, 0.0, 9)
    def testBetaMatrix3(self):
        #Исходные данные
        betasProcesses = [(["Pow1", "Pow2", "Pow5"],
                           "proc1")]
        energyPowersVarBetaNames       = ["Pow1",  "Pow2",  "Pow5"]
        processCoordinatesVarBetaNames = ["proc1", "proc1", "proc1"]
        Betas = [0.15, 0.33, 0.52]
        arrBetas = Betas
        
        #Эталонные значения
        etBetaValues = np.array([0.15, 0.33, 0.52], dtype=np.double)
        
        #Вызываем функцию
        fBetaMatrix = BetaMatrix(energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                                 processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                                 
                                 betasProcesses#Доли распределения некопенсированных теплот каждого процесса по энергетически степеня свободы
                                 )
        betaValues = fBetaMatrix(arrBetas)
        
        #Проверяем значения
        err = np.max(np.abs(betaValues - etBetaValues))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
