import numpy as np

from MathProtEnergyProc.HeatPowerValues import IntPotentialsOne

import unittest

#Модульные тесты
class TestIntPotentialsOne(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testIntPotentialsOne1(self):
        #Исходные данные
        stateCoordinatesNamesArr = ["x1","x2","x5","x7","x10","x20","x50"]#Имена координат состояния
        energyPowersNamesArr = ["e1","e2","e5"]#Имена энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1","x2","x5","x7","x10","x20","x50"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames     = ["e1","e1","e2","e1", "e2", "e5", "e5"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        jacSStateCoordinatesArr = np.array([2.1, -3.9, 4.5, 2.7, -8.1, 8.37, 3.81], dtype=np.double)#Якобиан энтропии по координатам состояния
        tempraturesArr = np.array([186.3, 315.3, 210.3], dtype=np.double)#Темперетауры по искомым потенциалам взаимодействия
        
        #Эталонный результат
        intPotentialEt = np.array([2.1*186.3, -3.9*186.3, 4.5*315.3, 2.7*186.3, -8.1*315.3, 8.37*210.3, 3.81*210.3], dtype=np.double)
        
        #Создаем функтор потенциала взаимодействия
        intPotentialsOneFun = IntPotentialsOne(stateCoordinatesNamesArr,#Имена координат состояния
                                               energyPowersNamesArr,#Имена энергетических степеней свободы
                                               
                                               stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                               energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                               
                                               isInvTemperatute=False#По температурам или по обратным температурам считаем
                                               )
        
        #Вызываем функцию
        intPotential = intPotentialsOneFun(jacSStateCoordinatesArr,#Якобиан энтропии по координатам состояния
                                           tempraturesArr#Темперетауры по искомым потенциалам взаимодействия
                                           )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
        
        #Создаем функтор потенциала взаимодействия
        intPotentialsOneFun = IntPotentialsOne(stateCoordinatesNamesArr,#Имена координат состояния
                                               energyPowersNamesArr,#Имена энергетических степеней свободы
                                               
                                               stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                               energyPowersVarPotentialsInterNames#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                               )
        
        #Вызываем функцию
        intPotential = intPotentialsOneFun(jacSStateCoordinatesArr,#Якобиан энтропии по координатам состояния
                                           tempraturesArr#Темперетауры по искомым потенциалам взаимодействия
                                           )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialsOne2(self):
        #Исходные данные
        stateCoordinatesNamesArr = ["x1","x2","x5","x23","x25","x3","x7","x10","x20","x15","x50"]#Имена координат состояния
        energyPowersNamesArr = ["e1","e2","e5","e7","e10","e20"]#Имена энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1","x2","x5","x23","x10","x25","x15","x3","x7","x20","x50"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames     = ["e1","e1","e7", "e1", "e5", "e1","e10","e7","e5", "e2","e20"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        jacSStateCoordinatesArr = np.array([-5.7, 3.69, 7.23, -1.5, 8.123, 3.81, 2.37, 1.329, 7.23, 1.59, 1.65], dtype=np.double)#Якобиан энтропии по координатам состояния
        tempraturesArr = np.array([186.3, 315.3, 219.3, 510.3, 645.63, 318.9], dtype=np.double)#Темперетауры по искомым потенциалам взаимодействия
        
        #Эталонный результат
        intPotentialEt = np.array([-5.7*186.3, 3.69*186.3, 7.23*510.3, -1.5*186.3, 1.329*219.3, 8.123*186.3, 1.59*645.63, 3.81*510.3, 2.37*219.3, 7.23*315.3, 1.65*318.9], dtype=np.double)
        
        #Создаем функтор потенциала взаимодействия
        intPotentialsOneFun = IntPotentialsOne(stateCoordinatesNamesArr,#Имена координат состояния
                                               energyPowersNamesArr,#Имена энергетических степеней свободы
                                               
                                               stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                               energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                               
                                               isInvTemperatute=False#По температурам или по обратным температурам считаем
                                               )
        
        #Вызываем функцию
        intPotential = intPotentialsOneFun(jacSStateCoordinatesArr,#Якобиан энтропии по координатам состояния
                                           tempraturesArr#Темперетауры по искомым потенциалам взаимодействия
                                           )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
        
        #Создаем функтор потенциала взаимодействия
        intPotentialsOneFun = IntPotentialsOne(stateCoordinatesNamesArr,#Имена координат состояния
                                               energyPowersNamesArr,#Имена энергетических степеней свободы
                                               
                                               stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                               energyPowersVarPotentialsInterNames#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                               )
        
        #Вызываем функцию
        intPotential = intPotentialsOneFun(jacSStateCoordinatesArr,#Якобиан энтропии по координатам состояния
                                           tempraturesArr#Темперетауры по искомым потенциалам взаимодействия
                                           )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialsOne3(self):
        #Исходные данные
        stateCoordinatesNamesArr = ["x1","x2","x5"]#Имена координат состояния
        energyPowersNamesArr = ["e1","e2","e5"]#Имена энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1","x2","x5"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames     = ["e1","e2","e5"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        jacSStateCoordinatesArr = np.array([2.7, 3.3, 6.9], dtype=np.double)#Якобиан энтропии по координатам состояния
        tempraturesArr = np.array([219.3, 375.9, 159.3], dtype=np.double)#Темперетауры по искомым потенциалам взаимодействия
        
        #Эталонный результат
        intPotentialEt = np.array([2.7*219.3, 3.3*375.9, 6.9*159.3], dtype=np.double)
        
        #Создаем функтор потенциала взаимодействия
        intPotentialsOneFun = IntPotentialsOne(stateCoordinatesNamesArr,#Имена координат состояния
                                               energyPowersNamesArr,#Имена энергетических степеней свободы
                                               
                                               stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                               energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                               
                                               isInvTemperatute=False#По температурам или по обратным температурам считаем
                                               )
        
        #Вызываем функцию
        intPotential = intPotentialsOneFun(jacSStateCoordinatesArr,#Якобиан энтропии по координатам состояния
                                           tempraturesArr#Темперетауры по искомым потенциалам взаимодействия
                                           )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
        
        #Создаем функтор потенциала взаимодействия
        intPotentialsOneFun = IntPotentialsOne(stateCoordinatesNamesArr,#Имена координат состояния
                                               energyPowersNamesArr,#Имена энергетических степеней свободы
                                               
                                               stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                               energyPowersVarPotentialsInterNames#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                               )
        
        #Вызываем функцию
        intPotential = intPotentialsOneFun(jacSStateCoordinatesArr,#Якобиан энтропии по координатам состояния
                                           tempraturesArr#Темперетауры по искомым потенциалам взаимодействия
                                           )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialsOne4(self):
        #Исходные данные
        stateCoordinatesNamesArr = ["x1","x2","x5","x7","x10","x20","x50"]#Имена координат состояния
        energyPowersNamesArr = ["e1","e2","e5"]#Имена энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1","x2","x5","x7","x10","x20","x50"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames     = ["e1","e1","e2","e1", "e2", "e5", "e5"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        jacSStateCoordinatesArr = np.array([2.1, -3.9, 4.5, 2.7, -8.1, 8.37, 3.81], dtype=np.double)#Якобиан энтропии по координатам состояния
        invTempraturesArr = np.array([1 / 186.3, 1 / 315.3, 1 / 210.3], dtype=np.double)#Темперетауры по искомым потенциалам взаимодействия
        
        #Эталонный результат
        intPotentialEt = np.array([2.1*186.3, -3.9*186.3, 4.5*315.3, 2.7*186.3, -8.1*315.3, 8.37*210.3, 3.81*210.3], dtype=np.double)
        
        #Создаем функтор потенциала взаимодействия
        intPotentialsOneFun = IntPotentialsOne(stateCoordinatesNamesArr,#Имена координат состояния
                                               energyPowersNamesArr,#Имена энергетических степеней свободы
                                               
                                               stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                               energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                               
                                               isInvTemperatute=True#По температурам или по обратным температурам считаем
                                               )
        
        #Вызываем функцию
        intPotential = intPotentialsOneFun(jacSStateCoordinatesArr,#Якобиан энтропии по координатам состояния
                                           invTempraturesArr#Темперетауры по искомым потенциалам взаимодействия
                                           )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialsOne5(self):
        #Исходные данные
        stateCoordinatesNamesArr = ["x1","x2","x5","x23","x25","x3","x7","x10","x20","x15","x50"]#Имена координат состояния
        energyPowersNamesArr = ["e1","e2","e5","e7","e10","e20"]#Имена энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1","x2","x5","x23","x10","x25","x15","x3","x7","x20","x50"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames     = ["e1","e1","e7", "e1", "e5", "e1","e10","e7","e5", "e2","e20"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        jacSStateCoordinatesArr = np.array([-5.7, 3.69, 7.23, -1.5, 8.123, 3.81, 2.37, 1.329, 7.23, 1.59, 1.65], dtype=np.double)#Якобиан энтропии по координатам состояния
        invTempraturesArr = np.array([1 / 186.3, 1 / 315.3, 1 / 219.3, 1 / 510.3, 1 / 645.63, 1 / 318.9], dtype=np.double)#Темперетауры по искомым потенциалам взаимодействия
        
        #Эталонный результат
        intPotentialEt = np.array([-5.7*186.3, 3.69*186.3, 7.23*510.3, -1.5*186.3, 1.329*219.3, 8.123*186.3, 1.59*645.63, 3.81*510.3, 2.37*219.3, 7.23*315.3, 1.65*318.9], dtype=np.double)
        
        #Создаем функтор потенциала взаимодействия
        intPotentialsOneFun = IntPotentialsOne(stateCoordinatesNamesArr,#Имена координат состояния
                                               energyPowersNamesArr,#Имена энергетических степеней свободы
                                               
                                               stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                               energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                               
                                               isInvTemperatute=True#По температурам или по обратным температурам считаем
                                               )
        
        #Вызываем функцию
        intPotential = intPotentialsOneFun(jacSStateCoordinatesArr,#Якобиан энтропии по координатам состояния
                                           invTempraturesArr#Темперетауры по искомым потенциалам взаимодействия
                                           )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialsOne6(self):
        #Исходные данные
        stateCoordinatesNamesArr = ["x1","x2","x5"]#Имена координат состояния
        energyPowersNamesArr = ["e1","e2","e5"]#Имена энергетических степеней свободы
        stateCoordinatesVarPotentialsInterNames = ["x1","x2","x5"]#Имена переменных потенциалов взаимодействия по координатам состояния
        energyPowersVarPotentialsInterNames     = ["e1","e2","e5"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
        jacSStateCoordinatesArr = np.array([2.7, 3.3, 6.9], dtype=np.double)#Якобиан энтропии по координатам состояния
        invTempraturesArr = np.array([1 / 219.3, 1 / 375.9, 1 / 159.3], dtype=np.double)#Темперетауры по искомым потенциалам взаимодействия
        
        #Эталонный результат
        intPotentialEt = np.array([2.7*219.3, 3.3*375.9, 6.9*159.3], dtype=np.double)
        
        #Создаем функтор потенциала взаимодействия
        intPotentialsOneFun = IntPotentialsOne(stateCoordinatesNamesArr,#Имена координат состояния
                                               energyPowersNamesArr,#Имена энергетических степеней свободы
                                               
                                               stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                                               energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                               
                                               isInvTemperatute=True#По температурам или по обратным температурам считаем
                                               )
        
        #Вызываем функцию
        intPotential = intPotentialsOneFun(jacSStateCoordinatesArr,#Якобиан энтропии по координатам состояния
                                           invTempraturesArr#Темперетауры по искомым потенциалам взаимодействия
                                           )
        
        #Проверяем значения
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
