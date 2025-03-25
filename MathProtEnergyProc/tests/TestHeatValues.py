import numpy as np

from MathProtEnergyProc.HeatPowerValues import HeatValues

import unittest

#Модульные тесты
class TestHeatValues(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testHeatValues1(self):
        #Исходные данные
        stateCoordinatesNames = ["x1","x2","x5","x7","x10"]#Имена координат состояния
        energyPowersNames = ["e1","e2","e5"]#Имена энергетических степеней свободы
        temperaturesEnergyPowersVarInvHeatCapacityNames = ["e2","e2","e1","e2","e5","e1","e5","e1","e5"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersVarInvHeatCapacityNames             = ["e5","e2","e1","e1","e1","e2","e2","e5","e5"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        temperaturesEnergyPowersVarHeatEffectNames = ["e5", "e2","e2", "e1","e1","e1","e1","e1","e2","e2","e2","e5","e5","e5", "e5"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        stateCoordinatesVarHeatEffectNames         = ["x2","x10","x5","x10","x1","x2","x5","x7","x1","x2","x7","x1","x5","x7","x10"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        jacSTemperatures = np.array([1.53, 0.21, 3.03], dtype=np.double)#Якобиан приведенной энтропии по температурам
        hesSTemperatures = np.array([[ 8.49, -0.81, 0.027],
                                     [-0.81,  6.93, 0.393],
                                     [0.027, 0.393,   9.3]], dtype=np.double)#Матрица Гесса приведенной энтропии по температурам
        hesSTemperaturesStateCoordinates = np.array([[8.49, -0.81, 0.027,  4.5,  -3.9],
                                                     [0.21,  6.93, 0.393, 6.51,   0.3],
                                                     [0.03, 0.063,   9.3, 2.31, -1.23]], dtype=np.double)#Матрица Гесса приведенной энтропии по температурам и координатам состояния
        temperatures = np.array([213.3, 369.3, 930.9], dtype=np.double)#Температуры
        
        #Эталонные значения
        heatCapacitiesEt = np.array([[2*1.53*213.3 + 8.49*213.3*213.3,               -0.81*213.3*213.3,              0.027*213.3*213.3],
                                     [              -0.81*369.3*369.3, 2*0.21*369.3 + 6.93*369.3*369.3,              0.393*369.3*369.3],
                                     [              0.027*930.9*930.9,               0.393*930.9*930.9, 2*3.03*930.9 + 9.3*930.9*930.9]], dtype=np.double)
        heatEffectsEt = np.array([[8.49*213.3*213.3, -0.81*213.3*213.3, 0.027*213.3*213.3,  4.5*213.3*213.3,  -3.9*213.3*213.3],
                                  [0.21*369.3*369.3,  6.93*369.3*369.3, 0.393*369.3*369.3, 6.51*369.3*369.3,   0.3*369.3*369.3],
                                  [0.03*930.9*930.9, 0.063*930.9*930.9,   9.3*930.9*930.9, 2.31*930.9*930.9, -1.23*930.9*930.9]], dtype=np.double)
        
        #Создаем функтор
        heatValues = HeatValues(stateCoordinatesNames,#Имена координат состояния
                                energyPowersNames,#Имена энергетических степеней свободы
                                
                                temperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                temperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                stateCoordinatesVarHeatEffectNames#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                )
        
        #Вызываем функцию
        (invHeatCapacitiesArr,
         redHeatEffectsArr) = heatValues(jacSTemperatures,#Якобиан приведенной энтропии по температурам
                                         hesSTemperatures,#Матрица Гесса приведенной энтропии по температурам
                                         hesSTemperaturesStateCoordinates,#Матрица Гесса приведенной энтропии по температурам и координатам состояния
                                         temperatures#Температуры
                                         )
        
        #Получаем матрицу теплоемкостей
        heatCapacities = heatValues.GetHeatCapacityMatrix()
        
        #Получаем матрицу тепловых эффектов
        heatEffects = heatValues.GetHeatEffectMatrix()
        
        #Приводим результат расчета обратной теплоемкости к матрице теплоемкостей
        invHeatCapacities = np.zeros_like(hesSTemperatures)
        invHeatCapacities[[1,1,0,1,2,0,2,0,2],
                          [2,1,0,0,0,1,1,2,2]] = invHeatCapacitiesArr
        
        #Приводим результат расчета приведенных эффектов к матрице приведенных эффектов
        redHeatEffects = np.zeros_like(hesSTemperaturesStateCoordinates)
        redHeatEffects[[2,1,1,0,0,0,0,0,1,1,1,2,2,2,2],
                       [1,4,2,4,0,1,2,3,0,1,3,0,2,3,4]] = redHeatEffectsArr
        
        #Проверяем значения
        err = np.max(np.abs(heatCapacities - heatCapacitiesEt))
        self.assertAlmostEqual(err, 0.0, 8)
        err = np.max(np.abs(heatEffects - heatEffectsEt))
        self.assertAlmostEqual(err, 0.0, 8)
        err = np.max(np.abs(np.dot(heatCapacities, invHeatCapacities) - np.eye(3)))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(np.dot(heatCapacities, redHeatEffects) + heatEffects))
        self.assertAlmostEqual(err, 0.0, 9)
    def testHeatValues2(self):
        #Исходные данные
        stateCoordinatesNames = ["x1","x2","x5","x7","x10","x15"]#Имена координат состояния
        energyPowersNames = ["e1","e2","e5","e7"]#Имена энергетических степеней свободы
        temperaturesEnergyPowersVarInvHeatCapacityNames = ["e5","e7","e5","e1","e2","e5","e1","e2","e7","e1","e2","e5","e7","e1","e2","e7"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersVarInvHeatCapacityNames             = ["e7","e1","e2","e1","e1","e1","e2","e2","e2","e5","e5","e5","e5","e7","e7","e7"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        temperaturesEnergyPowersVarHeatEffectNames = [ "e7","e1","e5","e2","e1","e1","e1", "e1", "e1","e2","e2","e2", "e2", "e2","e5","e5","e5", "e5", "e5","e7","e7","e7","e7", "e7"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        stateCoordinatesVarHeatEffectNames         = ["x10","x2","x7","x2","x5","x1","x7","x10","x15","x1","x5","x7","x10","x15","x1","x2","x5","x10","x15","x1","x2","x5","x7","x15"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        jacSTemperatures = np.array([0.33, 2.1, 5.37, 4.35], dtype=np.double)#Якобиан приведенной энтропии по температурам
        hesSTemperatures = np.array([[ 8.49,  0.51,  0.03, -1.65],
                                     [ 0.51,  6.93, 0.693,  6.51],
                                     [ 0.03, 0.693,   9.3, -1.23],
                                     [-1.65,  6.51, -1.23,  5.31]], dtype=np.double)#Матрица Гесса приведенной энтропии по температурам
        hesSTemperaturesStateCoordinates = np.array([[8.49, -0.81, 0.027,  4.5,  -3.9,  1.5],
                                                     [0.21,  6.93, 0.393, 6.51,   0.3, 9.51],
                                                     [0.03, 0.063,   9.3, 2.31, -1.23, 5.61],
                                                     [1.65, 2.061,  12.3, 5.31,  4.23, 8.61]], dtype=np.double)#Матрица Гесса приведенной энтропии по температурам и координатам состояния
        temperatures = np.array([243.3, 639.3, 960.3, 153.3], dtype=np.double)#Температуры
        
        #Эталонные значения
        heatCapacitiesEt = np.array([[2*0.33*243.3 + 8.49*243.3*243.3,               0.51*243.3*243.3,               0.03*243.3*243.3,               -1.65*243.3*243.3],
                                     [               0.51*639.3*639.3, 2*2.1*639.3 + 6.93*639.3*639.3,              0.693*639.3*639.3,                6.51*639.3*639.3],
                                     [               0.03*960.3*960.3,              0.693*960.3*960.3, 2*5.37*960.3 + 9.3*960.3*960.3,               -1.23*960.3*960.3],
                                     [              -1.65*153.3*153.3,               6.51*153.3*153.3,              -1.23*153.3*153.3, 2*4.35*153.3 + 5.31*153.3*153.3]], dtype=np.double)
        heatEffectsEt = np.array([[8.49*243.3*243.3, -0.81*243.3*243.3, 0.027*243.3*243.3,  4.5*243.3*243.3,  -3.9*243.3*243.3,  1.5*243.3*243.3],
                                  [0.21*639.3*639.3,  6.93*639.3*639.3, 0.393*639.3*639.3, 6.51*639.3*639.3,   0.3*639.3*639.3, 9.51*639.3*639.3],
                                  [0.03*960.3*960.3, 0.063*960.3*960.3,   9.3*960.3*960.3, 2.31*960.3*960.3, -1.23*960.3*960.3, 5.61*960.3*960.3],
                                  [1.65*153.3*153.3, 2.061*153.3*153.3,  12.3*153.3*153.3, 5.31*153.3*153.3,  4.23*153.3*153.3, 8.61*153.3*153.3]], dtype=np.double)
        
        #Создаем функтор
        heatValues = HeatValues(stateCoordinatesNames,#Имена координат состояния
                                energyPowersNames,#Имена энергетических степеней свободы
                                
                                temperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                temperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                stateCoordinatesVarHeatEffectNames#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                )
        
        #Вызываем функцию
        (invHeatCapacitiesArr,
         redHeatEffectsArr) = heatValues(jacSTemperatures,#Якобиан приведенной энтропии по температурам
                                         hesSTemperatures,#Матрица Гесса приведенной энтропии по температурам
                                         hesSTemperaturesStateCoordinates,#Матрица Гесса приведенной энтропии по температурам и координатам состояния
                                         temperatures#Температуры
                                         )
        
        #Получаем матрицу теплоемкостей
        heatCapacities = heatValues.GetHeatCapacityMatrix()
        
        #Получаем матрицу тепловых эффектов
        heatEffects = heatValues.GetHeatEffectMatrix()
        
        #Приводим результат расчета обратной теплоемкости к матрице теплоемкостей
        invHeatCapacities = np.zeros_like(hesSTemperatures)
        invHeatCapacities[[2,3,2,0,1,2,0,1,3,0,1,2,3,0,1,3],
                          [3,0,1,0,0,0,1,1,1,2,2,2,2,3,3,3]] = invHeatCapacitiesArr
        
        #Приводим результат расчета приведенных эффектов к матрице приведенных эффектов
        redHeatEffects = np.zeros_like(hesSTemperaturesStateCoordinates)
        redHeatEffects[[3,0,2,1,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3],
                       [4,1,3,1,2,0,3,4,5,0,2,3,4,5,0,1,2,4,5,0,1,2,3,5]] = redHeatEffectsArr
        
        #Проверяем значения
        err = np.max(np.abs(heatCapacities - heatCapacitiesEt))
        self.assertAlmostEqual(err, 0.0, 8)
        err = np.max(np.abs(heatEffects - heatEffectsEt))
        self.assertAlmostEqual(err, 0.0, 8)
        err = np.max(np.abs(np.dot(heatCapacities, invHeatCapacities) - np.eye(4)))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(np.dot(heatCapacities, redHeatEffects) + heatEffects))
        self.assertAlmostEqual(err, 0.0, 8)
    def testHeatValues3(self):
        #Исходные данные
        stateCoordinatesNames = ["x1","x2","x5"]#Имена координат состояния
        energyPowersNames = ["e1","e2"]#Имена энергетических степеней свободы
        temperaturesEnergyPowersVarInvHeatCapacityNames = ["e1","e2","e1","e2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
        energyPowersVarInvHeatCapacityNames             = ["e1","e1","e2","e2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        temperaturesEnergyPowersVarHeatEffectNames = ["e1","e1","e1","e2","e2","e2"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
        stateCoordinatesVarHeatEffectNames         = ["x1","x2","x5","x1","x2","x5"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
        jacSTemperatures = np.array([0.21, 3.9], dtype=np.double)#Якобиан приведенной энтропии по температурам
        hesSTemperatures = np.array([[ 8.49, -0.57],
                                     [-0.57,  6.93]], dtype=np.double)#Матрица Гесса приведенной энтропии по температурам
        hesSTemperaturesStateCoordinates = np.array([[8.49, -0.81, 0.027],
                                                     [0.21,  6.93, 0.393]], dtype=np.double)#Матрица Гесса приведенной энтропии по температурам и координатам состояния
        temperatures = np.array([123.3, 159.3], dtype=np.double)#Температуры
        
        #Эталонные значения
        heatCapacitiesEt = np.array([[2*0.21*123.3 + 8.49*123.3*123.3,              -0.57*123.3*123.3],
                                     [              -0.57*159.3*159.3, 2*3.9*159.3 + 6.93*159.3*159.3]], dtype=np.double)
        heatEffectsEt = np.array([[8.49*123.3*123.3, -0.81*123.3*123.3, 0.027*123.3*123.3],
                                  [0.21*159.3*159.3,  6.93*159.3*159.3, 0.393*159.3*159.3]], dtype=np.double)
        
        #Создаем функтор
        heatValues = HeatValues(stateCoordinatesNames,#Имена координат состояния
                                energyPowersNames,#Имена энергетических степеней свободы
                                
                                temperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                                energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                temperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                                stateCoordinatesVarHeatEffectNames#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                                )
        
        #Вызываем функцию
        (invHeatCapacitiesArr,
         redHeatEffectsArr) = heatValues(jacSTemperatures,#Якобиан приведенной энтропии по температурам
                                         hesSTemperatures,#Матрица Гесса приведенной энтропии по температурам
                                         hesSTemperaturesStateCoordinates,#Матрица Гесса приведенной энтропии по температурам и координатам состояния
                                         temperatures#Температуры
                                         )
        
        #Получаем матрицу теплоемкостей
        heatCapacities = heatValues.GetHeatCapacityMatrix()
        
        #Получаем матрицу тепловых эффектов
        heatEffects = heatValues.GetHeatEffectMatrix()
        
        #Приводим результат расчета обратной теплоемкости к матрице теплоемкостей
        invHeatCapacities = np.zeros_like(hesSTemperatures)
        invHeatCapacities[[0,1,0,1],
                          [0,0,1,1]] = invHeatCapacitiesArr
        
        #Приводим результат расчета приведенных эффектов к матрице приведенных эффектов
        redHeatEffects = np.zeros_like(hesSTemperaturesStateCoordinates)
        redHeatEffects[[0,0,0,1,1,1],
                       [0,1,2,0,1,2]] = redHeatEffectsArr
        
        #Проверяем значения
        err = np.max(np.abs(heatCapacities - heatCapacitiesEt))
        self.assertAlmostEqual(err, 0.0, 8)
        err = np.max(np.abs(heatEffects - heatEffectsEt))
        self.assertAlmostEqual(err, 0.0, 8)
        err = np.max(np.abs(np.dot(heatCapacities, invHeatCapacities) - np.eye(2)))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(np.dot(heatCapacities, redHeatEffects) + heatEffects))
        self.assertAlmostEqual(err, 0.0, 9)
                                                                                                                        
#Запустить тестирование                                                                                                 
if __name__ == "__main__":
    unittest.main()
