import numpy as np

from MathProtEnergyProc.HeatPowerValues.Base import IntPotentialIndepenCompFunT

import unittest

#Модульные тесты
class TestIntPotentialIndepenCompFunT(unittest.TestCase):
    def setUp(self):
        #Выполнить настройку тестов (если необходимо)
        pass
    
    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass
    
    #Модульные тесты
    def testIntPotentialIndepenCompFunT1(self):
        #Исходные данные
        intPotentialCubMatrix = np.array([[[1.11, -0.9],
                                           [-2.7, 0.03]], 
                                          [[ 0.81, 0.27],
                                           [-0.75, 0.33]]], dtype=np.double)
        tempratures = np.array([213.69, 453.39], dtype=np.double)#Темперетауры энергетических степеней свободы
        intPotentialKnowComp = np.array([4.5, -3.69, 0.15], dtype=np.double)#Известные составляющие потенциалов взаимодействия
        intPotentialCubMatrixForKnowComp = np.array([[[  0.3, 0.45],
                                                      [-0.21, -0.9]],
                                                     [[ 0.39,  0.483],
                                                      [-0.81, -0.951]],
                                                     [[0.513, 0.45],
                                                      [0.219, -0.9]]], dtype=np.double)#Базовая составляющая потенциалов взаимодействия
        jacSStateCoordinates = np.array([ 0.882 +  0.3*4.5 - 0.39*3.69 + 0.513*0.15, 0.675 + 0.45*4.5 - 0.483*3.69 + 0.45*0.15], dtype=np.double) / 213.69 + \
                               np.array([-0.315 - 0.21*4.5 + 0.81*3.69 + 0.219*0.15, 0.486 -  0.9*4.5 + 0.951*3.69 -  0.9*0.15], dtype=np.double) / 453.39
        
        #Эталонный результат
        allIntPotentialKnowCompEt = np.array([[  0.3*4.5 - 0.39*3.69 + 0.513*0.15, 0.45*4.5 - 0.483*3.69 + 0.45*0.15],
                                              [-0.21*4.5 + 0.81*3.69 + 0.219*0.15, -0.9*4.5 + 0.951*3.69 -  0.9*0.15]], dtype=np.double)#Базовая составляющая потенциала взаимодействия
        intPotentialCompEt = np.array([-0.3, 1.5], dtype=np.double)#Компоненты потенциалов взаимодействия
        intPotentialEt = np.array([[ 0.882 +  0.3*4.5 - 0.39*3.69 + 0.513*0.15, 0.675 + 0.45*4.5 - 0.483*3.69 + 0.45*0.15],
                                   [-0.315 - 0.21*4.5 + 0.81*3.69 + 0.219*0.15, 0.486 -  0.9*4.5 + 0.951*3.69 -  0.9*0.15]], dtype=np.double)#Эталонные потенциалы взаимодействия
        
        #Вызываем функцию
        ((intPotential,
          intPotentialComp),
          allIntPotentialKnowComp) = IntPotentialIndepenCompFunT(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                                                 tempratures,#Темперетауры энергетических степеней свободы
                                                                 intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                                                                 intPotentialKnowComp,#Независимые составляющие потенциалов взаимодействия
                                                                 intPotentialCubMatrixForKnowComp#Кубическая матрица по независимым составляющим потенциалов взаимодействия
                                                                 )
        
        #Проверяем значения
        err = np.max(np.abs(intPotentialComp - intPotentialCompEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(allIntPotentialKnowComp - allIntPotentialKnowCompEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialIndepenCompFunT2(self):
        #Исходные данные
        intPotentialCubMatrix = np.array([[[1.17, 0.69,  5.1],
                                           [-9.3, 0.75,  8.1],
                                           [9.33, 3.51, 11.7],
                                           [-6.3, 6.45, -8.7]], 
                                          [[-4.17,  3.69,   3.3],
                                           [-3.63,  2.25,   5.7],
                                           [ 9.63, -1.53,  15.3],
                                           [-4.53, -3.51, -11.7]], 
                                          [[7.17, 9.69,  8.1],
                                           [-6.3, 6.75, -2.1],
                                           [6.33, 6.51, 17.7],
                                           [-9.3, 3.45,  2.7]]], dtype=np.double)
        tempratures = np.array([243.69, 318.9, 753.39, 156.3], dtype=np.double)#Темперетауры энергетических степеней свободы
        intPotentialKnowComp = np.array([7.83, 3.09], dtype=np.double)#Известные составляющие потенциалов взаимодействия
        intPotentialCubMatrixForKnowComp = np.array([[[   0.3,   0.45, -0.27],
                                                      [ -0.21, -0.753, 0.333],
                                                      [  0.45,  0.915, 0.639],
                                                      [-0.369, -0.873, 0.573]],
                                                     [[  0.39,  0.483,   0.57],
                                                      [ -0.81,  0.927,  0.369],
                                                      [  0.51, -0.687, -0.663],
                                                      [-0.873, -0.351,  0.933]]], dtype=np.double)#Базовая составляющая потенциалов взаимодействия
        jacSStateCoordinates = np.array([ 47.1249 +   0.3*7.83 +  0.39*3.09, 18.2439 +  0.45*7.83 + 0.483*3.09,  16.245 -  0.27*7.83 +  0.57*3.09], dtype=np.double) / 243.69 + \
                               np.array([ -9.8991 -  0.21*7.83 -  0.81*3.09, 14.3775 - 0.753*7.83 + 0.927*3.09, -31.635 + 0.333*7.83 + 0.369*3.09], dtype=np.double) / 318.9 + \
                               np.array([-18.9567 +  0.45*7.83 +  0.51*3.09, 32.9913 + 0.915*7.83 - 0.687*3.09,  -3.321 + 0.639*7.83 - 0.663*3.09], dtype=np.double) / 753.39 + \
                               np.array([-15.2721 - 0.369*7.83 - 0.873*3.09, 32.5863 - 0.873*7.83 - 0.351*3.09,  62.559 + 0.573*7.83 + 0.933*3.09], dtype=np.double) / 156.3
        
        #Эталонный результат
        allIntPotentialKnowCompEt = np.array([[   0.3*7.83 +  0.39*3.09,   0.45*7.83 + 0.483*3.09, -0.27*7.83 +  0.57*3.09],
                                              [ -0.21*7.83 -  0.81*3.09, -0.753*7.83 + 0.927*3.09, 0.333*7.83 + 0.369*3.09],
                                              [  0.45*7.83 +  0.51*3.09,  0.915*7.83 - 0.687*3.09, 0.639*7.83 - 0.663*3.09],
                                              [-0.369*7.83 - 0.873*3.09, -0.873*7.83 - 0.351*3.09, 0.573*7.83 + 0.933*3.09]], dtype=np.double)#Базовая составляющая потенциала взаимодействия
        intPotentialCompEt = np.array([0.45, -4.83, 3.69], dtype=np.double)#Компоненты потенциалов взаимодействия
        intPotentialEt = np.array([[ 47.1249 +   0.3*7.83 +  0.39*3.09, 18.2439 +  0.45*7.83 + 0.483*3.09,  16.245 -  0.27*7.83 +  0.57*3.09],
                                   [ -9.8991 -  0.21*7.83 -  0.81*3.09, 14.3775 - 0.753*7.83 + 0.927*3.09, -31.635 + 0.333*7.83 + 0.369*3.09],
                                   [-18.9567 +  0.45*7.83 +  0.51*3.09, 32.9913 + 0.915*7.83 - 0.687*3.09,  -3.321 + 0.639*7.83 - 0.663*3.09],
                                   [-15.2721 - 0.369*7.83 - 0.873*3.09, 32.5863 - 0.873*7.83 - 0.351*3.09,  62.559 + 0.573*7.83 + 0.933*3.09]], dtype=np.double)#Эталонные потенциалы взаимодействия
        
        #Вызываем функцию
        ((intPotential,
          intPotentialComp),
          allIntPotentialKnowComp) = IntPotentialIndepenCompFunT(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                                                 tempratures,#Темперетауры энергетических степеней свободы
                                                                 intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                                                                 intPotentialKnowComp,#Независимые составляющие потенциалов взаимодействия
                                                                 intPotentialCubMatrixForKnowComp#Кубическая матрица по независимым составляющим потенциалов взаимодействия
                                                                 )
        
        #Проверяем значения
        err = np.max(np.abs(intPotentialComp - intPotentialCompEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(allIntPotentialKnowComp - allIntPotentialKnowCompEt))
        self.assertAlmostEqual(err, 0.0, 9)
    def testIntPotentialIndepenCompFunT3(self):
        #Исходные данные
        intPotentialCubMatrix = np.array([[[1.71,  -0.9, 2.1]], 
                                          [[3.21,  3.57, 0.3]], 
                                          [[6.21, -0.27, 9.3]]], dtype=np.double)
        tempratures = np.array([165.9], dtype=np.double)#Темперетауры энергетических степеней свободы
        intPotentialKnowComp = np.array([1.83, -9.09], dtype=np.double)#Известные составляющие потенциалов взаимодействия
        intPotentialCubMatrixForKnowComp = np.array([[[0.69,   0.75, -0.33]],
                                                     [[0.39, -0.483,  0.57]]], dtype=np.double)#Базовая составляющая потенциалов взаимодействия
        jacSStateCoordinates = np.array([36.72 + 0.69*1.83 - 0.39*9.09, 14.6151 + 0.75*1.83 + 0.483*9.09, 36.216 - 0.33*1.83 - 0.57*9.09], dtype=np.double) / 165.9
        
        #Эталонный результат
        allIntPotentialKnowCompEt = np.array([[0.69*1.83 - 0.39*9.09, 0.75*1.83 + 0.483*9.09, -0.33*1.83 - 0.57*9.09]], dtype=np.double)#Базовая составляющая потенциала взаимодействия
        intPotentialCompEt = np.array([-0.63, 4.23, 3.9], dtype=np.double)#Компоненты потенциалов взаимодействия
        intPotentialEt = np.array([[36.72 + 0.69*1.83 - 0.39*9.09, 14.6151 + 0.75*1.83 + 0.483*9.09, 36.216 - 0.33*1.83 - 0.57*9.09]], dtype=np.double)#Потенциалы взаимодействия
        
        #Вызываем функцию
        ((intPotential,
          intPotentialComp),
          allIntPotentialKnowComp) = IntPotentialIndepenCompFunT(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                                                 tempratures,#Темперетауры энергетических степеней свободы
                                                                 intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                                                                 intPotentialKnowComp,#Независимые составляющие потенциалов взаимодействия
                                                                 intPotentialCubMatrixForKnowComp#Кубическая матрица по независимым составляющим потенциалов взаимодействия
                                                                 )
        
        #Проверяем значения
        err = np.max(np.abs(intPotentialComp - intPotentialCompEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(intPotential - intPotentialEt))
        self.assertAlmostEqual(err, 0.0, 9)
        err = np.max(np.abs(allIntPotentialKnowComp - allIntPotentialKnowCompEt))
        self.assertAlmostEqual(err, 0.0, 9)
         
#Запустить тестирование
if __name__ == "__main__":
    unittest.main()
