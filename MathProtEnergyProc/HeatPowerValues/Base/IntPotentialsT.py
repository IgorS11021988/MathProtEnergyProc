import numpy as np
from scipy import linalg as lg

from .IntPotentialsBase import IntPotentialFromCompFun

#Функция матрицы Якоби энтропии по координатам состояния
def JacSStateCoordinatesFunT(tempratures,#Темперетауры энергетических степеней свободы
                             intPotential,#Матрица потенциалов взаимодействия
                             ):
    #Определяем и выводим матрицу Якоби энтропии по координатам процессов
    return np.sum(intPotential / tempratures.reshape(-1,1), axis=0)

#Функция матрицы Якоби энтропии по координатам состояния
def JacSStateCoordinatesIndepenCompFunT(tempratures,#Темперетауры энергетических степеней свободы
                                        intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                                        intPotentialComp,#Независимые составляющие потенциалов взаимодействия
                                        ):
    #составляющая потенциалов взаимодействия, обусловленная ее  известными составляющими
    intPotential = IntPotentialFromCompFun(intPotentialCubMatrix,
                                           intPotentialComp)
    
    #Определяем и выводим матрицу Якоби энтропии по координатам процессов
    return (JacSStateCoordinatesFunT(tempratures,
                                     intPotential),
            intPotential)

#Функция потенциалов взаимодействия 
def IntPotentialFunT(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                     tempratures,#Темперетауры энергетических степеней свободы
                     intPotentialCubMatrix#Кубическая матрица потенциалов взаимодействия
                     ):
    #Получаем матрицу для независимых составляющих потенциалов взаимодйствия
    intPotentialIndepenCompMatrix = np.sum(intPotentialCubMatrix / tempratures.reshape(1,-1,1),
                                           axis=1)
    
    #Решаем систему уравнений относительно независимых составляющих потенциалов взаимодействия
    intPotentialIndepenComp = lg.solve(intPotentialIndepenCompMatrix,
                                       jacSStateCoordinates,
                                       transposed=True)
    
    #Определяем и выводим потенциалы взаимодействия вместе с их независимыми составляющими
    return (IntPotentialFromCompFun(intPotentialCubMatrix,
                                    intPotentialIndepenComp),
            intPotentialIndepenComp)

#Функция потенциалов взаимодействия с учетом известных независимых составляющих
def IntPotentialBaseCompFunT(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                             tempratures,#Темперетауры энергетических степеней свободы
                             intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                             intPotentialBaseComp#Базовая составляющая потенциалов взаимодействия
                             ):
    #Оределяем составляющие потенциалов взаимодействия и матрицы Якоби энтропии по координатам процессов, обусловленные известными составляющими потенциалов взаимодействия
    jacSStateCoordinatesKnowComp = JacSStateCoordinatesFunT(tempratures,#Темперетауры энергетических степеней свободы
                                                            intPotentialBaseComp#Базовая составляющая потенциалов взаимодействия
                                                            )
    
    #Вычитаем из якобиана энтропии по координатам состояния составляющу, обусловленную известными составляющими потенциалов взаимодействия
    jacSStateCoordinatesNoKnowComp = jacSStateCoordinates - jacSStateCoordinatesKnowComp
    
    #Получаем и выводим потенциалвы взаимодействия и их независимые составляющие
    (intPotential,
     intPotentialIndepenComp) = IntPotentialFunT(jacSStateCoordinatesNoKnowComp,
                                                 tempratures,
                                                 intPotentialCubMatrix)
    
    #Выводим результат
    return (intPotential + intPotentialBaseComp,
            intPotentialIndepenComp)

#Функция потенциалов взаимодействия с учетом известных независимых составляющих
def IntPotentialIndepenCompFunT(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                tempratures,#Темперетауры энергетических степеней свободы
                                intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                                intPotentialKnowComp,#Независимые составляющие потенциалов взаимодействия
                                intPotentialCubMatrixForKnowComp#Кубическая матрица по независимым составляющим потенциалов взаимодействия
                                ):
    #Определяем базовую составляющую потенциалов взаимодействия
    allIntPotentialKnowComp = IntPotentialFromCompFun(intPotentialCubMatrixForKnowComp,
                                                      intPotentialKnowComp)
    
    #Получаем и выводим потенциалвы взаимодействия и их составляющую, обусловленную независимой составляющей
    return (IntPotentialBaseCompFunT(jacSStateCoordinates,
                                     tempratures,
                                     intPotentialCubMatrix,
                                     allIntPotentialKnowComp),
            allIntPotentialKnowComp)

#Функция потенциалов взаимодействия при условии единственности энергетических степеней свободы для соответсвующих потенциалов взаимодействия 
def IntPotentialFunTOne(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                        tempraturesIntPotentials#Темперетауры по искомым потенциалам взаимодействия
                        ):
    #Определяем и выводим потенциалы взаимодействия
    return jacSStateCoordinates*tempraturesIntPotentials
