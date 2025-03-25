import numpy as np
from scipy import linalg as lg

from MathProtEnergyProc.HeatPowerValues.Base.IntPotentialsBase import IntPotentialFromCompFun

#Функция матрицы Якоби энтропии по координатам состояния
def JacSStateCoordinatesFunSdU(jacSUPowerEnergies,#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы
                               intPotential,#Матрица потенциалов взаимодействия
                               ):
    #Определяем и выводим матрицу Якоби энтропии по координатам процессов
    return np.sum(intPotential*jacSUPowerEnergies.reshape(-1,1), axis=0)

#Функция матрицы Якоби энтропии по координатам состояния
def JacSStateCoordinatesIndepenCompFunSdU(jacSUPowerEnergies,#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы
                                          intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                                          intPotentialComp,#Независимые составляющие потенциалов взаимодействия
                                          ):
    #составляющая потенциалов взаимодействия, обусловленная ее  известными составляющими
    intPotential = IntPotentialFromCompFun(intPotentialCubMatrix,
                                           intPotentialComp)
    
    #Определяем и выводим матрицу Якоби энтропии по координатам процессов
    return (JacSStateCoordinatesFunSdU(jacSUPowerEnergies,
                                       intPotential),
            intPotential)

#Функция потенциалов взаимодействия 
def IntPotentialFunSdU(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                       jacSUPowerEnergies,#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы
                       intPotentialCubMatrix#Кубическая матрица потенциалов взаимодействия
                       ):
    #Получаем матрицу для независимых составляющих потенциалов взаимодйствия
    intPotentialIndepenCompMatrix = np.sum(intPotentialCubMatrix*jacSUPowerEnergies.reshape(1,-1,1),
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
def IntPotentialBaseCompFunSdU(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                               jacSUPowerEnergies,#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы
                               intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                               intPotentialBaseComp#Базовая составляющая потенциалов взаимодействия
                               ):
    #Оределяем составляющие потенциалов взаимодействия и матрицы Якоби энтропии по координатам процессов, обусловленные известными составляющими потенциалов взаимодействия
    jacSStateCoordinatesKnowComp = JacSStateCoordinatesFunSdU(jacSUPowerEnergies,#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы
                                                              intPotentialBaseComp#Базовая составляющая потенциалов взаимодействия
                                                              )
    
    #Вычитаем из якобиана энтропии по координатам состояния составляющу, обусловленную известными составляющими потенциалов взаимодействия
    jacSStateCoordinatesNoKnowComp = jacSStateCoordinates - jacSStateCoordinatesKnowComp
    
    #Получаем и выводим потенциалвы взаимодействия и их составляющую, обусловленную независимой составляющей
    (intPotential,
     intPotentialIndepenComp) = IntPotentialFunSdU(jacSStateCoordinatesNoKnowComp,
                                                   jacSUPowerEnergies,
                                                   intPotentialCubMatrix)
    
    #Выводим результат
    return (intPotential + intPotentialBaseComp,
            intPotentialIndepenComp)

#Функция потенциалов взаимодействия с учетом известных независимых составляющих
def IntPotentialIndepenCompFunSdU(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                                  jacSUPowerEnergies,#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы
                                  intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                                  intPotentialKnowComp,#Независимые составляющие потенциалов взаимодействия
                                  intPotentialCubMatrixForKnowComp#Кубическая матрица по независимым составляющим потенциалов взаимодействия
                                  ):
    #Определяем базовую составляющую потенциалов взаимодействия
    allIntPotentialKnowComp = IntPotentialFromCompFun(intPotentialCubMatrixForKnowComp,
                                                      intPotentialKnowComp)
    
    #Получаем и выводим потенциалвы взаимодействия и их составляющую, обусловленную независимой составляющей
    return (IntPotentialBaseCompFunSdU(jacSStateCoordinates,
                                       jacSUPowerEnergies,
                                       intPotentialCubMatrix,
                                       allIntPotentialKnowComp),
            allIntPotentialKnowComp)

#Функция потенциалов взаимодействия при условии единственности энергетических степеней свободы для соответсвующих потенциалов взаимодействия 
def IntPotentialFunSdUOne(jacSStateCoordinates,#Якобиан энтропии по координатам состояния
                          jacSUPowerEnergiesIntPotentials#Матрица Якоби энтропии по внутренним энергиям энергетических степеней свободы, соответствующая искомым потенциалам взаимодействия
                          ):
    #Определяем и выводим потенциалы взаимодействия
    return jacSStateCoordinates / jacSUPowerEnergiesIntPotentials
