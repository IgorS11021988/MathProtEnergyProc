from MathProtEnergyProc.CorrectionModel.HelpFunctions import CreateBlockMatrix, ConcatArrayValues
from MathProtEnergyProc.NonEqProcess.NonEqSystemBase import funKineticMatrixIndexes

#Функторы преобразования кинетических матриц в вектора значений элементов
class KineticMatrix(object):
    #Инициализация класса
    def __init__(self,
    
                 varKineticNames,#Имена сопряженностей между собой координат процессов
                 varKineticAffNames,#Имена сопряженностей между собой термодинамических сил
                 
                 arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов по кинетической матрице
                 ):
        #Формируем имена коордиат состояния кинетической матрицы
        kineticMatrixProcessCoordinatesNames = ConcatArrayValues(arrKineticMatrixProcessCoordinatesNames)
        
        #Формируем индексы величин
        (self.__ProcessCoordinatesVarKineticIndexes,
         self.__ProcessCoordinatesVarKineticAffIndexes
         ) = funKineticMatrixIndexes(varKineticNames,#Имена координат процессов по коэффициентам кинетической матрицы
                                     varKineticAffNames,#Имена термодинамических сил по координатам процессов по коэффициентам кинетической матрицы
                                     kineticMatrixProcessCoordinatesNames#Имена координат процессов по кинетической матрице
                                     )
    
    #Получение коэффициентов кинетической матрицы
    def __call__(self,
                
                 kineticSubMatrixes#Кинетические матрицы
                 ):
        #Формируем кинетическую матрицу из блоков
        kineticMatrix = CreateBlockMatrix(kineticSubMatrixes)
        
        #Выводим величины коэффициентов кинетической матрицы согласно их индексам
        return kineticMatrix[self.__ProcessCoordinatesVarKineticIndexes,
                             self.__ProcessCoordinatesVarKineticAffIndexes].copy()
    