from MathProtEnergyProc.CorrectionModel.HelpFunctions import CreateBlockMatrix, ConcatArrayValues
from MathProtEnergyProc.NonEqProcess.NonEqSystemQBase import funKineticMatrixQIndexes

#Функторы преобразования кинетических матриц в вектора значений элементов
class KineticMatrixQ(object):
    #Инициализация класса
    def __init__(self,
    
                 varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                 varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                 varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                 varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                 varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                 varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                 varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                 varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                 
                 arrKineticMatrixProcessCoordinatesNames#Массив массивов имен координат процессов (в том числе и перенесенных теплот) по кинетической матрице
                 ):
        #Формируем имена коордиат состояния кинетической матрицы
        kineticMatrixProcessCoordinatesNames = ConcatArrayValues(arrKineticMatrixProcessCoordinatesNames)#Имена координат процессов (в том числе и перенесенных теплот) по кинетической матрице
        
        #Формируем индексы величин
        (self.__ProcessCoordinatesVarKineticPCPCIndexes,
         self.__ProcessCoordinatesVarKineticPCPCAffIndexes,
         
         self.__ProcessCoordinatesVarKineticPCHeatIndexes,
         self.__ProcessCoordinatesVarKineticPCHeatAffIndexes,
         
         self.__ProcessCoordinatesVarKineticHeatPCIndexes,
         self.__ProcessCoordinatesVarKineticHeatPCAffIndexes,
         
         self.__ProcessCoordinatesVarKineticHeatHeatIndexes,
         self.__ProcessCoordinatesVarKineticHeatHeatAffIndexes
         ) = funKineticMatrixQIndexes(varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                                      varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                                      varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                                      varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                                      varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                                      varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                                      varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                                      varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                                      
                                      kineticMatrixProcessCoordinatesNames,#Имена координат процессов (в том числе и перенесенных теплот) по кинетической матрице
                                      kineticMatrixProcessCoordinatesNames#Имена координат процессов (в том числе и перенесенных теплот) по кинетической матрице
                                      )
    
    #Получение коэффициентов кинетической матрицы
    def __call__(self,
                
                 kineticSubMatrixes#Кинетические матрицы
                 ):
        #Формируем кинетическую матрицу из блоков
        kineticMatrix = CreateBlockMatrix(kineticSubMatrixes)
        
        #Выводим величины коэффициентов кинетической матрицы согласно их индексам
        return (kineticMatrix[self.__ProcessCoordinatesVarKineticPCPCIndexes,
                              self.__ProcessCoordinatesVarKineticPCPCAffIndexes].copy(),#Главный блок кинетической матрицы по координатам процесса (кроме переноса теплоты)
                kineticMatrix[self.__ProcessCoordinatesVarKineticPCHeatIndexes,
                              self.__ProcessCoordinatesVarKineticPCHeatAffIndexes].copy(),#Перекрестный блок кинетической матрицы по координатам процессов (кроме перенесенных теплот) и перенесенным теплотам
                kineticMatrix[self.__ProcessCoordinatesVarKineticHeatPCIndexes,
                              self.__ProcessCoordinatesVarKineticHeatPCAffIndexes].copy(),#Перекрестный блок кинетической матрицы по перенесенным теплотам и координатам процессов (кроме перенесенных теплот)
                kineticMatrix[self.__ProcessCoordinatesVarKineticHeatHeatIndexes,
                              self.__ProcessCoordinatesVarKineticHeatHeatAffIndexes].copy()#Главный блок кинетической матрицы по перенесенным теплотам
                )
    