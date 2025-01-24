import numpy as np
from scipy.linalg import solve

from MathProtEnergyProc.CorrectionModel.HelpFunctions import CreateBlockMatrix

#Функции кинетических матриц из их обратимых и необратимых составляющих
def KineticMatrixFromPosSubMatrix(posSubMatrixes,#Положительные определенные составляющие атрицы
                                  subMatrixesBalance#Податрицы баланса
                                  ):
    #Получаем блочно-диагональную кинетическую матрицу
    blockKineticMatrix = CreateBlockMatrix(posSubMatrixes)
    blockKineticMatrix = np.array(blockKineticMatrix, dtype=np.double)#Приводим к типу вещественных чисел повышенной точности
    
    #Получаем матрицу обратимых состааляющих
    MatrixBalance = np.hstack(subMatrixesBalance, dtype=np.double)
    
    #Выводим результат
    return np.dot(np.dot(MatrixBalance, blockKineticMatrix), MatrixBalance.transpose())
def KineticMatrixFromSymAsym(irRevCom,#Необратимые компоненты
                             symRevCom,#Симметричные обратимые компоненты
                             asymRevCom#Антисимметричные обратимые компоненты
                             ):
    #Выводим результат
    return KineticMatrixFromPosSubMatrix(irRevCom, symRevCom) + asymRevCom - asymRevCom.transpose()

#Функции кинетических матриц из их состаляющих увлечения потоков и эквивалентности термодинамических сил и их необратимых составляющих
def SymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                 cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                 kineticSubMatrix#Блок кинетической матрицы
                 ):#Получение симметричных обратимых составляющих из коэффициентов увлечения и эквивалентности
    #Решаем систему линейных уравнений относительно сиетричной обратимой составляющих и выводи результат
    A = kineticSubMatrix + kineticSubMatrix.transpose()
    b = np.dot(cfFacStream, kineticSubMatrix).transpose() + np.dot(kineticSubMatrix, cfEkvAffinities)
    return solve(A, b).transpose()#Выводим симетричную составляющую путе решения системы
def SymAndAsymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                        cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                        kineticSubMatrix#Блок кинетической матрицы
                        ):
    #Получаем антисиметричную необратимую составляющую
    symRevCom = SymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                             cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                             kineticSubMatrix#Блок кинетической матрицы
                             )
    
    #Получаем и выводим необратимую составляющую
    return (symRevCom,np.dot(cfFacStream - symRevCom, kineticSubMatrix))
def KineticMatrixFromFacStreamEkvAff(irRevCom,#Необратимая составляющая кинетической матрицы
                                     cfFacStream,#Коэффициенты увлечения потоков
                                     cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                     kineticSubMatrix#Блок кинетической матрицы
                                     ):
    #Определяем симметричную и антисимметричную обратимые составляющие
    (symRevCom, asymRevCom) = SymAndAsymRevMatrix(cfFacStream,#Коэффициенты увлечения потоков
                                                  cfEkvAffinities,#Коэффициенты эквивалетность термодинаических сил
                                                  kineticSubMatrix#Блок кинетической матрицы
                                                  )
    
    #Определяем матрицу симметричных обратимых составляющих
    (nStr,nSt) = symRevCom.shape#Число строк и столбцов
    balMatrix = np.block([[        np.eye(nStr),   symRevCom],
                          [np.zeros((nSt,nStr)), np.eye(nSt)]])
    
    #Определяем матрицу антисимметричных обратимых составляющих
    asymRevCom = np.block([[np.zeros((nStr,nStr)),          asymRevCom],
                           [ np.zeros((nSt,nStr)), np.zeros((nSt,nSt))]])
    
    #Определяем и выводим кинетическую матрицу
    return KineticMatrixFromSymAsym([irRevCom,kineticSubMatrix],#Необратимые компоненты
                                    [balMatrix],#Симметричные обратимые компоненты
                                    asymRevCom#Антисимметричные обратимые компоненты
                                    )
    