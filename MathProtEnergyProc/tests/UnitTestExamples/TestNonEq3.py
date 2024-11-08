import math as m
import numpy as np

#Функция состояния системысистемы
def CountState3(stateCoordinates,#Координаты состояния
                systemParameters#Параметры системы
                ):
    #Получаем числа молей
    x1_1 = stateCoordinates[0]
    x1_2 = stateCoordinates[1]
    x1_3 = stateCoordinates[2]
    x1_4 = stateCoordinates[3]
    x1_5 = stateCoordinates[4]
    x1_6 = stateCoordinates[5]
    
    #Базовые числа молей
    xs1_1 = systemParameters[0]
    xs1_2 = systemParameters[1]
    xs1_3 = systemParameters[2]
    xs1_4 = systemParameters[3]
    xs1_5 = systemParameters[4]
    xs1_6 = systemParameters[5]
    
    #Базовые химические потенциалы веществ
    mus1_1 = systemParameters[6]
    mus1_2 = systemParameters[7]
    mus1_3 = systemParameters[7]
    mus1_4 = systemParameters[8]
    mus1_5 = systemParameters[9]
    mus1_6 = systemParameters[10]
    mu0s1_1 = systemParameters[11]
    mu0s1_2 = systemParameters[12]
    mu0s1_3 = systemParameters[13]
    mu0s1_4 = systemParameters[14]
    mu0s1_5 = systemParameters[15]
    mu0s1_6 = systemParameters[16]
    
    #Базовые коэффициенты главной кинетической матрицы процессов
    A0sChem1_1_1 = systemParameters[17]
    A0sChem1_1_2 = systemParameters[18]
    A0sChem1_2_1 = systemParameters[19]
    A0sChem1_2_2 = systemParameters[20]
    A0sChem1_3_3 = systemParameters[21]
    AsChem1_1_1 = systemParameters[22]
    AsChem1_1_2 = systemParameters[23]
    AsChem1_2_1 = systemParameters[24]
    AsChem1_2_2 = systemParameters[25]
    AsChem1_3_3 = systemParameters[26]

    #Определяем химические потенциалы веществ
    mu1_1 = mus1_1*m.log(x1_1 / xs1_1) + mu0s1_1
    mu1_2 = mus1_2*m.log(x1_2 / xs1_2) + mu0s1_2
    mu1_3 = mus1_3*m.log(x1_3 / xs1_3) + mu0s1_3
    mu1_4 = mus1_4*m.log(x1_4 / xs1_4) + mu0s1_4
    mu1_5 = mus1_5*m.log(x1_5 / xs1_5) + mu0s1_5
    mu1_6 = mus1_6*m.log(x1_6 / xs1_6) + mu0s1_6
    
    #Определяем главную кинетическую матрицу процессов
    AChem1_1_1 = A0sChem1_1_1 + AsChem1_1_1 * m.fabs(x1_1 - xs1_1)
    AChem1_1_2 = A0sChem1_1_2 + AsChem1_1_2 * m.fabs(x1_2 - xs1_2)
    AChem1_2_1 = A0sChem1_2_1 + AsChem1_2_1 * m.fabs(x1_1 - xs1_1)
    AChem1_2_2 = A0sChem1_2_2 + AsChem1_2_2 * m.fabs(x1_4 - xs1_4)
    AChem1_3_3 = A0sChem1_3_3 + AsChem1_3_3 * m.fabs(x1_3 - xs1_3)
    
    #Выводим результат
    chemPot = np.array([mu1_1,mu1_2,mu1_3,mu1_4,mu1_5,mu1_6])#Массив потенциалов взаимодействия
    kineticMatrix = np.array([AChem1_1_1, AChem1_1_2, AChem1_2_1, AChem1_2_2, AChem1_3_3])
    balanceMatrix = []
    return (balanceMatrix,[],-chemPot,kineticMatrix)

#Функция расчета характеристик системы
def CountSystem3(stateCoordinates,#Координаты состояния
                 systemParameters,#Параметры системы
                 
                 nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                 nu1_8,nu1_9,nu1_10,nu1_11):
    #Получаем числа молей
    x1_1 = stateCoordinates[0]
    x1_2 = stateCoordinates[1]
    x1_3 = stateCoordinates[2]
    x1_4 = stateCoordinates[3]
    x1_5 = stateCoordinates[4]
    x1_6 = stateCoordinates[5]
        
    #Базовые числа молей
    xs1_1 = systemParameters[0]
    xs1_2 = systemParameters[1]
    xs1_3 = systemParameters[2]
    xs1_4 = systemParameters[3]
    xs1_5 = systemParameters[4]
    xs1_6 = systemParameters[5]
    
    #Базовые химические потенциалы веществ
    mus1_1 = systemParameters[6]
    mus1_2 = systemParameters[7]
    mus1_3 = systemParameters[7]
    mus1_4 = systemParameters[8]
    mus1_5 = systemParameters[9]
    mus1_6 = systemParameters[10]
    mu0s1_1 = systemParameters[11]
    mu0s1_2 = systemParameters[12]
    mu0s1_3 = systemParameters[13]
    mu0s1_4 = systemParameters[14]
    mu0s1_5 = systemParameters[15]
    mu0s1_6 = systemParameters[16]
    
    #Базовые коэффициенты главной кинетической матрицы процессов
    A0sChem1_1_1 = systemParameters[17]
    A0sChem1_1_2 = systemParameters[18]
    A0sChem1_2_1 = systemParameters[19]
    A0sChem1_2_2 = systemParameters[20]
    A0sChem1_3_3 = systemParameters[21]
    AsChem1_1_1 = systemParameters[22]
    AsChem1_1_2 = systemParameters[23]
    AsChem1_2_1 = systemParameters[24]
    AsChem1_2_2 = systemParameters[25]
    AsChem1_3_3 = systemParameters[26]
    
    #Определяем химические потенциалы веществ
    mu1_1 = mus1_1*m.log(x1_1 / xs1_1) + mu0s1_1
    mu1_2 = mus1_2*m.log(x1_2 / xs1_2) + mu0s1_2
    mu1_3 = mus1_3*m.log(x1_3 / xs1_3) + mu0s1_3
    mu1_4 = mus1_4*m.log(x1_4 / xs1_4) + mu0s1_4
    mu1_5 = mus1_5*m.log(x1_5 / xs1_5) + mu0s1_5
    mu1_6 = mus1_6*m.log(x1_6 / xs1_6) + mu0s1_6
    
    #Определяем главную кинетическую матрицу процессов
    AChem1_1_1 = A0sChem1_1_1 + AsChem1_1_1 * m.fabs(x1_1 - xs1_1)
    AChem1_1_2 = A0sChem1_1_2 + AsChem1_1_2 * m.fabs(x1_2 - xs1_2)
    AChem1_2_1 = A0sChem1_2_1 + AsChem1_2_1 * m.fabs(x1_1 - xs1_1)
    AChem1_2_2 = A0sChem1_2_2 + AsChem1_2_2 * m.fabs(x1_4 - xs1_4)
    AChem1_3_3 = A0sChem1_3_3 + AsChem1_3_3 * m.fabs(x1_3 - xs1_3)
    
    #Определяем термодинамические силы
    Xchem1_1 = nu1_1*mu1_1 + nu1_2*mu1_2 - nu1_3*mu1_4
    Xchem1_2 = nu1_4*mu1_3 + nu1_5*mu1_2 - nu1_6*mu1_6 - nu1_7*mu1_5
    Xchem1_3 = nu1_8*mu1_1 + nu1_9*mu1_3 + nu1_10*mu1_5 - nu1_11*mu1_2
    
    #Определяем скорости процессов
    vChem1_1 = AChem1_1_1*Xchem1_1 + AChem1_1_2*Xchem1_2
    vChem1_2 = AChem1_2_1*Xchem1_1 + AChem1_2_2*Xchem1_2
    vChem1_3 = AChem1_3_3*Xchem1_3
    
    #Определяем скорости изменения координат состояния
    vx1_1 = -nu1_1*vChem1_1 - nu1_8*vChem1_3
    vx1_2 = -nu1_2*vChem1_1 - nu1_5*vChem1_2 + nu1_11*vChem1_3
    vx1_3 = -nu1_4*vChem1_2 - nu1_9*vChem1_3
    vx1_4 = nu1_3*vChem1_1
    vx1_5 = nu1_7*vChem1_2 - nu1_10*vChem1_3
    vx1_6 = nu1_6*vChem1_2
    
    #Выводим результат
    chemPot = np.array([mu1_1,mu1_2,mu1_3,mu1_4,mu1_5,mu1_6])#Массив потенциалов взаимодействия
    Aff = np.array([Xchem1_1, Xchem1_2, Xchem1_3])#Термодинамические силы
    kineticMatrix = np.array([[AChem1_1_1, AChem1_1_2,        0.0],
                              [AChem1_2_1, AChem1_2_2,        0.0],
                              [       0.0,        0.0, AChem1_3_3]])
    vProcesses = np.array([vChem1_1, vChem1_2, vChem1_3])
    balanceMatrix = np.array([[-nu1_1,    0.0,  -nu1_8],
                              [-nu1_2, -nu1_5,  nu1_11],
                              [   0.0, -nu1_4,  -nu1_9],
                              [ nu1_3,    0.0,     0.0],
                              [   0.0,  nu1_7, -nu1_10],
                              [   0.0,  nu1_6,     0.0]])
    vx = np.array([vx1_1, vx1_2, vx1_3, vx1_4, vx1_5, vx1_6])
    return (-chemPot,Aff,kineticMatrix,vProcesses,balanceMatrix,vx)
    