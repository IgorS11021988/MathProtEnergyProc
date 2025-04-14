import math as m
import numpy as np


# Функция состояния системысистемы
def CountState1(stateCoordinates,  # Координаты состояния
                systemParameters  # Параметры системы
                ):
    # Получаем числа молей
    x1_1 = stateCoordinates[0]
    x1_2 = stateCoordinates[1]
    x1_3 = stateCoordinates[2]
    x1_4 = stateCoordinates[3]
    x1_5 = stateCoordinates[4]
    x1_6 = stateCoordinates[5]
    x2_1 = stateCoordinates[6]
    x2_4 = stateCoordinates[7]
    x2_7 = stateCoordinates[8]
    x2_8 = stateCoordinates[9]

    # Базовые числа молей
    xs1_1 = systemParameters[0]
    xs1_2 = systemParameters[1]
    xs1_3 = systemParameters[2]
    xs1_4 = systemParameters[3]
    xs1_5 = systemParameters[4]
    xs1_6 = systemParameters[5]
    xs2_1 = systemParameters[6]
    xs2_4 = systemParameters[7]
    xs2_7 = systemParameters[8]
    xs2_8 = systemParameters[9]

    # Базовые химические потенциалы веществ
    mus1_1 = systemParameters[10]
    mus1_2 = systemParameters[11]
    mus1_3 = systemParameters[12]
    mus1_4 = systemParameters[13]
    mus1_5 = systemParameters[14]
    mus1_6 = systemParameters[15]
    mus2_1 = systemParameters[16]
    mus2_4 = systemParameters[17]
    mus2_7 = systemParameters[18]
    mus2_8 = systemParameters[19]
    mu0s1_1 = systemParameters[20]
    mu0s1_2 = systemParameters[21]
    mu0s1_3 = systemParameters[22]
    mu0s1_4 = systemParameters[23]
    mu0s1_5 = systemParameters[24]
    mu0s1_6 = systemParameters[25]
    mu0s2_1 = systemParameters[26]
    mu0s2_4 = systemParameters[27]
    mu0s2_7 = systemParameters[28]
    mu0s2_8 = systemParameters[29]

    # Базовые коэффициенты кинетической матрицы процессов
    A0sChem1_1_1 = systemParameters[30]
    A0sChem1_1_2 = systemParameters[31]
    A0sChem1_2_1 = systemParameters[32]
    A0sChem1_2_2 = systemParameters[33]
    A0sChem1_3_3 = systemParameters[34]
    A0sChem2_1_1 = systemParameters[35]
    A0sChem2_2_2 = systemParameters[36]
    AsChem1_1_1 = systemParameters[37]
    AsChem1_1_2 = systemParameters[38]
    AsChem1_2_1 = systemParameters[39]
    AsChem1_2_2 = systemParameters[40]
    AsChem1_3_3 = systemParameters[41]
    AsChem2_1_1 = systemParameters[42]
    AsChem2_2_2 = systemParameters[43]
    Adiff2_2 = 4.53 + systemParameters[44]

    # Определяем химические потенциалы веществ
    mu1_1 = mus1_1 * m.log(x1_1 / xs1_1) + mu0s1_1
    mu1_2 = mus1_2 * m.log(x1_2 / xs1_2) + mu0s1_2
    mu1_3 = mus1_3 * m.log(x1_3 / xs1_3) + mu0s1_3
    mu1_4 = mus1_4 * m.log(x1_4 / xs1_4) + mu0s1_4
    mu1_5 = mus1_5 * m.log(x1_5 / xs1_5) + mu0s1_5
    mu1_6 = mus1_6 * m.log(x1_6 / xs1_6) + mu0s1_6
    mu2_1 = mus2_1 * m.log(x2_1 / xs2_1) + mu0s2_1
    mu2_4 = mus2_4 * m.log(x2_4 / xs2_4) + mu0s2_4
    mu2_7 = mus2_7 * m.log(x2_7 / xs2_7) + mu0s2_7
    mu2_8 = mus2_8 * m.log(x2_8 / xs2_8) + mu0s2_8

    # Определяем кинетическую матрицу процессов
    AChem1_1_1 = A0sChem1_1_1 + AsChem1_1_1 * m.fabs(x1_1 - xs1_1)
    AChem1_1_2 = A0sChem1_1_2 + AsChem1_1_2 * m.fabs(x1_2 - xs1_2)
    AChem1_2_1 = A0sChem1_2_1 + AsChem1_2_1 * m.fabs(x2_1 - xs2_1)
    AChem1_2_2 = A0sChem1_2_2 + AsChem1_2_2 * m.fabs(x1_4 - xs1_4)
    AChem1_3_3 = A0sChem1_3_3 + AsChem1_3_3 * m.fabs(x1_3 - xs1_3)
    AChem2_1_1 = A0sChem2_1_1 + AsChem2_1_1 * m.fabs(x2_1 - xs2_1)
    AChem2_2_2 = A0sChem2_2_2 + AsChem2_2_2 * m.fabs(x2_4 - xs2_4)

    # Выводим результат
    chemPot = np.array([mu1_1, mu1_2, mu1_3, mu1_4, mu1_5, mu1_6, mu2_1, mu2_4, mu2_7, mu2_8])  # Массив потенциалов взаимодействия
    kineticMatrix = np.array([AChem1_1_1, AChem1_1_2, AChem1_2_1, AChem1_2_2, AChem1_3_3, AChem2_1_1, AChem2_2_2, Adiff2_2])
    balanceMatrix = []
    return (balanceMatrix, [], -chemPot, kineticMatrix)


# Функция расчета характеристик системы
def CountSystem1(stateCoordinates,  # Координаты состояния
                 systemParameters,  # Параметры системы

                 nu1_1, nu1_2, nu1_3, nu1_4, nu1_5, nu1_6, nu1_7,
                 nu1_8, nu1_9, nu1_10, nu1_11, nu2_1, nu2_2, nu2_3,
                 nu2_4, nu2_5, nu2_6, Adiff1_1, Adiff1_2, Adiff2_1):
    # Получаем числа молей
    x1_1 = stateCoordinates[0]
    x1_2 = stateCoordinates[1]
    x1_3 = stateCoordinates[2]
    x1_4 = stateCoordinates[3]
    x1_5 = stateCoordinates[4]
    x1_6 = stateCoordinates[5]
    x2_1 = stateCoordinates[6]
    x2_4 = stateCoordinates[7]
    x2_7 = stateCoordinates[8]
    x2_8 = stateCoordinates[9]

    # Базовые числа молей
    xs1_1 = systemParameters[0]
    xs1_2 = systemParameters[1]
    xs1_3 = systemParameters[2]
    xs1_4 = systemParameters[3]
    xs1_5 = systemParameters[4]
    xs1_6 = systemParameters[5]
    xs2_1 = systemParameters[6]
    xs2_4 = systemParameters[7]
    xs2_7 = systemParameters[8]
    xs2_8 = systemParameters[9]

    # Базовые химические потенциалы веществ
    mus1_1 = systemParameters[10]
    mus1_2 = systemParameters[11]
    mus1_3 = systemParameters[12]
    mus1_4 = systemParameters[13]
    mus1_5 = systemParameters[14]
    mus1_6 = systemParameters[15]
    mus2_1 = systemParameters[16]
    mus2_4 = systemParameters[17]
    mus2_7 = systemParameters[18]
    mus2_8 = systemParameters[19]
    mu0s1_1 = systemParameters[20]
    mu0s1_2 = systemParameters[21]
    mu0s1_3 = systemParameters[22]
    mu0s1_4 = systemParameters[23]
    mu0s1_5 = systemParameters[24]
    mu0s1_6 = systemParameters[25]
    mu0s2_1 = systemParameters[26]
    mu0s2_4 = systemParameters[27]
    mu0s2_7 = systemParameters[28]
    mu0s2_8 = systemParameters[29]

    # Базовые коэффициенты кинетической матрицы процессов
    A0sChem1_1_1 = systemParameters[30]
    A0sChem1_1_2 = systemParameters[31]
    A0sChem1_2_1 = systemParameters[32]
    A0sChem1_2_2 = systemParameters[33]
    A0sChem1_3_3 = systemParameters[34]
    A0sChem2_1_1 = systemParameters[35]
    A0sChem2_2_2 = systemParameters[36]
    AsChem1_1_1 = systemParameters[37]
    AsChem1_1_2 = systemParameters[38]
    AsChem1_2_1 = systemParameters[39]
    AsChem1_2_2 = systemParameters[40]
    AsChem1_3_3 = systemParameters[41]
    AsChem2_1_1 = systemParameters[42]
    AsChem2_2_2 = systemParameters[43]
    Adiff2_2 = 4.53 + systemParameters[44]

    # Определяем химические потенциалы веществ
    mu1_1 = mus1_1 * m.log(x1_1 / xs1_1) + mu0s1_1
    mu1_2 = mus1_2 * m.log(x1_2 / xs1_2) + mu0s1_2
    mu1_3 = mus1_3 * m.log(x1_3 / xs1_3) + mu0s1_3
    mu1_4 = mus1_4 * m.log(x1_4 / xs1_4) + mu0s1_4
    mu1_5 = mus1_5 * m.log(x1_5 / xs1_5) + mu0s1_5
    mu1_6 = mus1_6 * m.log(x1_6 / xs1_6) + mu0s1_6
    mu2_1 = mus2_1 * m.log(x2_1 / xs2_1) + mu0s2_1
    mu2_4 = mus2_4 * m.log(x2_4 / xs2_4) + mu0s2_4
    mu2_7 = mus2_7 * m.log(x2_7 / xs2_7) + mu0s2_7
    mu2_8 = mus2_8 * m.log(x2_8 / xs2_8) + mu0s2_8

    # Определяем главную кинетическую матрицу процессов
    AChem1_1_1 = A0sChem1_1_1 + AsChem1_1_1 * m.fabs(x1_1 - xs1_1)
    AChem1_1_2 = A0sChem1_1_2 + AsChem1_1_2 * m.fabs(x1_2 - xs1_2)
    AChem1_2_1 = A0sChem1_2_1 + AsChem1_2_1 * m.fabs(x2_1 - xs2_1)
    AChem1_2_2 = A0sChem1_2_2 + AsChem1_2_2 * m.fabs(x1_4 - xs1_4)
    AChem1_3_3 = A0sChem1_3_3 + AsChem1_3_3 * m.fabs(x1_3 - xs1_3)
    AChem2_1_1 = A0sChem2_1_1 + AsChem2_1_1 * m.fabs(x2_1 - xs2_1)
    AChem2_2_2 = A0sChem2_2_2 + AsChem2_2_2 * m.fabs(x2_4 - xs2_4)

    # Определяем термодинамические силы
    Xchem1_1 = nu1_1 * mu1_1 + nu1_2 * mu1_2 - nu1_3 * mu1_4
    Xchem1_2 = nu1_4 * mu1_3 + nu1_5 * mu1_2 - nu1_6 * mu1_6 - nu1_7 * mu1_5
    Xchem1_3 = nu1_8 * mu1_1 + nu1_9 * mu1_3 + nu1_10 * mu1_5 - nu1_11 * mu1_2
    Xchem2_1 = nu2_1 * mu2_1 + nu2_2 * mu2_7 - nu2_3 * mu2_8
    Xchem2_2 = nu2_4 * mu2_4 + nu2_5 * mu2_1 - nu2_6 * mu2_7
    Xdiff1 = mu1_1 - mu2_1
    Xdiff2 = mu1_4 - mu2_4

    # Определяем скорости процессов
    vChem1_1 = AChem1_1_1 * Xchem1_1 + AChem1_1_2 * Xchem1_2
    vChem1_2 = AChem1_2_1 * Xchem1_1 + AChem1_2_2 * Xchem1_2
    vChem1_3 = AChem1_3_3 * Xchem1_3
    vChem2_1 = AChem2_1_1 * Xchem2_1
    vChem2_2 = AChem2_2_2 * Xchem2_2
    vDiff1 = Adiff1_1 * Xdiff1 + Adiff1_2 * Xdiff2
    vDiff2 = Adiff2_1 * Xdiff1 + Adiff2_2 * Xdiff2

    # Определяем скорости изменения координат состояния
    vx1_1 = -nu1_1 * vChem1_1 - nu1_8 * vChem1_3 - vDiff1
    vx1_2 = -nu1_2 * vChem1_1 - nu1_5 * vChem1_2 + nu1_11 * vChem1_3
    vx1_3 = -nu1_4 * vChem1_2 - nu1_9 * vChem1_3
    vx1_4 = nu1_3 * vChem1_1 - vDiff2
    vx1_5 = nu1_7 * vChem1_2 - nu1_10 * vChem1_3
    vx1_6 = nu1_6 * vChem1_2
    vx2_1 = -nu2_1 * vChem2_1 - nu2_5 * vChem2_2 + vDiff1
    vx2_4 = -nu2_4 * vChem2_2 + vDiff2
    vx2_7 = -nu2_2 * vChem2_1 + nu2_6 * vChem2_2
    vx2_8 = nu2_3 * vChem2_1

    # Выводим результат
    chemPot = np.array([mu1_1, mu1_2, mu1_3, mu1_4, mu1_5, mu1_6, mu2_1, mu2_4, mu2_7, mu2_8])  # Массив потенциалов взаимодействия
    Aff = np.array([Xchem1_1, Xchem1_2, Xchem1_3, Xchem2_1, Xchem2_2, Xdiff1, Xdiff2])  # Термодинамические силы
    kineticMatrix = np.array([[AChem1_1_1, AChem1_1_2,        0.0,        0.0,        0.0,      0.0,      0.0],
                              [AChem1_2_1, AChem1_2_2,        0.0,        0.0,        0.0,      0.0,      0.0],
                              [       0.0,        0.0, AChem1_3_3,        0.0,        0.0,      0.0,      0.0],
                              [       0.0,        0.0,        0.0, AChem2_1_1,        0.0,      0.0,      0.0],
                              [       0.0,        0.0,        0.0,        0.0, AChem2_2_2,      0.0,      0.0],
                              [       0.0,        0.0,        0.0,        0.0,        0.0, Adiff1_1, Adiff1_2],
                              [       0.0,        0.0,        0.0,        0.0,        0.0, Adiff2_1, Adiff2_2]])
    vProcesses = np.array([vChem1_1, vChem1_2, vChem1_3, vChem2_1, vChem2_2, vDiff1, vDiff2])
    balanceMatrix = np.array([[-nu1_1,    0.0,  -nu1_8,    0.0,    0.0, -1.0,  0.0],
                              [-nu1_2, -nu1_5,  nu1_11,    0.0,    0.0,  0.0,  0.0],
                              [   0.0, -nu1_4,  -nu1_9,    0.0,    0.0,  0.0,  0.0],
                              [ nu1_3,    0.0,     0.0,    0.0,    0.0,  0.0, -1.0],
                              [   0.0,  nu1_7, -nu1_10,    0.0,    0.0,  0.0,  0.0],
                              [   0.0,  nu1_6,     0.0,    0.0,    0.0,  0.0,  0.0],
                              [   0.0,    0.0,     0.0, -nu2_1, -nu2_5,  1.0,  0.0],
                              [   0.0,    0.0,     0.0,    0.0, -nu2_4,  0.0,  1.0],
                              [   0.0,    0.0,     0.0, -nu2_2,  nu2_6,  0.0,  0.0],
                              [   0.0,    0.0,     0.0,  nu2_3,    0.0,  0.0,  0.0]])
    vx = np.array([vx1_1, vx1_2, vx1_3, vx1_4, vx1_5, vx1_6, vx2_1, vx2_4, vx2_7, vx2_8])
    return (-chemPot, Aff, kineticMatrix, vProcesses, balanceMatrix, vx)
