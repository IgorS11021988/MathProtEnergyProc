import math as m
import numpy as np


# Функция состояния системысистемы
def CountStateQ3(stateCoordinates,  # Координаты состояния
                 reducedTemp,  # Приведенные температуры энергетических степеней свободы
                 systemParameters  # Параметры системы
                 ):
    # Получаем числа молей
    x1_1 = stateCoordinates[0]
    x1_2 = stateCoordinates[1]
    x1_3 = stateCoordinates[2]
    x1_4 = stateCoordinates[3]
    x1_5 = stateCoordinates[4]
    x1_6 = stateCoordinates[5]

    # Получаем температуры
    T1 = reducedTemp[0]

    # Базовые числа молей
    xs1_1 = systemParameters[0]
    xs1_2 = systemParameters[1]
    xs1_3 = systemParameters[2]
    xs1_4 = systemParameters[3]
    xs1_5 = systemParameters[4]
    xs1_6 = systemParameters[5]

    # Базовые химические потенциалы веществ
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

    # Базовые коэффициенты главной кинетической матрицы процессов
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

    # Приведенные тепловые эффекты энергетических степеней свободы
    H1_2 = 3.81 + systemParameters[27]
    H1_5 = 6.21 + systemParameters[28]

    # Определяем химические потенциалы веществ
    mu1_1 = mus1_1 * m.log(x1_1 / xs1_1) + mu0s1_1
    mu1_2 = mus1_2 * m.log(x1_2 / xs1_2) + mu0s1_2
    mu1_3 = mus1_3 * m.log(x1_3 / xs1_3) + mu0s1_3
    mu1_4 = mus1_4 * m.log(x1_4 / xs1_4) + mu0s1_4
    mu1_5 = mus1_5 * m.log(x1_5 / xs1_5) + mu0s1_5
    mu1_6 = mus1_6 * m.log(x1_6 / xs1_6) + mu0s1_6

    # Определяем главную кинетическую матрицу процессов
    AChem1_1_1 = A0sChem1_1_1 + AsChem1_1_1 * m.fabs(x1_1 - xs1_1)
    AChem1_1_2 = A0sChem1_1_2 + AsChem1_1_2 * m.fabs(x1_2 - xs1_2)
    AChem1_2_1 = A0sChem1_2_1 + AsChem1_2_1 * m.fabs(x1_1 - xs1_1)
    AChem1_2_2 = A0sChem1_2_2 + AsChem1_2_2 * m.fabs(x1_4 - xs1_4)
    AChem1_3_3 = A0sChem1_3_3 + AsChem1_3_3 * m.fabs(x1_3 - xs1_3)

    # Выводим результат
    chemPot = np.array([mu1_1, mu1_2, mu1_3, mu1_4, mu1_5, mu1_6])  # Массив потенциалов взаимодействия
    Temp = np.array([T1])
    Betas = np.array([])
    kineticMatrixCPCP = np.array([AChem1_1_1, AChem1_1_2, AChem1_2_1, AChem1_2_2, AChem1_3_3])
    kineticMatrixCPHeat = []
    kineticMatrixHeatCP = []
    kineticMatrixHeatHeat = []
    balanceMatrix = []
    invC = []
    powH = np.array([[H1_2, H1_5]])
    return (balanceMatrix, [], [], Temp, -chemPot, [], Betas,
            kineticMatrixCPCP, kineticMatrixCPHeat,
            kineticMatrixHeatCP, kineticMatrixHeatHeat,
            invC, powH)


# Функция расчета характеристик системы
def CountSystemQ3(stateCoordinates,  # Координаты состояния
                  reducedTemp,  # Приведенные температуры энергетических степеней свободы
                  systemParameters,  # Параметры системы

                  nu1_1, nu1_2, nu1_3, nu1_4, nu1_5, nu1_6, nu1_7,
                  nu1_8, nu1_9, nu1_10, nu1_11, invC1, H1_1):
    # Опорная температура
    Tbase = 293

    # Получаем числа молей
    x1_1 = stateCoordinates[0]
    x1_2 = stateCoordinates[1]
    x1_3 = stateCoordinates[2]
    x1_4 = stateCoordinates[3]
    x1_5 = stateCoordinates[4]
    x1_6 = stateCoordinates[5]

    # Получаем температуры
    T1 = reducedTemp[0]

    # Базовые числа молей
    xs1_1 = systemParameters[0]
    xs1_2 = systemParameters[1]
    xs1_3 = systemParameters[2]
    xs1_4 = systemParameters[3]
    xs1_5 = systemParameters[4]
    xs1_6 = systemParameters[5]

    # Базовые химические потенциалы веществ
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

    # Базовые коэффициенты главной кинетической матрицы процессов
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

    # Приведенные тепловые эффекты энергетических степеней свободы
    H1_2 = 3.81 + systemParameters[27]
    H1_5 = 6.21 + systemParameters[28]

    # Определяем химические потенциалы веществ
    mu1_1 = mus1_1 * m.log(x1_1 / xs1_1) + mu0s1_1
    mu1_2 = mus1_2 * m.log(x1_2 / xs1_2) + mu0s1_2
    mu1_3 = mus1_3 * m.log(x1_3 / xs1_3) + mu0s1_3
    mu1_4 = mus1_4 * m.log(x1_4 / xs1_4) + mu0s1_4
    mu1_5 = mus1_5 * m.log(x1_5 / xs1_5) + mu0s1_5
    mu1_6 = mus1_6 * m.log(x1_6 / xs1_6) + mu0s1_6

    # Определяем главную кинетическую матрицу процессов
    AChem1_1_1 = A0sChem1_1_1 + AsChem1_1_1 * m.fabs(x1_1 - xs1_1)
    AChem1_1_2 = A0sChem1_1_2 + AsChem1_1_2 * m.fabs(x1_2 - xs1_2)
    AChem1_2_1 = A0sChem1_2_1 + AsChem1_2_1 * m.fabs(x1_1 - xs1_1)
    AChem1_2_2 = A0sChem1_2_2 + AsChem1_2_2 * m.fabs(x1_4 - xs1_4)
    AChem1_3_3 = A0sChem1_3_3 + AsChem1_3_3 * m.fabs(x1_3 - xs1_3)

    # Определяем термодинамические силы
    Xchem1_1 = (nu1_1 * mu1_1 + nu1_2 * mu1_2 - nu1_3 * mu1_4) * Tbase / T1
    Xchem1_2 = (nu1_4 * mu1_3 + nu1_5 * mu1_2 - nu1_6 * mu1_6 - nu1_7 * mu1_5) * Tbase / T1
    Xchem1_3 = (nu1_8 * mu1_1 + nu1_9 * mu1_3 + nu1_10 * mu1_5 - nu1_11 * mu1_2) * Tbase / T1

    # Определяем скорости процессов
    vChem1_1 = AChem1_1_1 * Xchem1_1 + AChem1_1_2 * Xchem1_2
    vChem1_2 = AChem1_2_1 * Xchem1_1 + AChem1_2_2 * Xchem1_2
    vChem1_3 = AChem1_3_3 * Xchem1_3

    # Определяем некомпенсированные теплоты
    vHeatChem1_1 = (nu1_1 * mu1_1 + nu1_2 * mu1_2 - nu1_3 * mu1_4) * vChem1_1
    vHeatChem1_2 = (nu1_4 * mu1_3 + nu1_5 * mu1_2 - nu1_6 * mu1_6 - nu1_7 * mu1_5) * vChem1_2
    vHeatChem1_3 = (nu1_8 * mu1_1 + nu1_9 * mu1_3 + nu1_10 * mu1_5 - nu1_11 * mu1_2) * vChem1_3

    # Определяем теплоты, полученные энергетическими степенями свободы
    vQPow1 = vHeatChem1_1 + vHeatChem1_2 + vHeatChem1_3

    # Определяем скорости изменения координат состояния
    vx1_1 = -nu1_1 * vChem1_1 - nu1_8 * vChem1_3
    vx1_2 = -nu1_2 * vChem1_1 - nu1_5 * vChem1_2 + nu1_11 * vChem1_3
    vx1_3 = -nu1_4 * vChem1_2 - nu1_9 * vChem1_3
    vx1_4 = nu1_3 * vChem1_1
    vx1_5 = nu1_7 * vChem1_2 - nu1_10 * vChem1_3
    vx1_6 = nu1_6 * vChem1_2

    # Определяем скорости изменения внутренних энергеий энергетических степеней свободы
    vUPow1 = vQPow1 + mu1_1 * vx1_1 + mu1_2 * vx1_2 + mu1_3 * vx1_3 + mu1_4 * vx1_4 + mu1_5 * vx1_5 + mu1_6 * vx1_6

    # Определяем скорости изменения температур энергетических степеней свободы
    vT1 = invC1 * vUPow1 + H1_1 * vx1_1 + H1_2 * vx1_2 + H1_5 * vx1_5

    # Выводим результат
    chemPot = np.array([[mu1_1, mu1_2, mu1_3, mu1_4, mu1_5, mu1_6]])  # Массив потенциалов взаимодействия
    Temp = np.array([T1])
    Betas = np.array([[1],
                      [1],
                      [1]])
    Aff = np.array([Xchem1_1, Xchem1_2, Xchem1_3])  # Термодинамические силы
    AffHeat = np.array([])
    kineticMatrixCPCP = np.array([[AChem1_1_1, AChem1_1_2,        0.0],
                                  [AChem1_2_1, AChem1_2_2,        0.0],
                                  [       0.0,        0.0, AChem1_3_3]])
    kineticMatrixCPHeat = np.array([[],
                                    [],
                                    []])
    kineticMatrixHeatCP = np.array([])
    kineticMatrixHeatHeat = np.array([])
    vProcesses = np.array([vChem1_1, vChem1_2, vChem1_3])
    vQTransf = np.array([])
    vHeatProcesses = np.array([vHeatChem1_1, vHeatChem1_2, vHeatChem1_3])
    vQPows = np.array([vQPow1])
    heatTransferMatrix = np.array([[]])
    balanceMatrix = np.array([[-nu1_1,    0.0,  -nu1_8],
                              [-nu1_2, -nu1_5,  nu1_11],
                              [   0.0, -nu1_4,  -nu1_9],
                              [ nu1_3,    0.0,     0.0],
                              [   0.0,  nu1_7, -nu1_10],
                              [   0.0,  nu1_6,     0.0]])
    vx = np.array([vx1_1, vx1_2, vx1_3, vx1_4, vx1_5, vx1_6])
    vUPows = np.array([vUPow1])
    invC = np.array([[invC1]])
    powH = np.array([[H1_1, H1_2, 0.0, 0.0, H1_5, 0.0]])
    vT = np.array([vT1])
    return (-chemPot, np.array([]), Temp, Betas, Aff, AffHeat,
            kineticMatrixCPCP, kineticMatrixCPHeat,
            kineticMatrixHeatCP, kineticMatrixHeatHeat,
            vProcesses, vQTransf, vHeatProcesses,
            vQPows, heatTransferMatrix, balanceMatrix,
            vx, vUPows, invC, powH, vT)
