import math as m
import numpy as np

#Функция состояния системысистемы
def CountStateQ2(stateCoordinates,#Координаты состояния
                 reducedTemp,#Приведенные температуры энергетических степеней свободы
                 systemParameters#Параметры системы
                 ):
    #Получаем числа молей
    x1_1 = stateCoordinates[0]
    x1_2 = stateCoordinates[1]
    x1_3 = stateCoordinates[2]
    x1_4 = stateCoordinates[3]
    x1_5 = stateCoordinates[4]
    x1_6 = stateCoordinates[5]
    x1_7 = stateCoordinates[10]
    x2_1 = stateCoordinates[6]
    x2_4 = stateCoordinates[7]
    x2_7 = stateCoordinates[8]
    x2_8 = stateCoordinates[9]
    x2_9 = stateCoordinates[11]
    
    #Получаем температуры
    T1 = reducedTemp[0]
    T2 = reducedTemp[1]
    
    #Доли распределения некомпенсированных теплот по энергетических степеням свободы
    beta1_1 = systemParameters[0]
    beta1_2 = systemParameters[1]
    
    #Базовые числа молей
    xs1_1 = systemParameters[2]
    xs1_2 = systemParameters[3]
    xs1_3 = systemParameters[4]
    xs1_4 = systemParameters[5]
    xs1_5 = systemParameters[6]
    xs1_6 = systemParameters[7]
    xs1_7 = systemParameters[53]
    xs2_1 = systemParameters[8]
    xs2_4 = systemParameters[9]
    xs2_7 = systemParameters[10]
    xs2_8 = systemParameters[11]
    
    #Базовые химические потенциалы веществ
    mus1_1 = systemParameters[12]
    mus1_2 = systemParameters[13]
    mus1_3 = systemParameters[14]
    mus1_4 = systemParameters[15]
    mus1_5 = systemParameters[16]
    mus1_6 = systemParameters[17]
    mus1_7 = systemParameters[54]
    mus2_1 = systemParameters[18]
    mus2_4 = systemParameters[19]
    mus2_7 = systemParameters[20]
    mus2_8 = systemParameters[21]
    mu0s1_1 = systemParameters[22]
    mu0s1_2 = systemParameters[23]
    mu0s1_3 = systemParameters[24]
    mu0s1_4 = systemParameters[25]
    mu0s1_5 = systemParameters[26]
    mu0s1_6 = systemParameters[27]
    mu0s1_7 = systemParameters[55]
    mu0s2_1 = systemParameters[28]
    mu0s2_4 = systemParameters[29]
    mu0s2_7 = systemParameters[30]
    mu0s2_8 = systemParameters[31]
    
    #Базовые коэффициенты главной кинетической матрицы процессов
    A0sChem1_1_1 = systemParameters[32]
    A0sChem1_1_2 = systemParameters[33]
    A0sChem1_2_1 = systemParameters[34]
    A0sChem1_2_2 = systemParameters[35]
    A0sChem1_3_3 = systemParameters[36]
    A0sChem2_1_1 = systemParameters[37]
    A0sChem2_2_2 = systemParameters[38]
    AsChem1_1_1 = systemParameters[39]
    AsChem1_1_2 = systemParameters[40]
    AsChem1_2_1 = systemParameters[41]
    AsChem1_2_2 = systemParameters[42]
    AsChem1_3_3 = systemParameters[43]
    AsChem2_1_1 = systemParameters[44]
    AsChem2_2_2 = systemParameters[45]
    Adiff2_2 = 4.53 + systemParameters[46]
    
    #Базовые коэффициенты главной кинетической матрицы переноса теплоты
    AHeat = 11.1 + systemParameters[47]
    
    #Обратные теплоемкости энергетических степеней свободы
    invC2 = 3.39 + systemParameters[48]
    
    #Приведенные тепловые эффекты энергетических степеней свободы
    H1_2 = 3.81 + systemParameters[49]
    H1_5 = 6.21 + systemParameters[50]
    H2_1 = 0.51 + systemParameters[51]
    
    #Внешний поток теплоты
    ExtQEnPowOkr = systemParameters[52]
    
    #Внешние потоки веществ
    xExt1_2 = systemParameters[56]
    xExt2_8 = systemParameters[57]
        
    #Потенциалы взаимодействия между энергетическими степенями свободы
    muBet1_2 = systemParameters[58]
    
    #Перекрестные коэффициенты передачи теплоты
    ADiffHeat1 = systemParameters[59]
    AHeatDiff1 = systemParameters[60]
    
    #Определяем химические потенциалы веществ
    mu1_1 = mus1_1*m.log(x1_1 / xs1_1) + mu0s1_1
    mu1_2 = mus1_2*m.log(x1_2 / xs1_2) + mu0s1_2
    mu1_3 = mus1_3*m.log(x1_3 / xs1_3) + mu0s1_3
    mu1_4 = mus1_4*m.log(x1_4 / xs1_4) + mu0s1_4
    mu1_5 = mus1_5*m.log(x1_5 / xs1_5) + mu0s1_5
    mu1_6 = mus1_6*m.log(x1_6 / xs1_6) + mu0s1_6
    mu1_7 = mus1_7*m.log(x1_7 / xs1_7) + mu0s1_7
    mu2_1 = mus2_1*m.log(x2_1 / xs2_1) + mu0s2_1
    mu2_4 = mus2_4*m.log(x2_4 / xs2_4) + mu0s2_4
    mu2_7 = mus2_7*m.log(x2_7 / xs2_7) + mu0s2_7
    mu2_8 = mus2_8*m.log(x2_8 / xs2_8) + mu0s2_8
    
    #Определяем главную кинетическую матрицу процессов
    AChem1_1_1 = A0sChem1_1_1 + AsChem1_1_1 * m.fabs(x1_1 - xs1_1)
    AChem1_1_2 = A0sChem1_1_2 + AsChem1_1_2 * m.fabs(x1_2 - xs1_2)
    AChem1_2_1 = A0sChem1_2_1 + AsChem1_2_1 * m.fabs(x2_1 - xs2_1)
    AChem1_2_2 = A0sChem1_2_2 + AsChem1_2_2 * m.fabs(x1_4 - xs1_4)
    AChem1_3_3 = A0sChem1_3_3 + AsChem1_3_3 * m.fabs(x1_3 - xs1_3)
    AChem2_1_1 = A0sChem2_1_1 + AsChem2_1_1 * m.fabs(x2_1 - xs2_1)
    AChem2_2_2 = A0sChem2_2_2 + AsChem2_2_2 * m.fabs(x2_4 - xs2_4)
    
    #Выводим результат
    chemPot = np.array([mu1_1,mu1_2,mu1_3,mu1_4,mu1_5,mu1_6,mu1_7, mu2_1, mu2_4, mu2_7, mu2_8])#Массив потенциалов взаимодействия
    Temp = np.array([T1, T2])
    Betas = np.array([beta1_1, beta1_2])
    kineticMatrixCPCP = np.array([AChem1_1_1, AChem1_1_2, AChem1_2_1, AChem1_2_2, AChem1_3_3, AChem2_1_1, AChem2_2_2, Adiff2_2])
    kineticMatrixCPHeat = np.array([ADiffHeat1])
    kineticMatrixHeatCP = np.array([AHeatDiff1])
    kineticMatrixHeatHeat = np.array([[AHeat]])
    balanceMatrix = []
    invC = np.array([invC2])
    powH = np.array([[H1_2, H1_5, H2_1]])
    heatStreams = np.array([ExtQEnPowOkr])
    Streams = np.array([xExt1_2, xExt2_8])
    potBet = np.array([-muBet1_2])
    return (balanceMatrix,Streams,heatStreams,Temp,
            -chemPot,potBet,Betas,kineticMatrixCPCP,
            kineticMatrixCPHeat,kineticMatrixHeatCP,
            kineticMatrixHeatHeat,invC, powH)

#Функция расчета характеристик системы
def CountSystemQ2(stateCoordinates,#Координаты состояния
                  reducedTemp,#Приведенные температуры энергетических степеней свободы
                  systemParameters,#Параметры системы
                 
                  nu1_1,nu1_2,nu1_3,nu1_4,nu1_5,nu1_6,nu1_7,
                  nu1_8,nu1_9,nu1_10,nu1_11,nu1_12,nu1_13,nu1_14,
                  nu1_15,nu2_1,nu2_2,nu2_3,nu2_4,nu2_5,nu2_6,
                  nu2_7,nu2_8,nu2_9,beta2_1,beta2_2,Adiff1_1,
                  Adiff1_2,Adiff2_1,ADiffHeat2,AHeatDiff2,invC1,
                  H1_1,H2_7,Tokr,AHeatOkr,ExtQEnPow1,AChem1_4_4,
                  mu2_9,AChem2_3_3,AChem2_1_2,AChem2_2_1,Adiff3_3,
                  ADiffHeat3,AHeatDiff3,xExt1_6,muBet1_3,muBet1_5):
    #Опорная температура
    Tbase = 293
    
    #Получаем числа молей
    x1_1 = stateCoordinates[0]
    x1_2 = stateCoordinates[1]
    x1_3 = stateCoordinates[2]
    x1_4 = stateCoordinates[3]
    x1_5 = stateCoordinates[4]
    x1_6 = stateCoordinates[5]
    x1_7 = stateCoordinates[10]
    x2_1 = stateCoordinates[6]
    x2_4 = stateCoordinates[7]
    x2_7 = stateCoordinates[8]
    x2_8 = stateCoordinates[9]
    x2_9 = stateCoordinates[11]
    
    #Получаем температуры
    T1 = reducedTemp[0]
    T2 = reducedTemp[1]
    
    #Доли распределения некомпенсированных теплот по энергетических степеням свободы
    beta1_1 = systemParameters[0]
    beta1_2 = systemParameters[1]
    
    #Базовые числа молей
    xs1_1 = systemParameters[2]
    xs1_2 = systemParameters[3]
    xs1_3 = systemParameters[4]
    xs1_4 = systemParameters[5]
    xs1_5 = systemParameters[6]
    xs1_6 = systemParameters[7]
    xs1_7 = systemParameters[53]
    xs2_1 = systemParameters[8]
    xs2_4 = systemParameters[9]
    xs2_7 = systemParameters[10]
    xs2_8 = systemParameters[11]
    
    #Базовые химические потенциалы веществ
    mus1_1 = systemParameters[12]
    mus1_2 = systemParameters[13]
    mus1_3 = systemParameters[14]
    mus1_4 = systemParameters[15]
    mus1_5 = systemParameters[16]
    mus1_6 = systemParameters[17]
    mus1_7 = systemParameters[54]
    mus2_1 = systemParameters[18]
    mus2_4 = systemParameters[19]
    mus2_7 = systemParameters[20]
    mus2_8 = systemParameters[21]
    mu0s1_1 = systemParameters[22]
    mu0s1_2 = systemParameters[23]
    mu0s1_3 = systemParameters[24]
    mu0s1_4 = systemParameters[25]
    mu0s1_5 = systemParameters[26]
    mu0s1_6 = systemParameters[27]
    mu0s1_7 = systemParameters[55]
    mu0s2_1 = systemParameters[28]
    mu0s2_4 = systemParameters[29]
    mu0s2_7 = systemParameters[30]
    mu0s2_8 = systemParameters[31]
    
    #Базовые коэффициенты главной кинетической матрицы процессов
    A0sChem1_1_1 = systemParameters[32]
    A0sChem1_1_2 = systemParameters[33]
    A0sChem1_2_1 = systemParameters[34]
    A0sChem1_2_2 = systemParameters[35]
    A0sChem1_3_3 = systemParameters[36]
    A0sChem2_1_1 = systemParameters[37]
    A0sChem2_2_2 = systemParameters[38]
    AsChem1_1_1 = systemParameters[39]
    AsChem1_1_2 = systemParameters[40]
    AsChem1_2_1 = systemParameters[41]
    AsChem1_2_2 = systemParameters[42]
    AsChem1_3_3 = systemParameters[43]
    AsChem2_1_1 = systemParameters[44]
    AsChem2_2_2 = systemParameters[45]
    Adiff2_2 = 4.53 + systemParameters[46]
    
    #Базовые коэффициенты главной кинетической матрицы переноса теплоты
    AHeat = 11.1 + systemParameters[47]
    
    #Обратные теплоемкости энергетических степеней свободы
    invC2 = 3.39 + systemParameters[48]
    
    #Приведенные тепловые эффекты энергетических степеней свободы
    H1_2 = 3.81 + systemParameters[49]
    H1_5 = 6.21 + systemParameters[50]
    H2_1 = 0.51 + systemParameters[51]
    
    #Внешний поток теплоты
    ExtQEnPowOkr = systemParameters[52]
    
    #Внешние потоки веществ
    xExt1_2 = systemParameters[56]
    xExt2_8 = systemParameters[57]
    
    #Потенциалы взаимодействия между энергетическими степенями свободы
    muBet1_2 = systemParameters[58]
    
    #Перекрестные коэффициенты передачи теплоты
    ADiffHeat1 = systemParameters[59]
    AHeatDiff1 = systemParameters[60]
    
    #Определяем химические потенциалы веществ
    mu1_1 = mus1_1*m.log(x1_1 / xs1_1) + mu0s1_1
    mu1_2 = mus1_2*m.log(x1_2 / xs1_2) + mu0s1_2
    mu1_3 = mus1_3*m.log(x1_3 / xs1_3) + mu0s1_3
    mu1_4 = mus1_4*m.log(x1_4 / xs1_4) + mu0s1_4
    mu1_5 = mus1_5*m.log(x1_5 / xs1_5) + mu0s1_5
    mu1_6 = mus1_6*m.log(x1_6 / xs1_6) + mu0s1_6
    mu1_7 = mus1_7*m.log(x1_7 / xs1_7) + mu0s1_7
    mu2_1 = mus2_1*m.log(x2_1 / xs2_1) + mu0s2_1
    mu2_4 = mus2_4*m.log(x2_4 / xs2_4) + mu0s2_4
    mu2_7 = mus2_7*m.log(x2_7 / xs2_7) + mu0s2_7
    mu2_8 = mus2_8*m.log(x2_8 / xs2_8) + mu0s2_8
    
    #Определяем главную кинетическую матрицу процессов
    AChem1_1_1 = A0sChem1_1_1 + AsChem1_1_1 * m.fabs(x1_1 - xs1_1)
    AChem1_1_2 = A0sChem1_1_2 + AsChem1_1_2 * m.fabs(x1_2 - xs1_2)
    AChem1_2_1 = A0sChem1_2_1 + AsChem1_2_1 * m.fabs(x2_1 - xs2_1)
    AChem1_2_2 = A0sChem1_2_2 + AsChem1_2_2 * m.fabs(x1_4 - xs1_4)
    AChem1_3_3 = A0sChem1_3_3 + AsChem1_3_3 * m.fabs(x1_3 - xs1_3)
    AChem2_1_1 = A0sChem2_1_1 + AsChem2_1_1 * m.fabs(x2_1 - xs2_1)
    AChem2_2_2 = A0sChem2_2_2 + AsChem2_2_2 * m.fabs(x2_4 - xs2_4)
    
    #Определяем термодинамические силы
    Xchem1_1 = (nu1_1*mu1_1 + nu1_2*(mu1_2 + muBet1_2) - nu1_3*mu1_4)*Tbase / T1
    Xchem1_2 = (nu1_4*(mu1_3 + muBet1_3) + nu1_5*(mu1_2 + muBet1_2) - nu1_6*mu1_6 - nu1_7*(mu1_5 + muBet1_5))*Tbase / T1
    Xchem1_3 = (nu1_8*mu1_1 + nu1_9*(mu1_3 + muBet1_3) + nu1_10*(mu1_5 + muBet1_5) - nu1_11*(mu1_2 + muBet1_2))*Tbase / T1
    Xchem1_4 = (nu1_12*mu1_7 + nu1_13*(mu1_3 + muBet1_3) - nu1_14*mu1_1 - nu1_15*mu1_6)*Tbase / T1
    Xchem2_1 = (nu2_1*mu2_1 + nu2_2*mu2_7 - nu2_3*mu2_8)*Tbase / T2
    Xchem2_2 = (nu2_4*mu2_4 + nu2_5*mu2_1 - nu2_6*mu2_7)*Tbase / T2
    Xchem2_3 = (nu2_7*mu2_8 - nu2_8*mu2_9 - nu2_9*mu2_4)*Tbase / T2
    Xdiff1 = (beta1_1*Tbase / T1 + beta1_2*Tbase / T2)*(mu1_1 - mu2_1)
    Xdiff2 = (beta2_1*Tbase / T1 + beta2_2*Tbase / T2)*(mu1_4 - mu2_4)
    Xdiff3 = (beta2_1*Tbase / T1 + beta2_2*Tbase / T2)*(mu1_7 - mu2_7)
    XQ = Tbase / T2 - Tbase / T1
    XQokr = Tbase / Tokr - Tbase / T1
    
    #Определяем скорости процессов
    vChem1_1 = AChem1_1_1*Xchem1_1 + AChem1_1_2*Xchem1_2
    vChem1_2 = AChem1_2_1*Xchem1_1 + AChem1_2_2*Xchem1_2
    vChem1_3 = AChem1_3_3*Xchem1_3
    vChem1_4 = AChem1_4_4*Xchem1_4
    vChem2_1 = AChem2_1_1*Xchem2_1 + AChem2_1_2*Xchem2_2
    vChem2_2 = AChem2_2_1*Xchem2_1 + AChem2_2_2*Xchem2_2
    vChem2_3 = AChem2_3_3*Xchem2_3
    vDiff1 = Adiff1_1*Xdiff1 + Adiff1_2*Xdiff2 + ADiffHeat1*XQ
    vDiff2 = Adiff2_1*Xdiff1 + Adiff2_2*Xdiff2 + ADiffHeat2*XQ
    vDiff3 = Adiff3_3*Xdiff3 + ADiffHeat3*XQ
    
    #Определяем тепловой поток
    vQ = AHeat*XQ + AHeatDiff1*Xdiff1 + AHeatDiff2*Xdiff2 + AHeatDiff3*Xdiff3
    vQTrOkr = AHeatOkr*XQokr
    
    #Определяем некомпенсированные теплоты
    vHeatChem1_1 = (nu1_1*mu1_1 + nu1_2*(mu1_2 + muBet1_2) - nu1_3*mu1_4)*vChem1_1
    vHeatChem1_2 = (nu1_4*(mu1_3 + muBet1_3) + nu1_5*(mu1_2 + muBet1_2) - nu1_6*mu1_6 - nu1_7*(mu1_5 + muBet1_5))*vChem1_2
    vHeatChem1_3 = (nu1_8*mu1_1 + nu1_9*(mu1_3 + muBet1_3) + nu1_10*(mu1_5 + muBet1_5) - nu1_11*(mu1_2 + muBet1_2))*vChem1_3
    vHeatChem1_4 = (nu1_12*mu1_7 + nu1_13*(mu1_3 + muBet1_3) - nu1_14*mu1_1 - nu1_15*mu1_6)*vChem1_4
    vHeatChem2_1 = (nu2_1*mu2_1 + nu2_2*mu2_7 - nu2_3*mu2_8)*vChem2_1
    vHeatChem2_2 = (nu2_4*mu2_4 + nu2_5*mu2_1 - nu2_6*mu2_7)*vChem2_2
    vHeatChem2_3 = (nu2_7*mu2_8 - nu2_8*mu2_9 - nu2_9*mu2_4)*vChem2_3
    vHeatDiff1 = (mu1_1 - mu2_1)*vDiff1
    vHeatDiff2 = (mu1_4 - mu2_4)*vDiff2
    vHeatDiff3 = (mu1_7 - mu2_7)*vDiff3
    
    #Определяем теплоты, полученные энергетическими степенями свободы
    vQPow1 = vHeatChem1_1 + vHeatChem1_2 + vHeatChem1_3 + vHeatChem1_4 + beta1_1*vHeatDiff1 + beta2_1*vHeatDiff2 + beta2_1*vHeatDiff3 - vQ - vQTrOkr + ExtQEnPow1
    vQPow2 = vHeatChem2_1 + vHeatChem2_2 + vHeatChem2_3 + beta1_2*vHeatDiff1 + beta2_2*vHeatDiff2 + beta2_2*vHeatDiff3 + vQ
    vQokr = vQTrOkr + ExtQEnPowOkr
    
    #Определяем скорости изменения координат состояния
    vx1_1 = -nu1_1*vChem1_1 - nu1_8*vChem1_3 + nu1_14*vChem1_4 - vDiff1
    vx1_2 = -nu1_2*vChem1_1 - nu1_5*vChem1_2 + nu1_11*vChem1_3 + xExt1_2
    vx1_3 = -nu1_4*vChem1_2 - nu1_9*vChem1_3 - nu1_13*vChem1_4
    vx1_4 = nu1_3*vChem1_1 - vDiff2
    vx1_5 = nu1_7*vChem1_2 - nu1_10*vChem1_3
    vx1_6 = nu1_6*vChem1_2 + nu1_15*vChem1_4 + xExt1_6
    vx1_7 = -nu1_12*vChem1_4 - vDiff3
    vx2_1 = -nu2_1*vChem2_1 - nu2_5*vChem2_2 + vDiff1
    vx2_4 = -nu2_4*vChem2_2 + vDiff2 + nu2_9*vChem2_3
    vx2_7 = -nu2_2*vChem2_1 + nu2_6*vChem2_2 + vDiff3
    vx2_8 = nu2_3*vChem2_1 - nu2_7*vChem2_3 + xExt2_8
    vx2_9 = nu2_8*vChem2_3
    
    #Определяем скорости изменения внутренних энергеий энергетических степеней свободы
    vUPow1 = vQPow1 + mu1_1*vx1_1 + mu1_2*vx1_2 + mu1_3*vx1_3 + mu1_4*vx1_4 + mu1_5*vx1_5 + mu1_6*vx1_6 + mu1_7*vx1_7
    vUPow2 = vQPow2 + mu2_1*vx2_1 + mu2_4*vx2_4 + mu2_7*vx2_7 + mu2_8*vx2_8 + mu2_9*vx2_9
    
    #Определяем скорости изменения температур энергетических степеней свободы
    vT1 = invC1*vUPow1 + H1_1*vx1_1 + H1_2*vx1_2 + H1_5*vx1_5
    vT2 = invC2*vUPow2 + H2_1*vx2_1 + H2_7*vx2_7
    
    #Выводим результат
    chemPot = np.array([[mu1_1,mu1_2,mu1_3,mu1_4,mu1_5,mu1_6, mu1_7,     0,     0,     0,     0,     0],
                        [    0,    0,    0,    0,    0,    0,     0, mu2_1, mu2_4, mu2_7, mu2_8, mu2_9],
                        [    0,    0,    0,    0,    0,    0,     0,     0,     0,     0,     0,     0]])#Массив потенциалов взаимодействия
    Temp = np.array([T1, T2, Tokr])
    Betas = np.array([[      1,       0, 0],
                      [      1,       0, 0],
                      [      1,       0, 0],
                      [      1,       0, 0],
                      [      0,       1, 0],
                      [      0,       1, 0],
                      [      0,       1, 0],
                      [beta1_1, beta1_2, 0],
                      [beta2_1, beta2_2, 0],
                      [beta2_1, beta2_2, 0]])
    Aff = np.array([Xchem1_1, Xchem1_2, Xchem1_3, Xchem1_4, Xchem2_1, Xchem2_2, Xchem2_3, Xdiff1, Xdiff2, Xdiff3])#Термодинамические силы
    AffHeat = np.array([XQ, XQokr])
    kineticMatrixCPCP = np.array([[AChem1_1_1, AChem1_1_2,        0.0,        0.0,        0.0,        0.0,        0.0,      0.0,      0.0,      0.0],
                                  [AChem1_2_1, AChem1_2_2,        0.0,        0.0,        0.0,        0.0,        0.0,      0.0,      0.0,      0.0],
                                  [       0.0,        0.0, AChem1_3_3,        0.0,        0.0,        0.0,        0.0,      0.0,      0.0,      0.0],
                                  [       0.0,        0.0,        0.0, AChem1_4_4,        0.0,        0.0,        0.0,      0.0,      0.0,      0.0],
                                  [       0.0,        0.0,        0.0,        0.0, AChem2_1_1, AChem2_1_2,        0.0,      0.0,      0.0,      0.0],
                                  [       0.0,        0.0,        0.0,        0.0, AChem2_2_1, AChem2_2_2,        0.0,      0.0,      0.0,      0.0],
                                  [       0.0,        0.0,        0.0,        0.0,        0.0,        0.0, AChem2_3_3,      0.0,      0.0,      0.0],
                                  [       0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0, Adiff1_1, Adiff1_2,      0.0],
                                  [       0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0, Adiff2_1, Adiff2_2,      0.0],
                                  [       0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,      0.0,      0.0, Adiff3_3]])
    kineticMatrixCPHeat = np.array([       [0.0, 0.0],
                                           [0.0, 0.0],
                                           [0.0, 0.0],
                                           [0.0, 0.0],
                                           [0.0, 0.0],
                                           [0.0, 0.0],
                                           [0.0, 0.0],
                                    [ADiffHeat1, 0.0],
                                    [ADiffHeat2, 0.0],
                                    [ADiffHeat3, 0.0]])
    kineticMatrixHeatCP = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, AHeatDiff1, AHeatDiff2, AHeatDiff3],
                                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,        0.0,        0.0,        0.0]])
    kineticMatrixHeatHeat = np.array([[AHeat,      0.0],
                                      [  0.0, AHeatOkr]])
    vProcesses = np.array([vChem1_1, vChem1_2, vChem1_3, vChem1_4, vChem2_1, vChem2_2, vChem2_3, vDiff1, vDiff2, vDiff3])
    vQTransf = np.array([vQ, vQTrOkr])
    vHeatProcesses = np.array([vHeatChem1_1, vHeatChem1_2, vHeatChem1_3, vHeatChem1_4, vHeatChem2_1, vHeatChem2_2, vHeatChem2_3, vHeatDiff1, vHeatDiff2, vHeatDiff3])
    vQPows = np.array([vQPow1, vQPow2, vQokr])
    heatTransferMatrix = np.array([[-1, -1],
                                    [1,  0],
                                    [0,  1]])
    balanceMatrix = np.array([[-nu1_1,    0.0,  -nu1_8,  nu1_14,    0.0,    0.0,    0.0, -1.0,  0.0,  0.0],
                              [-nu1_2, -nu1_5,  nu1_11,     0.0,    0.0,    0.0,    0.0,  0.0,  0.0,  0.0],
                              [   0.0, -nu1_4,  -nu1_9, -nu1_13,    0.0,    0.0,    0.0,  0.0,  0.0,  0.0],
                              [ nu1_3,    0.0,     0.0,     0.0,    0.0,    0.0,    0.0,  0.0, -1.0,  0.0],
                              [   0.0,  nu1_7, -nu1_10,     0.0,    0.0,    0.0,    0.0,  0.0,  0.0,  0.0],
                              [   0.0,  nu1_6,     0.0,  nu1_15,    0.0,    0.0,    0.0,  0.0,  0.0,  0.0],
                              [   0.0,    0.0,     0.0, -nu1_12,    0.0,    0.0,    0.0,  0.0,  0.0, -1.0],
                              [   0.0,    0.0,     0.0,     0.0, -nu2_1, -nu2_5,    0.0,  1.0,  0.0,  0.0],
                              [   0.0,    0.0,     0.0,     0.0,    0.0, -nu2_4,  nu2_9,  0.0,  1.0,  0.0],
                              [   0.0,    0.0,     0.0,     0.0, -nu2_2,  nu2_6,    0.0,  0.0,  0.0,  1.0],
                              [   0.0,    0.0,     0.0,     0.0,  nu2_3,    0.0, -nu2_7,  0.0,  0.0,  0.0],
                              [   0.0,    0.0,     0.0,     0.0,    0.0,    0.0,  nu2_8,  0.0,  0.0,  0.0]])
    vx = np.array([vx1_1, vx1_2, vx1_3, vx1_4, vx1_5, vx1_6, vx1_7, vx2_1, vx2_4, vx2_7, vx2_8, vx2_9])
    vUPows = np.array([vUPow1, vUPow2, vQokr])
    invC = np.array([[invC1,  0.0, 0.0],
                     [  0.0,invC2, 0.0]])
    powH = np.array([[H1_1, H1_2, 0.0, 0.0, H1_5, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0],
                     [ 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0, H2_1, 0.0, H2_7, 0.0, 0.0]])
    vT = np.array([vT1, vT2])
    heatStreams = np.array([ExtQEnPow1, ExtQEnPowOkr])
    Streams = np.array([xExt1_2, xExt1_6, xExt2_8])
    potBet = np.array([[0.0,       0.0,       0.0, 0.0, -muBet1_5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                       [0.0, -muBet1_2,       0.0, 0.0,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                       [0.0,       0.0, -muBet1_3, 0.0,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                       [0.0,       0.0,       0.0, 0.0,       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
    return (-chemPot,potBet,Temp,Betas,Aff,AffHeat,
            kineticMatrixCPCP,kineticMatrixCPHeat,
            kineticMatrixHeatCP,kineticMatrixHeatHeat,
            vProcesses,vQTransf,vHeatProcesses,
            vQPows,heatTransferMatrix,balanceMatrix,
            vx,vUPows,invC,powH,vT,heatStreams,Streams)
    