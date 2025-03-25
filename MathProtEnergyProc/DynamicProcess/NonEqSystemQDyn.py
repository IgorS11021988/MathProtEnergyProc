import numpy as np

#Класс динамики системы, моделируемой методом математического прототипирования энергетических процессов
class NonEqSystemQDyn:
    #Инициализатор класса
    def __init__(self,
                 
                 nonEqSystemQ,#Система
                 
                 condFunction,#Функция условий протекания процессов
                 
                 extParametersFunction,#Функция внешних параметров
                 
                 integDynamic#Метод интегрирования дифференциальных уравнений
                 ):
        #Система
        self.__NEqSystemQ = nonEqSystemQ
        
        #Число координат состояния системы
        self.__NStateCoordinates = len(nonEqSystemQ.GetStateCoordinatesNames())
        
        #Функция условий протекания процессов
        self.__CondFunction = condFunction
        
        #Функция внешних параметров
        self.__ExtParametersFunction = extParametersFunction
        
        #Метод интегрирования дифференциальных уравнений
        self.__integDynamic = integDynamic
    
    #Функция правой части системы
    def fNonEqSystemQ(self, t,
                      stateCoordinates,
                      reducedTemp,
                      systemParameters):
        #Рассчитываем условия протекания процессов
        systemParameters = self.__CondFunction(t,systemParameters)
        
        #Рассчитываем состояние системы
        self.__NEqSystemQ.CountSystem(stateCoordinates,
                                      reducedTemp,
                                      systemParameters)
        
        #Скорость изменения координат состояния
        return np.hstack((self.__NEqSystemQ.GetVStateCoordinates(),
                          self.__NEqSystemQ.GetVReducedTemperaturesEnergyPowers())).copy()
    
    #Функция правой части динамики системы
    def fSystemDyn(self, t,
                   allStParameters,
                   systemParameters):
        #Выделяем координаты состояния и приведенные температуры
        (stateCoordinates,reducedTemp) = self.__funSplitAllStParameters(allStParameters)
        
        #Скорость изменения координат состояния
        return self.fNonEqSystemQ(t,
                                  stateCoordinates,
                                  reducedTemp,
                                  systemParameters)
        
    #Функция расчета внешних параметров
    def fExtParameters(self, t,
                       stateCoordinates,
                       reducedTemp,
                       systemParameters):
        #Вызов функции расчета внешних параметров
        return self.__ExtParametersFunction(t,
                                            stateCoordinates,
                                            reducedTemp,
                                            systemParameters)
    
    #Функция расчета внешних параметров
    def fExtParametersByAllParSt(self, t,
                                 allStParameters,
                                 systemParameters):
        #Выделяем координаты состояния и внешние температуры
        (stateCoordinates,reducedTemp) = self.__funSplitAllStParameters(allStParameters.transpose())
        
        #Вызов расчета внешних параметров
        return self.fExtParameters(t,
                                   stateCoordinates.transpose(),
                                   reducedTemp.transpose(),
                                   systemParameters)

    #Функция расчета динамики состояния системы
    def NonEqSystemDynamicQ(self,Tint,
                            stateCoordinates0,
                            reducedTemp0,
                            systemParameters,
                            t_eval = None):
        #Получаем динамику коордиат состояния
        (t,allStParameters) = self.__integDynamic(self.fSystemDyn,#Функция правой части динамики
                                                  Tint,#Время интегрирования
                                                  np.hstack((stateCoordinates0,reducedTemp0)),#Начальные координаты
                                                  systemParameters,#Параметры системы
                                                  dynPoints=t_eval#Точки динамики
                                                  )
        
        #Выводим результат
        return self.fExtParametersByAllParSt(t,
                                             allStParameters,
                                             systemParameters)
        
    #Функция выделения координат состояния и приведенных температур
    def __funSplitAllStParameters(self,
                                  allStParameters):
        #Выделяем координаты состояния и внешние температуры
        stateCoordinates = allStParameters[0:self.__NStateCoordinates]#Координаты состояния
        reducedTemp = allStParameters[self.__NStateCoordinates::]#Приведенные температуры энергетических степеней свободы
        
        #Выводим результат
        return (stateCoordinates,reducedTemp)

#Класс модели
class ModelQ(object):
    #Инициализатор
    def __init__(self,
                 systemParametersFunctor,#Функтор параметров системы (и ее начального состояния)
                 nonEqSystemDynQ#Класс динамики неравновесной системы
                 ):
        #Заполняем поля
        self.__SystemParametersFunctor = systemParametersFunctor#Функтор для расчета параметров системы
        self.__NonEqSystemDynQ = nonEqSystemDynQ#Класс динамики неравновесной системы
    
    #Получение динамики
    def CountDynamic(self,Tint,
                     measParameters,
                     t_eval = None):
        #Получаем параметры системы и ее начальное состояние
        (stateCoordinates0,
         reducedTemp0,
         systemParameters) = self.__SystemParametersFunctor.SystemParametersCount(measParameters)
        
        #Получаем динамику системы
        return self.__NonEqSystemDynQ.NonEqSystemDynamicQ(Tint,
                                                          stateCoordinates0,
                                                          reducedTemp0,
                                                          systemParameters,
                                                          t_eval = t_eval)
        