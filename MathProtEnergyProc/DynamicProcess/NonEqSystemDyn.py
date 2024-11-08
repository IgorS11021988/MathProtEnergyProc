import numpy as np

#Класс динамики системы, моделируемой методом математического прототипирования энергетических процессов
class NonEqSystemDyn(object):
    #Инициализатор класса
    def __init__(self,
                 
                 nonEqSystem,#Система
                 
                 condFunction,#Функция условий протекания процессов
                 
                 extParametersFunction,#Функция внешних параметров
                 
                 integDynamic#Метод интегрирования дифференциальных уравнений
                 ):
        #Система
        self.__NEqSystem = nonEqSystem
        
        #Функция условий протекания процессов
        self.__CondFunction = condFunction
        
        #Функция внешних параметров
        self.__ExtParametersFunction = extParametersFunction
        
        #Метод интегрирования дифференциальных уравнений
        self.__integDynamic = integDynamic
    
    #Функция правой части
    def fNonEqSystem(self, t,
                     stateCoordinates,
                     systemParameters):
        #Рассчитываем условия протекания процессов
        systemParameters = self.__CondFunction(t,systemParameters)
        
        #Рассчитываем состояние системы
        self.__NEqSystem.CountSystem(stateCoordinates,
                                     systemParameters)
        
        #Скорость протекания процессов
        return self.__NEqSystem.GetVStateCoordinates().copy()
        
    #Функция расчета внешних параметров
    def fExtParameters(self, t,
                       stateCoordinates,
                       systemParameters):
        #Вызов функции расчета внешних параметров
        return self.__ExtParametersFunction(t,
                                            stateCoordinates,
                                            systemParameters)

    #Функция расчета динамики состояния системы
    def NonEqSystemDynamic(self,Tint,
                           stateCoordinates0,
                           systemParameters,
                           t_eval = None):
        #Получаем динамику коордиат состояния
        (t,stateCoordinates) = self.__integDynamic(self.fNonEqSystem,#Функция правой части динамики
                                                   Tint,#Время интегрирования
                                                   np.array(stateCoordinates0, dtype=np.double),#Начальные координаты
                                                   systemParameters,#Параметры системы
                                                   dynPoints=t_eval#Точки динамики
                                                   )
        
        #Выводим результат
        return self.fExtParameters(t,
                                   stateCoordinates,
                                   systemParameters)

#Класс модели
class Model(object):
    #Инициализатор
    def __init__(self,
                 systemParametersFunctor,#Функтор параметров системы (и ее начального состояния)
                 nonEqSystemDyn#Класс динамики неравновесной системы
                 ):
        #Заполняем поля
        self.__SystemParametersFunctor = systemParametersFunctor#Функтор для расчета параметров системы
        self.__NonEqSystemDyn = nonEqSystemDyn#Класс динамики неравновесной системы
    
    #Получение динамики
    def CountDynamic(self,Tint,
                     measParameters,
                     t_eval = None):
        #Получаем параметры системы и ее начальное состояние
        (stateCoordinates0,
         systemParameters) = self.__SystemParametersFunctor.SystemParametersCount(measParameters)
        
        #Получаем динамику системы
        return self.__NonEqSystemDyn.NonEqSystemDynamic(Tint,
                                                        stateCoordinates0,
                                                        systemParameters,
                                                        t_eval = t_eval)
        