import numpy as np

from MathProtEnergyProc.RepeatLastRowsMatrix import RepeatLastRowsMatrix

#Класс расчета динамик
class CountDynamics:
    #Инициализатор класса
    def __init__(self,
                 nonEqSystemDyn,#Класс динамики неравновесной системы
                 compressingFunction#Функция сжатия
                 ):
        #Заполняем поля
        self.__NonEqSystemDyn = nonEqSystemDyn#Класс динамики неравновесной системы
        self.__CompressingFunction = compressingFunction#Функция сжатия
    
    #Расчет динамики со сжатием
    def ModelingDynamic(self, Tint,
                        stateCoordinates0,
                        systemParameters,
                        t_eval = None,
                        index = None):
        #Вызываем моделирование динамики
        dyn = self.__NonEqSystemDyn.NonEqSystemDynamic(Tint,
                                                       stateCoordinates0,
                                                       systemParameters,
                                                       t_eval = t_eval)
        
        #Сжимаем результаты моделирования
        return self.__CompressingFunction(dyn, index)

    #Моделирование динамик (вычислительный эксперимент)
    def ComputingExperiment(self, Tints,
                            stateCoordinates0,
                            systemParameters,
                            t_evals = None):
        #Приводим входные данные к массиву
        stateCoordinates0 = np.array(stateCoordinates0, dtype=np.double)#Начальные координаты состояния
        systemParameters = np.array(systemParameters, dtype=np.double)#Параметры системы
        
        #Получаем число динамик
        nStateCoordinates0 = stateCoordinates0.shape[0]#Число совокупностей координат состояния
        nSystemParameters = systemParameters.shape[0]#Число совокупностей параметров системы
        nDyn = np.max([nStateCoordinates0,
                       nSystemParameters])
        
        #Приводим числа строк в соответсвие
        Tints = RepeatLastRowsMatrix(Tints.reshape(-1,1),nDyn)
        stateCoordinates0 = RepeatLastRowsMatrix(stateCoordinates0,nDyn)
        systemParameters = RepeatLastRowsMatrix(systemParameters,nDyn)
        
        #Векторизуем расчеты со сжатием
        if not (t_evals is None):
            def vModDyn(Tint, ind):
                return self.ModelingDynamic(Tint,
                                            stateCoordinates0[ind],
                                            systemParameters[ind],
                                            t_eval = t_evals[ind],
                                            index = ind)
        else:
            def vModDyn(Tint, ind):
                return self.ModelingDynamic(Tint,
                                            stateCoordinates0[ind],
                                            systemParameters[ind],
                                            index = ind)
        vecModelingDynamic = np.vectorize(vModDyn)
        return vecModelingDynamic(Tints,np.arange(nDyn).reshape(-1,1))

#Класс векторизованной модели
class VectorModel(object):
    #Инициализатор
    def __init__(self,
                 systemParametersFunctor,#Функтор параметров системы (и ее начального состояния)
                 countDynamics#Класс параллельного расчета динамик
                 ):
        #Заполняем поля
        self.__SystemParametersFunctor = systemParametersFunctor#Функтор для расчета параметров системы
        self.__CountDynamics = countDynamics#Класс параллельного расчета динамик
    
    #Получение динамики
    def ModelingDynamics(self,Tints,
                         measParameters,
                         t_evals = None):
        #Получаем параметры системы и ее начальное состояние
        (stateCoordinates0,
         systemParameters) = self.__SystemParametersFunctor.SystemParametersCount(measParameters)
        
        #Получаем динамику системы
        return self.__CountDynamics.ComputingExperiment(Tints,
                                                        stateCoordinates0,
                                                        systemParameters,
                                                        t_evals = t_evals)        

#Класс обучения модели
class ModelLearning(object):
    #Инициализация модели
    def __init__(self,
                 modelLearningFunctor,#Функтор обучения модели параметров системы
                 countDynamics#Класс параллельного расчета динамик
                 ):
        #Заполняем поля
        self.__ModelLearningFunctor = modelLearningFunctor#Функтор обучения модели параметров системы
        self.__CountDynamics = countDynamics#Класс параллельного расчета динамик
    
    #Обучаем модель
    def ModelLearning(self,cTints,
                      cStateCoordinates0,
                      cSystemParameters,
                      ct_evals = None):
        #Выполняем вычислительный эксперимент
        comExpRez = self.__CountDynamics.ComputingExperiment(cTints,
                                                             cStateCoordinates0,
                                                             cSystemParameters,
                                                             t_evals = ct_evals)
    
        #Обучаем модель
        self.__ModelLearningFunctor.ModelSystemParametersCount(comExpRez,
                                                               cStateCoordinates0,
                                                               cSystemParameters)
    