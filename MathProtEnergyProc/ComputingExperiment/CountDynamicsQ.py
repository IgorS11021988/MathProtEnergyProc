import numpy as np

from MathProtEnergyProc.RepeatLastRowsMatrix import RepeatLastRowsMatrix


# Класс расчета динамик
class CountDynamicsQ:
    # Инициализатор класса
    def __init__(self,
                 nonEqSystemDynQ,  # Класс динамики неравновесной системы
                 compressingFunction  # Функция сжатия
                 ):
        # Заполняем поля
        self.__NonEqSystemDynQ = nonEqSystemDynQ  # Класс динамики неравновесной системы
        self.__CompressingFunction = compressingFunction  # Функция сжатия

    # Расчет динамики со сжатием
    def ModelingDynamicQ(self, Tint,
                         stateCoordinates0,
                         reducedTemp0,
                         systemParameters,
                         t_eval=None,
                         index=None):
        # Вызываем моделирование динамики
        dyn = self.__NonEqSystemDynQ.NonEqSystemDynamicQ(Tint,
                                                         stateCoordinates0,
                                                         reducedTemp0,
                                                         systemParameters,
                                                         t_eval=t_eval)

        # Сжимаем результаты моделирования
        return self.__CompressingFunction(dyn, index)

    # Моделирование динамик (вычислительный эксперимент)
    def ComputingExperimentQ(self, Tints,
                             stateCoordinates0,
                             reducedTemp0,
                             systemParameters,
                             t_evals=None):
        # Приводим входные данные к массиву
        stateCoordinates0 = np.array(stateCoordinates0, dtype=np.double)  # Начальные координаты состояния
        reducedTemp0 = np.array(reducedTemp0, dtype=np.double)  # Начальные приведенные температуры
        systemParameters = np.array(systemParameters, dtype=np.double)  # Параметры системы

        # Получаем число динамик
        nStateCoordinates0 = stateCoordinates0.shape[0]  # Число совокупностей координат состояния
        nReducedTemp0 = reducedTemp0.shape[0]  # Число совокупностей приведенных температур
        nSystemParameters = systemParameters.shape[0]  # Число совокупностей параметров системы
        nDyn = np.max([nStateCoordinates0,
                       nReducedTemp0,
                       nSystemParameters])

        # Приводим числа строк в соответсвие
        Tints = RepeatLastRowsMatrix(Tints.reshape(-1, 1), nDyn)
        stateCoordinates0 = RepeatLastRowsMatrix(stateCoordinates0, nDyn)
        reducedTemp0 = RepeatLastRowsMatrix(reducedTemp0, nDyn)
        systemParameters = RepeatLastRowsMatrix(systemParameters, nDyn)

        # Векторизуем расчеты со сжатием
        if not (t_evals is None):
            def vModDyn(Tint, ind):
                return self.ModelingDynamicQ(Tint,
                                             stateCoordinates0[ind],
                                             reducedTemp0[ind],
                                             systemParameters[ind],
                                             t_eval=t_evals[ind],
                                             index=ind)
        else:
            def vModDyn(Tint, ind):
                return self.ModelingDynamicQ(Tint,
                                             stateCoordinates0[ind],
                                             reducedTemp0[ind],
                                             systemParameters[ind],
                                             index=ind)
        vecModelingDynamic = np.vectorize(vModDyn)
        return vecModelingDynamic(Tints, np.arange(nDyn).reshape(-1, 1))


# Класс векторизованной модели
class VectorModelQ(object):
    # Инициализатор
    def __init__(self,
                 systemParametersFunctor,  # Функтор параметров системы (и ее начального состояния)
                 countDynamicsQ  # Класс динамики неравновесной системы
                 ):
        # Заполняем поля
        self.__SystemParametersFunctor = systemParametersFunctor  # Функтор для расчета параметров системы
        self.__CountDynamicsQ = countDynamicsQ  # Класс динамики неравновесной системы

    # Получение динамики
    def ModelingDynamicsQ(self, Tints,
                          measParameters,
                          t_evals=None):
        # Получаем параметры системы и ее начальное состояние
        (stateCoordinates0,
         reducedTemp0,
         systemParameters) = self.__SystemParametersFunctor.SystemParametersCount(measParameters)

        # Получаем динамику системы
        return self.__CountDynamicsQ.ComputingExperimentQ(Tints,
                                                          stateCoordinates0,
                                                          reducedTemp0,
                                                          systemParameters,
                                                          t_evals=t_evals)
