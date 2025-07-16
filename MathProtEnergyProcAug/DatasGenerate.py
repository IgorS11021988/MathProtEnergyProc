import numpy as np


# Гененрируем локально-равномерно распределенные величины
def GenerateRandomdatasInDiapasons(minValues,  # Минимальные значения величин
                                   maxValues,  # Максимальные значения величин
                                   nPoints  # Числа точек в соответствующих диапазонах
                                   ):  # Конкатенация с размножением двух матриц
    # Матрицы максимальных и минимальных величин приводим к массиву numpy
    aMinValues = np.array(minValues)  # Минимальные значения величин
    aMaxValues = np.array(maxValues)  # Максимальные значения величин

    # Получаем матрицу приращений
    deltaValues = aMaxValues - aMinValues

    # Приодим минимальные величины и матрицу приращений в соответсвие с индексами
    aMinValues = np.repeat(aMinValues, nPoints, axis=0)  # Минимальные величины
    deltaValues = np.repeat(deltaValues, nPoints, axis=0)  # Приращения величин

    # Получаем случайные величины
    (nRows, nColumns) = aMinValues.shape  # Число рядов и колонок
    return aMinValues + deltaValues * np.random.rand(nRows, nColumns)  # Случайные числа


# Приведение одномерного массива к матрице-столбцу
def ToMatrixOneColumn(matr):
    if len(matr.shape) == 1:
        return matr.reshape(-1, 1)
    else:
        return matr


# Приведение списка к массиву
def ToArrayNumPy(List):
    return ToMatrixOneColumn(np.array(List))


# Горизонтальная конкатенация матриц с размножением
def HStackMatrixRepeatTwo(values1,  # Величины 1
                          values2,  # Величины 2
                          ConcatValues=True,  # Нужно ли конкатеновать величины
                          concatValuesIndexesList=None  # Список окончательных индексов конкатенованных величин
                          ):
    # Приводим величины к матрицам numpy
    aValues1 = ToArrayNumPy(values1)
    aValues2 = ToArrayNumPy(values2)

    # Получаем число строк матрицы 1
    nRows1 = aValues1.shape[0]

    # Получаем размерность матрицы 2
    (nRows2, nColumns2) = aValues2.shape

    # Размножаем матрицы
    rez = (np.repeat(aValues1, nRows2, axis=0),  # Матрица 1
           np.full((nRows1, aValues2.size), aValues2.reshape(1, -1)).reshape(-1, nColumns2)  # Матрица 2
           )

    # Конкатенуем матрицы и выбираем индексы (при необходимости)
    isNeedConcatValuesIndexes = concatValuesIndexesList is not None
    if ConcatValues or isNeedConcatValuesIndexes:
        # Конкатенуем индексы
        rez = np.hstack(rez)

        # Выбираем индексы (при необходимости)
        if isNeedConcatValuesIndexes:
            _rez = []
            for indexes in concatValuesIndexesList:
                _rez.append(rez[:, indexes])
            rez = _rez

    # Выводим результат
    return rez


# Конкатенация произвольного числа матриц с размножением
def HStackMatrixRepeat(listValues,  # Список массивов данных величин
                       repeatValues,  # Нужно ли повторять величины
                       concatValuesIndexes=None  # Окончательные индексы конкатенованных величин
                       ):
    # Приводим к массивам numpy аргументы
    aRepeatValues = np.array(repeatValues)  # Приводим необходимость повторения величин

    # Получаем индексы
    inds = range(aRepeatValues.size)

    # Конкатенуем с размножением
    rezMatrix = ToArrayNumPy(listValues[0])  # Добавляем первый элемент в результирующую матрицу
    for ind in inds:
        if aRepeatValues[ind]:
            rezMatrix = HStackMatrixRepeatTwo(rezMatrix,  # Величины 1
                                              listValues[ind + 1]  # Величины 2
                                              )
        else:
            rezMatrix = np.hstack((rezMatrix,
                                   ToArrayNumPy(listValues[ind + 1])))

    # Делаем выборку в соответсвие с индексами (при необхолдимости)
    if concatValuesIndexes is not None:
        aConcatValuesIndexes = np.array(concatValuesIndexes).tolist()  # Приводим окончательные индексы конкатенованных величин
        return rezMatrix[:, aConcatValuesIndexes]
    else:
        return rezMatrix


# Конкатенация заданных массивов величин
def HVStackMatrixRepeat(listsConcatedValues,  # Списки массивов конкатенуемых данных величин
                        concatValuesIndexes=None  # Окончательные индексы конкатенованных величин
                        ):
    # Конкатенуем по вертикали
    ConcatValues = []
    for conVerValues in listsConcatedValues:
        # Получаем задание на конкатенацию величин
        (listValues, repeatValues) = conVerValues

        # Конкатенуем
        ConcatValues += HStackMatrixRepeat(listValues, repeatValues).tolist()
    ConcatValues = np.array(ConcatValues)  # Переводим в массив numpy

    # Берем индексы (при необходимости)
    if concatValuesIndexes is not None:
        aConcatValuesIndexes = np.array(concatValuesIndexes).tolist()  # Приводим окончательные индексы конкатенованных величин
        return ConcatValues[:, aConcatValuesIndexes]
    else:
        return ConcatValues


# Аугментируем данные для вычислительного эксперимента
def ComputingExperimentBegDatasAugmentation(settingParametersValues,  # Значения задаваемых параметров системы
                                            randomParametersValues,  # Значения случайно сгенерированных параметров системы
                                            parametersState0Indexes,  # Индексы начального состояния системы
                                            systemParametersIndexes  # Индексы параметров системы
                                            ):
    # Конкатенуем данные и выделяем параметры системы и ее начальное состояние
    return HStackMatrixRepeatTwo(settingParametersValues,
                                 randomParametersValues,
                                 ConcatValues=True,
                                 concatValuesIndexesList=[parametersState0Indexes,
                                                          systemParametersIndexes]
                                 )


def ComputingExperimentBegDatasAugmentationQ(settingParametersValues,  # Значения задаваемых параметров системы
                                             randomParametersValues,  # Значения случайно сгенерированных параметров системы
                                             parametersState0Indexes,  # Индексы начального состояния системы
                                             reducedTemperatures0Indexes,  # Индексы начальных приведенных температур системы
                                             systemParametersIndexes  # Индексы параметров системы
                                             ):
    # Конкатенуем данные и выделяем параметры системы и ее начальное состояние
    return HStackMatrixRepeatTwo(settingParametersValues,
                                 randomParametersValues,
                                 ConcatValues=True,
                                 concatValuesIndexesList=[parametersState0Indexes,
                                                          reducedTemperatures0Indexes,
                                                          systemParametersIndexes]
                                 )
