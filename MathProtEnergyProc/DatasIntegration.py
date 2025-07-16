import numpy as np


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
        # Конкатенуем по горизонтали размноженные матрицы
        rez = np.hstack(rez)

        # Выбираем индексы (при необходимости)
        if isNeedConcatValuesIndexes:
            if ConcatValues:
                # Формируем индексы
                _inds = []
                for indexes in concatValuesIndexesList:
                    _inds += indexes
                
                # Выбираем подматрицу в соответствие с индексами
                rez = rez[:, _inds]
            else:
                # Формируем массив выбранных подматриц
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
