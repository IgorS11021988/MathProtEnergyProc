import numpy as np

#Конкатенуем индексы с повторением
def ConcatRepeatIndexes(listIndexesDiapasons,#Список диапазонов индексов
                        selectIndexes = None#Индексы окончательного выбора 
                        ):
    #Число диапазонов индексов
    nRanges = len(listIndexesDiapasons)
    
    #Формируем список диапазонов индексов
    listArraysIndexesDiapasons = []
    for indexDiapason in listIndexesDiapasons:
        #Получаем минимум, максимум и шаг диапазона
        minValue = indexDiapason[0]#Минимум
        maxValue = indexDiapason[1]#Максимум
        stepValue = indexDiapason[2]#Шаг
        
        #Формируем и добавляем диапазон в список
        listArraysIndexesDiapasons.append(range(minValue,
                                                maxValue,
                                                stepValue))

    #Конкатенуем списки индексов с повторением
    rez = np.mgrid[listArraysIndexesDiapasons].reshape(nRanges,-1).transpose()
    
    #Выбираем при необходимоти индексы
    if not selectIndexes is None:
        return rez[:,selectIndexes]
    else:
        return rez
    
#Индексированные имена
def IndexedNamesFromIndexes(indexes,#Индексы
                            beginName,#Начало имени
                            endName = "",#Конец имени
                            sepName = "_"#Раздлитель имени
                            ):
    #Переводим индексы в строку
    aIndexes = np.array(indexes).astype("str").astype(object)
    
    #Проверяем число столбцов матрицы индексов
    if len(aIndexes.shape) > 1:
        #Конкатенуем строки матрицы индексов
        if aIndexes.shape[1] > 1:
            aIndexes = aIndexes[:,0] + np.sum(sepName + aIndexes[:,1::], axis=1)
        else:
            aIndexes.shape = (-1,)#Приводим к матрице-строке
    
    #Выводим имя
    return (beginName + aIndexes + endName).astype("str").tolist()
def IndexedNames(listIndexesDiapasons,#Список диапазонов индексов
                 beginName,#Начало имени
                 endName = "",#Конец имени
                 sepName = "_",#Раздлитель имени
                 selectIndexes = None#Индексы окончательного выбора 
                 ):
    #Получаем индексы
    Indexes = ConcatRepeatIndexes(listIndexesDiapasons,#Список диапазонов индексов
                                  selectIndexes = selectIndexes#Индексы окончательного выбора 
                                  )
    
    #Получаем индексированные имена
    return IndexedNamesFromIndexes(Indexes,#Индексы
                                   beginName,#Начало имени
                                   endName = endName,#Конец имени
                                   sepName = sepName#Раздлитель имени
                                   )
def AllIndexedNamesFromIndexes(indexes,#Индексы
                               partsNames#Части имен
                               ):
    #Формируем массив индексов
    aIndexes = np.array(indexes)
    
    #Формируем список индексированных имен
    listIndexedNames = []
    for partName in partsNames:
        #Получаем части имени и выбираемые индексы
        beginName = partName[0]#Начало имени
        endName = partName[1]#Конец имени
        sepName = partName[2]#Раздлитель имени
        selectIndexes = partName[3]#Индексы окончательного выбора
        
        #Формируем индексированные имена и добавляем в список индексированных имен
        listIndexedNames.append(IndexedNamesFromIndexes(aIndexes[:,selectIndexes],#Индексы
                                                        beginName,#Начало имени
                                                        endName = endName,#Конец имени
                                                        sepName = sepName#Раздлитель имени
                                                        ))

    #Выводим результат
    return listIndexedNames
def AllIndexedNames(listIndexesDiapasons,#Список диапазонов индексов
                    partsNames#Части имен
                    ):
    #Получаем индексы
    Indexes = ConcatRepeatIndexes(listIndexesDiapasons#Список диапазонов индексов
                                  )
    
    #Формируем имена
    return AllIndexedNamesFromIndexes(Indexes,#Индексы
                                      partsNames#Части имен
                                      )
