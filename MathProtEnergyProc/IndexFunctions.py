import numpy as np

#Получение индексов в массиве имен
def GetIndex(arrNamesSource,#Массив имен, в котором ищем индекс
             name#Имя, для которого ищем индекс
             ):
    #Выводим индекс
    return arrNamesSource.index(name)

def GetIndexes(arrNamesSource,#Массив имен, в котором ищем индекс
               arrNames#Массив имен, для которых ищем индекс
               ):
    #Формируем массив индексов
    arrIndexes = []#Искомые индексы
    for name in arrNames:
        #Получаем индекс в массиве имен источников
        index = GetIndex(arrNamesSource,name)
        
        #Добавляем индекс в массив
        arrIndexes.append(index)
    
    #Выводим массив индексов
    return arrIndexes.copy()

def Get2Index(arrNamesSource1,#Массив 1 имен, в котором ищем индекс
              arrNamesSource2,#Массив 2 имен, в котором ищем индекс
              name#Имя, для которого ищем индекс
              ):
    #Выводим индекс в массиве имен источников
    index1 = GetIndex(arrNamesSource1,name)
    return GetIndex(arrNamesSource2,index1)

def Get2Indexes(arrNamesSource1,#Массив 1 имен, в котором ищем индекс
                arrNamesSource2,#Массив 2 имен, в котором ищем индекс
                arrNames#Массив имен, для которых ищем индекс
                ):
    #Формируем массив индексов
    arrIndexes = []#Искомые индексы
    for name in arrNames:
        #Получаем индекс в массиве имен источников
        index = Get2Index(arrNamesSource1,
                          arrNamesSource2,
                          name)
        
        #Добавляем индекс в массив
        arrIndexes.append(index)
    
    #Выводим массив индексов
    return arrIndexes.copy()
