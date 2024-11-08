#Функция состояния для литий-ионного аккумулятора
def FunCharQ1(t,#Моменты времени
              stateCoordinates,#Координаты состояния
              reducedTemp,#Приведенные температуры
              systemParameters#Параметры системы  
              ):
    #Получаем координаты состояния
    Char1 = 0.33*stateCoordinates[:,0] + 0.63*stateCoordinates[:,1] + 0.93*stateCoordinates[:,4]
    Char2 = 0.21*stateCoordinates[:,2] + 0.99*stateCoordinates[:,1] + 3.27*stateCoordinates[:,5]
    
    #Учитываем приведенные температуры
    Char1 += 0.81*reducedTemp[:,0] ** 2
    Char2 += 0.51*reducedTemp[:,1] ** 3
    
    #Выводим результат
    return (t,Char1,Char2)
    