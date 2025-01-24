#Функция состояния для литий-ионного аккумулятора
def FunChar1(t,#Моменты времени
             stateCoordinates,#Координаты состояния
             systemParameters#Параметры системы  
             ):
    #Получаем координаты состояния
    Char1 = 0.33*stateCoordinates[:,0] + 0.63*stateCoordinates[:,1] + 0.93*stateCoordinates[:,4]
    Char2 = 0.21*stateCoordinates[:,2] + 0.99*stateCoordinates[:,1] + 3.27*stateCoordinates[:,5]
    
    #Выводим результат
    return (t,Char1.reshape(-1,1),Char2.reshape(-1,1))
    