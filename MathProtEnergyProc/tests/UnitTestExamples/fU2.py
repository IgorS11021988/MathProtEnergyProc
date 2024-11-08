import numpy as np

#Функция условий протекания процессов
def fU2(t,#Моменты времени
        systemParameters#Параметры системы  
        ):
    #Выводим результат
    return systemParameters*np.sin(1 + 0.1*t) + 0.3
    