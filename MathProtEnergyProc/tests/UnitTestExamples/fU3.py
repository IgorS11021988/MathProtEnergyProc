import numpy as np

#Функция условий протекания процессов
def fU3(t,#Моменты времени
        systemParameters#Параметры системы  
        ):
    #Выводим результат
    return systemParameters*np.cos(1 + 0.1*t) + 0.3
    