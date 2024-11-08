#Функция условий протекания процессов
def fU1(t,#Моменты времени
        systemParameters#Параметры системы  
        ):
    #Выводим результат
    return systemParameters*(1 + 0.1*t)
    