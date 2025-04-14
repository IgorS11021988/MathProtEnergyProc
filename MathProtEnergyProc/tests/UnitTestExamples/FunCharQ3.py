# Функция состояния для литий-ионного аккумулятора
def FunCharQ3(t,  # Моменты времени
              stateCoordinates,  # Координаты состояния
              reducedTemp,  # Приведенные температуры
              systemParameters  # Параметры системы
              ):
    # Получаем координаты состояния
    Char1 = 0.33 * stateCoordinates[:, 0]    +  0.63 * stateCoordinates[:, 1]    + 0.93 * stateCoordinates[:, 4]
    Char2 = 0.21 * stateCoordinates[:, 2]    +  0.99 * stateCoordinates[:, 0]**3 + 3.27 * stateCoordinates[:, 5]
    Char3 = 0.27 * stateCoordinates[:, 3]**2 + 0.963 * stateCoordinates[:, 4]    + 6.21 * stateCoordinates[:, 2]

    # Учитываем температуру
    Char1 += 3.15
    Char2 += 6.45 * reducedTemp[:, 0]
    Char3 += 9.75 * reducedTemp[:, 0] ** 2

    # Выводим результат
    return (t, Char1.reshape(-1, 1), Char2.reshape(-1, 1), Char3.reshape(-1, 1))
