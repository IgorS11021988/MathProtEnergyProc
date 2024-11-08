import numpy as np
from scipy import integrate as integ

#Функция интегрирования
def Integ(fDyn,#Функция правой части динамики
          Tint,#Время интегрирования
          Coordinates0,#Начальные координаты
          systemParameters,#Параметры системы
          dynPoints = None,#Точки динамики
             
          atol = 1e-6,#Абсолютная погрешность
          rtol = 1e-3,#Относительная погрешность
          
          max_step = np.inf,#Максимальный шаг
          
          method = "RK45"#Метод интегрирования дифференциальных уравнений
          ):
    #Получаем динамику коордиат состояния
    allStParameters = integ.solve_ivp(fDyn,
                                      (0.0, Tint),
                                      Coordinates0,#Начальное состояние
                                      t_eval=dynPoints,#Моменты времени
                                      atol=atol,#Абсолютная погрешность
                                      rtol=rtol,#Относительная погрешность
                                      max_step=max_step,#Максимальный шаг
                                      method=method,#Метод численного интегрирования
                                      args=(systemParameters,)
                                      )#Интегрирование дифференциальных уравнений
    t = allStParameters.t.reshape(-1,1)#Моменты времени
    allStParameters = allStParameters.y.transpose()#Координаты состояния
    
    #Выводим результат
    return (t,allStParameters)
    
#Интеграторы динамики системы
class standartIntegrateDyn(object):#Стандартный интегратор динамики
    #Инициализатор класса
    def __init__(self,
                 
                 atol = 1e-6,#Абсолютная погрешность
                 rtol = 1e-3,#Относительная погрешность
                 
                 max_step = np.inf,#Максимальный шаг
                 
                 method = "RK45"#Метод интегрирования дифференциальных уравнений
                 ):
        #Погрешность интегрирования
        self.__atol = atol#Абсолютная погрешность
        self.__rtol = rtol#Относительная погрешность
        
        #Максимальный шаг
        self.__maxStep = max_step
        
        #Метод интегрирования дифференциальных уравнений
        self.__method = method
    
    #Функция вызова
    def __call__(self,
                 fDyn,#Функция правой части динамики
                 Tint,#Время интегрирования
                 Coordinates0,#Начальные координаты
                 systemParameters,#Параметры системы
                 dynPoints = None#Точки динамики
                 ):
        #Выводим динамику 
        return Integ(fDyn,#Функция правой части динамики
                     Tint,#Время интегрирования
                     Coordinates0,#Начальные координаты
                     systemParameters,#Параметры системы
                     dynPoints = dynPoints,#Точки динамики
                        
                     atol = self.__atol,#Абсолютная погрешность
                     rtol = self.__rtol,#Относительная погрешность
                     
                     max_step = self.__maxStep,#Максимальный шаг
                     
                     method = self.__method#Метод интегрирования дифференциальных уравнений
                     )
                                          
class stepIntegrateDyn(object):#Шаговый интегратор динамики
    #Инициализатор класса
    def __init__(self,
                 
                 funIntegAttributes,#Функция аттрибутов интегрирования
                 
                 method = "RK45"#Метод интегрирования дифференциальных уравнений
                 ):
        #Функция аттрибутов интегрирования
        self.__funIntegAttributes = funIntegAttributes
        
        #Метод интегрирования дифференциальных уравнений
        self.__method = method

    #Функция вызова
    def __call__(self,
                 fDyn,#Функция правой части динамики
                 Tint,#Время интегрирования
                 Coordinates0,#Начальные координаты
                 systemParameters,#Параметры системы
                 dynPoints = None#Точки динамики
                 ):
        #Рассчитываем аттрибуты интегрирования
        (atol,rtol,max_step) = self.__funIntegAttributes(systemParameters,self.__method)
        
        #Выводим динамику 
        return Integ(fDyn,#Функция правой части динамики
                     Tint,#Время интегрирования
                     Coordinates0,#Начальные координаты
                     systemParameters,#Параметры системы
                     dynPoints = dynPoints,#Точки динамики
                        
                     atol = atol,#Абсолютная погрешность
                     rtol = rtol,#Относительная погрешность
                     
                     max_step = max_step,#Максимальный шаг
                     
                     method = self.__method#Метод интегрирования дифференциальных уравнений
                     )
