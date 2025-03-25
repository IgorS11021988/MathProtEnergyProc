import numpy as np
from scipy import linalg as lg

#Функция получения потенциалов взаимодействия по их независимым состааляющим
def IntPotentialFromCompFun(intPotentialCubMatrix,#Кубическая матрица потенциалов взаимодействия
                            intPotentialComp#Независимые составляющие потенциалов взаимодействия
                            ):
    #Определяем и выводим потенциалы взаимодействия
    return np.sum(intPotentialCubMatrix*intPotentialComp.reshape(-1,1,1), 
                  axis=0)
