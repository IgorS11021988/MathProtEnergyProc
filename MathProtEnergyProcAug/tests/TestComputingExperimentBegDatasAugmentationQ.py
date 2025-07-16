import numpy as np

from MathProtEnergyProc.DatasAugmentation import ComputingExperimentBegDatasAugmentationQ

import unittest


# Модульные тесты
class TestComputingExperimentBegDatasAugmentationQ(unittest.TestCase):
    def setUp(self):
        # Выполнить настройку тестов (если необходимо)
        pass

    def tearDown(self):
        # Выполнить завершающие действия (если необходимо)
        pass

    # Модульные тесты
    def testComputingExperimentBegDatasAugmentationQ1(self):
        # Исходные данные
        settingParametersValues = [[1.1, 2.2, 5.5, 7.7],
                                   [6.3, 2.3, 4.5, 5.5],
                                   [2.7, 2.1, 3.3, 9.3],
                                   [1.1, 5.4, 2.8, 1.1],
                                   [7.8, 6.9, 1.2, 2.3]]
        randomParametersValues = [[4.1, 3.2, 6.5],
                                  [8.3, 4.3, 8.5],
                                  [6.7, 5.1, 6.3]]
        parametersState0Indexes = [2, 3, 5]
        reducedTemperatures0Indexes = [0]
        systemParametersIndexes = [1, 4]

        # Преобразованные матрицы
        matrValuesRezEt = np.array([[1.1, 2.2, 5.5, 7.7, 4.1, 3.2, 6.5],
                                    [1.1, 2.2, 5.5, 7.7, 8.3, 4.3, 8.5],
                                    [1.1, 2.2, 5.5, 7.7, 6.7, 5.1, 6.3],
                                    [6.3, 2.3, 4.5, 5.5, 4.1, 3.2, 6.5],
                                    [6.3, 2.3, 4.5, 5.5, 8.3, 4.3, 8.5],
                                    [6.3, 2.3, 4.5, 5.5, 6.7, 5.1, 6.3],
                                    [2.7, 2.1, 3.3, 9.3, 4.1, 3.2, 6.5],
                                    [2.7, 2.1, 3.3, 9.3, 8.3, 4.3, 8.5],
                                    [2.7, 2.1, 3.3, 9.3, 6.7, 5.1, 6.3],
                                    [1.1, 5.4, 2.8, 1.1, 4.1, 3.2, 6.5],
                                    [1.1, 5.4, 2.8, 1.1, 8.3, 4.3, 8.5],
                                    [1.1, 5.4, 2.8, 1.1, 6.7, 5.1, 6.3],
                                    [7.8, 6.9, 1.2, 2.3, 4.1, 3.2, 6.5],
                                    [7.8, 6.9, 1.2, 2.3, 8.3, 4.3, 8.5],
                                    [7.8, 6.9, 1.2, 2.3, 6.7, 5.1, 6.3]])

        # Генерируем случайные данные
        matrValuesRez = ComputingExperimentBegDatasAugmentationQ(settingParametersValues,  # Значения задаваемых параметров системы
                                                                 randomParametersValues,  # Значения случайно сгенерированных параметров системы
                                                                 parametersState0Indexes,  # Индексы начального состояния системы
                                                                 reducedTemperatures0Indexes,  # Индексы начальных приведенных температур системы
                                                                 systemParametersIndexes  # Индексы параметров системы
                                                                 )

        # Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez[0] - matrValuesRezEt[:, parametersState0Indexes]))
        self.assertEqual(dMatrValuesRez, 0.0)
        dMatrValuesRez = np.max(np.abs(matrValuesRez[1] - matrValuesRezEt[:, reducedTemperatures0Indexes]))
        self.assertEqual(dMatrValuesRez, 0.0)
        dMatrValuesRez = np.max(np.abs(matrValuesRez[2] - matrValuesRezEt[:, systemParametersIndexes]))
        self.assertEqual(dMatrValuesRez, 0.0)

    def testComputingExperimentBegDatasAugmentationQ2(self):
        # Исходные данные
        settingParametersValues = [[1.1, 1.2, 4.5],
                                   [5.3, 1.3, 3.5],
                                   [1.7, 1.1, 1.3],
                                   [2.1, 3.4, 1.8],
                                   [4.8, 2.9, 0.2],
                                   [3.7, 2.5, 6.2]]
        randomParametersValues = [[3.1, 3.2, 6.5, 4.5, 5.6],
                                  [8.3, 2.3, 8.5, 1.1, 7.7]]
        parametersState0Indexes = [0, 2]
        reducedTemperatures0Indexes = [3, 6]
        systemParametersIndexes = [1, 5, 4]

        # Преобразованные матрицы
        matrValuesRezEt = np.array([[1.1, 1.2, 4.5, 3.1, 3.2, 6.5, 4.5, 5.6],
                                    [1.1, 1.2, 4.5, 8.3, 2.3, 8.5, 1.1, 7.7],
                                    [5.3, 1.3, 3.5, 3.1, 3.2, 6.5, 4.5, 5.6],
                                    [5.3, 1.3, 3.5, 8.3, 2.3, 8.5, 1.1, 7.7],
                                    [1.7, 1.1, 1.3, 3.1, 3.2, 6.5, 4.5, 5.6],
                                    [1.7, 1.1, 1.3, 8.3, 2.3, 8.5, 1.1, 7.7],
                                    [2.1, 3.4, 1.8, 3.1, 3.2, 6.5, 4.5, 5.6],
                                    [2.1, 3.4, 1.8, 8.3, 2.3, 8.5, 1.1, 7.7],
                                    [4.8, 2.9, 0.2, 3.1, 3.2, 6.5, 4.5, 5.6],
                                    [4.8, 2.9, 0.2, 8.3, 2.3, 8.5, 1.1, 7.7],
                                    [3.7, 2.5, 6.2, 3.1, 3.2, 6.5, 4.5, 5.6],
                                    [3.7, 2.5, 6.2, 8.3, 2.3, 8.5, 1.1, 7.7]])

        # Генерируем случайные данные
        matrValuesRez = ComputingExperimentBegDatasAugmentationQ(settingParametersValues,  # Значения задаваемых параметров системы
                                                                 randomParametersValues,  # Значения случайно сгенерированных параметров системы
                                                                 parametersState0Indexes,  # Индексы начального состояния системы
                                                                 reducedTemperatures0Indexes,  # Индексы начальных приведенных температур системы
                                                                 systemParametersIndexes  # Индексы параметров системы
                                                                 )

        # Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez[0] - matrValuesRezEt[:, parametersState0Indexes]))
        self.assertEqual(dMatrValuesRez, 0.0)
        dMatrValuesRez = np.max(np.abs(matrValuesRez[1] - matrValuesRezEt[:, reducedTemperatures0Indexes]))
        self.assertEqual(dMatrValuesRez, 0.0)
        dMatrValuesRez = np.max(np.abs(matrValuesRez[2] - matrValuesRezEt[:, systemParametersIndexes]))
        self.assertEqual(dMatrValuesRez, 0.0)

    def testComputingExperimentBegDatasAugmentationQ3(self):
        # Исходные данные
        settingParametersValues = [[3.1, 2.2, 5.5, 7.7, 10.1],
                                   [5.1, 2.3, 4.5, 4.5, 11.3],
                                   [3.5, 2.1, 3.3, 6.3, 13.5]]
        randomParametersValues = [[4.7, 3.2],
                                  [9.3, 4.3],
                                  [7.7, 5.1],
                                  [4.5, 3.4]]
        parametersState0Indexes = [3]
        reducedTemperatures0Indexes = [0, 6, 2]
        systemParametersIndexes = [1, 4, 5]

        # Преобразованные матрицы
        matrValuesRezEt = np.array([[3.1, 2.2, 5.5, 7.7, 10.1, 4.7, 3.2],
                                    [3.1, 2.2, 5.5, 7.7, 10.1, 9.3, 4.3],
                                    [3.1, 2.2, 5.5, 7.7, 10.1, 7.7, 5.1],
                                    [3.1, 2.2, 5.5, 7.7, 10.1, 4.5, 3.4],
                                    [5.1, 2.3, 4.5, 4.5, 11.3, 4.7, 3.2],
                                    [5.1, 2.3, 4.5, 4.5, 11.3, 9.3, 4.3],
                                    [5.1, 2.3, 4.5, 4.5, 11.3, 7.7, 5.1],
                                    [5.1, 2.3, 4.5, 4.5, 11.3, 4.5, 3.4],
                                    [3.5, 2.1, 3.3, 6.3, 13.5, 4.7, 3.2],
                                    [3.5, 2.1, 3.3, 6.3, 13.5, 9.3, 4.3],
                                    [3.5, 2.1, 3.3, 6.3, 13.5, 7.7, 5.1],
                                    [3.5, 2.1, 3.3, 6.3, 13.5, 4.5, 3.4]])

        # Генерируем случайные данные
        matrValuesRez = ComputingExperimentBegDatasAugmentationQ(settingParametersValues,  # Значения задаваемых параметров системы
                                                                 randomParametersValues,  # Значения случайно сгенерированных параметров системы
                                                                 parametersState0Indexes,  # Индексы начального состояния системы
                                                                 reducedTemperatures0Indexes,  # Индексы начальных приведенных температур системы
                                                                 systemParametersIndexes  # Индексы параметров системы
                                                                 )

        # Проверяем значения
        dMatrValuesRez = np.max(np.abs(matrValuesRez[0] - matrValuesRezEt[:, parametersState0Indexes]))
        self.assertEqual(dMatrValuesRez, 0.0)
        dMatrValuesRez = np.max(np.abs(matrValuesRez[1] - matrValuesRezEt[:, reducedTemperatures0Indexes]))
        self.assertEqual(dMatrValuesRez, 0.0)
        dMatrValuesRez = np.max(np.abs(matrValuesRez[2] - matrValuesRezEt[:, systemParametersIndexes]))
        self.assertEqual(dMatrValuesRez, 0.0)


# Запустить тестирование
if __name__ == "__main__":
    unittest.main()
