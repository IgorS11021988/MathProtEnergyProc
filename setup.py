from distutils.core import setup

# Функция установки
setup(name="MathProtEnergyProc",
      version="1.0",
      author="Igor Starostin",
      author_email="starostinigo@yandex.ru",
      description="System modeling by mathematical prototiping method",
      packages=["MathProtEnergyProc",
                "MathProtEnergyProc.CorrectionModel",
                "MathProtEnergyProc.HeatPowerValues",
                "MathProtEnergyProc.ComputingExperiment",
                "MathProtEnergyProc.DynamicProcess",
                "MathProtEnergyProc.NonEqProcess",
                "MathProtEnergyProc.tests",
                "MathProtEnergyProc.tests.UnitTestExamples"],
      scripts=["testMathProtEnergyProc.py"]
      )
