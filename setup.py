import subprocess as sp

#Устанавливаем библиотеки
command = ["conda","install","numpy","scipy"]
process = sp.Popen(command)
process.communicate()

#Устанавливаем библиотеку
command = ["pip","install","MathProtEnergyProc-1.0.tar.gz"]
sp.Popen(command)
