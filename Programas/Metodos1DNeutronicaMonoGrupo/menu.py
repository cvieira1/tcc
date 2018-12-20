import subprocess as sp

exe1 = sp.Popen(['mpicxx','-o','main','main.cpp'])
print("Escolha de que forma sera executado o problema:")
print("1 - Serial")
print("2 - Paralelo")
escolha = input()
if(escolha == '1'):
  exe = sp.Popen(['mpirun','-n','1','./main'])
if(escolha == '2'):
  exe = sp.Popen(['mpirun','-n','2','./main'])
    
