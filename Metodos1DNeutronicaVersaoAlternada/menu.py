import subprocess as sp

file = open("dadosEntradaMenezes/problemaModelo1.txt","r")
comentario = file.readline()
escolha = file.read(1)
if(escolha == '1'):
  sp.call(['g++','main.cpp'])
if(escolha == '2'):
  sp.call(['g++','-Xpreprocessor','-fopenmp','main.cpp','-lomp'])
tmp = sp.call("./a.out");
