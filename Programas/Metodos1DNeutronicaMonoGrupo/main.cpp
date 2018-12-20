#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include "Bibliotecas/Quadratura1D.h"
#include "Bibliotecas/DadosArquivo.h"
#include "Bibliotecas/Metodos1D.h"
#include "mpi.h"

#include <vector>

#include <boost/filesystem.hpp>

#if defined(_WIN32)
  #include <windows.h>
#elif defined(__linux__)
  #include <sstream>
  #include <unistd.h>
#elif defined(__APPLE__)
  #include <mach-o/dyld.h>
#endif

string find_executable(){
  unsigned int bufferSize = 512;
  std::vector<char> buffer(bufferSize + 1);

#if defined(_WIN32)
  ::GetModuleFileName(NULL, &buffer[0], bufferSize);

#elif defined(__linux__)
  // Get the process ID.
  int pid = getpid();

  // Construct a path to the symbolic link pointing to the process executable.
  // This is at /proc/<pid>/exe on Linux systems (we hope).
  std::ostringstream oss;
  oss << "/proc/" << pid << "/exe";
  std::string link = oss.str();

  // Read the contents of the link.
  int count = readlink(link.c_str(), &buffer[0], bufferSize);
  if(count == -1) throw std::runtime_error("Could not read symbolic link");
  buffer[count] = '\0';

#elif defined(__APPLE__)
  if(_NSGetExecutablePath(&buffer[0], &bufferSize))
  {
    buffer.resize(bufferSize);
    _NSGetExecutablePath(&buffer[0], &bufferSize);
  }

#else
  #error Cannot yet find the executable on this platform
#endif

  std::string s = &buffer[0];
  return s;
}

using namespace std;

int main(int argc,char *argv[])
{
    int mynode,totalnodes;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&totalnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
    string programPath = find_executable();
    programPath = programPath.substr(0,programPath.find("./main",0));
    Dados_Entrada EntradaCaio;
    Metodos1D metodos;
    lerDadosEntrada(programPath + "dadosEntrada.txt",EntradaCaio);
    //exibeDadosEntrada(EntradaCaio);
    metodos.MetodoDD(EntradaCaio,programPath + "dadosSaida.txt",mynode,totalnodes);
    return 0;
}
