#ifndef DADOSARQUIVO_H_INCLUDED
#define DADOSARQUIVO_H_INCLUDED

#include <cmath>
#include "Quadratura1D.h"
#include "gnuplot.h"

using namespace std;

void exibeDadosEntrada(Dados_Entrada EntradaCaio){
    cout << endl;
    cout << "Ordem da Quadratura: " << EntradaCaio.ordemQuad << endl;
    cout << "Numero de Regioes: " << EntradaCaio.numRegioes << endl;
    cout << "Numero de Zonas: " << EntradaCaio.numZonas << endl;
    for(int i = 0;i < EntradaCaio.numRegioes;i++){
      cout << "Tamanho da Regiao " << i+1 << ": " << EntradaCaio.tamanhoRegiao[i] << endl;
    }
    for(int i = 0;i < EntradaCaio.numRegioes;i++){
      cout << "Nodos da Regiao " << i+1 << ": " << EntradaCaio.nodosRegiao[i] << endl;
    }
    cout << "Periodicidade: " << EntradaCaio.periodicidade << endl;
    for(int i = 0;i < EntradaCaio.numRegioes;i++){
      cout << "Mapeamento da Regiao " << i+1 << ": " << EntradaCaio.mapeamento[i] << endl;
    }
    cout << "Criterio de Parada: " << EntradaCaio.cp << endl;
    cout << "Numero de Grupos de Energia: " << EntradaCaio.numGrupos << endl;
    cout << "Grau da Anisotropia: " << EntradaCaio.grauAnisotropia << endl;
    for(int i = 0;i < EntradaCaio.numZonas;i++){
      cout << "Dados Zona " << (i + 1) << endl;
      cout << "Sigma Total: " << endl;
      for(int j = 0;j < EntradaCaio.numGrupos;j++){
        cout << EntradaCaio.sigmaTot[i][j] << " ";
      }
      cout << endl;
      for(int j = 0;j <= EntradaCaio.grauAnisotropia;j++){
        cout << "Sigma Espalhamento " << j << ": " << endl;
        for(int l = 0;l < EntradaCaio.numGrupos;l++){
          for(int m = 0;m < EntradaCaio.numGrupos;m++){
            cout << EntradaCaio.sigmaEsp[i][j][l][m] << " ";
          }
          cout << endl;
        }
        cout << endl;
      }
    }
    cout << "Tipo das Condicao de Contorno: " << EntradaCaio.tipoCc[0] << " " << EntradaCaio.tipoCc[1] << endl;
    for(int i = 0;i < EntradaCaio.numGrupos;i++){
      cout << "Grupo de Energia " << i+1 << endl;
      cout << "Valor das Condicao de Contorno: " << EntradaCaio.valorCc[i][0] << " " << EntradaCaio.valorCc[i][1] << endl;
    }
    for(int i = 0;i < EntradaCaio.numGrupos;i++){
      cout << "Grupo de Energia " << i+1 << endl;
      for(int j = 0;j < EntradaCaio.numRegioes;j++){
        cout << "Valor da Fonte Fisica na Regiao " << j+1 << ": " << EntradaCaio.fonteFisica[i][j] << endl;
      }
    }
    for(int i = 0;i < EntradaCaio.numGrupos;i++){
      cout << "Grupo de Energia " << i+1 << endl;
      for(int j = 0;j < EntradaCaio.numRegioes;j++){
        cout << "Valor da Fonte Adjunta na Regiao " << j+1 << ": " << EntradaCaio.fonteAdjunta[i][j] << endl;
      }
    }
    cout << "Mapeamento dos Grupos de Energia que tem Detector: " << endl;
    for(int g = 0;g < EntradaCaio.numDetector;g++){
        cout << EntradaCaio.mapeaDetector[g] << endl;
    }
    cout << "Valores de Mi:" << endl;
    for(int i = 0;i < EntradaCaio.ordemQuad;i++){
      cout << EntradaCaio.MI[i] << endl;
    }
    cout << "Valores dos Pesos da Quadratura:" << endl;
    for(int i = 0;i < EntradaCaio.ordemQuad;i++){
      cout << EntradaCaio.wn[i] << endl;
    }
    cout << "Numero total de nodos no dominio: " << EntradaCaio.numNodos << endl;
    for(int i = 0;i < EntradaCaio.numRegioes;i++){
      cout << "Tamanho do nodo na Regiao " << i+1 << ": " << EntradaCaio.tamanhoNodo[i] << endl;
    }
    cout << "Tamanho total do dominio: " << EntradaCaio.tamanhoDominio << endl;
}


void lerDadosEntrada(string caminho,Dados_Entrada &EntradaCaio){
    int ldado,exp,counter;
    string linha;
    ifstream arq;
    bool primeiro = true;
    bool grupoDetector = false; // variavel que verifica se ja tem um detector em um grupo
    counter = 0;
    EntradaCaio.tipoCc = new int[2];
    EntradaCaio.numDetector = 0;
    arq.open(caminho);
    ldado = 0;
    if(arq.is_open()){
      while(!arq.eof()){
        getline(arq,linha);
        if(linha.find("//") == string::npos){ //Ira entrar se a linha nao for 1 comentario
          ldado++;
          switch(ldado){
            case 2:
              EntradaCaio.tipoProblema = stoi(linha);
              break;
            case 3:
              EntradaCaio.ordemQuad = stoi(linha);
              break;
            case 4:
              EntradaCaio.numRegioes = stoi(linha.substr(0,linha.find(" ")));
              EntradaCaio.numZonas = stoi(linha.substr((linha.find(" ")+1)));
              EntradaCaio.tamanhoRegiao = new double[EntradaCaio.numRegioes];
              EntradaCaio.nodosRegiao = new int[EntradaCaio.numRegioes];
              EntradaCaio.mapeamento = new int[EntradaCaio.numRegioes];
              EntradaCaio.tamanhoNodo = new double[EntradaCaio.numRegioes];
              EntradaCaio.somaNodosRegiao = new int[EntradaCaio.numRegioes];
              EntradaCaio.numNodos = 0;
              EntradaCaio.somaNodosRegiao[0] = -1;
              EntradaCaio.tamanhoDominio = 0;
              break;
            case 5:
              for(int i = 0; i < EntradaCaio.numRegioes;i++){
                EntradaCaio.tamanhoRegiao[i] = stod(linha.substr(0,linha.find(" ")));
                EntradaCaio.tamanhoDominio += EntradaCaio.tamanhoRegiao[i];
                linha = linha.substr((linha.find(" ")+1));
              }
              break;
            case 6:
              for(int i = 0; i < EntradaCaio.numRegioes;i++){
                EntradaCaio.nodosRegiao[i] = stoi(linha.substr(0,linha.find(" ")));
                EntradaCaio.numNodos += EntradaCaio.nodosRegiao[i];
                if(i < (EntradaCaio.numRegioes - 1)){
                    EntradaCaio.somaNodosRegiao[i + 1] = EntradaCaio.numNodos - 1;
                }
                EntradaCaio.tamanhoNodo[i] = EntradaCaio.tamanhoRegiao[i] / EntradaCaio.nodosRegiao[i];
                linha = linha.substr((linha.find(" ")+1));
              }
              break;
            case 7:
              EntradaCaio.periodicidade = stoi(linha);
              break;
            case 8:
              for(int i = 0; i < EntradaCaio.numRegioes;i++){
                EntradaCaio.mapeamento[i] = stoi(linha.substr(0,linha.find(" ")));
                linha = linha.substr((linha.find(" ")+1));
              }
              break;
            case 9:
              exp = (-1)*stoi(linha);
              EntradaCaio.cp = pow(10,exp);
              break;
            case 10:
              EntradaCaio.numGrupos = stoi(linha);
              EntradaCaio.fonteFisica = new double *[EntradaCaio.numGrupos];
              EntradaCaio.mapeaFonteFisica = new double *[EntradaCaio.numGrupos];
              EntradaCaio.fonteAdjunta = new double *[EntradaCaio.numGrupos];
              EntradaCaio.valorCc = new double *[EntradaCaio.numGrupos];
              for(int i = 0;i < EntradaCaio.numGrupos;i++){
                EntradaCaio.fonteFisica[i] = new double[EntradaCaio.numRegioes];
                EntradaCaio.mapeaFonteFisica[i] = new double[EntradaCaio.numRegioes];
                EntradaCaio.fonteAdjunta[i] = new double[EntradaCaio.numRegioes];
                EntradaCaio.valorCc[i] = new double[2];
              }
              EntradaCaio.sigmaTot = new double *[EntradaCaio.numZonas];
              for(int i = 0;i < EntradaCaio.numZonas;i++){
                EntradaCaio.sigmaTot[i] = new double [EntradaCaio.numGrupos];
              }
              break;
            case 11:
              EntradaCaio.grauAnisotropia = stoi(linha);
              EntradaCaio.sigmaEsp = new double ***[EntradaCaio.numZonas];
              for(int i = 0;i < EntradaCaio.numZonas;i++){
                EntradaCaio.sigmaEsp[i] = new double **[EntradaCaio.grauAnisotropia + 1];
                for(int j = 0;j <= EntradaCaio.grauAnisotropia;j++){
                    EntradaCaio.sigmaEsp[i][j] = new double *[EntradaCaio.numGrupos];
                    for(int l = 0;l < EntradaCaio.numGrupos;l++){
                        EntradaCaio.sigmaEsp[i][j][l] = new double [EntradaCaio.numGrupos];
                    }
                }
              }
              break;
            case 12:
              if(EntradaCaio.tipoProblema == 1){
                for(int i = 0; i < EntradaCaio.numZonas;i++){
                  primeiro = true;
                  if(linha.find("//") == string::npos){
                    for(int j = 0; j < EntradaCaio.numGrupos; j++){
                      EntradaCaio.sigmaTot[i][j] = stod(linha.substr(0,linha.find(" ")));
                      linha = linha.substr((linha.find(" ")+1));
                    }
                    getline(arq,linha);
                    for(int j2 = 0; j2 <= EntradaCaio.grauAnisotropia; j2++){
                      if(linha.find("//") != string::npos){
                        if(j2 > 0){
                          j2--;
                        }
                      }else{
                        if(primeiro){
                          j2 = 0;
                          primeiro = false;
                        }
                        for(int l = 0; l < EntradaCaio.numGrupos; l++){
                          for(int m = 0; m < EntradaCaio.numGrupos; m++){
                            EntradaCaio.sigmaEsp[i][j2][l][m] = stod(linha.substr(0,linha.find(" ")));
                            linha = linha.substr((linha.find(" ")+1));
                          }
                          getline(arq,linha);
                        }
                      }
                      if(j2 < EntradaCaio.grauAnisotropia){
                        getline(arq,linha);
                      }
                      if(EntradaCaio.grauAnisotropia == 0 && primeiro){
                        j2--;
                        getline(arq,linha);
                      }
                    }
                  }else{
                    i--;
                  }
                  if(i < EntradaCaio.numZonas - 1){
                    getline(arq,linha);
                  }
                }
              }else{
                for(int i = 0; i < EntradaCaio.numZonas;i++){
                  primeiro = true;
                  if(linha.find("//") == string::npos){
                    for(int j = 0; j < EntradaCaio.numGrupos; j++){
                      EntradaCaio.sigmaTot[i][j] = stod(linha.substr(0,linha.find(" ")));
                      linha = linha.substr((linha.find(" ")+1));
                    }
                    getline(arq,linha);
                    for(int j2 = 0; j2 <= EntradaCaio.grauAnisotropia; j2++){
                      if(linha.find("//") != string::npos){
                        if(j2 > 0){
                          j2--;
                        }
                      }else{
                        if(primeiro){
                          j2 = 0;
                          primeiro = false;
                        }
                        for(int l = 0; l < EntradaCaio.numGrupos; l++){
                          for(int m = 0; m < EntradaCaio.numGrupos; m++){
                            EntradaCaio.sigmaEsp[i][j2][m][l] = stod(linha.substr(0,linha.find(" ")));
                            linha = linha.substr((linha.find(" ")+1));
                          }
                          getline(arq,linha);
                        }
                      }
                      if(j2 < EntradaCaio.grauAnisotropia){
                        getline(arq,linha);
                      }
                      if(EntradaCaio.grauAnisotropia == 0 && primeiro){
                        j2--;
                        getline(arq,linha);
                      }
                    }
                  }else{
                    i--;
                  }
                  if(i < EntradaCaio.numZonas - 1){
                    getline(arq,linha);
                  }
                }
              }
              break;
            case 13:
              EntradaCaio.tipoCc[0] = stoi(linha.substr(0,linha.find(" ")));
              EntradaCaio.tipoCc[1] = stoi(linha.substr((linha.find(" ")+1)));
              break;
            case 14:
              for(int i = 0;i < EntradaCaio.numGrupos;i++){
                EntradaCaio.valorCc[i][0] = stod(linha.substr(0,linha.find(" ")));
                EntradaCaio.valorCc[i][1] = stod(linha.substr((linha.find(" ")+1)));
                getline(arq,linha);
              }
              break;
            case 15:
              int k;
              for(int i = 0;i < EntradaCaio.numGrupos;i++){
                k = 0;
                for(int j = 0;j < EntradaCaio.numRegioes;j++){
                  EntradaCaio.fonteFisica[i][j] = stod(linha.substr(0,linha.find(" ")));
                  if(EntradaCaio.fonteFisica[i][j] != 0){
                    EntradaCaio.mapeaFonteFisica[i][k] = j;
                    k++;
                  }
                  linha = linha.substr((linha.find(" ")+1));
                }
                for(int l = k;l < EntradaCaio.numRegioes;l++){
                    EntradaCaio.mapeaFonteFisica[i][l] = -1;
                }
                getline(arq,linha);
              }
              break;
            case 16:
              for(int i = 0;i < EntradaCaio.numGrupos;i++){
                for(int j = 0;j < EntradaCaio.numRegioes;j++){
                  EntradaCaio.fonteAdjunta[i][j] = stod(linha.substr(0,linha.find(" ")));
                  if(EntradaCaio.fonteAdjunta[i][j] != 0 && !grupoDetector){
                    EntradaCaio.numDetector++;
                    grupoDetector = true;
                  }
                  linha = linha.substr((linha.find(" ")+1));
                }
                grupoDetector = false;
                getline(arq,linha);
              }
              EntradaCaio.mapeaDetector = new double[EntradaCaio.numDetector];
              for(int g = 0;g < EntradaCaio.numGrupos;g++){
                 for(int i = 0;i < EntradaCaio.numRegioes;i++){
                   if(EntradaCaio.fonteAdjunta[g][i] != 0){
                     EntradaCaio.mapeaDetector[counter] = g;
                     counter++;
                     break;
                   }
                 }
              }
              break;
            default:
              break;
          }
        }
      }
    }
    Quadratura1D Quad1D_GSUS;
    Quad1D_GSUS.DadosQuadrat_GL(EntradaCaio);
    arq.close();
}

void gerarGraficosSaida(string caminho,Dados_Entrada EntradaCaio){
    ifstream arq;
    arq.open(caminho);
    gnuplot gp;
    int linhaCount = 7;
    int numLinhas = (EntradaCaio.numNodos / EntradaCaio.periodicidade) + 1;
    for(int i = 1;i <= EntradaCaio.numGrupos;i++){
      gp("set encoding iso_8859_1");
      gp("set terminal png");
      gp("set output 'saidaGrupo" + to_string(i) + ".png'");
      gp("set xlabel \"Posi\347\343o\"");
      gp("set ylabel \"Fluxo Escalar\"");
      gp("plot \"<(sed -n '" + to_string(7) + "," + to_string(linhaCount + numLinhas) + "p' dadosSaida.txt)\" using 1:"
          + to_string(i + 1) + " notitle with lines linestyle 1," + "\"<(sed -n '" + to_string(7) + "," + to_string(linhaCount + numLinhas) +
          "p' dadosSaida.txt)\" using 1:" + to_string(i + 1) +"title 'Solu\347\343o de Refer\352ncia'");
    }
}

void gerarGraficosSaida2(string caminho){
    ifstream arq;
    arq.open(caminho);
    gnuplot gp;
    gp("set encoding iso_8859_1");
    gp("set terminal png");
    gp("set output 'graficoTempoNCD.png'");
    gp("set xlabel \"N\372mero de Nodos\"");
//    gp("set xlabel \"Ordem da Quadratura\"");
    gp("set ylabel \"Tempo\"");
    gp("set style data linespoints");
    gp("plot \"Arquivos/dadosTempoNCD.txt\" title 'Serial',\"Arquivos/dadosTempoNCD2.txt\" title 'Paralelo com 2 threads',\"Arquivos/dadosTempoNCD3.txt\" title 'Paralelo com 4 threads'");
}

void gerarGraficosSaida3(string caminho){
    char comandoFinal[250];
    const char* comando1 = "plot \"Arquivos/semic2019/quad32/dadosTempoSerial.txt\" title 'S_{32} - Serial' ls 1,\"Arquivos/semic2019/quad32/dadosTempo2threads.txt\" title 'S_{32} - Paralelo (2 threads)' ls 2, ";
    const char* comando2 = "\"Arquivos/semic2019/quad32/dadosTempo4threads.txt\" title 'S_{32} - Paralelo (4 threads)' ls 3,";
    const char* comando3 = "\"Arquivos/semic2019/quad64/dadosTempoSerial.txt\" title 'S_{64} - Serial' ls 4,\"Arquivos/semic2019/quad64/dadosTempo2threads.txt\" title 'S_{64} - Paralelo (2 threads)' ls 5,";
    const char* comando4 = "\"Arquivos/semic2019/quad64/dadosTempo4threads.txt\" title 'S_{64} - Paralelo 4 threads' ls 6";
    strcpy(comandoFinal,comando1);
    strcat(comandoFinal,comando2);
    strcat(comandoFinal,comando3);
    strcat(comandoFinal,comando4);
    gnuplot gp;
    gp("set encoding iso_8859_1");
    gp("set terminal png");
    gp("set output 'graficoTemposSemic2019.png'");
    gp("set xlabel \"Discretiza\347\343o espacial (n\372mero de nodos)\"");
//    gp("set xlabel \"Ordem da Quadratura\"");
    gp("set ylabel \"Tempo (s)\"");
    gp("set style data linespoints");
    gp("set style line 1 lc rgb 'red' ps 2 pt 8");
    gp("set style line 2 lc rgb 'red' ps 2 pt 6");
    gp("set style line 3 lc rgb 'red' ps 2 pt 4");
    gp("set style line 4 lc rgb 'blue' ps 2 pt 9");
    gp("set style line 5 lc rgb 'blue' ps 2 pt 7");
    gp("set style line 6 lc rgb 'blue' ps 2 pt 5");
    gp("set key left top Left");
    gp(comandoFinal);
}

#endif // DADOSARQUIVO_H_INCLUDED
