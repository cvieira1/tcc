#ifndef DADOSARQUIVO_H_INCLUDED
#define DADOSARQUIVO_H_INCLUDED

#include <cmath>
#include "Quadratura1D.h"

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
        cout << "Valor da Fonte na Regiao " << j+1 << ": " << EntradaCaio.fonte[i][j] << endl;
      }
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
    int ldado,exp;
    string linha;
    ifstream arq;
    bool primeiro = true;
    EntradaCaio.tipoCc = new int[2];
    arq.open(caminho);
    ldado = 0;
    if(arq.is_open()){
      while(!arq.eof()){
        getline(arq,linha);
        if(linha.find("//") == string::npos){ //Ira entrar se a linha nao for 1 comentario
          ldado++;
          switch(ldado){
            case 2:
              EntradaCaio.ordemQuad = stoi(linha);
              break;
            case 3:
              EntradaCaio.numRegioes = stoi(linha.substr(0,linha.find(" ")));
              EntradaCaio.numZonas = stoi(linha.substr((linha.find(" ")+1)));
              EntradaCaio.tamanhoRegiao = new double[EntradaCaio.numRegioes];
              EntradaCaio.nodosRegiao = new int[EntradaCaio.numRegioes];
              EntradaCaio.mapeamento = new int[EntradaCaio.numRegioes];
              EntradaCaio.tamanhoNodo = new double[EntradaCaio.numRegioes];
              EntradaCaio.numNodos = 0;
              EntradaCaio.tamanhoDominio = 0;
              break;
            case 4:
              for(int i = 0; i < EntradaCaio.numRegioes;i++){
                EntradaCaio.tamanhoRegiao[i] = stod(linha.substr(0,linha.find(" ")));
                EntradaCaio.tamanhoDominio += EntradaCaio.tamanhoRegiao[i];
                linha = linha.substr((linha.find(" ")+1));
              }
              break;
            case 5:
              for(int i = 0; i < EntradaCaio.numRegioes;i++){
                EntradaCaio.nodosRegiao[i] = stoi(linha.substr(0,linha.find(" ")));
                EntradaCaio.numNodos += EntradaCaio.nodosRegiao[i];
                EntradaCaio.tamanhoNodo[i] = EntradaCaio.tamanhoRegiao[i] / EntradaCaio.nodosRegiao[i];
                linha = linha.substr((linha.find(" ")+1));
              }
              break;
            case 6:
              EntradaCaio.periodicidade = stoi(linha);
              break;
            case 7:
              for(int i = 0; i < EntradaCaio.numRegioes;i++){
                EntradaCaio.mapeamento[i] = stoi(linha.substr(0,linha.find(" ")));
                linha = linha.substr((linha.find(" ")+1));
              }
              break;
            case 8:
              exp = (-1)*stoi(linha);
              EntradaCaio.cp = pow(10,exp);
              break;
            case 9:
              EntradaCaio.numGrupos = stoi(linha);
              EntradaCaio.fonte = new double *[EntradaCaio.numGrupos];
              EntradaCaio.valorCc = new double *[EntradaCaio.numGrupos];
              for(int i = 0;i < EntradaCaio.numGrupos;i++){
                EntradaCaio.fonte[i] = new double[EntradaCaio.numRegioes];
                EntradaCaio.valorCc[i] = new double[2];
              }
              EntradaCaio.sigmaTot = new double *[EntradaCaio.numZonas];
              for(int i = 0;i < EntradaCaio.numZonas;i++){
                EntradaCaio.sigmaTot[i] = new double [EntradaCaio.numGrupos];
              }
              break;
            case 10:
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
            case 11:
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
              break;
            case 12:
              EntradaCaio.tipoCc[0] = stoi(linha.substr(0,linha.find(" ")));
              EntradaCaio.tipoCc[1] = stoi(linha.substr((linha.find(" ")+1)));
              break;
            case 13:
              for(int i = 0;i < EntradaCaio.numGrupos;i++){
                EntradaCaio.valorCc[i][0] = stod(linha.substr(0,linha.find(" ")));
                EntradaCaio.valorCc[i][1] = stod(linha.substr((linha.find(" ")+1)));
                getline(arq,linha);
              }
              break;
            case 14:
              for(int i = 0;i < EntradaCaio.numGrupos;i++){
                for(int j = 0;j < EntradaCaio.numRegioes;j++){
                  EntradaCaio.fonte[i][j] = stod(linha.substr(0,linha.find(" ")));
                  linha = linha.substr((linha.find(" ")+1));
                }
                getline(arq,linha);
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

#endif // DADOSARQUIVO_H_INCLUDED
