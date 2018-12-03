#ifndef DADOSARQUIVO_H_INCLUDED
#define DADOSARQUIVO_H_INCLUDED

#include <cmath>
#include "Quadratura1D.h"

using namespace std;

void exibeDadosEntrada(Dados_Entrada EntradaCaio){
    cout << "Ordem da Quadratura: " << EntradaCaio.ordemQuad << endl;
    cout << "Numero de Regioes: " << EntradaCaio.numRegioes << endl;
    cout << "Numero de Zonas: " << EntradaCaio.numZonas << endl;
    for(int i = 0;i<EntradaCaio.numRegioes;i++){
      cout << "Tamanho da Regiao " << i+1 << ": " << EntradaCaio.tamanhoRegiao[i] << endl;
    }
    for(int i = 0;i<EntradaCaio.numRegioes;i++){
      cout << "Nodos da Regiao " << i+1 << ": " << EntradaCaio.nodosRegiao[i] << endl;
    }
    cout << "Periodicidade: " << EntradaCaio.periodicidade << endl;
    for(int i = 0;i<EntradaCaio.numRegioes;i++){
      cout << "Mapeamento da Regiao " << i+1 << ": " << EntradaCaio.mapeamento[i] << endl;
    }
    cout << "Criterio de Parada: " << EntradaCaio.cp << endl;
    cout << "Grau da Anisotropia: " << EntradaCaio.grauAnisotropia << endl;
    for(int i = 0;i<EntradaCaio.numZonas;i++){
      cout << "Sigma Total da Zona " << i+1 << ": " << EntradaCaio.sigmaTotZona[i] << endl;
    }
    for(int i = 0;i<EntradaCaio.numZonas;i++){
      cout << "Sigma Espalhamento da Zona " << i+1 << ": " << EntradaCaio.sigmaEspZona[i] << endl;
    }
    cout << "Tipo das Condicao de Contorno: " << EntradaCaio.tipoCc[0] << " " << EntradaCaio.tipoCc[1] << endl;
    cout << "Valor das Condicao de Contorno: " << EntradaCaio.valorCc[0] << " " << EntradaCaio.valorCc[1] << endl;
    for(int i = 0;i<EntradaCaio.numRegioes;i++){
      cout << "Valor da Fonte na Regiao " << i+1 << ": " << EntradaCaio.fonte[i] << endl;
    }
    cout << "Valores de Mi:" << endl;
    for(int i = 0;i<EntradaCaio.ordemQuad;i++){
      cout << EntradaCaio.MI[i] << endl;
    }
    cout << "Valores dos Pesos da Quadratura:" << endl;
    for(int i = 0;i<EntradaCaio.ordemQuad;i++){
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
    EntradaCaio.tipoCc = new int[2];
    EntradaCaio.valorCc = new double[2];
    arq.open(caminho);
    ldado = 0;
    if(arq.is_open()){
      while(!arq.eof()){
        getline(arq,linha);
        if(linha.find("//") == string::npos){ //Ira entrar se a linha nao for 1 comentario
          ldado++;
          switch(ldado){
            case 1:
              EntradaCaio.ordemQuad = stoi(linha);
              break;
            case 2:
              EntradaCaio.numRegioes = stoi(linha.substr(0,linha.find(" ")));
              EntradaCaio.numZonas = stoi(linha.substr((linha.find(" ")+1)));
              EntradaCaio.tamanhoRegiao = new double[EntradaCaio.numRegioes];
              EntradaCaio.nodosRegiao = new int[EntradaCaio.numRegioes];
              EntradaCaio.mapeamento = new int[EntradaCaio.numRegioes];
              EntradaCaio.fonte = new double[EntradaCaio.numRegioes];
              EntradaCaio.sigmaTotZona = new double[EntradaCaio.numZonas];
              EntradaCaio.sigmaEspZona = new double[EntradaCaio.numZonas];
              EntradaCaio.tamanhoNodo = new double[EntradaCaio.numRegioes];
              EntradaCaio.numNodos = 0;
              EntradaCaio.tamanhoDominio = 0;
              break;
            case 3:
              for(int i = 0; i < EntradaCaio.numRegioes;i++){
                EntradaCaio.tamanhoRegiao[i] = stod(linha.substr(0,linha.find(" ")));
                EntradaCaio.tamanhoDominio += EntradaCaio.tamanhoRegiao[i];
                linha = linha.substr((linha.find(" ")+1));
              }
              break;
            case 4:
              for(int i = 0; i < EntradaCaio.numRegioes;i++){
                EntradaCaio.nodosRegiao[i] = stoi(linha.substr(0,linha.find(" ")));
                EntradaCaio.numNodos += EntradaCaio.nodosRegiao[i];
                EntradaCaio.tamanhoNodo[i] = EntradaCaio.tamanhoRegiao[i] / EntradaCaio.nodosRegiao[i];
                linha = linha.substr((linha.find(" ")+1));
              }
              break;
            case 5:
              EntradaCaio.periodicidade = stoi(linha);
              break;
            case 6:
              for(int i = 0; i < EntradaCaio.numRegioes;i++){
                EntradaCaio.mapeamento[i] = stoi(linha.substr(0,linha.find(" ")));
                linha = linha.substr((linha.find(" ")+1));
              }
              break;
            case 7:
              exp = (-1)*stoi(linha);
              EntradaCaio.cp = pow(10,exp);
              break;
            case 8:
              EntradaCaio.grauAnisotropia = stoi(linha);
              break;
            case 9:
              for(int i = 0; i < EntradaCaio.numZonas;i++){
                EntradaCaio.sigmaTotZona[i] = stod(linha.substr(0,linha.find(" ")));
                linha = linha.substr((linha.find(" ")+1));
              }
              break;
            case 10:
              for(int i = 0; i < EntradaCaio.numZonas;i++){
                EntradaCaio.sigmaEspZona[i] = stod(linha.substr(0,linha.find(" ")));
                linha = linha.substr((linha.find(" ")+1));
              }
              break;
            case 11:
              EntradaCaio.tipoCc[0] = stoi(linha.substr(0,linha.find(" ")));
              EntradaCaio.tipoCc[1] = stoi(linha.substr((linha.find(" ")+1)));
              break;
            case 12:
              EntradaCaio.valorCc[0] = stod(linha.substr(0,linha.find(" ")));
              EntradaCaio.valorCc[1] = stod(linha.substr((linha.find(" ")+1)));
              break;
            case 13:
              for(int i = 0; i < EntradaCaio.numRegioes;i++){
                EntradaCaio.fonte[i] = stod(linha.substr(0,linha.find(" ")));
                linha = linha.substr((linha.find(" ")+1));
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
