#ifndef METODOS1D_H_INCLUDED
#define METODOS1D_H_INCLUDED

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "Quadratura1D.h"
#include "mpi.h"
#include "legendre.h"

using namespace std;
using namespace Legendre;

class Metodos1D{
   public:
      void MetodoDD(Dados_Entrada &EntradaCaio,string caminhoSaida,string caminhoSaida2,int mynode,int totalnodes);
};

void Metodos1D::MetodoDD(Dados_Entrada &EntradaCaio,string caminhoSaida,string caminhoSaida2,int mynode,int totalnodes){
    /////***************************************
    /*Algumas variaveis descritas na funcao
       Sgim = Termo para calculo do fluxo angular(Baseado em seu grupo de energia,ordem de quadratura e posicao no dominio
       FESC = Fluxo Escalar (nao leva em consideracao o angulo)
       FESCmed = Fluxo Escalar medio no nodo
       FLUX = Fluxo angular (leva em consideracao o angulo)
       FLUXangmed = Fluxo angular medio no nodo
       arq = arquivo de saida que contem numero de iteracoes,erro,tempo de resolucao e fluxo escalar
       arq2 = arquivo de saida que contem a fuga nas regioes de contorno e a taxa de absorcao em todas as regioes
       tAbs = taxa de absorcao
    */
    int iter,startVal, endVal, accum,jfrontOriginal,jbackOriginal;
    double ***Sgmi,  **FESCold, maxval,  val, ***FLUX, **FESC, **FESCmed, **fuga, **tAbs,***Sgim, startTime, endTime;
    FILE* arq = fopen("dadosSaida.txt","wb+");
    FILE* arq2 = fopen("dadosSaida2.txt","wb+");
    startTime = MPI_Wtime();
    FESC = new double *[EntradaCaio.numGrupos];
    FESCold = new double *[EntradaCaio.numGrupos];
    FLUX = new double **[EntradaCaio.numGrupos];
    Sgim = new double **[EntradaCaio.numGrupos];
    FESCmed = new double *[EntradaCaio.numGrupos];
    fuga = new double *[EntradaCaio.numGrupos];
    tAbs = new double *[EntradaCaio.numGrupos];
    for(int g = 0;g < EntradaCaio.numGrupos;g++){
      FESC[g] = new double [EntradaCaio.numNodos + 1];
      FESCold[g] = new double [EntradaCaio.numNodos + 1];
      FLUX[g] = new double *[EntradaCaio.numNodos + 1];
      Sgim[g] = new double *[EntradaCaio.numNodos + 1];
      FESCmed[g] = new double [EntradaCaio.numNodos];
      fuga[g] = new double[2];
      tAbs[g] = new double [EntradaCaio.numRegioes];
      for (int i = 0;i < EntradaCaio.numNodos + 1;i++){
        FLUX[g][i] = new double [EntradaCaio.ordemQuad];
        Sgim[g][i] = new double [EntradaCaio.ordemQuad];
      }
    }
    maxval = 100;
    iter = 0;

    ///Inicializar FLUX e Sgim e FESCmed e FESC e FESCold
    for(int g = 0;g < EntradaCaio.numGrupos;g++){
      for (int i = 0;i < EntradaCaio.numNodos + 1;i++){
        FESC[g][i] = 0;
        FESCold[g][i] = 0;
        if(i < EntradaCaio.numNodos){
          FESCmed[g][i] = 0;
        }
        for (int m = 0;m < EntradaCaio.ordemQuad;m++){
            FLUX[g][i][m] = 0;
            Sgim[g][i][m] = 0;
        }
      }
    }

    ///Inicializar Cond Contorno
    for(int g = 0;g < EntradaCaio.numGrupos;g++){
      for (int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
         FLUX[g][0][m] = EntradaCaio.valorCc[g][0];
      }
    }

    for(int g = 0;g < EntradaCaio.numGrupos;g++){
      for (int m = EntradaCaio.ordemQuad / 2;m < EntradaCaio.ordemQuad;m++){
         FLUX[g][EntradaCaio.numNodos][m] = EntradaCaio.valorCc[g][1];
      }
    }

    ////////////////////////////
    int IZ, NA, jback, jfront;
    double H, XT, Q, DMI, NUM, DEN;

    ////////////////////////////
    double den, num;
    /////***********************************************///
    while (maxval > EntradaCaio.cp) {

 		maxval = 0.0;
 		///////////////////////////
    	///Varredura direita
        ///cout<<"Aquiiiiii Varredura direita"<<endl;
        if (EntradaCaio.tipoCc[0] == 2){
          for(int g = 0;g < EntradaCaio.numGrupos;g++){
            for (int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
                FLUX[g][0][m] = FLUX[g][0][m + EntradaCaio.ordemQuad / 2];
            }
          }
        }


        for(int g = 0;g < EntradaCaio.numGrupos;g++){
          jback = -1;
		  for (int piv = 0;piv < EntradaCaio.numRegioes;piv++) {
            IZ = EntradaCaio.mapeamento[piv];
            H = EntradaCaio.tamanhoNodo[piv];
            XT = EntradaCaio.sigmaTot[IZ - 1][g] * 0.5;
            Q = EntradaCaio.fonte[g][piv];
            NA = EntradaCaio.nodosRegiao[piv];
            for (int i = 0;i < NA;i++)  {
               jback = jback + 1;
                jfront = jback + 1;
                for (int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
                    DMI = EntradaCaio.MI[m] / H;
                    NUM = (DMI - XT) * FLUX[g][jback][m] + Sgim[g][jback][m] + Q;
                    DEN = DMI + XT;
                    FLUX[g][jfront][m] = NUM / DEN;
                }
            }
		  }
        }

 		///////////////////////////
        ///Varredura esquerda

          if (EntradaCaio.tipoCc[1] == 2){
            for(int g = 0;g < EntradaCaio.numGrupos;g++){
              for (int m = EntradaCaio.ordemQuad / 2;m < EntradaCaio.ordemQuad;m++){
                FLUX[g][EntradaCaio.numNodos][m] = FLUX[g][EntradaCaio.numNodos][m - EntradaCaio.ordemQuad / 2];
              }
            }
          }

        jfrontOriginal = jfront;
        jbackOriginal = jback;

        for(int g = 0;g < EntradaCaio.numGrupos;g++){
          jfront = jfrontOriginal;
          jback = jbackOriginal;
		  for (int piv = EntradaCaio.numRegioes - 1;piv > -1;piv--) {
            IZ = EntradaCaio.mapeamento[piv];
            H = EntradaCaio.tamanhoNodo[piv];
            XT = EntradaCaio.sigmaTot[IZ - 1][g] * 0.5;
            Q = EntradaCaio.fonte[g][piv];
            NA = EntradaCaio.nodosRegiao[piv];
            for (int i = 0;i < NA;i++) {
                 for (int m = EntradaCaio.ordemQuad / 2; m < EntradaCaio.ordemQuad; m++){
                     DMI = -EntradaCaio.MI[m] / H;
                     NUM = (DMI - XT) * FLUX[g][jfront][m] + Sgim[g][jback][m] + Q;
                     DEN = DMI + XT;
                     FLUX[g][jback][m] = NUM / DEN;
                 }
                 jfront = jfront - 1;
                 jback = jfront - 1;
            }
		  }
        }

        ///Calculo Sgim
        for(int g1 = 0;g1 < EntradaCaio.numGrupos;g1++){
          jback = -1;
          for (int piv = 0;piv < EntradaCaio.numRegioes;piv++){
             IZ = EntradaCaio.mapeamento[piv];
             NA = EntradaCaio.nodosRegiao[piv];
             for (int i = 0;i < NA;i++) {
                jback = jback + 1;
                jfront = jback + 1;
                for(int m = 0;m < EntradaCaio.ordemQuad;m++){
                  double soma = 0;
                  double soma2 = 0; //Varaiavel necessaria para calculo do FESCmed
                  double termoSoma2 = 0;
                  double somaFinal = 0;
                  for(int g2 = 0;g2 < EntradaCaio.numGrupos;g2++){
                    soma = 0;
                    soma2 = 0;
                    for (int n = 0;n < EntradaCaio.ordemQuad;n++){
                      double somaIntLocal = 0;
                      double somaIntGlobal = 0;
                      startVal = (EntradaCaio.grauAnisotropia + 1) * mynode / totalnodes;
                      endVal = (EntradaCaio.grauAnisotropia + 1) * (mynode + 1) / totalnodes;
                      for(int l = startVal;l < endVal;l++){
                        somaIntLocal += (2*l + 1) * 0.5 * EntradaCaio.sigmaEsp[IZ - 1][l][g1][g2] * Pn(l,EntradaCaio.MI[m]) * Pn(l,EntradaCaio.MI[n]);
                      }
                      MPI_Allreduce(&somaIntLocal,&somaIntGlobal,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                      termoSoma2 = EntradaCaio.wn[n] * (FLUX[g2][jfront][n] + FLUX[g2][jback][n]) * 0.5;
                      soma += somaIntGlobal * termoSoma2;
                      soma2 += termoSoma2;
                    }
                    somaFinal += soma;
                    FESCmed[g2][jback] = soma2;
                  }
                  Sgim[g1][jback][m] = somaFinal;
                }
            }
          }
        }

        ///Fluxo escalar
        for(int g = 0;g < EntradaCaio.numGrupos;g++){
          for (int i = 0;i < EntradaCaio.numNodos + 1;i++)  {
            double soma = 0;
            double soma_global = 0;
            startVal = EntradaCaio.ordemQuad * mynode / totalnodes;
            endVal = EntradaCaio.ordemQuad * (mynode + 1) / totalnodes;
            for (int m = startVal;m < endVal;m++){
                soma = soma + FLUX[g][i][m] * EntradaCaio.wn[m];
            }
            MPI_Allreduce(&soma,&soma_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            FESC[g][i] = 0.5 * soma_global;
          }
        }

		///Calculo da norma para o criterio de parada (convergencia do NBI) e
		for(int g = 0;g < EntradaCaio.numGrupos;g++){
		  for (int i = 0;i < EntradaCaio.numNodos + 1;i++)  {
            num = fabs(FESC[g][i] - FESCold[g][i]);
            den = FESC[g][i];
            val = num / den;           //
            if (maxval < val){
                maxval = val;
            }
			FESCold[g][i] = FESC[g][i];
		  }
		}
		iter++;
		cout << iter << endl;
		/////**********************************************
    } ///fecha o while (maxval>erro)

/// Calculo de fuga
    for(int g = 0;g < EntradaCaio.numGrupos;g++){
      double soma = 0;
      for(int m = EntradaCaio.ordemQuad / 2;m < EntradaCaio.ordemQuad;m++){
        soma += -EntradaCaio.MI[m] * FLUX[g][0][m] * EntradaCaio.wn[m];
      }
      fuga[g][0] = soma * 0.5;
    }

    for(int g = 0;g < EntradaCaio.numGrupos;g++){
      double soma = 0;
      for(int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
        soma += EntradaCaio.MI[m] * FLUX[g][EntradaCaio.numNodos][m] * EntradaCaio.wn[m];
      }
      fuga[g][1] = soma * 0.5;
    }

/// Calculo taxa de absorcao

    for(int g = 0;g < EntradaCaio.numGrupos;g++){
      int realCount = 0;
      for(int i = 0;i < EntradaCaio.numRegioes;i++){
        double soma = 0;
        double soma_global = 0;
        double sigmaEspTot = 0;
        IZ = EntradaCaio.mapeamento[i];
        for(int j = 0;j < EntradaCaio.nodosRegiao[i];j++){
          soma += EntradaCaio.tamanhoNodo[i] * FESCmed[g][j + realCount];
        }
        realCount += EntradaCaio.nodosRegiao[i];
        for(int j = 0;j < EntradaCaio.numGrupos;j++){
            sigmaEspTot += EntradaCaio.sigmaEsp[IZ - 1][0][j][g];
        }
        soma_global = soma;
        soma_global *= (EntradaCaio.sigmaTot[IZ - 1][g] - sigmaEspTot);
        tAbs[g][i] = soma_global;
      }
    }
    endTime = MPI_Wtime();
    double taxa = (EntradaCaio.periodicidade * EntradaCaio.tamanhoDominio / EntradaCaio.numNodos);
    if(mynode == 0){
        fprintf(arq,"Numero de iteracoes: %d\n",iter);
        fprintf(arq,"Erro: %.6g\n",maxval);
        fprintf(arq,"Tempo de Resolucao do Problema: %.6g segundos\n",endTime - startTime);
        fprintf(arq,"/////////////////////////////////////////////\n");
        fprintf(arq,"Posicao\tFluxo Escalar\n");
        for(int g = 0;g < EntradaCaio.numGrupos;g++){
            fprintf(arq,"\tGrupo de Energia %d",g + 1);
        }
        fprintf(arq,"\n");
        double realI = 0.0;
        for(int i = 0;i < EntradaCaio.numNodos + 1;i += EntradaCaio.periodicidade){
            fprintf(arq,"%.6g\t",realI);
            for(int g = 0;g < EntradaCaio.numGrupos;g++){
                fprintf(arq,"%.6g\t",FESC[g][i]);
            }
            fprintf(arq,"\n");
            realI += taxa;
        }
        fclose(arq);
        fprintf(arq2,"Fuga Esq: ");
        for(int g = 0;g < EntradaCaio.numGrupos;g++){
          fprintf(arq2,"%.6g  ",fuga[g][0]);
        }
        fprintf(arq2,"\n");
        fprintf(arq2,"Fuga Dir: ");
        for(int g = 0;g < EntradaCaio.numGrupos;g++){
          fprintf(arq2,"%.6g  ",fuga[g][1]);
        }
        fprintf(arq2,"\n");
        fprintf(arq2,"Taxa de Absorcao\n");
        fprintf(arq2,"\t\t");
        for(int i = 0;i < EntradaCaio.numRegioes;i++){
            fprintf(arq2,"\tRegiao %d",i+1);
        }
        fprintf(arq2,"\n");
        for(int g = 0;g < EntradaCaio.numGrupos;g++){
          fprintf(arq2,"Grupo de Energia %d",g + 1);
          for(int i = 0;i < EntradaCaio.numRegioes;i++){
            fprintf(arq2,"\t%.6g",tAbs[g][i]);
          }
          fprintf(arq2,"\n");
        }
        fclose(arq2);
    }else{
    }
    ///**********************************************
   }
#endif
