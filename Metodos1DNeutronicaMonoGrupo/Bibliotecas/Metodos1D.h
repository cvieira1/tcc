#ifndef METODOS1D_H_INCLUDED
#define METODOS1D_H_INCLUDED

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <stdio.h>
#include "Quadratura1D.h"
#include "mpi.h"

using namespace std;

class Metodos1D{
   public:
      void MetodoDD(Dados_Entrada &EntradaCaio,string caminhoSaida,string caminhoSaida2,int mynode,int totalnodes);
};

void Metodos1D::MetodoDD(Dados_Entrada &EntradaCaio,string caminhoSaida,string caminhoSaida2,int mynode,int totalnodes){
    /////***************************************
    /*Algumas variaveis descritas na funcao
       Smj = Termo para calculo do fluxo angular
       FESC = Fluxo Escalar (nao leva em consideracao o angulo)
       FESCmed = Fluxo Escalar medio no nodo
       FLUX = Fluxo angular (leva em consideracao o angulo)
       FLUXangmed = Fluxo angular medio no nodo
       arq = arquivo de saida que contem numero de iteracoes,erro,tempo de resolucao e fluxo escalar
       arq2 = arquivo de saida que contem a fuga nas regioes de contorno e a taxa de absorcao em todas as regioes
       tAbs = taxa de absorcao
    */
    int iter,startVal, endVal, accum;
    double *Sj,  *FESCold, maxval,  val, **FLUX, *FESC, *FESCmed, *fuga, *tAbs;
    MPI_Status status;
    FILE* arq = fopen("dadosSaida.txt","wb+");
    FILE* arq2 = fopen("dadosSaida2.txt","wb+");
    clock_t time;
    time = clock();
    FESC = new double[EntradaCaio.numNodos + 1];
    FESCold = new double[EntradaCaio.numNodos + 1];
    FLUX = new double *[EntradaCaio.numNodos + 1];
    Sj = new double [EntradaCaio.numNodos + 1];
    for (int i = 0;i < EntradaCaio.numNodos + 1;i++){
        FLUX[i] = new double [EntradaCaio.ordemQuad];
    }
    FESCmed = new double[EntradaCaio.numNodos];
    fuga = new double[2];
    tAbs = new double[EntradaCaio.numRegioes];
    maxval = 100;
    iter = 0;

    ///Inicializar FLUX
    for (int i = 0;i < EntradaCaio.numNodos + 1;i++){
        for (int j = 0;j < EntradaCaio.ordemQuad;j++){
            FLUX[i][j] = 0;
        }
        FESC[i] = 0;
    }

    ///Inicializar FESCmed
    for(int i = 0;i < EntradaCaio.numNodos;i++){
        FESCmed[i] = 0;
    }

    ///Inicializar Cond Contorno
    for (int j = 0;j < EntradaCaio.ordemQuad / 2;j++){
         FLUX[0][j] = EntradaCaio.valorCc[0];
    }

    for (int j = EntradaCaio.ordemQuad / 2;j < EntradaCaio.ordemQuad;j++){
         FLUX[EntradaCaio.numNodos][j] = EntradaCaio.valorCc[1];
    }

    ///Inicializar Sj
    for (int i = 0;i < EntradaCaio.numNodos + 1;i++){
        Sj[i] = 0;
    }

    ///Inicializar FESCold
    for (int i = 0;i < EntradaCaio.numNodos + 1;i++){
        FESCold[i] = 0.0;
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
		jback = -1;
        ///cout<<"Aquiiiiii Varredura direita"<<endl;

        if (EntradaCaio.tipoCc[0] == 2){
            for (int j = 0;j < EntradaCaio.ordemQuad / 2;j++){
                FLUX[0][j] = FLUX[0][j + EntradaCaio.ordemQuad / 2];
            }
        }

		for (int piv = 0;piv < EntradaCaio.numRegioes;piv++) {
            IZ = EntradaCaio.mapeamento[piv];
            H = EntradaCaio.tamanhoNodo[piv];
            XT = EntradaCaio.sigmaTotZona[IZ - 1] * 0.5;
            Q = EntradaCaio.fonte[piv];
            NA = EntradaCaio.nodosRegiao[piv];
            for (int i = 0;i < NA;i++)  {
                jback = jback + 1;
                jfront = jback + 1;
                for (int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
                    DMI = EntradaCaio.MI[m] / H;
                    NUM = (DMI - XT) * FLUX[jback][m] + Sj[jback] + Q;
                    DEN = DMI + XT;
                    FLUX[jfront][m] = NUM / DEN;
                }
            }
		}

 		///////////////////////////
        ///Varredura esquerda
        if (EntradaCaio.tipoCc[1] == 2){
            for (int j = EntradaCaio.ordemQuad / 2;j < EntradaCaio.ordemQuad;j++){
                FLUX[EntradaCaio.numNodos][j] = FLUX[EntradaCaio.numNodos][j - EntradaCaio.ordemQuad/2];
            }
        }

		for (int piv = EntradaCaio.numRegioes - 1;piv > -1;piv--) {
            IZ = EntradaCaio.mapeamento[piv];
            H = EntradaCaio.tamanhoNodo[piv];
            XT = EntradaCaio.sigmaTotZona[IZ - 1] * 0.5;
            Q = EntradaCaio.fonte[piv];
            NA = EntradaCaio.nodosRegiao[piv];
            for (int i = 0;i < NA;i++) {
                 for (int m = EntradaCaio.ordemQuad / 2; m < EntradaCaio.ordemQuad; m++){
                     DMI = -EntradaCaio.MI[m] / H;
                     NUM = (DMI - XT) * FLUX[jfront][m] + Sj[jback] + Q;
                     DEN = DMI + XT;
                     FLUX[jback][m] = NUM / DEN;
                 }
                 jfront = jfront - 1;
                 jback = jfront - 1;
            }
		}

        ///Calculo Sj
        jback = -1;
        for (int piv = 0;piv < EntradaCaio.numRegioes;piv++){
             IZ = EntradaCaio.mapeamento[piv];
             NA = EntradaCaio.nodosRegiao[piv];
             for (int i = 0;i < NA;i++) {
                jback = jback + 1;
                jfront = jback + 1;
                double soma = 0;
                double soma_global = 0;
                startVal = EntradaCaio.ordemQuad * mynode / totalnodes;
                endVal = EntradaCaio.ordemQuad * (mynode + 1) / totalnodes;
                for (int n = startVal;n < endVal;n++){
                    soma = soma + EntradaCaio.wn[n] * (FLUX[jfront][n] + FLUX[jback][n]) * 0.5;
                }
                FESCmed[jback] = soma;
                cout << soma << "teste" << endl;
                cin.get();
                MPI_Allreduce(&soma,&soma_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                soma_global = soma_global * 0.5 * EntradaCaio.sigmaEspZona[IZ - 1];
                Sj[jback] = soma_global;
            }
         }

        ///Fluxo escalar
        for (int i = 0;i < EntradaCaio.numNodos + 1;i++)  {
            double soma = 0;
            double soma_global = 0;
            startVal = EntradaCaio.ordemQuad * mynode / totalnodes;
            endVal = EntradaCaio.ordemQuad * (mynode + 1) / totalnodes;
            for (int m = startVal;m < endVal;m++){
                soma = soma + FLUX[i][m] * EntradaCaio.wn[m];
            }
            MPI_Allreduce(&soma,&soma_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            FESC[i] = 0.5 * soma_global;
        }

		///Calculo da norma para o criterio de parada (convergencia do NBI)
		for (int i = 0;i < EntradaCaio.numNodos + 1;i++)  {
            num = fabs(FESC[i] - FESCold[i]);
            den = FESC[i];
            val = num / den;           //
            if (maxval < val){
                maxval = val;
            }
		}
		for (int i = 0;i < EntradaCaio.numNodos + 1;i++){
			FESCold[i] = FESC[i];
		}
		iter++;
		/////**********************************************
    } ///fecha o while (maxval>erro)


/// Calculo de fuga
    double soma = 0;
    for(int n = EntradaCaio.ordemQuad / 2;n < EntradaCaio.ordemQuad;n++){
        soma += -EntradaCaio.MI[n] * FLUX[0][n] * EntradaCaio.wn[n];
    }
    fuga[0] = soma;

    soma = 0;
    for(int n = 0;n < EntradaCaio.ordemQuad / 2;n++){
        soma += EntradaCaio.MI[n] * FLUX[EntradaCaio.numNodos][n] * EntradaCaio.wn[n];
    }
    fuga[1] = soma;

/// Calculo taxa de absorcao
    int iter2 = 0;
    int iter3 = 0;
    for(int i = 0;i < EntradaCaio.numRegioes;i++){
       double soma = 0;
       double soma_global = 0;
       IZ = EntradaCaio.mapeamento[i];
       startVal = EntradaCaio.nodosRegiao[i] * mynode / totalnodes;
       endVal = EntradaCaio.nodosRegiao[i] * (mynode + 1) / totalnodes;
       iter2 = iter3 + startVal;
       for(int j = startVal;j < endVal;j++){
          soma += EntradaCaio.tamanhoNodo[i] * FESCmed[iter2];
          iter2++;
       }
       iter3 += EntradaCaio.nodosRegiao[i];
       MPI_Reduce(&soma,&soma_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
       soma_global *= (EntradaCaio.sigmaTotZona[IZ - 1] - EntradaCaio.sigmaEspZona[IZ - 1]);
       tAbs[i] = soma_global;
    }

    time = clock() - time;
    if(mynode == 0){
      fprintf(arq,"Numero de iteracoes: %d\n",iter);
      fprintf(arq,"Erro: %.6g\n",maxval);
      fprintf(arq,"Tempo de Resolucao do Problema: %.6g\n segundos",(time / (double) CLOCKS_PER_SEC));
      fprintf(arq,"/////////////////////////////////////////////\n");
      fprintf(arq,"Posicao Fluxo Escalar\n");
      for(int i = 0;i < EntradaCaio.numNodos + 1;i += EntradaCaio.periodicidade){
        fprintf(arq,"%d\t%.6g\n",i,FESC[i]);
      }
      fclose(arq);
      arq2 = fopen("dadosSaida2.txt","wb+");
      fprintf(arq2,"Fuga\n");
      fprintf(arq2,"Esquerda  Direita\n");
      fprintf(arq2,"%.6g  %.6g\n",fuga[0],fuga[1]);
      fprintf(arq2,"Taxa de Absorcao\n");
      for(int i = 0;i < EntradaCaio.numRegioes;i++){
        fprintf(arq2,"Regiao %d: %.6g\n",i,tAbs[i]);
      }
      fclose(arq2);
    }else{
    }
    ///**********************************************
   }
#endif
