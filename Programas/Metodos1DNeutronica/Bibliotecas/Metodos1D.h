#ifndef METODOS1D_H_INCLUDED
#define METODOS1D_H_INCLUDED

#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include "Quadratura1D.h"
#include "mpi.h"

class Metodos1D{
   public:
      void MetodoDD(Dados_Entrada &EntradaCaio,string caminhoSaida,int mynode,int totalnodes);
};

void Metodos1D::MetodoDD(Dados_Entrada &EntradaCaio,string caminhoSaida,int mynode,int totalnodes){
    ///***************************************
    /*Algumas variaveis descritas na funcao
       Smj = Termo para calculo do fluxo angular
       FESC = Fluxo Escalar (nao leva em consideracao o angulo)
       FESCmed = Fluxo Escalar medio no nodo
       FLUX = Fluxo angular (leva em consideracao o angulo)
       FLUXangmed = Fluxo angular medio no nodo
    */
    int i, j, piv, iter, startVal, endVal, accum;;
    double *Sj,  *FESCold, maxval,  val, **FLUX, *FESC, *FESCmed;
    ofstream  arq;
    clock_t time;
    MPI_Status status;
    if(mynode == 0){
      arq.open(caminhoSaida);
    }
    time = clock();
    FESC = new double[EntradaCaio.numNodos + 1];
    FESCold = new double[EntradaCaio.numNodos + 1];
    FLUX = new double *[EntradaCaio.numNodos + 1];
    Sj = new double [EntradaCaio.numNodos + 1];

    for (i = 0;i < EntradaCaio.numNodos + 1;i++){
        FLUX[i] = new double [EntradaCaio.ordemQuad];
    }

    FESCmed = new double[EntradaCaio.numNodos];
    maxval = 100;
    iter = 0;

    ///Inicializar FLUX
    for (i = 0;i < EntradaCaio.numNodos + 1;i++){
        for (j = 0;j < EntradaCaio.ordemQuad;j++){
            FLUX[i][j] = 0;
        }
        FESC[i] = 0;
    }

    ///Inicializar Cond Contorno
    for (j = 0;j < EntradaCaio.ordemQuad / 2;j++){
         FLUX[0][j] = EntradaCaio.valorCc[0];
    }

    for (j = EntradaCaio.ordemQuad / 2;j < EntradaCaio.ordemQuad;j++){
         FLUX[EntradaCaio.numNodos][j] = EntradaCaio.valorCc[1];
    }

    ///Inicializar Sj
    for (i = 0;i < EntradaCaio.numNodos + 1;i++){
            Sj[i] = 0;
    }

    ///Inicializar FESCold
    for (i = 0;i < EntradaCaio.numNodos + 1;i++){
        FESCold[i] = 0.0;
    }

    ////////////////////////////
    int IZ, NA, jback, jfront;
    double H, XT, Q, DMI, NUM, DEN;

    ////////////////////////////
    double den, num;
    ///***********************************************///
    while (maxval > EntradaCaio.cp) {

 		maxval = 0.0;
 		///////////////////////////
    	///Varredura direita
		jback = -1;
        ///cout<<"Aquiiiiii Varredura direita"<<endl;

        if (EntradaCaio.tipoCc[0] == 2){
            for (j = 0;j < EntradaCaio.ordemQuad / 2;j++){
                FLUX[0][j] = FLUX[0][j + EntradaCaio.ordemQuad / 2];
            }
        }

		for (piv = 0;piv < EntradaCaio.numRegioes;piv++) {
            IZ = EntradaCaio.mapeamento[piv];
            H = EntradaCaio.tamanhoNodo[piv];
            XT = EntradaCaio.sigmaTotZona[IZ - 1] * 0.5;
            Q = EntradaCaio.fonte[piv];
            NA = EntradaCaio.nodosRegiao[piv];
            for (i = 0;i < NA;i++)  {
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
            for (j = EntradaCaio.ordemQuad / 2;j < EntradaCaio.ordemQuad;j++){
                FLUX[EntradaCaio.numNodos][j] = FLUX[EntradaCaio.numNodos][j - EntradaCaio.ordemQuad/2];
            }
        }

		for (piv = EntradaCaio.numRegioes - 1;piv > -1;piv--) {
            IZ = EntradaCaio.mapeamento[piv];
            H = EntradaCaio.tamanhoNodo[piv];
            XT = EntradaCaio.sigmaTotZona[IZ - 1] * 0.5;
            Q = EntradaCaio.fonte[piv];
            NA = EntradaCaio.nodosRegiao[piv];
            for (i = 0;i < NA;i++) {
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
        for (piv = 0;piv < EntradaCaio.numRegioes;piv++){
             IZ = EntradaCaio.mapeamento[piv];
             NA = EntradaCaio.nodosRegiao[piv];
             for (i = 0;i < NA;i++) {
                jback = jback + 1;
                jfront = jback + 1;
                double soma = 0;
                for (int n = 0;n < EntradaCaio.ordemQuad;n++){
                    soma = soma + EntradaCaio.wn[n] * (FLUX[jfront][n] + FLUX[jback][n]) * 0.5;
                }
                soma = soma * 0.5 * EntradaCaio.sigmaEspZona[IZ - 1];
                Sj[jback] = soma;
            }
         }
         ///FLUXo escalar
         for (i = 0;i < EntradaCaio.numNodos + 1;i++)  {
            double soma = 0;
            for (int m = 0;m < EntradaCaio.ordemQuad;m++)
                soma = soma + FLUX[i][m] * EntradaCaio.wn[m];
            FESC[i] = 0.5 * soma;
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
		cout << iter << endl;
		///**********************************************
    } ///fecha o while (maxval>erro)
    time = clock() - time;
    arq << "Numero de iteracoes: " << iter << endl;
    arq << "Erro: " << maxval << endl;
    arq << "Tempo de Resolucao do Problema: " << ((double)time) / CLOCKS_PER_SEC << " segundos" << endl;
    arq << "/////////////////////////////////////////////" << endl;
    arq << "Posicao Fluxo Escalar" << endl;
    for(int i = 0;i < EntradaCaio.numNodos + 1;i++){
        arq << i << "\t" << FESC[i] << endl;
    }
    arq.close();
    ///**********************************************
   }
#endif
