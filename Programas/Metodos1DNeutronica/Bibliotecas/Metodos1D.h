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
      void MetodoDD(Dados_Entrada &EntradaCaio,string caminhoSaida,string caminhoSaida2,int mynode,int totalnodes);
};

void Metodos1D::MetodoDD(Dados_Entrada &EntradaCaio,string caminhoSaida,string caminhoSaida2,int mynode,int totalnodes){
    ///***************************************
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
    //ofstream  arq,arq2;
    MPI_File arq,arq2;
    clock_t time;
    MPI_Status status;
    //if(mynode == 0){
       //arq.open(caminhoSaida);
       //arq2.open(caminhoSaida2);
    //}
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
    ///***********************************************///
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
                startVal = EntradaCaio.ordemQuad * mynode / totalnodes;
                endVal = EntradaCaio.ordemQuad * (mynode + 1) / totalnodes;
                for (int n = startVal;n < endVal;n++){
                    soma = soma + EntradaCaio.wn[n] * (FLUX[jfront][n] + FLUX[jback][n]) * 0.5;
                }
                FESCmed[jback] = soma;
                if(mynode != 0){
                  MPI_Send(&soma,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
                }else{
                  for(int j=1;j<totalnodes;j=j+1){
                    MPI_Recv(&accum,1,MPI_DOUBLE,j,1,MPI_COMM_WORLD, &status);
                    soma += soma;
                  }
                }
                soma = soma * 0.5 * EntradaCaio.sigmaEspZona[IZ - 1];
                Sj[jback] = soma;
            }
         }
        ///Fluxo escalar
        for (int i = 0;i < EntradaCaio.numNodos + 1;i++)  {
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
		///**********************************************
    } ///fecha o while (maxval>erro)
    time = clock() - time;

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
    for(int i = 0;i < EntradaCaio.numRegioes;i++){
       double soma = 0;
       IZ = EntradaCaio.mapeamento[i];
       for(int j = 0;j < EntradaCaio.nodosRegiao[i];j++){
          soma += EntradaCaio.tamanhoNodo[i] * FESCmed[iter2];
          iter2++;
       }
       soma *= (EntradaCaio.sigmaTotZona[IZ - 1] - EntradaCaio.sigmaEspZona[IZ - 1]);
       tAbs[i] = soma;
    }
    if(mynode == 0){
        cout << "teste2" << endl;
    }
      if(mynode == 1){
    cout << "teste" << endl;
      int rc = MPI_File_open(MPI_COMM_SELF, "dadosSaida4.txt",MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&arq);
       if(rc)
{
    cout << "teste" << endl;
}
            //fprintf(f,"%d \n",i);
            int *buf;
buf = (int *)malloc( 1 * sizeof(int) );
      buf[0] = mynode;
            MPI_File_write(arq,buf,1, MPI_INT,&status);
           MPI_File_close(&arq);}
/*    arq << "Numero de iteracoes: " << iter << endl;
    arq << "Erro: " << maxval << endl;
    arq << "Tempo de Resolucao do Problema: " << ((double)time) / CLOCKS_PER_SEC << " segundos" << endl;
    arq << "/////////////////////////////////////////////" << endl;
    arq << "Posicao Fluxo Escalar" << endl;
    for(int i = 0;i < EntradaCaio.numNodos + 1;i += EntradaCaio.periodicidade){
        arq << i << "\t" << FESC[i] << endl;
    }
    arq.close();
    arq2 << "Fuga" << endl;
    arq2 << "Esquerda " << " Direita" << endl;
    arq2 << fuga[0] << "  " << fuga[1] << endl;
    arq2 << "Taxa de Absorcao" << endl;
    for(int i = 0;i < EntradaCaio.numRegioes;i++){
        arq2 << "Regiao "<< i << ": " << tAbs[i] << endl;
    }
    arq2.close();*/
    ///**********************************************
   }
#endif
