#ifndef METODOS1D_H_INCLUDED
#define METODOS1D_H_INCLUDED

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <ctime>
#include <omp.h>
#include "Quadratura1D.h"
#include "legendre.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace std;
using namespace Legendre;
using Eigen::MatrixXd;
using Eigen::EigenSolver;

class Metodos1D{
   public:
      void MetodoDD(Dados_Entrada &EntradaCaio,string caminhoSaida,string caminhoSaida2);
      void metodoMatrizResposta(Dados_Entrada &EntradaCaio);
      MatrixXd gerarMatrizResposta(Dados_Entrada &EntradaCaio,int piv);
};

void Metodos1D::MetodoDD(Dados_Entrada &EntradaCaio,string caminhoSaida,string caminhoSaida2){
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
    int iter;
    double ***Sgmi, maxval,  val, **FLUX, **FESC, **FESCmed, **fuga, **tAbs,**Sgim;
    double startTime,  endTime;
    clock_t startTime2, endTime2;
    FILE* arq = fopen("dadosSaida.txt","wb+");
    FILE* arq2 = fopen("dadosSaida2.txt","wb+");
    #ifdef _OPENMP
        startTime = omp_get_wtime();
    #else
        startTime2 = clock();
    #endif
    FESC = new double *[EntradaCaio.numGrupos];
    FLUX = new double *[EntradaCaio.numNodos + 1];
    Sgim = new double *[EntradaCaio.numNodos];
    FESCmed = new double *[EntradaCaio.numGrupos];
    fuga = new double *[EntradaCaio.numGrupos];
    tAbs = new double *[EntradaCaio.numGrupos];
    for(int g = 0;g < EntradaCaio.numGrupos;g++){
        FESC[g] = new double [EntradaCaio.numNodos + 1];
        FESCmed[g] = new double [EntradaCaio.numNodos];
        fuga[g] = new double[2];
        tAbs[g] = new double [EntradaCaio.numRegioes];
    }
    for(int i = 0;i < EntradaCaio.numNodos + 1;i++){
        FLUX[i] = new double [EntradaCaio.ordemQuad * EntradaCaio.numGrupos];
        if(i < EntradaCaio.numNodos){
            Sgim[i] = new double [EntradaCaio.ordemQuad * EntradaCaio.numGrupos];
        }
    }

    maxval = 100;
    iter = 0;


    ///Inicializar FLUX e Sgim e FESCmed e FESC
    int igm;
    for(int g = 0;g < EntradaCaio.numGrupos;g++){
        for(int i = 0;i < EntradaCaio.numNodos + 1;i++){
            igm = g * EntradaCaio.ordemQuad - 1;
            FESC[g][i] = 0;
            FLUX[i][igm] = 0;
            if(i < EntradaCaio.numNodos){
                FESCmed[g][i] = 0;
                Sgim[i][igm] = 0;
            }
        }
    }

    ///Inicializar Cond Contorno
    for(int g = 0;g < EntradaCaio.numGrupos;g++){
        igm = g * EntradaCaio.ordemQuad - 1;
        for(int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
            igm++;
            FLUX[0][igm] = EntradaCaio.valorCc[g][0];
        }
    }

    for(int g = 0;g < EntradaCaio.numGrupos;g++){
        igm = (g * EntradaCaio.ordemQuad) - ((EntradaCaio.ordemQuad / 2) + 1);
        for(int m = EntradaCaio.ordemQuad / 2;m < EntradaCaio.ordemQuad;m++){
            igm++;
            FLUX[EntradaCaio.numNodos][igm] = EntradaCaio.valorCc[g][1];
        }
    }

    ////////////////////////////
    omp_set_num_threads(2);
    while(maxval > EntradaCaio.cp) {
 		maxval = 0.0;

        #pragma omp parallel sections
        {
 		///////////////////////////
    	///Varredura direita
        #pragma omp section
        {
        int igm;
        if(EntradaCaio.tipoCc[0] == 2){
            for(int g = 0;g < EntradaCaio.numGrupos;g++){
                igm = g * EntradaCaio.ordemQuad - 1;
                for(int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
                    igm++;
                    FLUX[0][igm] = FLUX[0][igm + EntradaCaio.ordemQuad / 2];
                }
            }
        }

        int NA, IZ, jfront, igm1, ig1n;
        double H, XT, Q, soma1, somaFinal, somaIntGlobal, newFLUX, DMI, NUM, DEN, den, num;
        int jback = -1;
        for(int piv = 0;piv < EntradaCaio.numRegioes;piv++) {
            IZ = EntradaCaio.mapeamento[piv];
            H = EntradaCaio.tamanhoNodo[piv];
            NA = EntradaCaio.nodosRegiao[piv];
            for(int i = 0;i < NA;i++)  {
                jback = jback + 1;
                jfront = jback + 1;
                for(int g1 = 0;g1 < EntradaCaio.numGrupos;g1++){
                    XT = EntradaCaio.sigmaTot[IZ - 1][g1] * 0.5;
                    Q = EntradaCaio.fonte[g1][piv];
                    igm1 = g1 * EntradaCaio.ordemQuad - 1;
                    for(int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
                        igm1++;
                        soma1 = 0;
                        somaFinal = 0;
                        DMI = EntradaCaio.MI[m] / H;
                        NUM = (DMI - XT) * FLUX[jback][igm1] + Sgim[jback][igm1] + Q;
                        DEN = DMI + XT;
                        newFLUX = NUM / DEN;
                        num = fabs(newFLUX - FLUX[jfront][igm1]);
                        den = fabs(newFLUX);
                        val = num / den;
                        if (maxval < val){
                            maxval = val;
                        }
                        FLUX[jfront][igm1] = newFLUX;
                        for(int g2 = 0;g2 < EntradaCaio.numGrupos;g2++){
                            soma1 = 0;
                            ig1n = g2 * EntradaCaio.ordemQuad - 1;
                            for(int n = 0;n < EntradaCaio.ordemQuad;n++){
                                somaIntGlobal = 0;
                                ig1n++;
                                for(int l = 0;l <= EntradaCaio.grauAnisotropia;l++){
                                    somaIntGlobal += (2*l + 1) * 0.5 * EntradaCaio.sigmaEsp[IZ - 1][l][g1][g2] * Pn(l,EntradaCaio.MI[m]) * Pn(l,EntradaCaio.MI[n]);
                                }
                                soma1 += somaIntGlobal * EntradaCaio.wn[n] * (FLUX[jfront][ig1n] + FLUX[jback][ig1n]) * 0.5;
                            }
                            somaFinal += soma1;
                        }
                        Sgim[jback][igm1] = somaFinal;
                    }
                }
            }
        }


        }

 		///////////////////////////
        ///Varredura esquerda
        #pragma omp section
        {
        int igm, NA, IZ, jfront, igm1, ig1n;
        double H, XT, Q, soma1, somaFinal, somaIntGlobal, newFLUX, DMI, NUM, DEN, den, num;
        int jback = EntradaCaio.numNodos - 1;
        jfront = EntradaCaio.numNodos;
        if(EntradaCaio.tipoCc[1] == 2){
            for(int g = 0;g < EntradaCaio.numGrupos;g++){
                igm = (g * EntradaCaio.ordemQuad) + (EntradaCaio.ordemQuad / 2) - 1;
                for(int m = EntradaCaio.ordemQuad / 2;m < EntradaCaio.ordemQuad;m++){
                    igm++;
                    FLUX[EntradaCaio.numNodos][igm] = FLUX[EntradaCaio.numNodos][igm - EntradaCaio.ordemQuad / 2];
                }
            }
        }

        for(int piv = EntradaCaio.numRegioes - 1;piv > -1;piv--) {
            IZ = EntradaCaio.mapeamento[piv];
            H = EntradaCaio.tamanhoNodo[piv];
            NA = EntradaCaio.nodosRegiao[piv];
            for(int i = 0;i < NA;i++) {
                for(int g1 = 0;g1 < EntradaCaio.numGrupos;g1++){
                    XT = EntradaCaio.sigmaTot[IZ - 1][g1] * 0.5;
                    Q = EntradaCaio.fonte[g1][piv];
                    igm1 = (g1 * EntradaCaio.ordemQuad) + (EntradaCaio.ordemQuad / 2) - 1;
                    for(int m = EntradaCaio.ordemQuad / 2; m < EntradaCaio.ordemQuad; m++){
                        igm1++;
                        soma1 = 0;
                        somaFinal = 0;
                        DMI = -EntradaCaio.MI[m] / H;
                        NUM = (DMI - XT) * FLUX[jfront][igm1] + Sgim[jback][igm1] + Q;
                        DEN = DMI + XT;
                        newFLUX = NUM / DEN;
                        num = fabs(newFLUX - FLUX[jback][igm1]);
                        den = fabs(newFLUX);
                        val = num / den;
                        if (maxval < val){
                            maxval = val;
                        }
                        FLUX[jback][igm1] = newFLUX;
                        for(int g2 = 0;g2 < EntradaCaio.numGrupos;g2++){
                            soma1 = 0;
                            ig1n = g2 * EntradaCaio.ordemQuad - 1;
                            for(int n = 0;n < EntradaCaio.ordemQuad;n++){
                                somaIntGlobal = 0;
                                ig1n++;
                                for(int l = 0;l <= EntradaCaio.grauAnisotropia;l++){
                                    somaIntGlobal += (2*l + 1) * 0.5 * EntradaCaio.sigmaEsp[IZ - 1][l][g1][g2] * Pn(l,EntradaCaio.MI[m]) * Pn(l,EntradaCaio.MI[n]);
                                }
                                soma1 += somaIntGlobal * EntradaCaio.wn[n] * (FLUX[jfront][ig1n] + FLUX[jback][ig1n]) * 0.5;
                            }
                            somaFinal += soma1;
                        }
                        Sgim[jback][igm1] = somaFinal;
                    }
                }
                jfront = jfront - 1;
                jback = jfront - 1;
            }
        }

        }

        }

		iter++;
		cout << iter << endl;
    }

///Fluxo escalar e Fluxo escalar medio
    double soma;
    for(int i = 0;i < EntradaCaio.numNodos;i++){
        for(int g = 0;g < EntradaCaio.numGrupos;g++){
            igm = g * EntradaCaio.ordemQuad - 1;
            soma = 0;
            for(int m = 0;m < EntradaCaio.ordemQuad;m++){
                igm++;
                soma += FLUX[i][igm] * EntradaCaio.wn[m];
            }
            FESC[g][i] = 0.5 * soma;
            FESCmed[g][i] = FESC[g][i] + FESC[g][i + 1];
        }
    }


/// Calculo de fuga
    for(int g = 0;g < EntradaCaio.numGrupos;g++){
        double soma = 0;
        igm = (g * EntradaCaio.ordemQuad) + (EntradaCaio.ordemQuad / 2) - 1;
        for(int m = EntradaCaio.ordemQuad / 2;m < EntradaCaio.ordemQuad;m++){
            igm++;
            soma += -EntradaCaio.MI[m] * FLUX[0][igm] * EntradaCaio.wn[m];
        }
        fuga[g][0] = soma * 0.5;
    }

    for(int g = 0;g < EntradaCaio.numGrupos;g++){
        double soma = 0;
        igm = g * EntradaCaio.ordemQuad - 1;
        for(int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
            igm++;
            soma += EntradaCaio.MI[m] * FLUX[EntradaCaio.numNodos][igm] * EntradaCaio.wn[m];
        }
        fuga[g][1] = soma * 0.5;
    }

/// Calculo taxa de absorcao
    int IZ;
    for(int g1 = 0;g1 < EntradaCaio.numGrupos;g1++){
        int realCount = 0;
        for(int i = 0;i < EntradaCaio.numRegioes;i++){
            double soma = 0;
            double sigmaEspTot = 0;
            IZ = EntradaCaio.mapeamento[i];
            for(int j = 0;j < EntradaCaio.nodosRegiao[i];j++){
                soma += EntradaCaio.tamanhoNodo[i] * FESCmed[g1][j + realCount];
            }
            realCount += EntradaCaio.nodosRegiao[i];
            for(int g2 = 0;g2 < EntradaCaio.numGrupos;g2++){
                sigmaEspTot += EntradaCaio.sigmaEsp[IZ - 1][0][g2][g1];
            }
            soma *= (EntradaCaio.sigmaTot[IZ - 1][g1] - sigmaEspTot);
            tAbs[g1][i] = soma;
        }
    }
    #ifdef _OPENMP
        endTime = omp_get_wtime() - startTime;
    #else
        endTime2 = clock();
        endTime = (endTime2 - startTime2) / (double) CLOCKS_PER_SEC;
    #endif
    cout << "Tempo: " << endTime << endl;
    double taxa = (EntradaCaio.periodicidade * EntradaCaio.tamanhoDominio / EntradaCaio.numNodos);
    fprintf(arq,"Numero de iteracoes: %d\n",iter);
    fprintf(arq,"Erro: %.6g\n",maxval);
    fprintf(arq,"Tempo de Resolucao do Problema: %.6g segundos\n",endTime);
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
    ///**********************************************
   }

void Metodos1D::metodoMatrizResposta(Dados_Entrada &EntradaCaio){
    double ordemMatriz = EntradaCaio.numGrupos * EntradaCaio.ordemQuad;
    MatrixXd AA((int)ordemMatriz, (int)ordemMatriz);
    for(int piv = 0;piv < EntradaCaio.numZonas;piv++){
        AA = gerarMatrizResposta(EntradaCaio, piv);
        EigenSolver<MatrixXd> eS(AA);
        cout << "Autovalores de AA:" << endl << eS.eigenvalues() << endl;
        cout << "Matriz de AutoVetores de AA: " << endl << eS.eigenvectors() << endl;
    }
}

MatrixXd Metodos1D::gerarMatrizResposta(Dados_Entrada &EntradaCaio,int piv){
    int aux;
    int IZ = EntradaCaio.mapeamento[piv];
    int igm = -1;
    double ordemMatriz = EntradaCaio.numGrupos * EntradaCaio.ordemQuad;
    MatrixXd AA((int)ordemMatriz, (int)ordemMatriz);
    cout << "piv:" << piv << endl;
    for(int g = 0;g < EntradaCaio.numGrupos;g++){
        for(int m = 0;m < EntradaCaio.ordemQuad;m++){
            igm++;
            int jg1n = -1;
            for(int g1 = 0;g1 < EntradaCaio.numGrupos;g1++){
                for(int n = 0;n < EntradaCaio.ordemQuad;n++){
                    jg1n++;
                    if((g == g1) && (m == n)){
                        aux = 1;
                    }else{
                        aux = 0;
                    }
                    double soma = 0;
                    for(int l = 0;l <= EntradaCaio.grauAnisotropia;l++){
                        soma += (2 * l + 1) * 0.5 * EntradaCaio.sigmaEsp[IZ - 1][l][g][g1] / EntradaCaio.sigmaTot[IZ - 1][g] *
                        Pn(l, EntradaCaio.MI[m]) * Pn(l, EntradaCaio.MI[n]) * EntradaCaio.wn[n];
                    }
                    AA(igm, jg1n) = EntradaCaio.sigmaTot[IZ - 1][g] / EntradaCaio.MI[m] * (aux - soma);
                    cout << AA(igm, jg1n) << " ";
                }
            }
            cout << endl;
        }
    }
    return AA;
}
#endif

