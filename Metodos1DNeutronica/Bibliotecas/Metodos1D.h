#ifndef METODOS1D_H_INCLUDED
#define METODOS1D_H_INCLUDED

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <ctime>
#include <omp.h>
#include <algorithm>
#include <iterator>
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
      void MetodoDD(Dados_Entrada &EntradaCaio, int gAdjunto);
      void metodoMatrizResposta(Dados_Entrada &EntradaCaio);
      MatrixXd gerarMatrizResposta(Dados_Entrada &EntradaCaio,int piv);
};

void Metodos1D::MetodoDD(Dados_Entrada &EntradaCaio,int gAdjunto){
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
    double ***Sgmi,  **FESCold, maxval,  val, **FLUX, **FESC, **FESCmed, **fuga, **tAbs, **Sgim, **Q;
    double startTime,  endTime;
    clock_t startTime2, endTime2;
    string nomeArqString;
    string nomeArq2String;
    if(EntradaCaio.tipoProblema == 1){
      nomeArqString = "dadosSaidaFisico.txt";
      nomeArq2String = "dadosSaida2Fisico.txt";
    }else{
      nomeArqString = "dadosSaidaAdjuntoGrupo"+to_string(gAdjunto)+".txt";
      nomeArq2String = "dadosSaida2Adjunto.txt";
    }
    const char *nomeArq, *nomeArq2;
    nomeArq = nomeArqString.c_str();
    nomeArq2 = nomeArq2String.c_str();
    FILE* arq = fopen(nomeArq,"wb+");
    FILE* arq2 = fopen(nomeArq2,"a+");
    #ifdef _OPENMP
        startTime = omp_get_wtime();
    #else
        startTime2 = clock();
    #endif
    FESC = new double *[EntradaCaio.numGrupos];
    FESCold = new double *[EntradaCaio.numGrupos];
    FLUX = new double *[EntradaCaio.numNodos +1];
    Sgim = new double *[EntradaCaio.numNodos];
    FESCmed = new double *[EntradaCaio.numGrupos];
    fuga = new double *[EntradaCaio.numGrupos];
    tAbs = new double *[EntradaCaio.numGrupos];
    Q = new double *[EntradaCaio.numGrupos];
    for(int g = 0;g < EntradaCaio.numGrupos;g++){
        FESC[g] = new double [EntradaCaio.numNodos + 1];
        FESCold[g] = new double [EntradaCaio.numNodos + 1];
        FESCmed[g] = new double [EntradaCaio.numNodos];
        fuga[g] = new double[2];
        tAbs[g] = new double [EntradaCaio.numRegioes];
        for(int i = 0;i < EntradaCaio.numRegioes;i++){
            tAbs[g][i] = 0;
        }
        Q[g] = new double [EntradaCaio.numRegioes];
    }
    for(int i = 0;i < EntradaCaio.numNodos + 1;i++){
        FLUX[i] = new double [EntradaCaio.ordemQuad * EntradaCaio.numGrupos];
        if(i < EntradaCaio.numNodos){
            Sgim[i] = new double [EntradaCaio.ordemQuad * EntradaCaio.numGrupos];
        }
    }

    maxval = 100;
    iter = 0;

    ///Inicializar FLUX e Sgim e FESCmed e FESC e FESCold
    int igm;
    for(int g = 0;g < EntradaCaio.numGrupos;g++){
        for(int i = 0;i < EntradaCaio.numNodos + 1;i++){
            igm = g * EntradaCaio.ordemQuad - 1;
            FESC[g][i] = 0;
            FESCold[g][i] = 0;
            FLUX[i][igm] = 0;
            if(i < EntradaCaio.numNodos){
                FESCmed[g][i] = 0;
                Sgim[i][igm] = 0;
            }
        }
    }
    ///Inicializar Fonte Q
    if(EntradaCaio.tipoProblema == 1){
      for(int g = 0;g < EntradaCaio.numGrupos;g++){
        for(int i = 0;i < EntradaCaio.numRegioes;i++){
          Q[g][i] = EntradaCaio.fonteFisica[g][i];
        }
      }
    }else{
      for(int i = 0;i < EntradaCaio.numRegioes;i++){
        Q[gAdjunto][i] = EntradaCaio.fonteAdjunta[gAdjunto][i];
      }
    }

    ///Inicializar Cond Contorno
    if(EntradaCaio.tipoProblema == 1){
      for(int g = 0;g < EntradaCaio.numGrupos;g++){
        igm = g * EntradaCaio.ordemQuad - 1;
        for(int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
            igm++;
              FLUX[0][igm] = EntradaCaio.valorCc[g][0];
        }
      }
    }

    if(EntradaCaio.tipoProblema == 1){
      for(int g = 0;g < EntradaCaio.numGrupos;g++){
        igm = (g * EntradaCaio.ordemQuad) - ((EntradaCaio.ordemQuad / 2) + 1);
        for(int m = EntradaCaio.ordemQuad / 2;m < EntradaCaio.ordemQuad;m++){
            igm++;
            FLUX[EntradaCaio.numNodos][igm] = EntradaCaio.valorCc[g][1];
        }
      }
    }

    ////////////////////////////
    int IZ, NA, jback, jfront;
    double H, XT, DMI, NUM, DEN;

    ////////////////////////////
    double den, num;
    /////***********************************************///

    while(maxval > EntradaCaio.cp) {
 		maxval = 0.0;
 		///////////////////////////
    	///Varredura direita
        ///cout<<"Aquiiiiii Varredura direita"<<endl;
        if(EntradaCaio.tipoCc[0] == 2){
            for(int g = 0;g < EntradaCaio.numGrupos;g++){
                igm = g * EntradaCaio.ordemQuad - 1;
                for(int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
                    igm++;
                    FLUX[0][igm] = FLUX[0][igm + EntradaCaio.ordemQuad / 2];
                }
            }
        }

        jback = -1;
        for(int piv = 0;piv < EntradaCaio.numRegioes;piv++) {
            IZ = EntradaCaio.mapeamento[piv];
            H = EntradaCaio.tamanhoNodo[piv];
            NA = EntradaCaio.nodosRegiao[piv];
            for(int i = 0;i < NA;i++)  {
                jback = jback + 1;
                jfront = jback + 1;
                for(int g = 0;g < EntradaCaio.numGrupos;g++){
                    XT = EntradaCaio.sigmaTot[IZ - 1][g] * 0.5;
                    igm = g * EntradaCaio.ordemQuad - 1;
                    for(int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
                        igm++;
                        DMI = EntradaCaio.MI[m] / H;
                        NUM = (DMI - XT) * FLUX[jback][igm] + Sgim[jback][igm] + Q[g][piv];
                        DEN = DMI + XT;
                        FLUX[jfront][igm] = NUM / DEN;
                    }
                }
            }
        }

 		///////////////////////////
        ///Varredura esquerda
        if(EntradaCaio.tipoCc[1] == 2){
            for(int g = 0;g < EntradaCaio.numGrupos;g++){
                igm = (g * EntradaCaio.ordemQuad) + (EntradaCaio.ordemQuad / 2) - 1;
                for(int m = EntradaCaio.ordemQuad / 2;m < EntradaCaio.ordemQuad;m++){
                    igm++;
                    FLUX[EntradaCaio.numNodos][igm] = FLUX[EntradaCaio.numNodos][igm - EntradaCaio.ordemQuad / 2];
                }
            }
        }

        jfrontOriginal = jfront;
        jbackOriginal = jback;

        jfront = jfrontOriginal;
        jback = jbackOriginal;
        for(int piv = EntradaCaio.numRegioes - 1;piv > -1;piv--) {
            IZ = EntradaCaio.mapeamento[piv];
            H = EntradaCaio.tamanhoNodo[piv];
            NA = EntradaCaio.nodosRegiao[piv];
            for(int i = 0;i < NA;i++) {
                for(int g = 0;g < EntradaCaio.numGrupos;g++){
                    XT = EntradaCaio.sigmaTot[IZ - 1][g] * 0.5;
                    igm = (g * EntradaCaio.ordemQuad) + (EntradaCaio.ordemQuad / 2) - 1;
                    for(int m = EntradaCaio.ordemQuad / 2; m < EntradaCaio.ordemQuad; m++){
                        igm++;
                        DMI = -EntradaCaio.MI[m] / H;
                        NUM = (DMI - XT) * FLUX[jfront][igm] + Sgim[jback][igm] + Q[g][piv];
                        DEN = DMI + XT;
                        FLUX[jback][igm] = NUM / DEN;
                    }
                }
                jfront = jfront - 1;
                jback = jfront - 1;
            }
        }


        ///Calculo Sgim
        double soma1,soma2,somaFinal,somaIntGlobal;
        int i,g1,igm1,g2,l,n,ig1n,m;
        jback = -1;
        #pragma omp parallel for firstprivate(IZ,NA,igm1,ig1n,jfront,jback,soma1,soma2,somaFinal,somaIntGlobal) num_threads(4)
        for(int piv = 0;piv < EntradaCaio.numRegioes;piv++){
            IZ = EntradaCaio.mapeamento[piv];
            NA = EntradaCaio.nodosRegiao[piv];
            jback = EntradaCaio.somaNodosRegiao[piv];
            for(int i = 0;i < NA;i++) {
                jback = jback + 1;
                jfront = jback + 1;
                for(int g1 = 0;g1 < EntradaCaio.numGrupos;g1++){
                    igm1 = g1 * EntradaCaio.ordemQuad - 1;
                    soma2 = 0;
                    for(int m = 0;m < EntradaCaio.ordemQuad;m++){
                        igm1++;
                        soma1 = 0;
                        somaFinal = 0;
                        soma2 += FLUX[jback][igm1] * EntradaCaio.wn[m];
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
                    FESC[g1][jback] = 0.5 * soma2;
                    FESCmed[g1][jback] = FESC[g1][jback] + FESC[g1][jfront];
                    num = fabs(FESC[g1][jback] - FESCold[g1][jback]);
                    den = FESC[g1][jback];
                    val = num / den;
                    if(maxval < val){
                      maxval = val;
                    }
                    FESCold[g1][jback] = FESC[g1][jback];
                }
            }
        }
        double soma;
        #pragma omp parallel for firstprivate(igm,soma,num,den) num_threads(4)
        for(int g3 = 0;g3 < EntradaCaio.numGrupos;g3++){
            igm = g3 * EntradaCaio.ordemQuad - 1;
            soma = 0;
            for(int m1 = 0;m1 < EntradaCaio.ordemQuad;m1++){
                igm++;
                soma += FLUX[jfront][igm] * EntradaCaio.wn[m1];
            }
            FESC[g3][jfront] = 0.5 * soma;
            num = fabs(FESC[g3][jfront] - FESCold[g3][jfront]);
            den = FESC[g3][jfront];
            val = num / den;
            if (maxval < val){
                maxval = val;
            }
            FESCold[g3][jfront] = FESC[g3][jfront];
        }
		iter++;
		cout << iter << endl;
		/////**********************************************
    } ///fecha o while (maxval>erro)

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

/// Calculo Taxa de absorcao e Leitura do Detector (Caso Adjunto)
    double *lDetectorFonte;
    double lDetectorContorno;
    int *regioesAdicionadas;
    int *p;
    regioesAdicionadas = new int[EntradaCaio.numRegioes];
    int counterFonteFisica = 0;
    for(int g = 0;g < EntradaCaio.numGrupos;g++){
      for(int i = 0;i < EntradaCaio.numRegioes;i++){
        p = find(regioesAdicionadas,regioesAdicionadas + EntradaCaio.numRegioes,i);
        if(EntradaCaio.fonteFisica[g][i] != 0 && p == (regioesAdicionadas + EntradaCaio.numRegioes)){
            regioesAdicionadas[counterFonteFisica] = i;
            counterFonteFisica++;
        }
      }
    }
    lDetectorFonte = new double[counterFonteFisica];
    for(int g = 0;g < counterFonteFisica;g++){
      lDetectorFonte[g] = 0;
    }
    lDetectorContorno = 0;
    if(EntradaCaio.tipoProblema == 1){
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
            soma_global = soma;
            for(int j = 0;j < EntradaCaio.numGrupos;j++){
                sigmaEspTot += EntradaCaio.sigmaEsp[IZ - 1][0][j][g];
            }
            soma_global *= (EntradaCaio.sigmaTot[IZ - 1][g] - sigmaEspTot);
            tAbs[g][i] = soma_global;
        }
      }
    }else{
      int newG;
      double soma;
      for(int g1 = 0;g1 < counterFonteFisica;g1++){
        int realCount = 0;
        newG = regioesAdicionadas[g1];
        soma = 0;
        for(int g2 = 0;g2 < EntradaCaio.numGrupos;g2++){
            realCount = EntradaCaio.somaNodosRegiao[newG] + 1;
            for(int j = 0;j < EntradaCaio.nodosRegiao[newG];j++){
                soma += EntradaCaio.tamanhoNodo[newG] * FESCmed[g2][j + realCount];
            }
            soma *= EntradaCaio.fonteFisica[g2][newG];
            lDetectorFonte[g1] += soma;
        }
      }
      for(int g = 0;g < EntradaCaio.numGrupos;g++){
        soma = 0;
        igm = g * EntradaCaio.ordemQuad - 1;
        for(int m = 0;m < EntradaCaio.ordemQuad / 2;m++){
            igm++;
            soma += -EntradaCaio.MI[m + EntradaCaio.ordemQuad / 2] * FLUX[EntradaCaio.numNodos][igm] * EntradaCaio.wn[m] * EntradaCaio.valorCc[g][1];
        }
        lDetectorContorno += soma;
        soma = 0;
        igm = (g * EntradaCaio.ordemQuad) + (EntradaCaio.ordemQuad / 2) - 1;
        for(int m = EntradaCaio.ordemQuad / 2;m < EntradaCaio.ordemQuad;m++){
            igm++;
            soma += EntradaCaio.MI[m - EntradaCaio.ordemQuad / 2] * FLUX[0][igm] * EntradaCaio.wn[m] * EntradaCaio.valorCc[g][0];
        }
        lDetectorContorno += soma;
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
    if(EntradaCaio.tipoProblema == 1){
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
    }else{
      double lDetectorTotal = 0;
      double *somaTamanhoRegiao = 0;
      double tamanhoMin,tamanhoMax = 0;
      somaTamanhoRegiao = new double[EntradaCaio.numRegioes];
      somaTamanhoRegiao[0] = 0;
      for(int i = 1;i < EntradaCaio.numRegioes;i++){
         somaTamanhoRegiao[i] = somaTamanhoRegiao[i - 1] + EntradaCaio.tamanhoRegiao[i - 1];
      }
      fprintf(arq2,"Grupo de Energia %d\n",gAdjunto);
      for(int i = 0;i < counterFonteFisica;i++){
        tamanhoMin = somaTamanhoRegiao[regioesAdicionadas[i]];
        tamanhoMax = somaTamanhoRegiao[regioesAdicionadas[i]] + EntradaCaio.tamanhoRegiao[regioesAdicionadas[i]];
        fprintf(arq2,"Leitura Detector Fonte %d(%.1fcm - %.1fcm): %.6g\n",i,tamanhoMin,tamanhoMax,lDetectorFonte[i]);
        lDetectorTotal += lDetectorFonte[i];
      }
      fprintf(arq2,"Leitura Detector Contorno: %.6g\n",lDetectorContorno);
      lDetectorTotal += lDetectorContorno;
      fprintf(arq2,"Leitura Detector Total: %.6g\n\n",lDetectorTotal);
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

