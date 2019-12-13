#ifndef QUADRATURA1D_H
#define QUADRATURA1D_H

#include <fstream>
#include <math.h>
#include <stdio.h>      //printf, fopen
#include <stdlib.h>
#include <iostream>
#include <string>

using namespace std;

/*
   DESCRICAO DOS ATRIBUTOS DOS DADOS DE ENTRADA
   ordemQuad = Ordem da Quadratura
   numRegioes = Numero de Nodos no Dominio
   numZonas = Numero de Zonas(materiais) no Dominio
   periodicidade = Periodicidade com que os resultados serao exibidos
   mapeamento = Mapeamento das Regios com suas Zonas
   cp = Criterio de Parada(10ö-?)
   grauAnisotropia = Grau da Anisotropia
   tipoCc = Tipo de cada Condicao de Contorno
   tamanhoRegiao = Comprimento de cada Regiao
   sigmaTot = Sigma Total para cada Grupo de Energia em cada Zona
   sigmaEspZona = Sigma de Espalhamento dos Graus de Anisotropia para cada Grupo de Energia em cada Zona
   valorCc = Valor de cada Condicao de Contorno
   fonte = Valor da Fonte em cada Regiao
   MI = Valor dos Mi de acordo com a Ordem da Quadratura
   wn = Valor dos Pesos de acordo com a Ordem da Quadratura
   numNodos = Numero total de nodos no dominio
   tamanhoNodo = Comprimento do Nodo de cada Regiao
   tamanhoDominio = Comprimento total do Dominio
*/

struct Dados_Entrada {
    int ordemQuad, numRegioes, numZonas, *nodosRegiao, periodicidade, *mapeamento, grauAnisotropia, *tipoCc, numNodos, numGrupos, tipoProblema, *somaNodosRegiao;
    int numDetector, *numFonteFisica;
    double *tamanhoRegiao, **sigmaTot, ****sigmaEsp, **valorCc, **fonteFisica, *tamanhoNodo, tamanhoDominio, cp, **fonteAdjunta, *mapeaDetector, **mapeaFonteFisica;
    double *MI, *wn;
};

class Quadratura1D {

    public:
        void DadosQuadrat_GL(Dados_Entrada &EntradaCaio);

};

///Método para atribuir os valores de MI e wn em problemas 1D
void Quadratura1D::DadosQuadrat_GL(Dados_Entrada &EntradaCaio) {
    EntradaCaio.MI=new double[EntradaCaio.ordemQuad]; //MIs da Quadratura
    EntradaCaio.wn=new double[EntradaCaio.ordemQuad]; //OMEGAs da Quadratura
    ////////////////////////////////////////////////
    //// Seleção do ordem da quadratura (mi e wn)
    //////////////////////////////////////////////
    //Calculo de los coeficientes de Gauss
    int n=EntradaCaio.ordemQuad;
    float E(10.0e-8), fracn, pin,z, pz, p1, p0, z1,dpz;
    int m=(n+1)/2;
    fracn=1.0-(1.0-1.0/n)/(8.0*n*n);
    pin=M_PI/(n+0.5);
    //////////////////////////
    for(int i=0; i<=m-1; i++) {//Pra calcular so os zeros no negativos do polinomio
        z=fracn*cos(((i+1)-0.25)*pin);
        int iter=0;
        while (1)  {
            iter++;
            pz=1.0;
            p1=0.0;
            for (int j=0; j<=n-1; j++) {
                p0=p1;
                p1=pz;
                pz=((2.0*j+1.0)*z*p1-j*p0)/(j+1.0);
            }
            dpz=(n*(p1-z*pz))/(1-z*z);
            z1=z;
            z=z1-pz/dpz;

            if (fabs(z-z1) <= E)
                break;
        }
        if (fabs(z-z1)<=E) {
            EntradaCaio.MI[i]=z;
            EntradaCaio.MI[n-i-1]=-z;
            EntradaCaio.wn[i]=2.0/((1.0-z*z)*dpz*dpz);
            EntradaCaio.wn[n-i-1]=EntradaCaio.wn[i];
        }
        else {
            EntradaCaio.MI[i]=0;
            EntradaCaio.MI[n-i-1]=0;
            EntradaCaio.wn[i]=0;
            EntradaCaio.wn[n-i-1]=0;
        }
        for (int i = 0; i < n/2; i++)  {
            EntradaCaio.MI[i]=-EntradaCaio.MI[i+n/2];
            EntradaCaio.wn[i]=EntradaCaio.wn[i+n/2];
        }

    }
}

#endif // QUADRATURA1D_H
