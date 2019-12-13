#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <math.h>

#include "legendre.h"

using namespace std ;
using namespace Legendre ;

string archivo_dados;

int main()
{
    int NZ, NR, Li, G, tipo_met, tipo_prob;

    archivo_dados="ArchivoGerarCrossSectionsMultigrupo1D.txt";
    ifstream arch_entrada;
    string comentario;
    //cout << archivo << endl;
    arch_entrada.open(archivo_dados.c_str());
    if (arch_entrada==NULL)
//            cout <<"Arquivo nao encontrado."<<endl;
    abort();
    else {
        getline(arch_entrada,comentario); //toma la linea del archivo de entrada
        arch_entrada >> NR;   ///Numero de zonas
        arch_entrada >> NZ;   ///Numero de zonas
        getline(arch_entrada,comentario); //toma la linea del archivo de entrada
        getline(arch_entrada,comentario); //toma la linea del archivo de entrada
        arch_entrada >> tipo_prob;   ///Tipo de Problema
        getline(arch_entrada,comentario); //toma la linea del archivo de entrada
        getline(arch_entrada,comentario); //toma la linea del archivo de entrada
        arch_entrada >> tipo_met;   ///Tipo de Metodo
        getline(arch_entrada,comentario); //toma la linea del archivo de entrada
        getline(arch_entrada,comentario); //toma la linea del archivo de entrada
        arch_entrada >> G;   ///Numero de Grupos
        getline(arch_entrada,comentario); //toma la linea del archivo de entrada
        getline(arch_entrada,comentario); //toma la linea del archivo de entrada
        arch_entrada >> Li;   ///Grau de Anisotropia (0-Isotropico, 1-Linearmente Anisotropico)
    }
    arch_entrada.close();

    ///**************************************///
    ofstream arch_saida1;  // Para arquivos de saída
    string   nomedoarquivo1;
    ostringstream  auxiliar1;
    if (tipo_prob==1) {
        auxiliar1 << "FisicoSGF_locDet_XSA_MultiGrupo.txt"; // cria o novo nome parte vari?vel mais a parte fixa entre aspas
        nomedoarquivo1 = auxiliar1.str();
        arch_saida1.open(nomedoarquivo1.c_str());
        ///////////////////////////////////////////
        arch_saida1<<"////Numero de Detectores"<<endl;
        arch_saida1<<NR<<endl;
        arch_saida1<<"////Seção de choque de absorção do detector (Por Grupo de Energia, g=1..G) - Posição Inicial - Posição Final"<<endl;
        double XSSl, Hgg1, d5, d10, XT, loc(0);
        for (int piv=1; piv<NR+1; piv++){
            for (int g=1; g<G+1; g++) {
                    double somaXS0=0;
                    for (int g1=1; g1<G+1; g1++) {
                        Hgg1=0.7-((g+g1)/200.0);
                        //system("pause");
                        if (g>g1)
                            XSSl=0;     ///XSSl=(piv+20)/21.0*(g/(100.0*(g1-g+1)))*pow(Hgg1,l);
                        else
                            XSSl=(piv+20)/21.0*(g/(100.0*(g1-g+1)))*pow(Hgg1,0);  ///(2*l+1)*
        //                cout<<Hgg1<<"  "<<l<<"  "<<pow(Hgg1,l) <<endl;
        //                cout<<(i+20)/21.0<<"  "<<100.0*(g-g1+1)<<"  "<<g1/(100.0*(g-g1+1)) <<endl;
                        somaXS0+=XSSl;
                    }
                    if (g==5)
                    d5=1;
                    else d5=0;
                    if (g==10)
                        d10=1;
                    else d10=0;
                    XT=pow((piv+20.0)/21.0, 5)*(g/10.0-0.15*d5-0.15*d10);
                    arch_saida1<<XT-somaXS0<<" ";
            }
            arch_saida1<<loc<<" "<<loc+piv+1<<endl;
            loc+=piv+1;
        }
    }
    else if (tipo_prob==2) {
        auxiliar1 << "FontesLeitDetAdjunto1D_MultiGrupo.txt"; // cria o novo nome parte vari?vel mais a parte fixa entre aspas
        nomedoarquivo1 = auxiliar1.str();
        arch_saida1.open(nomedoarquivo1.c_str());
        ///////////////////////////////////////////
        arch_saida1<<"////Fontes Fisicas para a Leitura do Detector Problema Adjunto SGF"<<endl;
        for (int g=1; g<G+1; g++){
            arch_saida1<<"//g= "<<g<<endl;
            for (int i=1; i<NR+1; i++)
                arch_saida1<<"0"<<" ";
            arch_saida1<<endl;
        }
    }
    ///**************************************///
    ofstream arch_saida2;  // Para arquivos de saída
    string   nomedoarquivo2;
    ostringstream  auxiliar2;
    auxiliar2 << "Dimensoes_1D2D_NumGrupos.txt"; // cria o novo nome parte vari?vel mais a parte fixa entre aspas
    nomedoarquivo2 = auxiliar2.str();
    arch_saida2.open(nomedoarquivo2.c_str());
    ///////////////////////////////////////////
    arch_saida2<<"//Dimensoes (1-1D e 2-2D)"<<endl;
    arch_saida2<<"1"<<endl;
    arch_saida2<<"//Número de Grupos de Energia"<<endl;
    arch_saida2<<G<<endl;
    arch_saida2<<"//Grau da Anisotropia do Espalhamento (L)"<<endl;
    arch_saida2<<Li<<endl;
    ///**************************************///
    ///**************************************///
    ofstream arch_saida;  // Para arquivos de saída
    string   nomedoarquivo;
    ostringstream  auxiliar;
    auxiliar << "Dados_1D_MultiGrupo.txt"; // cria o novo nome parte vari?vel mais a parte fixa entre aspas
    nomedoarquivo = auxiliar.str();
    arch_saida.open(nomedoarquivo.c_str());
    ///*****************************************///
    arch_saida<<"//Imprimir Matrizes Autovalores, Autovalores/Autovetores, Mat Theta, Mat GG (1-Sim, 0-Nao)"<<endl;
    arch_saida<<"0 0 0 0"<<endl;
    arch_saida<<"//Tipo de problema (1-Físico e 2-Adjunto)"<<endl;
    arch_saida<<tipo_prob<<endl;
    arch_saida<<"//Método(1-DD e 2-SGF) Deslocamento (1-Sem desloc)e 2-Com desloc)"<<endl;
    if (tipo_met==1)
        arch_saida<<"1 2"<<endl;
    else if (tipo_met==2)
        arch_saida<<"2 2"<<endl;
    arch_saida<<"//Ordem da Quadratura GL"<<endl;
    arch_saida<<"16"<<endl;
    arch_saida<<"//NR e NZ"<<endl;
    arch_saida<<NR<<" "<<NZ<<endl;
    arch_saida<<"//Dados por Regiao"<<endl;
    arch_saida<<"//Comprimento"<<endl;
    for (int i=1; i<NR+1; i++)
        arch_saida<<i+1<<" ";
    arch_saida<<endl;
    arch_saida<<"//Nodos por Regiao"<<endl;
    if (tipo_met==1) {
        for (int i=1; i<NR+1; i++)
            arch_saida<<(i+1)*10<<" ";
        arch_saida<<endl;
        arch_saida<<"//Periodicidade"<<endl;
        arch_saida<<"1"<<endl;
    }
    else if (tipo_met==2) {
        for (int i=1; i<NR+1; i++)
            arch_saida<<"1"<<" ";
        arch_saida<<endl;
        arch_saida<<"//Periodicidade"<<endl;
        arch_saida<<"100"<<endl;
    }
    ///Mapeamento (aqui fixamos o mapeamento igual a 1 em cada regiao, tem q ser mudado para rodar problema)
    arch_saida<<"//Mapeamento"<<endl;
    for (int i=1; i<NR+1; i++)
        arch_saida<<i<<" ";
    arch_saida<<endl;
    arch_saida<<"//CP e-?"<<endl;
    arch_saida<<"6"<<endl;
    arch_saida<<"//Grau da Anisotropia do Espalhamento (L)"<<endl;
    arch_saida<<Li<<endl;
    ///Dados Por Zonas
    double XT, XSSl, Hgg1, d5, d10;
    //arch_saida<<"*****************************"<<endl;
    arch_saida<<"//Dados por Zona"<<endl;
    int loc=0;
    for (int piv=1; piv<NZ+1; piv++) {
        arch_saida<<"///Zona "<<piv<<endl;
        arch_saida<<"//Sigma Total g (g=1.."<<G<<")"<<endl;
        for (int g=1; g<G+1; g++){
            if (g==5)
                d5=1;
            else d5=0;
            if (g==10)
                d10=1;
            else d10=0;
            XT=pow((piv+20.0)/21.0, 5)*(g/10.0-0.15*d5-0.15*d10);
            arch_saida<<XT<<" ";
        }
        arch_saida<<endl;
        //////////////////////////////////
        for (int l=0; l<Li+1; l++){
            arch_saida<<"//Sigma Espalhamento "<<l<<" g'g (g linha; g' coluna)"<<endl;
            for (int g=1; g<G+1; g++){
                for (int g1=1; g1<G+1; g1++) {
                    Hgg1=0.7-((g+g1)/200.0);
                    if (g1>g)
                        XSSl=0;     ///XSSl=(piv+20)/21.0*(g/(100.0*(g1-g+1)))*pow(Hgg1,l);
                    else
                        XSSl=(piv+20)/21.0*(g1/(100.0*(g-g1+1)))*pow(Hgg1,l);  ///(2*l+1)*
                    arch_saida<<XSSl<<" ";
                }
                arch_saida<<endl;
            }
        }
        loc+=piv+1;
    }
    arch_saida<<"//Tipo de Condicoes de Contorno (Esq, Dir) (1-Prescrita, 2-Reflexiva)"<<endl;
    arch_saida<<"1 1"<<endl;
    ///Fontes Físicas ou Sigma do Detector
    arch_saida<<"//Valor da condicao de contorno prescrita (Esq, Dir) "<<endl;
     for (int g=1; g<G+1; g++){
        arch_saida<<"//g= "<<g<<endl;
        if (g==1)
            arch_saida<<"1"<<" "<<"0"<<endl;
        else
            arch_saida<<"0"<<" "<<"0"<<endl;
    }
    ///Fontes Físicas ou Sigma do Detector
    //arch_saida<<"*****************************"<<endl;
    if (tipo_prob==1) {
        arch_saida<<"//Fonte Fisica"<<endl;
         for (int g=1; g<G+1; g++){
            arch_saida<<"//g= "<<g<<endl;
            for (int i=1; i<NR+1; i++)
                arch_saida<<"0"<<" ";
            arch_saida<<endl;
        }
    }
    else if (tipo_prob==2) {
        arch_saida<<"//Fonte Adjunta (Sigma de Absorcao do Detector)"<<endl;
         for (int g=1; g<G+1; g++){
            arch_saida<<"//g= "<<g<<endl;
            for (int piv=1; piv<NR+1; piv++) {
                    double somaXS0=0;
                    for (int g1=1; g1<G+1; g1++) {
                        Hgg1=0.7-((g+g1)/200.0);
                        //system("pause");
                        if (g>g1)
                            XSSl=0;     ///XSSl=(piv+20)/21.0*(g/(100.0*(g1-g+1)))*pow(Hgg1,l);
                        else
                            XSSl=(piv+20)/21.0*(g/(100.0*(g1-g+1)))*pow(Hgg1,0);  ///(2*l+1)*
        //                cout<<Hgg1<<"  "<<l<<"  "<<pow(Hgg1,l) <<endl;
        //                cout<<(i+20)/21.0<<"  "<<100.0*(g-g1+1)<<"  "<<g1/(100.0*(g-g1+1)) <<endl;
                        somaXS0+=XSSl;
                    }
                    if (g==5)
                    d5=1;
                    else d5=0;
                    if (g==10)
                        d10=1;
                    else d10=0;
                    XT=pow((piv+20.0)/21.0, 5)*(g/10.0-0.15*d5-0.15*d10);
                    arch_saida<<XT-somaXS0<<" ";
            }
            arch_saida<<endl;
        }
    }

    //cout<<NZ<<"\t"<<Li<<"\t"<<G<<"\t"<<endl;
    arch_saida.close();
    arch_saida1.close();
    arch_saida2.close();

    ///**************************************************///
    ///Evaluando Polinomios de Legendre TESTE
//    double PolLeg;
//    int n=8;
//    double x=0.4;
//    PolLeg=Pn(n, x) ;
//    cout<<"P"<<n<<"("<<x<<")= "<<PolLeg<<endl;

    ///**************************************************///
    ///**************************************************///






    return 0;
}
