#include <iostream>
#include <cmath>
#include <ctime>
#include <omp.h>

using namespace std;

// constants

const int loopsize = 20; // 20 million

const short int loopcount = 10;

// variables

int a[loopsize],b[loopsize],c[loopsize];


void doCalcs() {

#pragma omp parallel for

for ( int i=0; i < loopsize;i++) {

b[i] = (rand() % loopsize);

c[i] = rand() % loopsize;

}

#pragma omp parallel for

for (int i=0; i< loopsize;i++) {

a[i] = b[i]*c[i];

}

}

void subdomain(float *x, int iStart, int iPoints){
    int i;
    for(i = 0;i < iPoints;i++){
        x[iStart + i] = 123.456;
    }
}

void sub(float *x, int nPoints){
    int iam,nt,iPoints,iStart;
    #pragma omp parallel default(shared) private(iam,nt,iPoints,iStart)
    {
        iam = omp_get_thread_num();
        nt = omp_get_num_threads();
        iPoints = nPoints / nt;
        iStart = iam * iPoints;
        if(iam == (nt - 1)){
            iPoints = nPoints - iStart;
        }
        subdomain(x, iStart, iPoints);
    }
}
int main()
{
    #pragma omp parallel
    {
        #pragma omp single
        printf("Beginning work 1");
        #pragma omp single
        printf("Finishing work 1");
        #pragma omp single nowait
        printf("finished work 1");
    }
//    double startTime,endTime;
//    clock_t startTime2,endTime2;
    /*const int size = 256;
    double sinTable[size];
    #ifdef _OPENMP
    startTime = omp_get_wtime();
    #endif
    #pragma omp parallel for
    for(int n = 0;n < size;n++){
        sinTable[n] = sin(2 * M_PI * n / size);
    }
    doCalcs();
    #ifdef _OPENMP
    endTime = omp_get_wtime();
    #endif
    cout << (endTime - startTime) << endl;*/
/*omp_set_num_threads(2);

double startTime;
    #ifdef _OPENMP
    startTime = omp_get_wtime();
    #endif


double times = 0;

double endTime;

#pragma omp parallel for

for (int loop=0;loop < loopcount;loop++) {

    #ifdef _OPENMP
    startTime = omp_get_wtime();
    #endif

doCalcs();

    #ifdef _OPENMP
    endTime = omp_get_wtime();
    #endif

if (loop>0) {

std::cout << (endTime - startTime) << std::endl;
times += (endTime - startTime);
}

}
  cout << "Total: " << times << endl;*/
    return 0;

}
