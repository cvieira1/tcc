#ifndef GNUPLOT_H_INCLUDED
#define GNUPLOT_H_INCLUDED
#include <iostream>
#include <string>

using namespace std;

class gnuplot{
public:
    gnuplot();
    ~gnuplot();
    void operator()(const string & command);
protected:
    FILE *gnuplotpipe;
};

gnuplot::gnuplot(){
    gnuplotpipe = popen("gnuplot -persistent","w");
    if(!gnuplotpipe){
        cerr << ("Gnuplot not found !");
    }
}

gnuplot::~gnuplot(){
    fprintf(gnuplotpipe,"exit\n");
    pclose(gnuplotpipe);
}

void gnuplot::operator()(const string & command){
    fprintf(gnuplotpipe,"%s\n",command.c_str());
    fflush(gnuplotpipe);
}


#endif // GNUPLOT_H_INCLUDED
