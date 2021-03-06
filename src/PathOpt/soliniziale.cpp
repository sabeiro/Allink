#include "stopwatch.hpp"
#include "matrix.h"
#include <iostream>
#include <vector>
#include <stdlib.h>

int main (){
    std::cout << "Hello World!";
    //leggo i dati
    //determino il numero delle infermiere, dei servizi, le distanze e i tempi
    unsigned ns=7;
    std::vector<unsigned> tempi={3,4,2,6,7,8,9};
    unsigned L=30; //costo massimo per tour
    //matrice delle distanze
    Matrix d(ns+1,ns+1,0);
    std::cout << "d";
    for (unsigned i = 0; i < ns; ++i){
	for (unsigned j = 0; j < ns; ++j){
	    d[i][j]=(rand() % 11);
	    if(i==j){d[i][j]=0;}
	}
    }
    stopwatch sw;
    sw.reset();
    sw.start();
    //calcolo la matrice dei costi
    Matrix costi(ns,ns,0);
    costi.riempi_cost(d,tempi);
    std::cout << "costi";
    Matrix tours(ns+1,ns,0);
    tours.sol_iniziale(costi,L);
    sw.stop();
    std::cout<<"Tempo: "<< sw.elapsed_ms()<<" ms "<<std::endl;
    return 0;
}

