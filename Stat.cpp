//
// Created by Сергей Титов on 26/06/15.
//

#include "Stat.h"
#include "iostream"


void Stat::collect(double* _np, double* _cp, int m_dim) {



    while (_cp != _np){

        _cp += m_dim;
        s_cycle[n_cycle]+=1;
    }


    n_cycle +=1; // счетчик кол - ва циклов
}

void Stat::save() {

}