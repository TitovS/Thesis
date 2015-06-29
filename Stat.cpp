//
// Created by Сергей Титов on 26/06/15.
//

#include "Stat.h"
#include "iostream"
#include "math.h"

void Stat::collect(double* _np, double* _cp, int m_dim) {

double x = 0;

    while (_cp != _np){



        for (int i = 0; i < m_dim; ++i) {

          x += (_cp[i] - ((_cp+1)[i]))*(_cp[i] - ((_cp+1)[i]));

        }

        l_cycle[n_cycle]+= sqrt(x);
        s_cycle[n_cycle]+= 1;
        x = 0;
        _cp += m_dim;

    }


    n_cycle +=1; // счетчик кол - ва циклов
}

void Stat::save() {

}