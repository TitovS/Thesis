//
// Created by Сергей Титов on 26/06/15.
//

#include "Stat.h"
#include "iostream"
#include "math.h"
#include <fstream>

void Stat::collect(double* _np, double* _cp, int m_dim) {


    double x = 0; // длинна между точками цикла
    n_cycle +=1; // счетчик кол - ва циклов

    std::ofstream out;
    std::string s = std::to_string(n_cycle);
    s += ".txt";
    char const *pchar = s.c_str();
    out.open(pchar, std::ofstream::out | std::ofstream::trunc);

    while (_cp != _np){

        for (int i = 0; i < m_dim; ++i) {

            x += (_cp[i] - _cp[i+m_dim])*(_cp[i] - _cp[i+m_dim]);
            out << _cp[i] << " ";

        }
        out  << "\n";
        l_cycle[n_cycle] += sqrt(x);
        s_cycle[n_cycle] += 1;
        x = 0;


        _cp +=m_dim;

    }

    out.close();

}
