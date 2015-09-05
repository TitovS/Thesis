//
// Created by Сергей Титов on 26/06/15.
//

#include "Stat.h"
#include "iostream"
#include "math.h"
#include <fstream>

void Stat::collect(double* _np, double* _cp, int m_dim, int num_cycle) {



    double x = 0; // длинна между точками цикла


    //Создаем файл с номером цикла
    std::ofstream out;
    std::string s = std::to_string(u_cycle);
    s += ".txt";
    char const *pchar = s.c_str();
    out.open(pchar, std::ofstream::out | std::ofstream::trunc);

    //Проходим по массиву траекторий
    for (int j = 0; j < num_cycle; ++j){


        for (int i = 0; i < m_dim; ++i) {

            x += (_cp[i] - _cp[i + m_dim]) * (_cp[i] - _cp[i + m_dim]);
            out << _cp[i] << " ";

        }
        out << "\n";
        l_cycle[u_cycle] += sqrt(x);
        s_cycle[u_cycle] += 1;
        x = 0;

        _cp += m_dim;
    }
    u_cycle += 1; // счетчик кол - ва циклов

}






