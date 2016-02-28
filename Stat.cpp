//
// Created by Сергей Титов on 26/06/15.
//

#include "Stat.h"
#include "iostream"
#include "math.h"
#include <fstream>


Stat::Stat(void) {
    n_cycle = 0;
    u_cycle = 0;
    l_cycle = new double [100000]; //TODO Оформить динмаичесские массивы
    s_cycle = new int [100000];
}

void Stat::collect(double* _cp, int m_dim, int num_cycle, double a_index, bool mode) {

    if (mode == false) {
        double x = 0; // длинна между точками цикла
        //Создаем файл с номером цикла
        std::ofstream out;
        std::string a = std::to_string(a_index);
        std::string s = std::to_string(u_cycle);
        s = a + "_" + s + ".txt";
        char const *pchar = s.c_str();
        out.open(pchar, std::ofstream::out | std::ofstream::trunc);

        //Проходим по массиву траекторий
        for (int j = 0; j < num_cycle + 1; ++j) {

            for (int i = 0; i < m_dim; ++i) {

                x += (_cp[i] - _cp[i + m_dim]) * (_cp[i] - _cp[i + m_dim]);
                out << _cp[i] << " "; //записываем в файл

            }
            out << "\n";

            l_cycle[u_cycle] += sqrt(x); //длина
            s_cycle[u_cycle] += 1; //кол-во точек в цикле
            x = 0;

            _cp += m_dim; // переход к следующей точке
        }
        out << "\n";
    }

    u_cycle += 1; // счетчик кол - ва циклов


}

void Stat::Reset() {

    u_cycle = 0;
    n_cycle = 0;


}

Stat::~Stat(){

    delete[] l_cycle;
    delete[] s_cycle;

}
