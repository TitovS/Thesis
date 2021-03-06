//
// Created by Сергей Титов on 26/06/15.
//

#include "Stat.h"
#include "iostream"
#include "math.h"
#include <fstream>
#include <iomanip>


Stat::Stat(void) {
    n_cycle = 0;
    u_cycle = 0;
    l_cycle = 0;
    s_cycle = 0;

    char const *pchar = "MainTable.txt";
    std::ofstream out;
    out.open(pchar, std::ofstream::out | std::ofstream::trunc);
    out.close();
}

void Stat::collect(double* _cp, int m_dim, int num_cycle, double a_index, bool mode) {


    if (!mode) {
        l_cycle = 0;
        s_cycle = 0;

        //std::ofstream out;
        //std::string a = std::to_string(a_index);
        //std::string s = std::to_string(u_cycle);
        //s = a + "_" + s + ".txt";
        //char const *pchar = s.c_str();
        //out.open(pchar, std::ofstream::out | std::ofstream::trunc);

        double x = 0; // длинна между точками цикла

        //Проходим по массиву траекторий
        for (int j = 0; j < num_cycle + 1; ++j) {

            for (int i = 0; i < m_dim; ++i) {

                x += (_cp[i] - _cp[i + m_dim]) * (_cp[i] - _cp[i + m_dim]);
                //out << _cp[i] << " "; //записываем в файл
            }
            //out << "\n";

            l_cycle += sqrt(x); //длина
            s_cycle += 1; //кол-во точек в цикле
            x = 0;

            _cp += m_dim; // переход к следующей точке


        }
        //out << "\n";

        //out.close();

        //Запись в файл

        char const *pchar2 = "MainTable.txt";
        std::ofstream out2;
        out2.open(pchar2, std::ofstream::out | std::ofstream::app);
        out2 << std::setprecision(9) << a_index << ' ' << s_cycle <<' ' << l_cycle <<"\n";
        out2.close();

    }


    u_cycle += 1; // счетчик кол - ва циклов


}

void Stat::Reset() {

    u_cycle = 0;
    n_cycle = 0;


}

Stat::~Stat(){



}
