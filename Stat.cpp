//
// Created by Сергей Титов on 26/06/15.
//

#include "Stat.h"
#include "iostream"
#include "math.h"
#include <fstream>

void Stat::collect(double* _np, double* _cp, int m_dim) {

    int check = 0; //
    n_cycle += 1; // счетчик кол - ва циклов


    //Запись первого цикла
    if (firsts.size() == 0)
    {
        for (int i = 0; i < m_dim; ++i) {
            firsts.push_back(_np[i]);
        }
    }


    //создание и помещение итератора на первую позицию
    std::list<double>::iterator it;
    it = firsts.begin();

    while (it != firsts.end()) {

        for (int i = 0; i < m_dim; ++i) { //сравниваем первые 3 значения

            if (*it != _np[i]) { //если координаты не свпадают то
                check += 1; //отмечаем это
                std::advance(it,m_dim-i); //перемешаем итератор к следующей тройке
                break; //заканчиваем  с этой тройкой
            }
            else {
                it++; //переход к следующей координате
            }

        }


    if (check == firsts.size()/m_dim) //если координаты не совпали ни с одним из циклов то сохраняем его.
    {
        for (int i = 0; i < m_dim; ++i) {
            firsts.push_back(_np[i]);
        }

        double x = 0; // длинна между точками цикла
        u_cycle += 1;

        //Создаем файл с номером цикла
        std::ofstream out;
        std::string s = std::to_string(u_cycle);
        s += ".txt";
        char const *pchar = s.c_str();
        out.open(pchar, std::ofstream::out | std::ofstream::trunc);

        //Проходим по массиву траекторий
        while (_cp != _np) {


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

        out.close();
        break;
    }

    }
}






