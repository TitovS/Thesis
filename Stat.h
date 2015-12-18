//
// Created by Сергей Титов on 26/06/15.
//

#ifndef THESIS_STAT_H
#define THESIS_STAT_H


#include "list"

class Stat {
public:
//Статистики циклов

    //Массив длин циклов
    double* l_cycle;

    //Массив шагов циклов
    int* s_cycle;

    //Количество циклов
    int n_cycle;

    //Уникальных циклов
    int u_cycle;

//Функции

    void SInit();

    //Сбор статистики о цикле
    void collect(double*, int, int, int);

    void Reset();

    ~Stat();
};


#endif //THESIS_STAT_H
