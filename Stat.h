//
// Created by Сергей Титов on 26/06/15.
//

#ifndef THESIS_STAT_H
#define THESIS_STAT_H


class Stat {
public:
//Статистики циклов

    //Массив длин циклов
    double* l_cycle;

    //Массив шагов циклов
    int* s_cycle;


//Функции

    //Сбор статистики о цикле
    void collect(double*, double*, int);

    //Сохранение статистики в файле
    void save();


    //Количество циклов
    int n_cycle;
};


#endif //THESIS_STAT_H
