//
// Created by Сергей Титов on 10/05/15.
//


#ifndef THESIS_LORENZ_H
#define THESIS_LORENZ_H

#include "map"
#include "list"
#include "Stat.h"

class Lorenz
{
public:

    // конструктор
    Lorenz (void);

    // задание параметров системы
    // Параметры:
    // double - b
    // double - sigma
    // double - r
    void SetParam (double, double, double);

    // вычисление n точек фазовой траектории
    // Параметры:
    // double* - массив начальных условий
    // double - шаг вычисления
    // int - количество вычисляемых точек
    void GetTr (double*, double, int, double);

    // сохранение траектории в файл
    // Параметры:
    // int - начальная точка
    // int - интервал пропуска точек
    // int - конечная точка
    void Save (char*, int, int, int);

    // деконструктор
    ~Lorenz ();

private:
//____________________________________

    //Система Лоренца

        //параметры
        double m_b;
        double m_sigma;
        double m_r;

        // размерность системы
        int m_dim;

        // Статистика системы

        Stat S;

//____________________________________

    // Метод Рунге-Кутта

        // Функция
        // double*& - текущая точка
        // double*& - следующая точка
        void Runge_Cutta (double*& , double*&);

        // массив точек фазовой траектории
        double* m_tr;

        // количество точек в траектории
        int m_n;

        // шаг алгоритма
        double m_dt;

//____________________________________

   // Вычисление фазовой Скорости в заданной точке

        // Функция:
        // double*& - текущая точка
        void PhVelocity (double*&);

        // вектор фазовой скорости
        double* m_v;

//____________________________________

    // Метод пространтсвенно-временой дискретизации

        // Функция
        // double*& - текущая точка
        // double шаг решетки
        void GridTr (double*&);

        // Шаг решетки
        double a;

//____________________________________

    // Хэш-функция и проверка циклов

        // Параметры
        // double*& - полученная точка
        //
        void HashFun (double*&);

        // хэш-массив
        std::map<int,std::list<double*>> m_hash;

        // Размер хэш-функции
        int h_n;

        // Хэш значения
        int p1;
        int p2;
        int p3;

//____________________________________

};



#endif //CLASSES_CPOINT_H