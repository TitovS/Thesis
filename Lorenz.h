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

    // Задание параметров системы
        // Параметры:
        // double - b
        // double - sigma
        // double - r
        void SetParam (double, double, double);

    // Вычисление n точек фазовой траектории
        // Параметры:
        // double* - массив начальных условий
        // double - шаг вычисления
        // int - количество вычисляемых точек
        // double - шаг рештки
        // int - шаг на решетке
        // bool - получение траектории (true) или получение циклов(false)
        void GetTr (double*, double, double, int, bool = false);


    // Сохранение траектории в файл
        // Параметры:
        // int - начальная точка
        // int - интервал пропуска точек
        // int - конечная точка
        void Save (int, int, int);

    // Функция составления уравнения прямой для свдига
        void GetLine();

    // Динамическое увеличение массива
        void DynamicMemory(int, int&, double*&, double*&);

    // Деконструктор
        ~Lorenz ();

    // Статистика системы
        Stat S;

    // Пареметры уравнения свдига
        double* vector; //направляющий вектор
        double* dot; //точка на прямой

    // размерность системы
        int m_dim;


private:
//____________________________________

    //Система Лоренца

        //параметры
        double m_b;
        double m_sigma;
        double m_r;

        //показатели перехода по циклам
        int Cycles_done; //переходов совершено
        int num_main; //сквозная нумерация
        int num_first; //номер начальной точки вычисления
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

        // Функция ЦПВД
        // double*& - текущая точка
        // double шаг решетки
        void GridTrCPVD (double*&,double*& );

        // Фунция обынчая
        void GridTr (double*&);

        // Размер решетки
        double a;

        // Шаг на решетке
        int b;

//____________________________________

    // Хэш-функция циклов

        int HashFun(double *);

    // Проверка циклов

        // Параметры
        // double* - полученная точка
        void CycleCheck(double *);

        // хэш-массив
        std::list<std::pair<double*,int>>* m_hash;

        // Размер хэш-функции
        int h_n;

        // Хэш значения
        int p1;
        int p2;
        int p3;

//____________________________________



};



#endif //CLASSES_CPOINT_H