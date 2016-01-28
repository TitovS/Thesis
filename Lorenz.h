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

    // Вычисление n точек фазовой траектории
        // Параметры:
        // double* - массив начальных условий
        // double - шаг решетки
        // int - шаг на решетке
        void GetTr (double*, double, int, Stat*);

    // Вычисление всех циклов
        // Stat* - статискика
        // double - шаг решетки
        void GetCycles (Stat*,double);

    // Поикс точки разрыва
        // Параметры
        // Stat* - класс статистики
        // double - левая точка
        // double - правая точка точка
        // double - допустимая погрешность
        void GetBreak (Stat*, double, double, double, int);
        bool BreakpointSerachMode;
        double* a_breakpoints;
    // Сохранение траектории в файл
        // Параметры:
        // Stat* - класс статистики
        //
        void Save (Stat*);

    // Функция составления уравнения прямой для свдига
        void GetLine();

    // Динамическое увеличение массива
        void DynamicMemory(int, int&, double*&, double*&);

    // Обнуление перменных
        void Reset();

    // Деконструктор
        ~Lorenz ();

    // Пареметры уравнения свдига
        double* vector; //направляющий вектор
        double* dot; //точка на прямой

    // Размерность системы
        int const m_dim;

    // Нумерация
        int num_main; //сквозная нумерация
        int num_first;//номер начальной точки вычисления

private:
//____________________________________

    //Система Лоренца

        //параметры
        double const m_b;
        double const m_sigma;
        double const m_r;

        //показатели перехода по циклам
        int Cycles_done; //переходов совершено

        // массив точек фазовой траектории
        double* m_tr;
        double* m_tr1;
        double* m_tr2;
        int k;

        // количество точек в траектории
        int m_n;

//____________________________________

    // Метод Рунге-Кутта

        // Функция
        // double*& - текущая точка
        // double*& - следующая точка
        void Runge_Cutta (double*& , double*&);

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
        //void GridTrCPVD (double*&,double*& );

        // Фунция обынчая
        void GridTr (double*&);

        // Размер решетки
        double a;

//____________________________________

    // Хэш-функция циклов

        int HashFun(double *);

    // Проверка циклов

        // Параметры
        // double* - полученная точка
        // int - счетчик свдига по а
        void CycleCheck(double *, Stat*);

        // хэш-массив
        std::list<int>* m_hash;

        // Размер хэш-функции
        int h_n;

        // Хэш значения
        int const p1;
        int const p2;
        int const p3;

//____________________________________
};



#endif //CLASSES_CPOINT_H