//
// Created by Сергей Титов on 10/05/15.
//


#ifndef THESIS_LORENZ_H
#define THESIS_LORENZ_H

#include "map"
#include "list"
#include "Stat.h"
#include "Grid.h"

class Lorenz
{
public:

    // конструктор
    Lorenz (void);

    // Вычисление всех циклов
        // Stat* - класс статистики
        // double - шаг решетки
        void GetCycles (Stat*,double);

    // Поиск точки разрыва
        // Параметры
        // Stat* - класс статистики
        // Grid* - класс решетки
        void GetBreak (Stat*, Grid*);

    // Сохранение траектории в файл
        void Save ();

    // Обнуление перменных
        void Reset();

    // Деконструктор
        ~Lorenz ();


private:
//_____________Система________________

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

        // Пареметры уравнения свдига
        double* vector; //направляющий вектор
        double* dot; //точка на прямой

        // Размерность системы
        int const m_dim;

        // Нумерация
        int num_main; //сквозная нумерация
        int num_first;//номер начальной точки вычисления

        // Маркер поиска разрывов
        bool BreakpointSearchMode;


//_____________Вычисление________________

    // Вычисление n точек фазовой траектории

        // Параметры:
        // double* - массив начальных условий
        // double - шаг решетки
        // int - шаг на решетке
        // Stat* - класс статистики
        void GetTr (double*, double, int, Stat*);

    // Метод Рунге-Кутта

        // Функция
        // double*& - текущая точка
        // double*& - следующая точка
        void Runge_Cutta (double*& , double*&);

        // шаг алгоритма
        double m_dt;

   // Вычисление фазовой Скорости в заданной точке

        // Функция:
        // double*& - текущая точка
        void PhVelocity (double*&);

        // вектор фазовой скорости
        double* m_v;

//_____________Дискретизация________________

    // Метод пространтсвенно-временой дискретизации

        // Функция ЦПВД
        // double*& - текущая точка
        // double шаг решетки
        //void GridTrCPVD (double*&,double*& );

        // Фунция обынчая
        void GridTr (double*&);

        // Размер решетки
        double a;

//_____________Хэш-функция________________

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


//_____________Дополнительное________________

    // Динамическое увеличение массива
    void DynamicMemory(int, int&, double*&, double*&);

    // Функция составления уравнения прямой для свдига
    void GetLine();


};



#endif //CLASSES_CPOINT_H