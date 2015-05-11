//
// Created by Сергей Титов on 10/05/15.
//
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
    void GetTr (double*, double, int);

    // сохранение траектории в файл
    // Параметры:
    // int - начальная точка
    // int - интервал пропуска точек
    // int - конечная точка
    void Save (char*, int, int, int);

    // деконструктор
    ~Lorenz ();

private:

    // метод Эйлера
    // Параметры:
    // double*& - текущая точка
    // double*& - следующая точка
    void Euler (double*&, double*&);


    // Метод Рунге-Кутта
    // Параметры
    // double*& - текущая точка
    // double*& - следующая точка
    void Runge_Cutta (double*& , double*&, double*&);

    // вычисление фазовой скорости в заданной точке
    // Параметр:
    // double*& - текущая точка
    void PhVelocity (double*&);

    // Метод пространтсвенно-временой дискретизации
    // Параметры
    // double*& - текущая точка
    // double шаг решетки
    void GridTr (double*&, double);


    // массив точек фазовой траектории
    double* m_tr;

    // вектор фазовой скорости
    double* m_v;

    // количество точек в траектории
    int m_n;

    // размерность системы Лоренца
    int m_dim;

    // шаг алгоритма
    double m_dt;

    // параметры системы Лоренца
    double m_b;
    double m_sigma;
    double m_r;

    // текущая точка
    int m_k;
};
