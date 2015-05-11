#include "Lorenz.h"
#include <iostream>
#include <fstream>
#include <math.h>

// конструктор
Lorenz::Lorenz (void)
{
    // размерность системы
    m_dim = 3;

    // вектор фазовой скорости
    m_v = new double [m_dim];

    // количество точек в траектории
    m_n = 0;
}


// задание параметров системы
void Lorenz::SetParam (double _p1, double _p2, double _p3)
{
    m_b = _p1;
    m_sigma = _p2;
    m_r = _p3;
}


// вычисление фазовой траектории
void Lorenz::GetTr (double* _inc, double _dt, int _n)
{
    // шаг метода
    m_dt = _dt;

    // удаление старого массива точек траектории
    if (m_n != 0)
        delete [] m_tr;

    // количество вычисляемых точек
    m_n = _n;

    // создание массива точек фазовой траектории
    double* m_tr = new (std::nothrow) double [m_dim * m_n];
    double* tp;//Я не знаю в чем дело, но работает только так.
    // запись начальных условий
    for (int i = 0; i < m_dim; i++)
        m_tr [i] = _inc [i];

    // вычисление траектории
    for (int i = 0; i < m_n - 1; i++)
    {
        // указатель на текущую точку
        double* cp = m_tr + i * m_dim;
        // указатель на следующую точку
        double* np = m_tr + (i + 1) * m_dim;

        // вычисление следующей точки методом Рунге-Кутта

        //Runge_Cutta(np, cp, tp);
        Euler(np,cp);
    }
}


// вычисление следующей точки фазовой траектории по текущей Эйлер
void Lorenz::Euler (double*& _np, double*& _cp)
{
    // вычисление фазовой скорости
    PhVelocity (_cp);

    // метод Эйлера
    for (int i = 0; i < m_dim; i++)
        _np [i] = _cp [i] + m_v [i] * m_dt;

}

// вычисление следующей точки фазовой траектории по текущей Рунге-Кутта
void Lorenz::Runge_Cutta (double*& _np, double*& _cp, double*& _tp)
{
    PhVelocity (_cp);

    double n[4*m_dim];
    // n1
    for (int i=0; i< m_dim; i++) {
        n[i]=m_v[i];
    }

    // n2
    for (int i=0; i< m_dim; i++) {  //
        _tp[i] =_cp[i] + n[i]*m_dt/2;
    }
    PhVelocity(_tp);
    for (int i=0; i< m_dim; i++) {
        n[i*2] = m_v[i];
    }

    // n3
    for (int i=0; i< m_dim; i++) {  //
        _tp[i] =_cp[i] + n[i * 2]*m_dt/2;
    }
    PhVelocity(_tp);
    for (int i=0; i< m_dim; i++) {
        n[i*3] = m_v[i];
    }

    // n4
    for (int i=0; i< m_dim; i++) {  //
        _tp[i] =_cp[i] + n[i * 3]*m_dt;
    }
    PhVelocity(_tp);
    for (int i=0; i< m_dim; i++) {
        n[i*4] = m_v[i];
    }



    for (int i=0; i< m_dim; i++) {
        _np[i]=_cp[i]+m_dt/6*(n[i]+2*n[i*2]+2*n[i*3]+n[i*4])  ;
    }

}

// вычисление фазовой скорости в заданной точке
void Lorenz::PhVelocity (double*& _cp)

{

    m_v [0] = m_sigma * (_cp [1] - _cp [0]);
    m_v [1] = (m_r - _cp [2]) * _cp [0] - _cp [1];
    m_v [2] = _cp [0] * _cp [1] - m_b * _cp [2];

}

// дикретизация полученной траектории
void Lorenz::GridTr (double*& _cp, double a)
{
    for (int i = 0; i<m_dim; i++) {
        _cp[i] = round(_cp[i]/a)*a +a/2;
    }


}




// сохранение траектории в файл
void Lorenz::Save (char* _fln,int _be, int _dn, int _en)
{
    // создание файлового потока вывода
    std::ofstream out;
    // связывание потока с файлом
    out.open (_fln);

    // вывод результата
    out.precision (15); // задание количества значащих цифр в выводимых числах
    for (int i = _be; i < _en; i += _dn * m_dim)
    {
        for (int j = 0; j < m_dim; j++)
        { std::cout << m_tr [i * m_dim + j] << " ";
        }
        std::cout << "\n";
    }

    // закрытие файла
    out.close ();
}


// деструктор
Lorenz::~Lorenz ()
{
    delete [] m_tr;

    delete [] m_v;
}
