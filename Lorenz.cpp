#include "Lorenz.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <set>

// конструктор
Lorenz::Lorenz (void)
{
    // размерность системы
    m_dim = 3;

    // вектор фазовой скорости
    m_v = new double [m_dim];

    // количество точек в траектории
    m_n = 0;

    m_tr=0;

    // хэш-значения
    p1=73856093;
    p2=19349663;
    p3=83492791;

    h_n=10000;
}


// задание параметров системы
void Lorenz::SetParam (double _p1, double _p2, double _p3)
{
    m_b = _p1;
    m_sigma = _p2;
    m_r = _p3;
}


// вычисление фазовой траектории
void Lorenz::GetTr (double* _init, double _dt, int _n, double _a)
{
    // шаг метода
    m_dt = _dt;

    // количество вычисляемых точек
    m_n = _n;

    // Шаг решетки
    a=_a;

    // создание массива точек фазовой траектории
    m_tr = new (std::nothrow) double [m_dim * m_n];

    // создание хэш массива


    // запись начальных условий
    for (int i = 0; i < m_dim; i++)
        m_tr [i] = _init [i];

    // вычисление траектории
    for (int i = 0; i < m_n - 1; i++)
    {
        // указатель на текущую точку
        double* cp = m_tr + i * m_dim;
        // указатель на следующую точку
        double* np = m_tr + (i + 1) * m_dim;

        // вычисление следующей точки методом Рунге-Кутта
        Runge_Cutta(np, cp);

        //Дискретизация полученной траектории
        GridTr(np);

        //Хэширование полученных траекторий
        HashFun(np);
    }

}


// вычисление следующей точки фазовой траектории по текущей - Эйлер
void Lorenz::Euler (double*& _np, double*& _cp)
{
    // вычисление фазовой скорости
    PhVelocity (_cp);

    // метод Эйлера
    for (int i = 0; i < m_dim; i++)
        _np [i] = _cp[i] + m_v [i] * m_dt;

}

// вычисление следующей точки фазовой траектории по текущей - Рунге-Кутта
void Lorenz::Runge_Cutta (double*& _np, double*& _cp )
{
    PhVelocity (_cp);

    double* tp = new double [m_dim];
    double n[12];
    // n1
    for (int i=0; i< m_dim; i++) {
        n[i]=m_v[i];
    }

    // n2
    for (int i=0; i< m_dim; i++) {
        tp[i] = _cp[i] + n[i]*m_dt/2;

    }
    PhVelocity(tp);
    for (int i=0; i< m_dim; i++) {
        n[i+m_dim] = m_v[i];
    }

    // n3
    for (int i=0; i< m_dim; i++) {  //
        tp[i] =_cp[i] + n[i+m_dim]*m_dt/2;
    }
    PhVelocity(tp);
    for (int i=0; i< m_dim; i++) {
        n[i + 2*m_dim] = m_v[i];
    }

    // n4
    for (int i=0; i< m_dim; i++) {  //
        tp[i] =_cp[i] + n[i + 2*m_dim]*m_dt;
    }
    PhVelocity(tp);
    for (int i=0; i< m_dim; i++) {
        n[i+ 3*m_dim] = m_v[i];
    }

    for (int i=0; i< m_dim; i++) {
        _np[i]=_cp[i]+(m_dt/6)*(n[i]+2*n[i+m_dim]+2*n[i+2*m_dim]+n[i+3*m_dim])  ;
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
void Lorenz::GridTr (double*& _np)
{
    for (int i = 0; i<m_dim; i++) {
        _np[i] = (round(_np[i]/a))*a +a/2;
    }
}

// хэш функция
void Lorenz::HashFun (double*& _np){

    //Хэширование полученной точки
    int _hp0=_np[0]*pow(a,-1); //переходы в целочисленные значения для проведения побитовых операций
    int _hp1=_np[1]*pow(a,-1);
    int _hp2=_np[2]*pow(a,-1);

    int _hp =abs(((_hp0*p1)^(_hp1*p2)^(_hp2*p3))%h_n); //при мелкой сетке выходит за границы - пока проблема не решена


    std::list<double>::iterator it;  //Созадние итератора
    int n_n=0; //дополнительная переменная для подсчетка кол-ва совпадений координат между двух точек


    //Проверка коодинат
    if (m_hash[_hp].size() == 0 ) {  //если в листе ничего нет, то мы записываем туда три координаты
        for (int i = 0; i < m_dim; ++i) {
            m_hash[_hp].push_back(_np[i]);
        }

    }

    else { //если в листе есть что-то то мы проверям совпадения
        for (it = m_hash[_hp].begin(); //будем проверять каждый элемент листа
             it != m_hash[_hp].end();) {

            for (int k = 0; k < m_dim; ++k) {

                if (*it == _np[k]) {  //сравниваем содержимое
                    n_n = n_n + 1; //если равны, то прибавляем
                    it++; //переходим к следующему элементу
                }
                else { it++; } //переходим к следующему элементу
            }
        }
        //после проверки записываем новые элементы
        for (int i = 0; i < m_dim; ++i) {
            m_hash[_hp].push_back(_np[i]);
        }

        //Если все три координаты совпали, то отмечаем наличие цикла
        if (n_n == 3) { stat_ncycle = stat_ncycle + 1; }
    }
}




// сохранение траектории в файл
void Lorenz::Save (char* _fln,int _be, int _dn, int _en)
{
    // создание файлового потока вывода
   /* std::ofstream out;

    // связывание потока с файлом
    out.open (_fln);

    // вывод результата
    out.precision (15); // задание количества значащих цифр в выводимых числах
    for (int i = _be; i < _en; i += _dn * m_dim)
    {
        for (int j = 0; j < m_dim; j++)
        { std::cout << m_tr [i * m_dim + j] << " ";
        }
        std::cout  << "\n";
    }

 */
    std::cout  << stat_ncycle;
    // закрытие файла
    //out.close ();
}


// деструктор
Lorenz::~Lorenz ()
{
    if (m_tr != 0){
    delete [] m_tr;}

    delete [] m_v;
}
