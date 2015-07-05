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


    // количество точек в траектории
    m_n = 0;

    m_tr=0;

    // хэш-значения
    p1=73856093;
    p2=19349663;
    p3=83492791;

    h_n=10000;

    S.n_cycle=0;
    S.l_cycle =  new double [100000];
    S.s_cycle =  new int [100000];

    c_old =0;
}


// задание параметров системы
void Lorenz::SetParam (double _p1, double _p2, double _p3)
{
    m_b = _p1;
    m_sigma = _p2;
    m_r = _p3;
}


// вычисление фазовой траектории
void Lorenz::GetTr (double* _init, double _dt, int _n, double _a, int _b)
{
    // шаг метода
    m_dt = _dt;

    // количество вычисляемых точек
    m_n = _n;

    // Шаг решетки
    a=_a;
    b=_b;

    // создание массива точек фазовой траектории
    m_tr = new (std::nothrow) double [m_dim * m_n];
    m_v = new double [m_dim];

    // запись начальных условий
    for (int i = 0; i < m_dim; i++)
        m_tr [i] = _init [i];

    // вычисление траектории
    for (int i = 0; i < m_n-1 ; i++)
    {
        // указатель на текущую точку
        double* cp = m_tr + i * m_dim;
        // указатель на следующую точку
        double* np = m_tr + (i + 1) * m_dim;

        // вычисление следующей точки методом Рунге-Кутта и Дискретизация полученной траектории
        Runge_Cutta(np, cp);

        //дискретизация при шаге на решетке
        if (i%b == 0) {
            GridTr(np);
            HashFun(np);}
        //Если обнаржуен цикл, то вычисление заканчиывается.
        if (S.n_cycle > c_old) {

            c_old = S.n_cycle;
            std::cout  << " " << S.s_cycle[S.n_cycle]<<  " " << S.l_cycle[S.n_cycle]<<"\n ";
            break;
            }

    }

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


void Lorenz::GridTr (double*& _np)
{
    for (int i = 0; i<m_dim; i++) {
        _np[i] = (round(_np[i]/a))*a +a/2;
    }
}

// дикретизация полученной траектории методом Цпвд
void Lorenz::GridTrCPVD (double*& _np,double*& _cp)
{

    double* x = new double [m_dim]; //Доп. переменная для создания точек
    double* xr = new double [m_dim]; //Доп. переменная для результата точек

    for (int i = 0; i < m_dim; ++i) {
        x[i]=_cp[i];
    }
    // Создания точек на решетке вокруг искомой
    for (int l = -1; l < 2; ++l) {
        for (int i = -1; i <2; ++i) {
            for (int j = -1; j <2; ++j) {


                x[0] = _cp[0] +l*a;
                x[1] = _cp[1] +i*a;
                x[2] = _cp[2] +j*a;

                Runge_Cutta(x, _cp); //вычисление траекторый

                for (int k = 0; k < m_dim; ++k) {
                    xr[k] += x[k];  // сумма по каждой координате
                }
            }
        }
    }

    for (int j = 0; j < m_dim; ++j) {
        _np[j] = round ((xr[j]/(m_dim*m_dim*m_dim))/a)*a +a/2;
    }


}

// Хэш функция
void Lorenz::HashFun (double* _np){

    //Хэширование полученной точки
    int h[m_dim];
    for (int i = 0; i < m_dim; ++i) {
        h[i]=_np[i]*1/a; //переходы в целочисленные значения для проведения побитовых операций
    }

    int hp =abs(((h[0]*p1)^(h[1]*p2)^(h[2]*p3))%h_n); //при мелкой сетке выходит за границы - пока проблема не решена


    std::list<double*>::iterator it;  //Созадние итератора
    int n_n=0; //дополнительная переменная для подсчетка кол-ва совпадений координат между двух точек
    int io = 0;
    //Проверка коодинат
    if (m_hash[hp].size() == 0 ) {  //если в листе ничего нет, то мы записываем туда три координаты
        m_hash[hp].push_back(_np);
    }
    else { //если в листе есть что-то то мы проверям совпадения
        for (it = m_hash[hp].begin(); it != m_hash[hp].end(); it++) {

            n_n=memcmp(*it, _np, m_dim*sizeof(_np));
            if (n_n == 0) {
                S.collect(_np, *it, m_dim);
            }
        }
        if (S.s_cycle[S.n_cycle] == 0)
        {m_hash[hp].push_back(_np);}
    }
}




// сохранение траектории в файл
void Lorenz::Save (int _be, int _dn, int _en) {
    //S.save;
    //создание файлового потока вывода
    std::ofstream out;

    // связывание потока с файлом
    out.open("tr1.txt");

    // вывод результата
    //out.precision (15); // задание количества значащих цифр в выводимых числах
    for (int i = _be; i < _en; i += _dn * m_dim) {
        for (int j = 0; j < m_dim; j++) {
            out << m_tr[i * m_dim + j] << " ";
        }
        out << "\n";
    }
// закрытие файла
    out.close();

//    for (int k = 1; k < S.n_cycle+1; ++k) {
//        std::cout << k << " " << S.s_cycle[k]<<  " " << S.l_cycle[k]<<"\n ";
//    }

    m_tr = 0;
    m_v = 0;
    m_hash.clear();

//    typedef std::map<int, std::list<double*>>::const_iterator MapIterator;
//    for (MapIterator iter = m_hash.begin(); iter != m_hash.end(); iter++)
//    {
//        std::cout << "Key: " << iter->first << "   " << "Values:" << " ";
//        typedef std::list<double*>::const_iterator ListIterator;
//        for (ListIterator list_iter = iter->second.begin(); list_iter != iter->second.end(); list_iter++)
//        {std::cout << " " << **list_iter ;}
//        std::cout << "\n ";
//}
}


// деструктор
Lorenz::~Lorenz ()
{
//    if (m_tr != 0){
//        delete [] m_tr;}
  //  delete [] m_v;
    m_hash.clear();
}
