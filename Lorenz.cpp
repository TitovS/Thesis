#include "Lorenz.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>

// конструктор
Lorenz::Lorenz (void)
{
    // размерность системы
    m_dim = 3;

    // количество точек в траектории
    m_n = 8000*m_dim;

    // хэш-значения
    p1=73856093;
    p2=19349663;
    p3=83492791;

    // Размер хэш таблицы
    h_n=100000;
    m_hash = new std::list<int> [h_n];

    // Инициализация класса статистики
    S.n_cycle=0;
    S.u_cycle=0;
    S.l_cycle =  new double [100000];
    S.s_cycle =  new int [100000];

    //Показатель перехода на следуюшую тракторию
    Cycles_done = 0;

    //количество вычисленных точек
    num_main = 0;

    vector =  new double [3];
    dot =  new double [3];
}


// задание параметров системы
void Lorenz::SetParam (double _p1, double _p2, double _p3)
{
    m_b = _p1;
    m_sigma = _p2;
    m_r = _p3;
}


// вычисление фазовой траектории
void Lorenz::GetTr (double* _init, double _dt, double _a, int _b, bool baseline)
{
    // Шаг вычислительного метода
    m_dt = _dt;

    // Параметры решетки
    a=_a;
    b=_b;

    num_first = num_main; // присвоение начальной точки данной траектории

    // Создание массивов точек фазовой траектории
    double* m_tr1 = new (std::nothrow) double [m_n];
    double* m_tr2 = 0;// = new (std::nothrow) double [m_n];
    // Указатель на акутальный масссив
    m_tr=m_tr1;

    // Массив для вчисления точек
    m_v = new double [m_dim];

    // Запись начальных условий
    memcpy (m_tr, _init, m_dim * sizeof (double));


    bool cycle = false; // Инвариант цикла - вычисление до нахождения цикла

    int k = 0; // Иднидактор динмачиеского массива 0 - m_tr1 / 1 - m_tr2 | Смотри DynamicMemory

    while (!cycle)

    {
        DynamicMemory(num_main,k,m_tr1,m_tr2); // Проверка и увеличение массива в случае необходимости

        // Указатель на текущую точку
        double* cp = m_tr + num_main * m_dim;
        // Указатель на следующую точку
        double* np = m_tr + (num_main + 1) * m_dim;

        // Вычисление следующей точки методом Рунге-Кутта
        Runge_Cutta(np, cp);

        //дискретизация при шаге на решетке
        if (num_main%b == 0 ) {
            GridTr(np);
            CycleCheck(np);}

        //Если обнаржуен цикл, то вычисление заканчиывается.
        if (S.n_cycle > Cycles_done) {
            Cycles_done = S.n_cycle;
            cycle = true;
        }

        num_main+=1; // сквозная нумерация точки
    }

    if (!baseline){
        m_tr = 0;
        m_v = 0;
    }
    if (k == 0) {
        delete[] m_tr1;
    } else {
        delete[] m_tr2;
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
        _np[i]=_cp[i]+(m_dt/6)*(n[i]+2*n[i+m_dim]+2*n[i+2*m_dim]+n[i+3*m_dim]) ;
    }

    delete [] tp;
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
int Lorenz::HashFun(double *_np) {

    //Хэширование полученной точки
    int h[m_dim];

    for (int i = 0; i < m_dim; ++i) {
        h[i]=_np[i]*1/a; //переходы в целочисленные значения для проведения побитовых операций
    }

    return abs(((h[0]*p1)^(h[1]*p2)^(h[2]*p3))%h_n); //TODO Он временами выдает отрицательные результаты, вероятно выходит за границы
    
}





void Lorenz::CycleCheck(double *_np){

    int hp = HashFun(_np); // Присвоение

    std::list<int>::iterator it;  //Созадние итератора TODO Пары можно убрать, теперь мы можем просто сохарнять номер


    if (*_np != m_tr[(num_main+1)*m_dim]) { std::cout << "error"; }
    //Проверка коодинат

    if (m_hash[hp].size() != 0 ) {   //если в листе есть что-то то мы проверям совпадения
        //идем по списку из пар которые назодятся по ключу
        for (it = m_hash[hp].begin(); it != m_hash[hp].end(); it++) {
            //сравниваем значения (гет 0)

            if ((m_tr[*it*m_dim] == *_np)) {

                if (m_tr[*it*m_dim + 1] == *(_np + 1)) {

                    if (m_tr[*it*m_dim + 2] == *(_np + 2)) { //TODO подумать над другим методом сравнения

                        if (*it > num_first) {

                            S.collect(_np, &m_tr[*it], m_dim, num_main - *it);
                        }

                        S.n_cycle += 1;
                    }
                }
            }
        }
    }
    //Запись координаты
    m_hash[hp].push_back(num_main+1);
}




// сохранение траектории в файл
void Lorenz::Save (int _be, int _dn, int _en) {
    //S.save;
    //создание файлового потока вывода
    std::ofstream out;

    // связывание потока с файлом
    out.open("tr1.txt");

    // вывод результата
    for (int i = _be; i < _en; i += _dn * m_dim) {
        for (int j = 0; j < m_dim; j++) {
          //  out << m_tr[i * m_dim + j] << " ";
        }
    //    out << "\n";
    }
// закрытие файла
    out.close();

    for (int j = 0; j <S.u_cycle; ++j) {
        std::cout  << S.s_cycle[j]<<  " " << S.l_cycle[j]<<"\n";
    }

}



void Lorenz::GetLine() {

    dot[0] = sqrt((m_r-1)*m_b) ;
    dot[1] =dot[0];
    dot[2] =m_r-1;

    vector[0] = -1*(dot[0] + dot[0]);
    vector[1] = vector[0];
    vector[2] = 0;
}


void Lorenz::DynamicMemory(int i, int& k, double*& m_tr1, double*& m_tr2 ) {

    if (i * m_dim >= m_n - m_dim) { //

        if (k == 0) {
            //delete[] m_tr2;
            m_n = m_n * 2; //увеличение размера массива в двое
            m_tr2 = new(std::nothrow) double[m_n];
            memcpy(m_tr2, m_tr1, m_n / 2 * sizeof(double));
            m_tr = m_tr2;
            delete[] m_tr1;
            m_tr1 = NULL;
            k=1;
        }

        else{
            //delete[] m_tr1;
            m_n = m_n * 2; //увеличение размера массива в двое
            m_tr1 = new(std::nothrow) double[m_n];
            memcpy(m_tr1, m_tr2, m_n / 2 * sizeof(double));
            m_tr = m_tr1;
            delete[] m_tr2;
            m_tr2 = NULL;
            k=0;
        }

    }



}
// деструктор
Lorenz::~Lorenz ()
{
    if (m_tr != 0){
        delete [] m_tr;}

    delete [] m_v;
    delete [] vector;
    delete [] dot;
    delete [] m_hash;

    m_tr = 0;
    m_v = 0;
    m_n = 8000*m_dim;
}