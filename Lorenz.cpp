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

    // вектор фазовой скорости


    // количество точек в траектории
    m_n = 8000*m_dim;

    // хэш-значения
    p1=73856093;
    p2=19349663;
    p3=83492791;

    // Размер хэш таблицы
    h_n=100000;
    m_hash = new std::list<std::pair<double*,int>> [h_n];

    // Инициализация класса статистики
    S.n_cycle=0;
    S.u_cycle=0;
    S.l_cycle =  new double [100000];
    S.s_cycle =  new int [100000];

    //Показатель перехода на следуюшую тракторию
    Cycle_cheсk = 0;

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
    // шаг метода
    m_dt = _dt;

    // Шаг решетки
    a=_a;
    b=_b;

    num_first = num_main; // начальная точка вычисления

    // создание массива точек фазовой траектории
    double* m_tr1 = new (std::nothrow) double [m_n]; // Поменять способ задавать массив
    double* m_tr2 = new (std::nothrow) double [m_n];
    m_tr=m_tr1;

    m_v = new double [m_dim];

    // запись начальных условий
    memcpy (m_tr, _init, m_dim * sizeof (double));


    // вычисление траектории
    //for (int i = 0; i < m_n-1 ; i++)
    bool unique = true;
    int i = 0;
    int k = 0;

    while (unique)

    {

        if (i * m_dim >= m_n - m_dim) { //

            if (k == 0) {
                m_n = m_n * 2; //увеличение размера массива в двое
                m_tr2 = new(std::nothrow) double[m_n];
                memcpy(m_tr2, m_tr1, m_n / 2 * sizeof(double));
                m_tr = m_tr2;
                delete[] m_tr1;
                k=1;
            }

            else{
                m_n = m_n * 2; //увеличение размера массива в двое
                m_tr1 = new(std::nothrow) double[m_n];
                memcpy(m_tr1, m_tr2, m_n / 2 * sizeof(double));
                m_tr = m_tr1;
                delete[] m_tr2;
                k=0;
            }

        }

        // указатель на текущую точку
        double* cp = m_tr + i * m_dim;
        // указатель на следующую точку
        double* np = m_tr + (i + 1) * m_dim;

        // вычисление следующей точки методом Рунге-Кутта и Дискретизация полученной траектории
        Runge_Cutta(np, cp);

        //дискретизация при шаге на решетке
        if (i%b == 0 ) {
            GridTr(np);
            HashFun(np);}

        //Если обнаржуен цикл, то вычисление заканчиывается.
        if (S.n_cycle > Cycle_cheсk) {
            Cycle_cheсk = S.n_cycle;
            unique = false;
        }

        num_main+=1; // сквозная нумерация точки
        i+=1; // счетчик по массиву


        //Проверка переполнености массива



    }

    if (!baseline){
        m_tr = 0;
        m_v = 0;
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
void Lorenz::HashFun (double* _np){

    //Хэширование полученной точки
    long h[m_dim];

    for (int i = 0; i < m_dim; ++i) {
        h[i]=_np[i]*1/a; //переходы в целочисленные значения для проведения побитовых операций
    }

    long hp =abs(((h[0]*p1)^(h[1]*p2)^(h[2]*p3))%h_n); //TODO Он временами выдает отрицательные результаты, веростно выходит за границы

    std::list<std::pair<double*,int>>::iterator it;  //Созадние итератора

    //Проверка коодинат
    if (m_hash[hp].size() != 0 ) {   //если в листе есть что-то то мы проверям совпадения
        //идем по списку из пар которые назодятся по ключу
        for (it = m_hash[hp].begin(); it != m_hash[hp].end(); it++) {
            //сравниваем значения (гет 0)
            if (memcmp(std::get<0>(*it), _np, m_dim*sizeof(_np)) == 0)
            {
                //проверяем, что это новый цикл
                if(std::get<1>(*it) > num_first){

                    S.collect(_np, &m_tr[std::get<1>(*it)], m_dim, num_main - std::get<1>(*it));
                }

                S.n_cycle+=1;
            }
        }
    }
    //Запись координаты
    m_hash[hp].push_back(std::make_pair(_np,num_main));
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

    double x[] = {1,1,1};

    GetTr(x, 0.065, 0.065 ,10, true);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> distribution(1, num_main);
    int r = distribution(generator);


    memcpy (dot, &m_tr[(num_main- 1) * m_dim],m_dim * sizeof (double));

    for (int i = 0; i < m_dim; ++i) {
    //    dot[i] =m_tr[(num_main- 1) * m_dim + i];
        vector[i] = m_tr[(num_main- 1) * m_dim + i - r] - m_tr[(num_main- 1) * m_dim + i];
    }

    delete [] m_hash;
    num_main = 0;
    m_hash = new std::list<std::pair<double*,int>> [h_n];

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