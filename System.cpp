#include "System.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <string>
#include <cstring> 



    // конструктор
    System::System(void)

    {   std::ifstream fin("init/SysInit.txt"); // opens the text file
        std::string name;

         fin >> name >> steps;

        //Размерность системы
        fin >> name >> m_dim;

        //Характеристики системы
        fin >> name >> A;
        fin >> name >> B;

        //Хэш Значения
        fin >> name >> p1;
        fin >> name >> p2;
        fin >> name >> p3;

        //Начальное количество точек в траектории
        fin >> name >> m_n;

        // Размер хэш таблицы
        fin >> name >> h_n;


        fin.close();

        m_n = m_n*m_dim;


        m_hash = new std::list<int> [h_n];

        //количество вычисленных точек
        num_main = 0;
        num_first = 0;

        //Создание динмаического массива траекторий см. DynamicMemory
        m_tr1 = new (std::nothrow) double [m_n];
        m_tr2 = 0;
        k = 0;
        m_tr = m_tr1;

        //Прямая для вычисления всех значений
        vector =  new double [m_dim];
        dot =  new double [m_dim];

        //Моды
        BreakpointSearchMode = false;
    }



// Вычисление циклов

    // Получение координат для построения прямой
    void System::GetLine() {


        //Lorenz
        // dot[0] = sqrt((m_r-1)*m_b) ;
        // dot[1] =dot[0];
        // dot[2] =m_r-1;

        // vector[0] = -1*(dot[0] + dot[0]);
        // vector[1] = vector[0];
        // vector[2] = 0;



        dot[0] = 0;
        dot[1] = sqrt(B);
        dot[2] = B;

        vector[0] = 0;
        vector[1] = -2*sqrt(B);
        vector[2] = 0;
    }

    // Вычисление всех циклов при заданной дискретизации
    void System::GetCycles(Stat *S, double a) {

        //Переменные времени и начальной точки
        double t = 0.;
        double *x = new double[m_dim];

        // Получаем одну траекторию и создаем прямую, для находждения всех циклов
        GetLine();
        memcpy(x, dot, m_dim * sizeof(double));

        double k_max;

        k_max = round(3.25/a);

        for (int k = 0; k < k_max; ++k) { // 200 шагов - при данном t повзовлет пройти всю прямую

            GetTr(x, a, 10, S); // Начальные данные, величина решетки, шаг на решетке, класс статистики

            for (int j = 0; j < m_dim; ++j) { // Следующий шаг
                x[j] = dot[j] + vector[j] * t;
            }

            t+=a;
        }

        //Запись идет только в том случае, если не поиск разрыва
        if (!BreakpointSearchMode) {
            //Save();
            //std::cout << '=' << std::flush;
        }

        delete[] x;
    };

    // Вычисление фазовой траектории
    void System::GetTr (double* _init, double _a, int _b, Stat* S) {
        // Шаг вычислительного метода
        m_dt = (_a*_b);

        // Параметры решетки
        a=_a;


        num_first = num_main; // присвоение начальной точки данной траектории
        Cycles_done = S -> n_cycle; // запись количества выполенных циклов статистику

        // Массив для вчисления точек
        m_v = new double [m_dim];

        // Запись начальных условий
        GridTr(_init);
        memcpy (&m_tr[num_first*m_dim], _init, m_dim * sizeof (double));

        bool cycle = false; // Инвариант цикла - вычисление до нахождения цикла


        while (!cycle)

        {
            // Проверка массива
            DynamicMemory(num_main, k, m_tr1, m_tr2);

            // Указатель на текущую точку
            double* cp = m_tr + num_main * m_dim;

            // Указатель на следующую точку
            double* np = m_tr + (num_main + 1) * m_dim;

            // Вычисление следующей точки методом Рунге-Кутта
            Runge_Cutta(np, cp);

            //дискретизация при шаге на решетке
            if (num_main%1 == 0 ) { //TODO каждые сколько шагов проверять?
                GridTr(np);
                CycleCheck(np, S);}

            //Если обнаржуен цикл, то вычисление заканчивается
            if (S -> n_cycle > Cycles_done) {
                Cycles_done = S -> n_cycle;
                cycle = true;
            }
            // сквозная нумерация точки
            num_main+=1;

        }

        delete[] m_v;
    }

    // Поиск разрыва
    void System::GetBreak(Stat *S, Grid *A) {
        //TODO сделать проверку на новые циклы при поиске разрыва
        // Маркер поиска разрыва

        // Отчистка старой статистики
        S->Reset();

        // Создание локальных перменных - кол-во циклов для середин и крайней точки
        double amid = (A->a_right-A->a_left)/2;
        int num_cycles_right=0;
        int num_cycles_middle=0;

        // Сохранение старого шага
        double oldgrid = A->a_left;


        // Цикл поиска точки разрыва с точностью eps
        while (amid > A->eps){

            BreakpointSearchMode = true;

            //Поиск правой части
            GetCycles(S,A->a_right);
            num_cycles_right = S->u_cycle;
            S->Reset();

            BreakpointSearchMode = false;

            //Поиск центральной части
            GetCycles(S,(A->a_left + amid));
            num_cycles_middle = S->u_cycle;
            S->Reset();

            //Сравнение
            if (num_cycles_right == num_cycles_middle) A->a_right -= amid;

                //Сдвиг левой части на центральную
            else {
                A->a_left += amid;
                // Если центральная часть отличается от старых реультатов, то сохраняется
                //if (num_cycles_middle != A->grid_results[A->grid_num-1])
                    //A->Save(num_cycles_middle);


            }

            amid = (A->a_right-A->a_left)/2;
        }

        // A->Save(num_cycles_middle);
        BreakpointSearchMode = false;

        A->a_left = oldgrid;

    }

// Вычисление траекторий

    // Вычисление следующей точки фазовой траектории по текущей - Рунге-Кутта
    void System::Runge_Cutta (double*& _np, double*& _cp ) {
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

    // Вычисление фазовой скорости в заданной точке
    void System::PhVelocity (double*& _cp) {


        m_v [0] = -1*A*_cp [0] + B*_cp [1] - (_cp [2]*_cp [1]);
        m_v [1] = _cp [0];
        m_v [2] = -1*_cp [2]+(_cp[1]*_cp[1]);

    //Lorenz
        //m_v [0] = m_sigma * (_cp [1] - _cp [0]);
        //m_v [1] = (m_r - _cp [2]) * _cp [0] - _cp [1];
        //m_v [2] = _cp [0] * _cp [1] - m_b * _cp [2];

    //Duffing
        // m_v [0] = _cp [1];
        // m_v [1] = _cp [0]*_cp [0]*_cp [0] - m_r*_cp [1]+m_b*sin(_cp[2]);
        // m_v [2] = 1;

        double X = sqrt(m_v[0]*m_v[0] + m_v[1]*m_v[1] +m_v[2]*m_v[2]);

        m_v [0] = m_v[0]/X;
        m_v [1] = m_v[1]/X;
        m_v [2] = m_v[2]/X;
    }

// Проверка цикла

    // Дискретизация  методом Рану
    void System::GridTr (double*& _np) {
        for (int i = 0; i< m_dim; i++) {
            _np[i] = (floor(_np[i]/a))*a +a/2;
        }
    }

    // Хэш функция
    int System::HashFun(double *_np) {

        //Хэширование полученной точки
        int h[m_dim];

        for (int i = 0; i < m_dim; ++i) {
            h[i]=_np[i]*1/a; //переходы в целочисленные значения для проведения побитовых операций
        }

        return abs(((h[0]*p1)^(h[1]*p2)^(h[2]*p3))%h_n); //TODO Он временами выдает отрицательные результаты, вероятно выходит за границы

    }

    // Проверка на обнаружение цикла
    void System::CycleCheck(double *_np, Stat* S){

        int hp = HashFun(_np); // Присвоение

        std::list<int>::iterator it;  // Созадние итератора

        // Проверка коодинат

        if (m_hash[hp].size() != 0 ) {   // Если в листе есть что-то, то мы проверям совпадения
            // Идем по листу, который назодятся по ключу
            for (it = m_hash[hp].begin(); it != m_hash[hp].end(); it++) {
                //double s = m_hash[hp].size();
                //std::cout << s;
                if (m_tr[*it*m_dim] == *_np) {
                    //std::cout << m_tr[*it*m_dim] << ' ' << *_np <<' ' << '1' <<"\n";
                    if (m_tr[*it*m_dim + 1] == *(_np + 1)) {
                      //  std::cout << m_tr[*it*m_dim+1] << ' ' << *(_np+1)<<' ' << '2' << "\n";
                        if (m_tr[*it*m_dim + 2] == *(_np + 2)) {
                        //    std::cout << m_tr[*it*m_dim+2] << ' ' << *(_np+2) <<' ' << '3'<< "\n";
                            S -> n_cycle += 1; // Цикл обнуаржуен

                            if (*it > num_first) { // Проверка на новизну цикла
                                S -> collect(&m_tr[(*it)*m_dim], m_dim, num_main + 1 - *it, a, BreakpointSearchMode); // Ссылка на начальную точку цикла, размерность, размер цикла, размер решетки
                            }
                        }
                    }
                }
            }
        }

        //Запись начальной координаты цикла в лист (+1 - нумерация акутальной точки выставляется в конце gettr)
        m_hash[hp].push_back(num_main+1);
    }

// Технические функции

    // Реалзация динамического массива для траектории
    void System::DynamicMemory(int i, int& k, double*& m_tr1, double*& m_tr2 ) {

        if (i * m_dim >= m_n - m_dim) { //

            if (k == 0) {
                //delete[] m_tr2;
                m_n = m_n * 2; //увеличение размера массива в двое
                m_tr2 = new(std::nothrow) double[m_n];
                memcpy(m_tr2, m_tr1, m_n / 2 * sizeof(double));
                m_tr = m_tr2;
                delete[] m_tr1;
                m_tr1 = NULL;
                k = 1;
            }

            else {
                //delete[] m_tr1;
                m_n = m_n * 2; //увеличение размера массива в двое
                m_tr1 = new(std::nothrow) double[m_n];
                memcpy(m_tr1, m_tr2, m_n / 2 * sizeof(double));
                m_tr = m_tr1;
                delete[] m_tr2;
                m_tr2 = NULL;
                k = 0;
            }
        }
    }

    // Обнуление параметров
    void System::Reset(){

        // Размер хэш таблицы
        h_n=100000;  //TODO динмаический хэш.
        delete[] m_hash;
        m_hash = new std::list<int> [h_n];

        num_main = 0;

        if (k == 0 ){

            memset(m_tr1 ,0, m_n* sizeof(double));
            m_tr = m_tr1;
            delete[] m_tr2;
            m_tr2 = NULL;
        }

        else {

            m_tr1 = new(std::nothrow) double[m_n];
            memset(m_tr1 ,0, m_n* sizeof(double));
            m_tr = m_tr1;
            delete[] m_tr2;
            m_tr2 = NULL;
            k = 0;

        }

        memset(vector,0,3* sizeof(double));
        memset(dot,0,3* sizeof(double));

    }

    // Сохранение траектории в файл
    void System::Save() {
        //S.save;
        //создание файлового потока вывода
        std::ofstream out;

        // связывание потока с файлом
        out.open("tr1.txt");

        // вывод результата
        for (int i = 0; i < num_main; i +=  m_dim) {
            for (int j = 0; j < m_dim; j++) {
                out << m_tr[i * m_dim + j] << " ";
            }
            out << "\n";
        }
        // закрытие файла
        out.close();

        /*for (int j = 0; j < S -> u_cycle; ++j) { //TODO разбораться с тем как закрывать
            if (S->s_cycle[j] > 20) {
                std::cout << S->s_cycle[j] << " " << S->l_cycle[j] << "\n";
            }
        }*/
    }

    // деструктор
    System::~System(){

        if (m_tr1 != 0){
            delete[] m_tr1;}

        if (m_tr2 != 0){
            delete[] m_tr2;}

        delete[] vector;
        delete[] dot;
        delete [] m_hash;

    }

    // дикретизация полученной траектории методом Цпвд
    /*void System::GridTrCPVD (double*& _np,double*& _cp)
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


    }*/
