//
// Created by Сергей Титов on 14/02/16.
//

#ifndef THESIS_GRID_H
#define THESIS_GRID_H


class Grid {

public:

    double*grid_results; // Результирующий массив

    int grid_num; // Количество точек всего

    double a_left; // Актуальный размер решетки

    double a_right; // Прошлый размер решетки

    double step; // Шаг изменения решетки

    double eps; // Предельная точность поиска разрыва

//___________________________________

    Grid(void); // Инициализациия

    void Save(int); // Сохранение актуальной точки

    void Grid_step_est(); // Подсчет оптимального шага решетки

    void Grid_make_step(); // Уменьшение решетки на один шаг

    void Save_in_file(); // Сохранение результатов в файл

    ~Grid(); // Деструктор

};


#endif //THESIS_GRID_H
