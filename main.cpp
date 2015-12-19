//
//  main.cpp
//  Thesis
//
//  Created by Сергей Титов on 04/05/15.
//  Copyright (c) 2015 Сергей Титов. All rights reserved.
//

#include "Lorenz.h"
#include "iostream"


int main (void) {

    Stat* S = new Stat; // Создание класса статистики
    Lorenz* L = new Lorenz; // Создание класса аттрактора

    double a = 0.065; //Величина решетки
    double eps = 0.001;
    double *x = new double[L->m_dim];



    for (int i = 0; i < 50; ++i) {

        double t = 0.; // Время

        L->GetLine(); // Получаем одну траекторию и создаем пряму, для находждения всех циклов
        memcpy(x, L->dot, L->m_dim * sizeof(double));

        for (int k = 0; k < 200; ++k) { // 200 шагов - при данном t повзовлет пройти всю прямую

            L->GetTr(x, a, 10, S); // Начальные данные, величина решетки, шаг на решетке, класс статистики

            for (int j = 0; j < L->m_dim; ++j) { // Следующий шаг
                x[j] = L->dot[j] + L->vector[j] * t;
            }

            t += 0.065;
        }

        L->Save(S, i);// Сохранение траектории и вывод результатов
        std::cout << "\n";

        L->Reset();// Обновление класса аттрактора
        S->Reset();// Обновление класса статистики
        a-=0.001;
    }

    return (0);
}
