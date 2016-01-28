//
//  main.cpp
//  Thesis
//
//  Created by Сергей Титов on 04/05/15.
//  Copyright (c) 2015 Сергей Титов. All rights reserved.
//

#include "Lorenz.h"
#include "iostream"
#include "Stat.h"


int main (void) {

    Stat* S = new Stat; // Создание класса статистики
    Lorenz* L = new Lorenz; // Создание класса аттрактора

    double a_right = 0.065; //Величина решетки
    double a_left=0;
    double eps = 0.000000001;

    int n_last_cycles = 0;
    int num_breakpoits = 0; //для точного количества надо делить на два - в массиве храняться пары.


    L->GetCycles(S,a_right);
    n_last_cycles = S->u_cycle;
    S->Reset();// Обновление класса статистики

    a_left = a_right;
    a_right -= 0.00001;

    num_breakpoits += 2;


    for (int i = 1; i < 10; ++i) {

        L->GetCycles(S,a_right);


        if (n_last_cycles != S->u_cycle) L->GetBreak(S, a_right, a_left, eps, num_breakpoits);

        n_last_cycles = S->u_cycle;
        S->Reset();// Обновление класса статистики

        a_left = a_right;
        a_right -= 0.00001;

        num_breakpoits += 2;
    }

    for (int j = 0; j < 20; ++j) {
        std::cout << L->a_breakpoints[j] << ' ' << L->a_breakpoints[j+1] <<"\n";
        j+=1;

    }


    return (0);
}
