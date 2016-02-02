//
//  main.cpp
//  Thesis
//
//  Created by Сергей Титов on 04/05/15.
//  Copyright (c) 2015 Сергей Титов. All rights reserved.
//

#include "Lorenz.h"
#include "iostream"
#include "fstream"
#include "Stat.h"


int main (void) {

    Stat* S = new Stat; // Создание класса статистики
    Lorenz* L = new Lorenz; // Создание класса аттрактора

    double a_right = 0.065; //Величина решетки
    double a_left=0;
    double eps = 0.000000001;

    int n_last_cycles = 0;



    L->GetCycles(S,a_right);
    n_last_cycles = S->u_cycle;
    S->Reset();// Обновление класса статистики

    a_left = a_right;
    a_right -= 0.00001;

    //L->num_breakpoints += 4;


    for (int i = 1; i < 10; ++i) {

        L->GetCycles(S,a_right);


        if (n_last_cycles != S->u_cycle) {

            L->GetBreak(S, a_right, a_left, eps);
            L->num_breakpoints += 4;
        }


        n_last_cycles = S->u_cycle;


        a_left = a_right;
        a_right -= 0.00001;

        L -> a_breakpoints[L->num_breakpoints] = a_right;
        L -> a_breakpoints[L->num_breakpoints] = n_last_cycles;
        L->  num_breakpoints +=2;
        S->Reset();// Обновление класса статистики
    }

    std::ofstream out;
    char const *pchar = "MainTable";
    out.open(pchar, std::ofstream::out | std::ofstream::app);



    for (int j = 0; j < L-> num_breakpoints/2; ++j) {
        std::cout << L->a_breakpoints[j] << ' ' << L->a_breakpoints[j+1] <<"\n";
        out << L->a_breakpoints[j]<< ' ' << L->a_breakpoints[j+1] <<"\n";
    }

    out.close();
    return (0);
}
