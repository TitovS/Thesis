//
//  main.cpp
//  Thesis
//
//  Created by Сергей Титов on 04/05/15.
//  Copyright (c) 2015 Сергей Титов. All rights reserved.
//

#include "Lorenz.h"
//#include "iostream"
//#include "fstream"
//#include "Stat.h"
//#include "Grid.h"


int main (void) {

    Stat* S = new Stat; // Создание класса статистики
    Lorenz* L = new Lorenz; // Создание класса аттрактора
    Grid* A= new Grid; // Создание класса решетки

    int n_last_cycles = 0; // Количество циклов в прошлой итерации


    L->GetCycles(S, A->a_left); // Количество циклов при нынешней
    n_last_cycles = S->u_cycle; // Запись
    A->Save(n_last_cycles);// Сохранение результатов

    S->Reset(); // Обнуление класса статистики
    A->Grid_make_step(); // Переход на следующую решетку

    for (int i = 1; i < 10; ++i) {

        L->GetCycles(S,A->a_left); // Поиск всех циклов
        A->Save(S->u_cycle); // Запись реузльтатов
        if (n_last_cycles != S->u_cycle) L->GetBreak(S,A); // Если количество циклов не совпадает, то ищем точку разрыва

        n_last_cycles = S->u_cycle; // Возвращаемся на движение по прямой

        A->Grid_make_step(); // Переход на следующую решеткуПереход на следующую решетку
        S->Reset();// Обнуление класса статистики
    }

    A->Save_in_file();

    delete A;
    delete L;
    delete S;

    return (0);
}
