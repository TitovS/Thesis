//
// Created by Сергей Титов on 14/02/16.
//

#include "Grid.h"
#include <iostream>
#include <fstream>
#include <iomanip>

Grid::Grid() {

    grid_results = new (std::nothrow) double[10000000]; //TODO может имеет смысл сделать вектором для динмаичегого размера

    grid_num=0;

    a_right = 0;

    a_left=0.065;

    eps = 0.000000001;

    step = 0.00001;


}

void Grid::Grid_make_step() {

    a_right=a_left;
    a_left-=step;

}

void Grid::Save(int num_cycles){

    grid_results [grid_num] = a_left; // Размер решетки
    grid_results [grid_num+1] = num_cycles; // Соответсвующее количество циклов

    grid_num+=2;

}

void Grid::Save_in_file() {

    std::ofstream out;
    char const *pchar = "MainTable";
    out.open(pchar, std::ofstream::out | std::ofstream::app);


    for (int j = 0; j < grid_num/2; ++j) {
        std::cout << std::setprecision(10)<<  grid_results[j] << ' ' <<  grid_results[j+1] <<"\n";
        out <<  grid_results[j]<< ' ' <<  grid_results[j+1] <<"\n";

        j+=1;
    }

    out.close();
}

void Grid::Grid_step_est() {



}


Grid::~Grid() {

    delete[]grid_results;

}