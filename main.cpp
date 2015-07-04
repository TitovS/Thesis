//
//  main.cpp
//  Thesis
//
//  Created by Сергей Титов on 04/05/15.
//  Copyright (c) 2015 Сергей Титов. All rights reserved.
//

#include "Lorenz.h"


int main (void) {
    Lorenz L;

    L.SetParam(8. / 3., 10, 28); // b, sigma, r

    double x[] = {1, 1, 1}; // начальные условия

    L.GetTr(x, 0.07, 265000, 0.065 , 10); // шаг метода, количество точек , величина решетки. шаг на решетке

    char* t= "tr1"; //Имя файла

    L.Save(0, 3, 100000); //  начальная точка, интервал пропуска, конечная точка

    return (0);
}
