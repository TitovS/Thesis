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

    double x[] = {2.2425, 3.7375, 15.6975}; // начальные условия

    for (int i = 0; i < 10; ++i) {

       L.GetTr(x, 0.065, 265000, 0.065 , 10); // шаг метода, количество точек , величина решетки. шаг на решетке

        L.Save(0, 3, 100000); //  начальная точка, интервал пропуска, конечная точка

        x[0]-= 0.01;

    }
    return (0);
}
