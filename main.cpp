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

    double x[] = {3.2445, 3.7375, 15.6975}; // начальные условия
    double t = 0.01;

    for (int i = 0; i < 200; ++i) {

        L.GetTr(x, 0.01, 265000, 0.1 ,10); // шаг метода, количество точек , величина решетки. шаг на решетке

        x[0]-= 23.2*t - 6.55;
        x[1]-= 30.8*t - 10.45;
        x[2]-= 1.78*t + 1.675;

        t+=0.01;

    }

    L.Save(0, 3, 100000);

    return (0);
}
