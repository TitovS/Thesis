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
    L.GetLine();

    double* x = new double[L.m_dim];

    memcpy (x, L.dot, L.m_dim * sizeof (double));

    double t = 0.01;

    for (int i = 0; i < 1000; ++i) {

        L.GetTr(x, 0.0065, 0.0065 ,10); // шаг метода, количество точек , величина решетки. шаг на решетке

        for (int j = 0; j < L.m_dim; ++j) {
            x[j] = L.vector[j] * t - L.dot[j];
        }

        t+=0.0065;

    }

    L.Save(0, 3, 100000);

    return (0);
}
