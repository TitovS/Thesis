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
    Lorenz L;

    L.SetParam(8. / 3., 10, 28); // b, sigma, r
    L.GetLine();

    double* x = new double[L.m_dim];

    memcpy (x, L.dot, L.m_dim * sizeof (double));

    double t = 0.;

    for (int i = 0; i < 1500; ++i) {

        L.GetTr(x, 0.065, 0.065 ,10); // шаг метода, количество точек , величина решетки. шаг на решетке

        for (int j = 0; j < L.m_dim; ++j) {
            x[j] = L.dot[j] + L.vector[j] * t;

        }

        std::cout << x[0] << " "<< std::endl;
        t+=0.00065;

    }

    L.Save(0, 3, 100000);

    return (0);
}
